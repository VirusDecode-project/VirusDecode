package VirusDecode.backend.service;

import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.repository.UserRepository;
import jakarta.transaction.Transactional;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.scheduling.annotation.Scheduled;
import org.springframework.security.crypto.password.PasswordEncoder;
import org.springframework.stereotype.Service;

import java.time.LocalDateTime;
import java.time.temporal.ChronoUnit;
import java.util.List;
import java.util.Optional;
@Log4j2
@Service
public class UserService {
    private final JsonDataService jsonDataService;
    private final HistoryService historyService;
    private final UserRepository userRepository;
    private final PasswordEncoder passwordEncoder;

    @Autowired
    public UserService(JsonDataService jsonDataService, HistoryService historyService,
                       UserRepository userRepository, PasswordEncoder passwordEncoder) {
        this.jsonDataService = jsonDataService;
        this.historyService = historyService;
        this.userRepository = userRepository;
        this.passwordEncoder = passwordEncoder;
    }


    public User findUserByLoginId(String loginId) {
        return userRepository.findByLoginId(loginId);
    }

    public User findUserByUserId(Long userId) {
        return userRepository.findById(userId).orElse(null);
    }

    @Transactional
    public User createUser(SignUpDto signUpDto, String role) {
        User newUser = new User();
        newUser.setFirstName(signUpDto.getFirstName());
        newUser.setLastName(signUpDto.getLastName());
        newUser.setLoginId(signUpDto.getLoginId());
        newUser.setPassword(passwordEncoder.encode(signUpDto.getPassword()));
        newUser.setRole(role);

        return userRepository.save(newUser);
    }

    public boolean checkPassword(User user, String password) {
        return passwordEncoder.matches(password, user.getPassword());
    }

    // userId로 유저 객체를 반환
    public Optional<User> getUserById(Long userId) {
        return userRepository.findById(userId);
    }

    public Long getUserIdByLoginId(String loginId) {
        User user = userRepository.findByLoginId(loginId);
        return user != null ? user.getId() : null;  // user가 있으면 userId 반환, 없으면 null 반환
    }


    /*
        ┌───────────── 초 (0 - 59)
        │ ┌───────────── 분 (0 - 59)
        │ │ ┌───────────── 시간 (0 - 23)
        │ │ │ ┌───────────── 일 (1 - 31)
        │ │ │ │ ┌───────────── 월 (1 - 12)
        │ │ │ │ │ ┌───────────── 요일 (0 - 7) (일요일=0 또는 7)
        │ │ │ │ │ │
        * * * * * *
     */
//    @Scheduled(cron = "0 */5 * * * ?")  // 5분마다 실행
    @Scheduled(cron = "0 0 3 * * ?")  // 매일 오전 3시에 실행
    @Transactional
    public void deleteGuestUsers() {
        List<User> guestUsers = userRepository.findUsersByRole("GUEST");
        int guest_cnt = 0;
        for (User user : guestUsers) {
            if (user.getCreatedAt().isBefore(LocalDateTime.now().minus(24, ChronoUnit.HOURS))) {
                guest_cnt++;
                Long userId = user.getId();
                List<String> historyList = historyService.getHistoryNamesByUserId(userId);
                for(String historyName : historyList){
                    History history = historyService.getHistory(historyName, userId);
                    jsonDataService.deleteJsonData(history);
                    historyService.deleteHistory(historyName, userId);
                }
                userRepository.deleteUserById(user.getId());
            }
        }

        if (guest_cnt > 0) {
            log.info("Deleted {} GUEST users who were created more than 24 hours ago.", guest_cnt);
        }
    }

    public void copySampleHistoriesToNewUser(User newUser) {
        Long guestUserId = getUserIdByLoginId("Guest");
        List<String> guestHistoryNames = historyService.getHistoryNamesByUserId(guestUserId);

        for (String historyName : guestHistoryNames) {
            History history = historyService.getHistory(historyName, guestUserId);
            JsonData originalJsonData = jsonDataService.getJsonData(history);
            if (originalJsonData != null) {
                History newHistory = new History();
                newHistory.setHistoryName(historyName);
                newHistory.setUser(newUser);
                historyService.createHistory(newHistory);

                JsonData newJsonData = new JsonData();
                newJsonData.setReferenceId(originalJsonData.getReferenceId());
                newJsonData.setAlignment(originalJsonData.getAlignment());
                newJsonData.setLinearDesign(originalJsonData.getLinearDesign());
                newJsonData.setPdb(originalJsonData.getPdb());
                newJsonData.setHistory(newHistory);
                jsonDataService.saveJsonData(newJsonData);
            }
        }
    }
}
