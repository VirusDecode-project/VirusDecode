package VirusDecode.backend.controller;
import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.dto.UserLoginDto;
import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.service.HistoryService;
import VirusDecode.backend.service.JsonDataService;
import VirusDecode.backend.service.UserService;
import com.google.gson.Gson;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import jakarta.servlet.http.HttpSession;

import java.util.HashMap;
import java.util.List;
import java.util.Map;


@RestController
@RequestMapping("/api/auth")
public class UserController {
    private final UserService userService;
    private final HistoryService historyService;
    private final JsonDataService jsonDataService;

    @Autowired
    public UserController(UserService userService, HistoryService historyService, JsonDataService jsonDataService){
        this.userService = userService;
        this.historyService = historyService;
        this.jsonDataService = jsonDataService;
    }

    @PostMapping("/login")
    public ResponseEntity<String> login(@RequestBody UserLoginDto loginDto, HttpSession session) {
        User user = userService.findUserByLoginId(loginDto.getLoginId());
        session.setMaxInactiveInterval(3600);

        if (user == null) {
            return ResponseEntity.status(401).body("유효하지 않은 ID 입니다.");
        }

        if (userService.checkPassword(user, loginDto.getPassword())) {
            session.setAttribute("userId", user.getId());
            return ResponseEntity.ok("User logged in successfully.");
        } else {
            return ResponseEntity.status(401).body("비밀번호가 틀렸습니다.");
        }
    }
    @PostMapping("/signup")
    public ResponseEntity<String> signup(@RequestBody SignUpDto signupDto) {
        if (userService.findUserByLoginId(signupDto.getLoginId()) != null) {
            return ResponseEntity.status(400).body("이미 존재하는 ID 입니다.");
        }

        User newUser = userService.createUser(signupDto, "USER");
        return ResponseEntity.ok("User created successfully with ID: " + newUser.getId());
    }

    @PostMapping("/guest-login")
    public ResponseEntity<String> guestLogin(HttpSession session) {
        session.setMaxInactiveInterval(3600);
        String sessionId = session.getId();
        String uniqueLoginId = "Guest_" + sessionId.substring(0, 6);
        User existingUser = userService.findUserByLoginId(uniqueLoginId);
        if (existingUser != null) {
            // 이미 존재하는 유저가 있으면 그 유저로 로그인 처리
            session.setAttribute("userId", existingUser.getId());
            return ResponseEntity.ok("Existing user logged in with ID: " + uniqueLoginId);
        }

        // signupDto 생성 및 LoginId에 고유한 값 설정
        SignUpDto signupDto = new SignUpDto();
        signupDto.setLoginId(uniqueLoginId);  // 세션 ID를 포함한 고유한 ID
        signupDto.setPassword("default_password");  // 기본 비밀번호 설정
        signupDto.setFirstName("Guest");
        signupDto.setLastName("Guest");

        User newUser = userService.createUser(signupDto, "GUEST");

        Long guestUserId = userService.getUserIdByLoginId("Guest");
        List<String> guestHistoryNames = historyService.getHistoryNamesByUserId(guestUserId);

        if (guestHistoryNames.isEmpty()) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history for Guest user");
        }

        // 각 history_name에 대해 JsonData를 복사하여 새 유저로 저장
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
        session.setAttribute("userId", newUser.getId());
        return ResponseEntity.ok("New temporary user created and logged in with ID: " + uniqueLoginId);
    }

    @PostMapping("/userinfo")
    public ResponseEntity<String> getUserInfo(HttpSession session) {
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }

        User user = userService.findUserByUserId(userId);
        if(user==null){
            return ResponseEntity.status(400).body("유저 이름을 찾을 수 없습니다.");
        }


        Map<String, String> combinedJson = new HashMap<>();
        combinedJson.put("userName", user.getFirstName());
        combinedJson.put("loginId", user.getLoginId());

        String userInfo = new Gson().toJson(combinedJson);

        return ResponseEntity.ok(userInfo);
    }

    @PostMapping("/logout")
    public ResponseEntity<String> logout(HttpSession session) {
        session.invalidate();  // 세션 무효화
        return ResponseEntity.ok("User logged out successfully.");
    }


}
