package VirusDecode.backend.service;

import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.entity.User;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;

import jakarta.servlet.http.HttpSession;
import java.util.List;

@Service
public class GuestLoginService {
    private final UserService userService;
    private final HistoryService historyService;
    private final JsonDataService jsonDataService;

    @Autowired
    public GuestLoginService(UserService userService, HistoryService historyService, JsonDataService jsonDataService) {
        this.userService = userService;
        this.historyService = historyService;
        this.jsonDataService = jsonDataService;
    }

    public ResponseEntity<String> loginAsGuest(HttpSession session) {
        session.setMaxInactiveInterval(3600);
        String sessionId = session.getId();
        String uniqueLoginId = "Guest_" + sessionId.substring(0, 6);
        User existingUser = userService.findUserByLoginId(uniqueLoginId);
        if (existingUser != null) {
            session.setAttribute("userId", existingUser.getId());
            return ResponseEntity.ok("Existing user logged in with ID: " + uniqueLoginId);
        }

        SignUpDto signupDto = new SignUpDto();
        signupDto.setLoginId(uniqueLoginId);
        signupDto.setPassword("default_password");
        signupDto.setFirstName("Guest");
        signupDto.setLastName("Guest");

        User newUser = userService.createUser(signupDto, "GUEST");

        Long guestUserId = userService.getUserIdByLoginId("Guest");
        List<String> guestHistoryNames = historyService.getHistoryNamesByUserId(guestUserId);

        if (guestHistoryNames.isEmpty()) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history for Guest user");
        }

        copyGuestHistoriesToNewUser(guestHistoryNames, guestUserId, newUser);

        session.setAttribute("userId", newUser.getId());
        return ResponseEntity.ok("New temporary user created and logged in with ID: " + uniqueLoginId);
    }

    private void copyGuestHistoriesToNewUser(List<String> guestHistoryNames, Long guestUserId, User newUser) {
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
