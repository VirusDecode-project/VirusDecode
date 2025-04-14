package virusdecode.backend.user.controller;

import virusdecode.backend.user.dto.SignUpDto;
import virusdecode.backend.user.dto.UserInfoDto;
import virusdecode.backend.user.dto.UserLoginDto;
import virusdecode.backend.user.entity.User;
import virusdecode.backend.user.service.GuestLoginService;
import virusdecode.backend.user.service.UserService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import jakarta.servlet.http.HttpSession;

import static virusdecode.backend.util.UserSessionUtil.getAuthenticatedUserId;

@RestController
@RequestMapping("/api/auth")
public class UserController {
    private final UserService userService;
    private final GuestLoginService guestLoginService;

    @Autowired
    public UserController(UserService userService, GuestLoginService guestLoginService){
        this.userService = userService;
        this.guestLoginService = guestLoginService;
    }

    @PostMapping("/login")
    public ResponseEntity<?> login(@RequestBody UserLoginDto loginDto, HttpSession session) {
        session.setMaxInactiveInterval(3600);
        UserInfoDto userInfo = userService.login(loginDto.getLoginId(), loginDto.getPassword());
        Long userId = userService.getUserIdByLoginId(userInfo.getLoginId());
        session.setAttribute("userId", userId);
        return ResponseEntity.ok(userInfo);
    }

    @PostMapping("/signup")
    public ResponseEntity<String> signup(@RequestBody SignUpDto signupDto) {
        User newUser = userService.createUser(signupDto, "USER");
        return ResponseEntity.ok("User created successfully with ID: " + newUser.getId());
    }

    @PostMapping("/guest-login")
    public ResponseEntity<String> guestLogin(HttpSession session) {
        session.setMaxInactiveInterval(3600);

        String sessionId = session.getId();
        String shortSessionId = sessionId.length() >= 6 ? sessionId.substring(0, 6) : sessionId;
        String uniqueLoginId = "Guest_" + shortSessionId;

        // 새로운 게스트 유저십니다.
        User user = guestLoginService.loginOrCreateGuest(uniqueLoginId);
        session.setAttribute("userId", user.getId());

        return ResponseEntity.ok("New temporary user created and logged in with ID: " + uniqueLoginId);
    }

    @PostMapping("/userinfo")
    public ResponseEntity<?> getUserInfo(HttpSession session) {
        Long userId = getAuthenticatedUserId(session);

        // 세션 아이디로 유저를 찾아냅니다.
        UserInfoDto userInfo = userService.fetchUserInfo(userId);

        return ResponseEntity.ok(userInfo);
    }

    @PostMapping("/logout")
    public ResponseEntity<String> logout(HttpSession session) {
        session.invalidate();  // 세션 무효화
        return ResponseEntity.ok("User logged out successfully.");
    }
}
