package VirusDecode.backend.User.controller;
import VirusDecode.backend.User.dto.SignUpDto;
import VirusDecode.backend.User.dto.UserInfoDto;
import VirusDecode.backend.User.dto.UserLoginDto;
import VirusDecode.backend.User.entity.User;
import VirusDecode.backend.User.service.GuestLoginService;
import VirusDecode.backend.User.service.UserService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import jakarta.servlet.http.HttpSession;



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
        if (userInfo == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("유효하지 않는 회원 정보입니다.");
        }

        Long userId = userService.getUserIdByLoginId(userInfo.getLoginId());
        session.setAttribute("userId", userId);

        return ResponseEntity.ok(userInfo);
    }
    @PostMapping("/signup")
    public ResponseEntity<String> signup(@RequestBody SignUpDto signupDto) {
        if (userService.findUserByLoginId(signupDto.getLoginId()) != null) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body("이미 존재하는 ID 입니다.");
        }

        User newUser = userService.createUser(signupDto, "USER");

        return ResponseEntity.ok("User created successfully with ID: " + newUser.getId());
    }

    @PostMapping("/guest-login")
    public ResponseEntity<String> guestLogin(HttpSession session) {
        session.setMaxInactiveInterval(3600);

        // 이미 존재하는 분인가요?
        String sessionId = session.getId();
        String shortSessionId = sessionId.length() >= 6 ? sessionId.substring(0, 6) : sessionId;
        String uniqueLoginId = "Guest_" + shortSessionId;


        User existingUser = userService.findUserByLoginId(uniqueLoginId);
        if (existingUser != null) {
            session.setAttribute("userId", existingUser.getId());
            return ResponseEntity.ok("Existing user logged in with ID: " + uniqueLoginId);
        }

        // 새로운 게스트 유저십니다.
        User newUser = guestLoginService.createGuestUser(uniqueLoginId);

        session.setAttribute("userId", newUser.getId());
        return ResponseEntity.ok("New temporary user created and logged in with ID: " + uniqueLoginId);
    }

    @PostMapping("/userinfo")
    public ResponseEntity<?> getUserInfo(HttpSession session) {
        Long userId = (Long) session.getAttribute("userId");

        // 세션에 아이디가 없습니다.
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }

        // 세션 아이디로 유저를 찾아냅니다.
        UserInfoDto userInfo = userService.fetchUserInfo(userId);

        // 유저가 없군요?
        if(userInfo == null){
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body("유저 이름을 찾을 수 없습니다.");
        }

        // 세션에 저장된 유저입니다.
        return ResponseEntity.ok(userInfo);
    }

    @PostMapping("/logout")
    public ResponseEntity<String> logout(HttpSession session) {
        session.invalidate();  // 세션 무효화
        return ResponseEntity.ok("User logged out successfully.");
    }


}
