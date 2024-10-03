package VirusDecode.backend.controller;
import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.dto.UserLoginDto;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.service.UserService;
import com.google.gson.Gson;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import jakarta.servlet.http.HttpSession;

import java.util.HashMap;
import java.util.Map;


@RestController
@RequestMapping("/api/auth")
public class UserController {
    private final UserService userService;

    @Autowired
    public UserController(UserService userService){
        this.userService = userService;
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

        User newUser = userService.createUser(signupDto);
        return ResponseEntity.ok("User created successfully with ID: " + newUser.getId());
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
