package VirusDecode.backend.controller;
import VirusDecode.backend.dto.UserLoginDto;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.service.UserService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import jakarta.servlet.http.HttpSession;


@RestController
@RequestMapping("/user")
public class UserController {
    private final UserService userService;

    @Autowired
    public UserController(UserService userService){
        this.userService = userService;
    }

    @PostMapping("/login")
    public ResponseEntity<String> login(@RequestBody UserLoginDto loginDto, HttpSession session) {
        User user = userService.findUserByUsername(loginDto.getUsername());
        session.setMaxInactiveInterval(3600);

        if (user == null) {
            user = userService.createUser(loginDto);
            session.setAttribute("userId", user.getId());
            return ResponseEntity.ok("New user created and logged in.");
        }

        // 기존 유저 로그인 처리
        if (userService.checkPassword(user, loginDto.getPassword())) {
            session.setAttribute("userId", user.getId());
            return ResponseEntity.ok("User logged in successfully.");
        } else {
            return ResponseEntity.status(401).body("Invalid password.");
        }
    }
}
