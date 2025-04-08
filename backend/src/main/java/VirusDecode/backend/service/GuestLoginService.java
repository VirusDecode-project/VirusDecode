package VirusDecode.backend.service;

import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.entity.User;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

@Service
public class GuestLoginService {
    private final UserService userService;

    @Autowired
    public GuestLoginService(UserService userService) {
        this.userService = userService;
    }

    public User createGuestUser(String uniqueLoginId) {
        SignUpDto signupDto = new SignUpDto("Guest", "Guest", uniqueLoginId, "default_password");
        return userService.createUser(signupDto, "GUEST");
    }
}
