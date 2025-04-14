package virusdecode.backend.user.service;

import virusdecode.backend.user.dto.SignUpDto;
import virusdecode.backend.user.entity.User;
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

    public User loginOrCreateGuest(String uniqueLoginId){
        User existingUser = userService.findUserByLoginId(uniqueLoginId);
        if (existingUser != null) {
            return existingUser;
        }
        return createGuestUser(uniqueLoginId);
    }
}
