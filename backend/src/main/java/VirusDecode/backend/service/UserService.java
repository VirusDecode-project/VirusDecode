package VirusDecode.backend.service;

import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.repository.UserRepository;
import jakarta.transaction.Transactional;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.crypto.password.PasswordEncoder;
import org.springframework.stereotype.Service;

import java.util.Optional;
@Service
public class UserService {

    @Autowired
    private UserRepository userRepository;

    @Autowired
    private PasswordEncoder passwordEncoder;


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

    public String getUsernameByUserId(Long userId) {
        Optional<User> user = userRepository.findById(userId);
        return user.map(User::getFirstName).orElse(null);
    }

}
