package VirusDecode.backend.service;

import VirusDecode.backend.dto.UserLoginDto;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.repository.UserRepository;
import jakarta.transaction.Transactional;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.Optional;
@Service
public class UserService {

    @Autowired
    private UserRepository userRepository;

    public User findUserByUsername(String username) {
        return userRepository.findByUsername(username);
    }

    @Transactional
    public User createUser(UserLoginDto loginDto) {
        User newUser = new User();
        newUser.setUsername(loginDto.getUsername());
        newUser.setPassword(loginDto.getPassword());
        // 추가적인 사용자 정보 설정
        return userRepository.save(newUser);
    }

    public boolean checkPassword(User user, String password) {
        // 비밀번호 확인 로직 구현
        return user.getPassword().equals(password); // 단순 비교 예시
    }

    // userId로 유저 객체를 반환
    public Optional<User> getUserById(Long userId) {
        return userRepository.findById(userId);
    }
}
