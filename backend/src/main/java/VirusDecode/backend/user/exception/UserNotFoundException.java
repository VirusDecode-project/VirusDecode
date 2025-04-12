package VirusDecode.backend.user.exception;

// 유저 정보를 찾을 수 없음
public class UserNotFoundException extends RuntimeException {
    public UserNotFoundException(String message) {
        super(message);
    }
}