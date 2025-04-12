package VirusDecode.backend.user.exception;

// 로그인 실패
public class InvalidLoginException extends RuntimeException {
    public InvalidLoginException(String message) {
        super(message);
    }
}
