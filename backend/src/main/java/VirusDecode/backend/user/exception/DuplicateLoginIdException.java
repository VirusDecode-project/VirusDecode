package VirusDecode.backend.user.exception;

// 이미 존재하는 로그인 ID 예외
public class DuplicateLoginIdException extends RuntimeException {
    public DuplicateLoginIdException(String message) {
        super(message);
    }
}