package virusdecode.backend.common.exception;

// 인증되지 않은 사용자
public class UnauthenticatedUserException extends RuntimeException {
    public UnauthenticatedUserException(String message) {
        super(message);
    }
}