package VirusDecode.backend.history.exception;

// 히스토리를 찾을 수 없음
public class HistoryNotFoundException extends RuntimeException{
    public HistoryNotFoundException(String message) {
        super(message);
    }
}
