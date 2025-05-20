package virusdecode.backend.analysis.exception;

// 분석 정보를 찾을 수 없음
public class AnalysisNotFoundException extends RuntimeException{
    public AnalysisNotFoundException(String message) {
        super(message);
    }
}
