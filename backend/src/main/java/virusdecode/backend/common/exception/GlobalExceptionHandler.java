package virusdecode.backend.common.exception;

import virusdecode.backend.analysis.exception.*;
import virusdecode.backend.analysis.exception.*;
import virusdecode.backend.bioinput.exception.FastaFileSaveFailException;
import virusdecode.backend.history.exception.HistoryNotFoundException;
import virusdecode.backend.user.exception.DuplicateLoginIdException;
import virusdecode.backend.user.exception.InvalidLoginException;
import virusdecode.backend.user.exception.UserNotFoundException;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.jdbc.support.MetaDataAccessException;
import org.springframework.web.bind.annotation.ExceptionHandler;
import org.springframework.web.bind.annotation.RestControllerAdvice;

@RestControllerAdvice
public class GlobalExceptionHandler {

    // User 예외 처리
    @ExceptionHandler(InvalidLoginException.class)
    public ResponseEntity<String> handleInvalidLogin(InvalidLoginException ex) {
        return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body(ex.getMessage());
    }

    @ExceptionHandler(DuplicateLoginIdException.class)
    public ResponseEntity<String> handleDuplicateLogin(DuplicateLoginIdException ex) {
        return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(ex.getMessage());
    }

    @ExceptionHandler(UnauthenticatedUserException.class)
    public ResponseEntity<String> handleUnauthenticated(UnauthenticatedUserException ex) {
        return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body(ex.getMessage());
    }

    @ExceptionHandler(UserNotFoundException.class)
    public ResponseEntity<String> handleUserNotFound(UserNotFoundException ex) {
        return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(ex.getMessage());
    }

    // History 예외 처리
    @ExceptionHandler(HistoryNotFoundException.class)
    public ResponseEntity<String> handleHistoryNotFound(HistoryNotFoundException ex) {
        ex.printStackTrace();
        return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(ex.getMessage());
    }

    // Analysis 예외 처리
    @ExceptionHandler(AnalysisNotFoundException.class)
    public ResponseEntity<String> handleAnalysisNotFound(AnalysisNotFoundException ex) {
        return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(ex.getMessage());
    }
    @ExceptionHandler(LinearDesignFailException.class)
    public ResponseEntity<String> handleLinearDesignFail(LinearDesignFailException ex) {
        return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(ex.getMessage());
    }
    @ExceptionHandler(PdbFailException.class)
    public ResponseEntity<String> handlePdbFail(PdbFailException ex) {
        return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(ex.getMessage());
    }
    @ExceptionHandler(InvalidSequenceRangeException.class)
    public ResponseEntity<String> handleInvalidRange(InvalidSequenceRangeException ex) {
        return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(ex.getMessage());
    }
    @ExceptionHandler(InvalidAnalysisDataException.class)
    public ResponseEntity<String> handleInvalidAnalysisDataException(InvalidAnalysisDataException ex){
        return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body(ex.getMessage());
    }


    // 바이오 데이터 입력 에러
    @ExceptionHandler(FastaFileSaveFailException.class)
    public ResponseEntity<String> handleFastaFileSaveFail(FastaFileSaveFailException ex){
        return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body(ex.getMessage());
    }

    @ExceptionHandler(MetaDataAccessException.class)
    public ResponseEntity<String> handleMetaDataAccess(MetaDataAccessException ex){
        return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body(ex.getMessage());
    }

    // 바이오 파이썬 에러
    @ExceptionHandler(BioPythonException.class)
    public ResponseEntity<String> handleBioPython(Exception ex){
        return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body(ex.getMessage());
    }

    // 예외 종류별로 추가 가능
    @ExceptionHandler(Exception.class)
    public ResponseEntity<String> handleGeneralException(Exception ex) {
        ex.printStackTrace(); // 로깅용
        return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("서버 오류가 발생했습니다.");
    }
}
