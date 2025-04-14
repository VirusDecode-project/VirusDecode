package virusdecode.backend.analysis.exception;

// Pdb 실패
public class PdbFailException extends RuntimeException{
    public PdbFailException(String message) {
        super(message);
    }
}
