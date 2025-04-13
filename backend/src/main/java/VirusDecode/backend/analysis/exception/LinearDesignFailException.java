package VirusDecode.backend.analysis.exception;

// LinearDesign 실패
public class LinearDesignFailException extends RuntimeException{
    public LinearDesignFailException(String message) {
        super(message);
    }
}
