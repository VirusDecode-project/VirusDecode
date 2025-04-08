package VirusDecode.backend.common;

public class ApiResponse<T> {
    private boolean success;
    private int status;
    private String message;
    private T data;

    public ApiResponse(boolean success, int status, String message, T data) {
        this.success = success;
        this.status = status;
        this.message = message;
        this.data = data;
    }
}
