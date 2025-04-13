package VirusDecode.backend.util;

import VirusDecode.backend.common.exception.UnauthenticatedUserException;
import jakarta.servlet.http.HttpSession;

public class UserSessionUtil {
    public static Long getAuthenticatedUserId(HttpSession session) {
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            throw new UnauthenticatedUserException("User not authenticated");
        }
        return userId;
    }
}
