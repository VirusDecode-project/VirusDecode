package VirusDecode.backend.controller;

import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.dto.UserLoginDto;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.service.GuestLoginService;
import VirusDecode.backend.service.UserService;
import com.google.gson.Gson;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;

import jakarta.servlet.http.HttpSession;

import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.*;

class UserControllerTest {

    @Mock
    private UserService userService;

    @Mock
    private GuestLoginService guestLoginService;

    @Mock
    private HttpSession session;

    @InjectMocks
    private UserController userController;

    @BeforeEach
    void setup() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    void testLogin_Success() {
        UserLoginDto loginDto = new UserLoginDto();
        loginDto.setLoginId("testUser");
        loginDto.setPassword("password");

        User mockUser = new User();
        mockUser.setId(1L);
        mockUser.setLoginId("testUser");

        when(userService.findUserByLoginId("testUser")).thenReturn(mockUser);
        when(userService.checkPassword(mockUser, "password")).thenReturn(true);

        ResponseEntity<String> response = userController.login(loginDto, session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("User logged in successfully.", response.getBody());
        verify(session).setAttribute("userId", mockUser.getId());
    }

    @Test
    void testLogin_InvalidUser() {
        UserLoginDto loginDto = new UserLoginDto();
        loginDto.setLoginId("nonExistentUser");
        loginDto.setPassword("password");

        when(userService.findUserByLoginId("nonExistentUser")).thenReturn(null);

        ResponseEntity<String> response = userController.login(loginDto, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("유효하지 않은 ID 입니다.", response.getBody());
    }

    @Test
    void testLogin_InvalidPassword() {
        UserLoginDto loginDto = new UserLoginDto();
        loginDto.setLoginId("testUser");
        loginDto.setPassword("wrongPassword");

        User mockUser = new User();
        mockUser.setId(1L);
        mockUser.setLoginId("testUser");

        when(userService.findUserByLoginId("testUser")).thenReturn(mockUser);
        when(userService.checkPassword(mockUser, "wrongPassword")).thenReturn(false);

        ResponseEntity<String> response = userController.login(loginDto, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("비밀번호가 틀렸습니다.", response.getBody());
    }

    @Test
    void testSignup_Success() {
        SignUpDto signUpDto = new SignUpDto();
        signUpDto.setLoginId("newUser");
        signUpDto.setPassword("password");

        when(userService.findUserByLoginId("newUser")).thenReturn(null);

        User newUser = new User();
        newUser.setId(1L);
        newUser.setLoginId("newUser");

        when(userService.createUser(signUpDto, "USER")).thenReturn(newUser);

        ResponseEntity<String> response = userController.signup(signUpDto);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("User created successfully with ID: " + newUser.getId(), response.getBody());
    }

    @Test
    void testSignup_UserAlreadyExists() {
        SignUpDto signUpDto = new SignUpDto();
        signUpDto.setLoginId("existingUser");

        when(userService.findUserByLoginId("existingUser")).thenReturn(new User());

        ResponseEntity<String> response = userController.signup(signUpDto);

        assertEquals(HttpStatus.BAD_REQUEST, response.getStatusCode());
        assertEquals("이미 존재하는 ID 입니다.", response.getBody());
    }

    @Test
    void testGuestLogin_Success() {
        when(guestLoginService.loginAsGuest(session)).thenReturn(ResponseEntity.ok("New temporary user created and logged in"));

        ResponseEntity<String> response = userController.guestLogin(session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("New temporary user created and logged in", response.getBody());
    }

    @Test
    void testGetUserInfo_Success() {
        Long userId = 1L;
        User mockUser = new User();
        mockUser.setId(userId);
        mockUser.setFirstName("First");
        mockUser.setLoginId("testUser");

        when(session.getAttribute("userId")).thenReturn(userId);
        when(userService.findUserByUserId(userId)).thenReturn(mockUser);

        ResponseEntity<String> response = userController.getUserInfo(session);

        Map<String, String> expectedJson = new HashMap<>();
        expectedJson.put("userName", "First");
        expectedJson.put("loginId", "testUser");

        String expectedResponse = new Gson().toJson(expectedJson);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals(expectedResponse, response.getBody());
    }
    @Test
    void testGetUserInfo_UserNotFound() {
        // Given
        Long userId = 1L;

        // 세션에 userId가 존재하는 경우
        when(session.getAttribute("userId")).thenReturn(userId);
        // 유저를 찾을 수 없는 경우 null 반환
        when(userService.findUserByUserId(userId)).thenReturn(null);

        // When
        ResponseEntity<String> response = userController.getUserInfo(session);

        // Then
        assertEquals(HttpStatus.BAD_REQUEST, response.getStatusCode());
        assertEquals("유저 이름을 찾을 수 없습니다.", response.getBody());
    }

    @Test
    void testGetUserInfo_NotAuthenticated() {
        when(session.getAttribute("userId")).thenReturn(null);

        ResponseEntity<String> response = userController.getUserInfo(session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("User not authenticated", response.getBody());
    }

    @Test
    void testLogout() {
        ResponseEntity<String> response = userController.logout(session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("User logged out successfully.", response.getBody());
        verify(session).invalidate();
    }
}
