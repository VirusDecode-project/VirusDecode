package VirusDecode.backend.controller;

import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.dto.UserInfoDto;
import VirusDecode.backend.dto.UserLoginDto;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.service.GuestLoginService;
import VirusDecode.backend.service.UserService;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;

import jakarta.servlet.http.HttpSession;

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

        UserInfoDto mockUserInfo = new UserInfoDto("testUser", "testName");

        when(userService.login("testUser", "password")).thenReturn(mockUserInfo);
        when(userService.getUserIdByLoginId("testUser")).thenReturn(1L);

        ResponseEntity<?> response = userController.login(loginDto, session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals(mockUserInfo, response.getBody());
        verify(session).setAttribute("userId", mockUser.getId());
    }

    @Test
    void testLogin_InvalidUser() {
        UserLoginDto loginDto = new UserLoginDto();
        loginDto.setLoginId("nonExistentUser");
        loginDto.setPassword("password");

        when(userService.login("testUser", "password")).thenReturn(null);
        when(userService.findUserByLoginId("nonExistentUser")).thenReturn(null);

        ResponseEntity<?> response = userController.login(loginDto, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("유효하지 않는 회원 정보입니다.", response.getBody());
    }

//    @Test
//    void testLogin_InvalidPassword() {
//        UserLoginDto loginDto = new UserLoginDto();
//        loginDto.setLoginId("testUser");
//        loginDto.setPassword("wrongPassword");
//
//        User mockUser = new User();
//        mockUser.setId(1L);
//        mockUser.setLoginId("testUser");
//
//        when(userService.login("testUser", "password")).thenReturn(null);
////        when(userService.findUserByLoginId("testUser")).thenReturn(mockUser);
////        when(userService.checkPassword(mockUser, "wrongPassword")).thenReturn(false);
//
//        ResponseEntity<?> response = userController.login(loginDto, session);
//
//        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
//        assertEquals("비밀번호가 틀렸습니다.", response.getBody());
//    }

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
    void testGetUserInfo_Success() {
        // given
        Long userId = 1L;
        UserInfoDto userInfo = new UserInfoDto("testUser", "홍길동");

        when(session.getAttribute("userId")).thenReturn(userId);
        when(userService.fetchUserInfo(userId)).thenReturn(userInfo);

        // when
        ResponseEntity<?> response = userController.getUserInfo(session);

        // then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals(userInfo, response.getBody());
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
        ResponseEntity<?> response = userController.getUserInfo(session);

        // Then
        assertEquals(HttpStatus.BAD_REQUEST, response.getStatusCode());
        assertEquals("유저 이름을 찾을 수 없습니다.", response.getBody());
    }

    @Test
    void testGetUserInfo_NotAuthenticated() {
        when(session.getAttribute("userId")).thenReturn(null);

        ResponseEntity<?> response = userController.getUserInfo(session);

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
