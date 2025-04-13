package VirusDecode.backend.user.service;

import VirusDecode.backend.user.dto.SignUpDto;
import VirusDecode.backend.user.entity.User;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.*;
import org.mockito.junit.jupiter.MockitoExtension;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.*;

@ExtendWith(MockitoExtension.class)
class GuestLoginServiceTest {

    @InjectMocks
    private GuestLoginService guestLoginService;

    @Mock
    private UserService userService;

    @BeforeEach
    void setUp() {
//        MockitoAnnotations.openMocks(this);
    }

    @Test
    @DisplayName("기존 게스트 유저가 존재하면 그대로 반환")
    void loginOrCreateGuest_ShouldReturnExistingUser() {
        User mockUser = new User();
        mockUser.setLoginId("Guest_123456");

        when(userService.findUserByLoginId("Guest_123456")).thenReturn(mockUser);

        User result = guestLoginService.loginOrCreateGuest("Guest_123456");

        assertEquals("Guest_123456", result.getLoginId());
        verify(userService, never()).createUser(any(SignUpDto.class), anyString());
    }

    @Test
    @DisplayName("게스트 유저가 없으면 새로 생성")
    void loginOrCreateGuest_ShouldCreateNewGuestUser() {
        User newUser = new User();
        newUser.setLoginId("Guest_654321");

        when(userService.findUserByLoginId("Guest_654321")).thenReturn(null);
        when(userService.createUser(any(SignUpDto.class), eq("GUEST"))).thenReturn(newUser);

        User result = guestLoginService.loginOrCreateGuest("Guest_654321");

        assertEquals("Guest_654321", result.getLoginId());
        verify(userService).createUser(any(SignUpDto.class), eq("GUEST"));
    }

    @Test
    @DisplayName("createGuestUser는 createUser를 GUEST 권한으로 호출")
    void createGuestUser_ShouldCallCreateUserWithGuestRole() {
        User newUser = new User();
        newUser.setLoginId("Guest_777");

        when(userService.createUser(any(SignUpDto.class), eq("GUEST"))).thenReturn(newUser);

        User result = guestLoginService.createGuestUser("Guest_777");

        assertEquals("Guest_777", result.getLoginId());
        verify(userService).createUser(argThat(dto ->
                dto.getLoginId().equals("Guest_777") &&
                        dto.getFirstName().equals("Guest") &&
                        dto.getLastName().equals("Guest") &&
                        dto.getPassword().equals("default_password")
        ), eq("GUEST"));
    }
}
