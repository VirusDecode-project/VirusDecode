package VirusDecode.backend.user.controller;

import VirusDecode.backend.user.dto.SignUpDto;
import VirusDecode.backend.user.dto.UserInfoDto;
import VirusDecode.backend.user.dto.UserLoginDto;
import VirusDecode.backend.user.entity.User;
import VirusDecode.backend.user.service.GuestLoginService;
import VirusDecode.backend.user.service.UserService;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.http.MediaType;
import org.springframework.test.web.servlet.MockMvc;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.*;

@WebMvcTest(UserController.class)
public class UserControllerTest {

    @Autowired
    private MockMvc mockMvc;

    @MockBean
    private UserService userService;

    @MockBean
    private GuestLoginService guestLoginService;

    private final ObjectMapper objectMapper = new ObjectMapper();

    @Test
    @DisplayName("로그인 성공")
    void login_Success() throws Exception {
        UserLoginDto loginDto = new UserLoginDto("testuser", "password123");
        UserInfoDto userInfoDto = new UserInfoDto("testuser", "John");

        Mockito.when(userService.login("testuser", "password123")).thenReturn(userInfoDto);
        Mockito.when(userService.getUserIdByLoginId("testuser")).thenReturn(1L);

        mockMvc.perform(post("/api/auth/login")
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(objectMapper.writeValueAsString(loginDto)))
                .andExpect(status().isOk())
                .andExpect(jsonPath("$.loginId").value("testuser"))
                .andExpect(jsonPath("$.userName").value("John"));
    }

    @Test
    @DisplayName("회원가입 성공")
    void signup_Success() throws Exception {
        SignUpDto signUpDto = new SignUpDto("John", "Doe", "newuser", "pass1234");
        User createdUser = new User();
        createdUser.setId(10L);

        Mockito.when(userService.createUser(Mockito.any(SignUpDto.class), Mockito.eq("USER")))
                .thenReturn(createdUser);

        mockMvc.perform(post("/api/auth/signup")
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(objectMapper.writeValueAsString(signUpDto)))
                .andExpect(status().isOk())
                .andExpect(content().string("User created successfully with ID: 10"));
    }

    @Test
    @DisplayName("게스트 로그인 성공")
    void guestLogin_Success() throws Exception {
        User guestUser = new User();
        guestUser.setId(100L);
        guestUser.setLoginId("Guest_123456");

        Mockito.when(guestLoginService.loginOrCreateGuest(Mockito.anyString())).thenReturn(guestUser);

        mockMvc.perform(post("/api/auth/guest-login"))
                .andExpect(status().isOk())
                .andExpect(content().string(org.hamcrest.Matchers.containsString("Guest_")));
    }

    @Test
    @DisplayName("유저 세션 정보 조회 성공")
    void getUserInfo_Authenticated() throws Exception {
        UserInfoDto userInfoDto = new UserInfoDto("testuser", "John");

        Mockito.when(userService.fetchUserInfo(1L)).thenReturn(userInfoDto);

        mockMvc.perform(post("/api/auth/userinfo")
                        .sessionAttr("userId", 1L))
                .andExpect(status().isOk())
                .andExpect(jsonPath("$.loginId").value("testuser"))
                .andExpect(jsonPath("$.userName").value("John"));
    }

    @Test
    @DisplayName("유저 세션 정보 조회 실패")
    void getUserInfo_Unauthenticated() throws Exception {
        mockMvc.perform(post("/api/auth/userinfo"))
                .andExpect(status().isUnauthorized())
                .andExpect(content().string("User not authenticated"));
    }

    @Test
    @DisplayName("로그아웃 성공")
    void logout_Success() throws Exception {
        mockMvc.perform(post("/api/auth/logout")
                        .sessionAttr("userId", 1L))
                .andExpect(status().isOk())
                .andExpect(content().string("User logged out successfully."));
    }
}
