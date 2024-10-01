//package VirusDecode.backend.controller;
//
//import VirusDecode.backend.dto.SignUpDto;
//import VirusDecode.backend.dto.UserLoginDto;
//import VirusDecode.backend.entity.User;
//import VirusDecode.backend.service.UserService;
//import com.google.gson.Gson;
//import org.junit.jupiter.api.Assertions;
//import org.junit.jupiter.api.BeforeEach;
//import org.junit.jupiter.api.Test;
//import org.mockito.InjectMocks;
//import org.mockito.Mock;
//import org.mockito.MockitoAnnotations;
//import org.springframework.http.MediaType;
//import org.springframework.mock.web.MockHttpSession;
//import org.springframework.test.web.servlet.MockMvc;
//import org.springframework.test.web.servlet.setup.MockMvcBuilders;
//
//import static org.junit.jupiter.api.Assertions.assertNull;
//import static org.mockito.ArgumentMatchers.any;
//import static org.mockito.Mockito.*;
//import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
//import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.*;
//
//class UserControllerTest {
//
//    private MockMvc mockMvc;
//
//    @Mock
//    private UserService userService;
//
//    @InjectMocks
//    private UserController userController;
//
//    private MockHttpSession mockSession;
//
//    @BeforeEach
//    void setUp() {
//        MockitoAnnotations.openMocks(this);
//        mockMvc = MockMvcBuilders.standaloneSetup(userController).build();
//        mockSession = new MockHttpSession();
//    }
//
//    @Test
//    void testLogin_Successful() throws Exception {
//        // Given
//        UserLoginDto loginDto = new UserLoginDto();
//        loginDto.setLoginId("user123");
//        loginDto.setPassword("securePassword");
//
//        User user = new User();
//        user.setId(1L);
//        user.setLoginId("user123");
//        user.setPassword("hashedPassword");
//
//        when(userService.findUserByLoginId("user123")).thenReturn(user);
//        when(userService.checkPassword(user, "securePassword")).thenReturn(true);
//
//        // When & Then
//        mockMvc.perform(post("/auth/login")
//                        .session(mockSession)
//                        .contentType(MediaType.APPLICATION_JSON)
//                        .content(new Gson().toJson(loginDto)))
//                .andExpect(status().isOk())
//                .andExpect(content().string("User logged in successfully."));
//
//        verify(userService, times(1)).findUserByLoginId("user123");
//        verify(userService, times(1)).checkPassword(user, "securePassword");
//        Assertions.assertEquals(1L, mockSession.getAttribute("userId"));
//    }
//
//    @Test
//    void testLogin_InvalidLoginId() throws Exception {
//        // Given
//        UserLoginDto loginDto = new UserLoginDto();
//        loginDto.setLoginId("invalidUser");
//        loginDto.setPassword("somePassword");
//
//        when(userService.findUserByLoginId("invalidUser")).thenReturn(null);
//
//        // When & Then
//        mockMvc.perform(post("/auth/login")
//                        .session(mockSession)
//                        .contentType(MediaType.APPLICATION_JSON)
//                        .content(new Gson().toJson(loginDto)))
//                .andExpect(status().isUnauthorized())
//                .andExpect(content().string("Invalid login ID."));
//
//        verify(userService, times(1)).findUserByLoginId("invalidUser");
//        verify(userService, never()).checkPassword(any(), any());
//        assertNull(mockSession.getAttribute("userId"));
//    }
//
//    @Test
//    void testLogin_InvalidPassword() throws Exception {
//        // Given
//        UserLoginDto loginDto = new UserLoginDto();
//        loginDto.setLoginId("user123");
//        loginDto.setPassword("wrongPassword");
//
//        User user = new User();
//        user.setId(1L);
//        user.setLoginId("user123");
//        user.setPassword("hashedPassword");
//
//        when(userService.findUserByLoginId("user123")).thenReturn(user);
//        when(userService.checkPassword(user, "wrongPassword")).thenReturn(false);
//
//        // When & Then
//        mockMvc.perform(post("/auth/login")
//                        .session(mockSession)
//                        .contentType(MediaType.APPLICATION_JSON)
//                        .content(new Gson().toJson(loginDto)))
//                .andExpect(status().isUnauthorized())
//                .andExpect(content().string("Invalid password."));
//
//        verify(userService, times(1)).findUserByLoginId("user123");
//        verify(userService, times(1)).checkPassword(user, "wrongPassword");
//        assertNull(mockSession.getAttribute("userId"));
//    }
//
//    @Test
//    void testSignup_UserAlreadyExists() throws Exception {
//        // Given
//        SignUpDto signupDto = new SignUpDto();
//        signupDto.setLoginId("user123");
//        signupDto.setPassword("securePassword");
//        signupDto.setFirstName("John");
//        signupDto.setLastName("Doe");
//
//        when(userService.findUserByLoginId("user123")).thenReturn(new User());
//
//        // When & Then
//        mockMvc.perform(post("/auth/signup")
//                        .contentType(MediaType.APPLICATION_JSON)
//                        .content(new Gson().toJson(signupDto)))
//                .andExpect(status().isBadRequest())
//                .andExpect(content().string("Login ID already exists."));
//
//        verify(userService, times(1)).findUserByLoginId("user123");
//        verify(userService, never()).createUser(any(SignUpDto.class));
//    }
//
//    @Test
//    void testSignup_Successful() throws Exception {
//        // Given
//        SignUpDto signupDto = new SignUpDto();
//        signupDto.setLoginId("newUser123");
//        signupDto.setPassword("securePassword");
//        signupDto.setFirstName("John");
//        signupDto.setLastName("Doe");
//
//        User newUser = new User();
//        newUser.setId(1L);
//        newUser.setLoginId("newUser123");
//
//        when(userService.findUserByLoginId("newUser123")).thenReturn(null);
//        when(userService.createUser(any(SignUpDto.class))).thenReturn(newUser);
//
//        // When & Then
//        mockMvc.perform(post("/auth/signup")
//                        .contentType(MediaType.APPLICATION_JSON)
//                        .content(new Gson().toJson(signupDto)))
//                .andExpect(status().isOk())
//                .andExpect(content().string("User created successfully with ID: 1"));
//
//        verify(userService, times(1)).findUserByLoginId("newUser123");
//        verify(userService, times(1)).createUser(any(SignUpDto.class));
//    }
//}
