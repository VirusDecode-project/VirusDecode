package virusdecode.backend.user.service;

import virusdecode.backend.analysis.service.AnalysisService;
import virusdecode.backend.history.entity.History;
import virusdecode.backend.history.service.HistoryService;
import virusdecode.backend.user.dto.SignUpDto;
import virusdecode.backend.user.dto.UserInfoDto;
import virusdecode.backend.user.entity.User;
import virusdecode.backend.user.exception.DuplicateLoginIdException;
import virusdecode.backend.user.exception.InvalidLoginException;
import virusdecode.backend.user.exception.UserNotFoundException;
import virusdecode.backend.user.repository.UserRepository;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.*;
import org.mockito.junit.jupiter.MockitoExtension;
import org.springframework.security.crypto.password.PasswordEncoder;
import org.springframework.test.util.ReflectionTestUtils;

import java.time.LocalDateTime;
import java.util.List;
import java.util.Optional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.*;

@ExtendWith(MockitoExtension.class)
class UserServiceTest {

    @InjectMocks
    private UserService userService;

    @Mock
    private UserRepository userRepository;
    @Mock
    private HistoryService historyService;
    @Mock
    private AnalysisService analysisService;
    @Mock
    private PasswordEncoder passwordEncoder;

    private User user;

    @BeforeEach
    void setup() {
        user = new User("John", "Doe", "testuser", "encodedpassword", "USER");
    }

    @Test
    @DisplayName("로그인 성공")
    void login_ShouldReturnUserInfo_WhenCredentialsValid() {
        when(userRepository.findByLoginId("testuser")).thenReturn(user);
        when(passwordEncoder.matches("password123", "encodedpassword")).thenReturn(true);

        UserInfoDto result = userService.login("testuser", "password123");

        assertEquals("testuser", result.getLoginId());
        assertEquals("John", result.getUserName());
    }

    @Test
    @DisplayName("로그인 실패 - 유저 정보 불일치")
    void login_ShouldThrowInvalidLoginException_WhenPasswordIncorrect() {
        when(userRepository.findByLoginId("testuser")).thenReturn(user);
        when(passwordEncoder.matches(anyString(), anyString())).thenReturn(false);

        assertThrows(InvalidLoginException.class, () -> {
            userService.login("testuser", "wrongpass");
        });
    }


    @Test
    @DisplayName("회원가입 성공")
    void createUser_ShouldSaveAndReturnUser_WhenLoginIdIsUnique() {
        // given
        SignUpDto dto = new SignUpDto("John", "Doe", "newuser", "pass1234");

        when(userRepository.findByLoginId("newuser")).thenReturn(null);
        when(passwordEncoder.encode("pass1234")).thenReturn("encodedpass");

        // 실제 저장될 User 객체를 미리 만들어서 리턴값으로 지정
        User savedUser = new User("John", "Doe", "newuser", "encodedpass", "USER");

        when(userRepository.save(any(User.class))).thenReturn(savedUser);

        // when
        User createdUser = userService.createUser(dto, "USER");

        // then
        assertEquals("newuser", createdUser.getLoginId());
        assertEquals("encodedpass", createdUser.getPassword());
    }


    @Test
    @DisplayName("회원가입 실패 - 중복 아이디")
    void createUser_ShouldThrowException_WhenLoginIdAlreadyExists() {
        when(userRepository.findByLoginId("testuser")).thenReturn(user);

        SignUpDto dto = new SignUpDto("John", "Doe", "testuser", "pass1234");

        assertThrows(DuplicateLoginIdException.class, () -> {
            userService.createUser(dto, "USER");
        });
    }

    @Test
    @DisplayName("유저 정보 조회 성공")
    void fetchUserInfo_ShouldReturnCorrectInfo() {
        when(userRepository.findById(1L)).thenReturn(Optional.of(user));

        UserInfoDto result = userService.fetchUserInfo(1L);

        assertEquals("testuser", result.getLoginId());
        assertEquals("John", result.getUserName());
    }

    @Test
    @DisplayName("유저 조회 실패")
    void getUserById_ShouldThrowException_WhenUserNotFound() {
        when(userRepository.findById(99L)).thenReturn(Optional.empty());

        assertThrows(UserNotFoundException.class, () -> userService.getUserById(99L));
    }

    @Test
    @DisplayName("게스트 정보 삭제 성공")
    void deleteGuestUsers_ShouldDeleteOldGuests() {
        User guest1 = new User();
        guest1.setRole("GUEST");
        guest1.setLoginId("Guest_123");
        guest1.setCreatedAt(LocalDateTime.now().minusHours(25));
        ReflectionTestUtils.setField(guest1, "id", 2L);

        when(userRepository.findUsersByRole("GUEST")).thenReturn(List.of(guest1));
        when(historyService.getHistoryNamesByUserId(2L)).thenReturn(List.of("History1"));
        History mockHistory = new History();
        when(historyService.getHistory("History1", 2L)).thenReturn(mockHistory);

        userService.deleteGuestUsers();

        verify(analysisService).deleteAnalysisData(mockHistory);
        verify(historyService).deleteHistory("History1", 2L);
        verify(userRepository).deleteUserById(2L);
    }

}