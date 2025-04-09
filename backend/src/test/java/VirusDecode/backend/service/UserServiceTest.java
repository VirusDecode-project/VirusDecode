package VirusDecode.backend.service;

import VirusDecode.backend.analysis.service.AnalysisService;
import VirusDecode.backend.user.dto.SignUpDto;
import VirusDecode.backend.user.service.UserService;
import VirusDecode.backend.history.entity.History;
import VirusDecode.backend.analysis.entity.Analysis;
import VirusDecode.backend.user.entity.User;
import VirusDecode.backend.user.repository.UserRepository;
import VirusDecode.backend.history.service.HistoryService;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.security.crypto.password.PasswordEncoder;

import java.time.LocalDateTime;
import java.util.List;
import java.util.Optional;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;

class UserServiceTest {

    @Mock
    private AnalysisService analysisService;

    @Mock
    private HistoryService historyService;

    @Mock
    private UserRepository userRepository;

    @Mock
    private PasswordEncoder passwordEncoder;

    @InjectMocks
    private UserService userService;

    @BeforeEach
    void setUp() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    void testFindUserByLoginId() {
        // Given
        String loginId = "testUser";
        User mockUser = new User();
        mockUser.setLoginId(loginId);
        when(userRepository.findByLoginId(loginId)).thenReturn(mockUser);

        // When
        User user = userService.findUserByLoginId(loginId);

        // Then
        assertNotNull(user);
        assertEquals(loginId, user.getLoginId());
        verify(userRepository, times(1)).findByLoginId(loginId);
    }

    @Test
    void testFindUserByUserId() {
        // Given
        Long userId = 1L;
        User mockUser = new User();
        mockUser.setId(userId);
        when(userRepository.findById(userId)).thenReturn(Optional.of(mockUser));

        // When
        User user = userService.findUserByUserId(userId);

        // Then
        assertNotNull(user);
        assertEquals(userId, user.getId());
        verify(userRepository, times(1)).findById(userId);
    }

    @Test
    void testCreateUser() {
        // Given
        SignUpDto signUpDto = new SignUpDto();
        signUpDto.setFirstName("First");
        signUpDto.setLastName("Last");
        signUpDto.setLoginId("testUser");
        signUpDto.setPassword("password");
        String role = "USER";

        // Mock the password encoding and UserRepository save method
        when(passwordEncoder.encode(signUpDto.getPassword())).thenReturn("encodedPassword");

        User mockUser = new User();
        mockUser.setId(1L);
        mockUser.setLoginId(signUpDto.getLoginId());
        mockUser.setPassword("encodedPassword");
        mockUser.setRole(role);

        when(userRepository.save(any(User.class))).thenReturn(mockUser);

        // When
        User newUser = userService.createUser(signUpDto, role);

        // Then
        assertNotNull(newUser);
        assertEquals("encodedPassword", newUser.getPassword());
        assertEquals(signUpDto.getLoginId(), newUser.getLoginId());
        assertEquals(role, newUser.getRole());
        verify(userRepository, times(1)).save(any(User.class));
    }


    @Test
    void testCheckPassword() {
        // Given
        User user = new User();
        user.setPassword("encodedPassword");
        when(passwordEncoder.matches("password", "encodedPassword")).thenReturn(true);

        // When
        boolean isMatch = userService.checkPassword(user, "password");

        // Then
        assertTrue(isMatch);
        verify(passwordEncoder, times(1)).matches("password", "encodedPassword");
    }

    @Test
    void testGetUserById() {
        // Given
        Long userId = 1L;
        User mockUser = new User();
        mockUser.setId(userId);
        when(userRepository.findById(userId)).thenReturn(Optional.of(mockUser));

        // When
        Optional<User> user = userService.getUserById(userId);

        // Then
        assertTrue(user.isPresent());
        assertEquals(userId, user.get().getId());
        verify(userRepository, times(1)).findById(userId);
    }

    @Test
    void testGetUserIdByLoginId() {
        // Given
        String loginId = "testUser";
        User mockUser = new User();
        mockUser.setId(1L);
        when(userRepository.findByLoginId(loginId)).thenReturn(mockUser);

        // When
        Long userId = userService.getUserIdByLoginId(loginId);

        // Then
        assertEquals(1L, userId);
        verify(userRepository, times(1)).findByLoginId(loginId);
    }
    @Test
    void testGetUserIdByLoginId_Null() {
        // Given
        String loginId = "nonExistentUser";
        when(userRepository.findByLoginId(loginId)).thenReturn(null);

        // When
        Long userId = userService.getUserIdByLoginId(loginId);

        // Then
        assertNull(userId);
        verify(userRepository, times(1)).findByLoginId(loginId);
    }

    @Test
    void testDeleteGuestUsers() {
        // Given
        User guestUser = new User();
        guestUser.setId(1L);
        guestUser.setLoginId("Guest_123");
        guestUser.setRole("GUEST");
        guestUser.setCreatedAt(LocalDateTime.now().minusDays(1));

        when(userRepository.findUsersByRole("GUEST")).thenReturn(List.of(guestUser));
        when(historyService.getHistoryNamesByUserId(guestUser.getId())).thenReturn(List.of("History1"));

        History mockHistory = new History();
        when(historyService.getHistory("History1", guestUser.getId())).thenReturn(mockHistory);

        // When
        userService.deleteGuestUsers();

        // Then
        verify(analysisService, times(1)).deleteAnalysisData(mockHistory);
        verify(historyService, times(1)).deleteHistory("History1", guestUser.getId());
        verify(userRepository, times(1)).deleteUserById(guestUser.getId());
    }

    @Test
    void testDeleteGuestUsers_UserNotOldEnough() {
        // Given
        User guestUser = new User();
        guestUser.setId(1L);
        guestUser.setLoginId("Guest_123");
        guestUser.setRole("GUEST");
        guestUser.setCreatedAt(LocalDateTime.now().minusHours(23)); // 23시간 전, 조건을 만족하지 않음

        when(userRepository.findUsersByRole("GUEST")).thenReturn(List.of(guestUser));

        // When
        userService.deleteGuestUsers();

        // Then
        // jsonDataService 및 historyService가 호출되지 않는지 검증
        verify(analysisService, never()).deleteAnalysisData(any());
        verify(historyService, never()).deleteHistory(anyString(), anyLong());
        verify(userRepository, never()).deleteUserById(guestUser.getId());
    }

    @Test
    void testCopySampleHistoriesToNewUser() {
        // Given
        String guestLoginId = "Guest";
        Long guestUserId = 1L;
        User guestUser = new User();
        guestUser.setId(guestUserId);
        guestUser.setLoginId(guestLoginId);

        // New user for whom histories will be copied
        User newUser = new User();
        newUser.setId(2L);

        // Set up mock return values
        when(userRepository.findByLoginId(guestLoginId)).thenReturn(guestUser); // User 객체 반환
        when(historyService.getHistoryNamesByUserId(guestUserId)).thenReturn(List.of("SampleHistory"));

        History guestHistory = new History();
        guestHistory.setHistoryName("SampleHistory");
        guestHistory.setUser(guestUser);

        Analysis guestAnalysis = new Analysis();
        guestAnalysis.setReferenceId("Ref123");
        guestAnalysis.setAlignment("AlignmentData");
        guestAnalysis.setLinearDesign("LinearDesignData");
        guestAnalysis.setPdb("PdbData");

        when(historyService.getHistory("SampleHistory", guestUserId)).thenReturn(guestHistory);
        when(analysisService.getAnalysisData(guestHistory)).thenReturn(guestAnalysis);

        // When
        userService.copySampleHistoriesToNewUser(newUser);

        // Then
        verify(historyService, times(1)).createHistory(any(History.class));
        verify(analysisService, times(1)).saveAnalysisData(any(Analysis.class));
    }
    @Test
    void testCopySampleHistoriesToNewUser_OriginalJsonDataIsNull() {
        // Given
        String guestLoginId = "Guest";
        Long guestUserId = 1L;
        User guestUser = new User();
        guestUser.setId(guestUserId);
        guestUser.setLoginId(guestLoginId);

        // New user for whom histories will be copied
        User newUser = new User();
        newUser.setId(2L);

        // Set up mock return values
        when(userRepository.findByLoginId(guestLoginId)).thenReturn(guestUser);
        when(historyService.getHistoryNamesByUserId(guestUserId)).thenReturn(List.of("SampleHistory"));

        History guestHistory = new History();
        guestHistory.setHistoryName("SampleHistory");
        guestHistory.setUser(guestUser);

        // Here, originalJsonData will be null
        when(historyService.getHistory("SampleHistory", guestUserId)).thenReturn(guestHistory);
        when(analysisService.getAnalysisData(guestHistory)).thenReturn(null); // Simulate null condition

        // When
        userService.copySampleHistoriesToNewUser(newUser);

        // Then
        verify(historyService, never()).createHistory(any(History.class));
        verify(analysisService, never()).saveAnalysisData(any(Analysis.class));
    }

    @Test
    void testCopySampleHistoriesToNewUser_HistoryIsNull() {
        // Given
        String guestLoginId = "Guest";
        Long guestUserId = 1L;
        User guestUser = new User();
        guestUser.setId(guestUserId);
        guestUser.setLoginId(guestLoginId);

        // New user for whom histories will be copied
        User newUser = new User();
        newUser.setId(2L);

        // Set up mock return values
        when(userRepository.findByLoginId(guestLoginId)).thenReturn(guestUser);
        when(historyService.getHistoryNamesByUserId(guestUserId)).thenReturn(List.of("SampleHistory"));

        // Here, guestHistory will be null
        when(historyService.getHistory("SampleHistory", guestUserId)).thenReturn(null); // Simulate null condition

        // When
        userService.copySampleHistoriesToNewUser(newUser);

        // Then
        verify(analysisService, never()).getAnalysisData(any()); // getJsonData should not be called
        verify(historyService, never()).createHistory(any(History.class));
        verify(analysisService, never()).saveAnalysisData(any(Analysis.class));
    }


}
