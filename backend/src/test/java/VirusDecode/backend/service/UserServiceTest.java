package VirusDecode.backend.service;

import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.repository.UserRepository;
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
    private JsonDataService jsonDataService;

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
        verify(jsonDataService, times(1)).deleteJsonData(mockHistory);
        verify(historyService, times(1)).deleteHistory("History1", guestUser.getId());
        verify(userRepository, times(1)).deleteUserById(guestUser.getId());
    }
}
