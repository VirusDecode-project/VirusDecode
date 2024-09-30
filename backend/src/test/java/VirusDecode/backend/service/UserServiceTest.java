package VirusDecode.backend.service;

import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.repository.UserRepository;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import java.util.Optional;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;

class UserServiceTest {

    @Mock
    private UserRepository userRepository;

    @InjectMocks
    private UserService userService;

    @BeforeEach
    void setUp() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    void testFindUserByLoginId_UserExists() {
        // Given
        User user = new User();
        user.setLoginId("testUser");
        when(userRepository.findByLoginId("testUser")).thenReturn(user);

        // When
        User result = userService.findUserByLoginId("testUser");

        // Then
        assertNotNull(result);
        assertEquals("testUser", result.getLoginId());
        verify(userRepository, times(1)).findByLoginId("testUser");
    }

    @Test
    void testFindUserByLoginId_UserDoesNotExist() {
        // Given
        when(userRepository.findByLoginId("nonExistentUser")).thenReturn(null);

        // When
        User result = userService.findUserByLoginId("nonExistentUser");

        // Then
        assertNull(result);
        verify(userRepository, times(1)).findByLoginId("nonExistentUser");
    }

//    @Test
//    void testCreateUser_Success() {
//        // Given
//        SignUpDto signUpDto = new SignUpDto();
//        signUpDto.setFirstName("John");
//        signUpDto.setLastName("Doe");
//        signUpDto.setLoginId("johndoe123");
//        signUpDto.setPassword("securePassword");
//
//        User user = new User();
//        user.setFirstName("John");
//        user.setLastName("Doe");
//        user.setLoginId("johndoe123");
//        user.setPassword("securePassword");
//
//        when(userRepository.save(any(User.class))).thenReturn(user);
//
//        // When
//        User result = userService.createUser(signUpDto);
//
//        // Then
//        assertNotNull(result);
//        assertEquals("johndoe123", result.getLoginId());
//        assertEquals("securePassword", result.getPassword());
//        verify(userRepository, times(1)).save(any(User.class));
//    }
//
//    @Test
//    void testCheckPassword_ValidPassword() {
//        // Given
//        User user = new User();
//        user.setPassword("correctPassword");
//
//        // When
//        boolean result = userService.checkPassword(user, "correctPassword");
//
//        // Then
//        assertTrue(result);
//    }

//    @Test
//    void testCheckPassword_InvalidPassword() {
//        // Given
//        User user = new User();
//        user.setPassword("correctPassword");
//
//        // When
//        boolean result = userService.checkPassword(user, "wrongPassword");
//
//        // Then
//        assertFalse(result);
//    }

    @Test
    void testGetUserById_UserExists() {
        // Given
        User user = new User();
        user.setId(1L);
        when(userRepository.findById(1L)).thenReturn(Optional.of(user));

        // When
        Optional<User> result = userService.getUserById(1L);

        // Then
        assertTrue(result.isPresent());
        assertEquals(1L, result.get().getId());
        verify(userRepository, times(1)).findById(1L);
    }

    @Test
    void testGetUserById_UserDoesNotExist() {
        // Given
        when(userRepository.findById(1L)).thenReturn(Optional.empty());

        // When
        Optional<User> result = userService.getUserById(1L);

        // Then
        assertFalse(result.isPresent());
        verify(userRepository, times(1)).findById(1L);
    }
}
