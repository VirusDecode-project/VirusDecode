package VirusDecode.backend.service;

import VirusDecode.backend.dto.UserLoginDto;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.repository.UserRepository;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import java.util.Optional;

import static org.junit.jupiter.api.Assertions.*;
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
    void testFindUserByUsername() {
        String username = "testUser";
        User user = new User();
        user.setUsername(username);
        when(userRepository.findByUsername(username)).thenReturn(user);

        User result = userService.findUserByUsername(username);

        assertNotNull(result);
        assertEquals(username, result.getUsername());
        verify(userRepository, times(1)).findByUsername(username);
    }

    @Test
    void testCreateUser() {
        UserLoginDto loginDto = new UserLoginDto();
        loginDto.setUsername("newUser");
        loginDto.setPassword("password");

        User user = new User();
        user.setUsername("newUser");
        user.setPassword("password");

        when(userRepository.save(any(User.class))).thenReturn(user);

        User result = userService.createUser(loginDto);

        assertNotNull(result);
        assertEquals("newUser", result.getUsername());
        assertEquals("password", result.getPassword());
        verify(userRepository, times(1)).save(any(User.class));
    }

    @Test
    void testCheckPassword() {
        User user = new User();
        user.setPassword("password");

        boolean result = userService.checkPassword(user, "password");

        assertTrue(result);

        result = userService.checkPassword(user, "wrongPassword");

        assertFalse(result);
    }

    @Test
    void testGetUserById() {
        Long userId = 1L;
        User user = new User();
        when(userRepository.findById(userId)).thenReturn(Optional.of(user));

        Optional<User> result = userService.getUserById(userId);

        assertTrue(result.isPresent());
        assertEquals(user, result.get());
        verify(userRepository, times(1)).findById(userId);
    }
}
