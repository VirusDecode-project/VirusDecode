package VirusDecode.backend.service;

import VirusDecode.backend.dto.SignUpDto;
import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.entity.User;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;

import jakarta.servlet.http.HttpSession;

import java.util.List;
import java.util.Optional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.anyLong;
import static org.mockito.Mockito.*;

public class GuestLoginServiceTest {

    @Mock
    private UserService userService;

    @Mock
    private HistoryService historyService;

    @Mock
    private JsonDataService jsonDataService;

    @Mock
    private HttpSession session;

    @InjectMocks
    private GuestLoginService guestLoginService;

    @BeforeEach
    public void setup() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    public void testLoginAsGuest_ExistingUser() {
        // Given
        String sessionId = "abcdef123456";
        String uniqueLoginId = "Guest_abcdef";
        User existingUser = new User();
        existingUser.setId(1L);
        existingUser.setLoginId(uniqueLoginId);

        when(session.getId()).thenReturn(sessionId);
        when(userService.findUserByLoginId(uniqueLoginId)).thenReturn(existingUser);

        // When
        ResponseEntity<String> response = guestLoginService.loginAsGuest(session);

        // Then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("Existing user logged in with ID: " + uniqueLoginId, response.getBody());
        verify(session, times(1)).setAttribute("userId", existingUser.getId());
        verify(userService, never()).createUser(any(SignUpDto.class), any());
    }

    @Test
    public void testLoginAsGuest_NewUser_NoHistory() {
        // Given
        String sessionId = "abcdef123456";
        String uniqueLoginId = "Guest_abcdef";
        User newUser = new User();
        newUser.setId(2L);
        newUser.setLoginId(uniqueLoginId);

        when(session.getId()).thenReturn(sessionId);
        when(userService.findUserByLoginId(uniqueLoginId)).thenReturn(null);
        when(userService.createUser(any(SignUpDto.class), eq("GUEST"))).thenReturn(newUser);
        when(userService.getUserIdByLoginId("Guest")).thenReturn(1L);
        when(historyService.getHistoryNamesByUserId(1L)).thenReturn(List.of());

        // When
        ResponseEntity<String> response = guestLoginService.loginAsGuest(session);

        // Then
        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("There is no history for Guest user", response.getBody());
        verify(session, never()).setAttribute("userId", newUser.getId());
    }

    @Test
    public void testLoginAsGuest_NewUser_WithHistory() {
// Given
        String sessionId = "abcdef123456";
        String uniqueLoginId = "Guest_abcdef";
        User newUser = new User();
        newUser.setId(2L);
        newUser.setLoginId(uniqueLoginId);
        Long guestUserId = 1L;
        List<String> guestHistoryNames = List.of("history1", "history2");

        when(session.getId()).thenReturn(sessionId);
        when(userService.findUserByLoginId(uniqueLoginId)).thenReturn(null);
        when(userService.createUser(any(SignUpDto.class), eq("GUEST"))).thenReturn(newUser);
        when(userService.getUserIdByLoginId("Guest")).thenReturn(guestUserId);
        when(historyService.getHistoryNamesByUserId(guestUserId)).thenReturn(guestHistoryNames);

        History mockHistory = new History();
        when(historyService.getHistory("history1", guestUserId)).thenReturn(mockHistory);
        when(historyService.getHistory("history2", guestUserId)).thenReturn(mockHistory);
        when(jsonDataService.getJsonData(mockHistory)).thenReturn(new JsonData());

// Here we use when().thenReturn() instead of doNothing()
        when(historyService.createHistory(any(History.class))).thenReturn(new History());
        when(jsonDataService.saveJsonData(any(JsonData.class))).thenReturn(null);

// When
        ResponseEntity<String> response = guestLoginService.loginAsGuest(session);

// Then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("New temporary user created and logged in with ID: " + uniqueLoginId, response.getBody());
        verify(session, times(1)).setAttribute("userId", newUser.getId());
        verify(historyService, times(2)).createHistory(any(History.class));
        verify(jsonDataService, times(2)).saveJsonData(any(JsonData.class));
    }


    @Test
    public void testLoginAsGuest_WithNullOriginalJsonData() {
        // Given
        String sessionId = "session123456";
        Long guestUserId = 2L;
        Long newUserId = 3L;
        String uniqueLoginId = "Guest_sessio";
        when(session.getId()).thenReturn(sessionId);

        User mockGuestUser = new User();
        mockGuestUser.setId(guestUserId);
        when(userService.findUserByLoginId(uniqueLoginId)).thenReturn(null);

        SignUpDto signUpDto = new SignUpDto();
        signUpDto.setLoginId(uniqueLoginId);
        signUpDto.setPassword("default_password");
        signUpDto.setFirstName("Guest");
        signUpDto.setLastName("Guest");

        User mockNewUser = new User();
        mockNewUser.setId(newUserId);
        when(userService.createUser(any(SignUpDto.class), eq("GUEST"))).thenReturn(mockNewUser);

        when(userService.getUserIdByLoginId("Guest")).thenReturn(guestUserId);
        List<String> guestHistoryNames = List.of("history1");
        when(historyService.getHistoryNamesByUserId(guestUserId)).thenReturn(guestHistoryNames);

        History mockHistory = new History();
        when(historyService.getHistory("history1", guestUserId)).thenReturn(mockHistory);

        // Return `null` for `originalJsonData`
        when(jsonDataService.getJsonData(mockHistory)).thenReturn(null);

        // When
        ResponseEntity<String> response = guestLoginService.loginAsGuest(session);

        // Then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("New temporary user created and logged in with ID: " + uniqueLoginId, response.getBody());

        // Verify that saveJsonData is never called because originalJsonData is `null`
        verify(jsonDataService, never()).saveJsonData(any(JsonData.class));
        // Verify that createHistory is called for the new user's history creation
        verify(historyService, never()).createHistory(any(History.class));
        verify(session, times(1)).setAttribute("userId", newUserId);
    }
}
