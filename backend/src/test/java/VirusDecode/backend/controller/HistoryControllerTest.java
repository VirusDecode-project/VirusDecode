package VirusDecode.backend.controller;

import VirusDecode.backend.dto.HistoryDto;
import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.service.HistoryService;
import VirusDecode.backend.service.JsonDataService;
import com.google.gson.Gson;
import jakarta.servlet.http.HttpSession;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;

import java.util.*;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;

class HistoryControllerTest {

    @Mock
    private JsonDataService jsonDataService;

    @Mock
    private HistoryService historyService;

    @Mock
    private HttpSession session;

    @InjectMocks
    private HistoryController historyController;

    @BeforeEach
    void setup() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    void testGetListHistory_Success() {
        Long userId = 1L;
        List<String> mockHistoryList = Arrays.asList("history1", "history2");

        when(session.getAttribute("userId")).thenReturn(userId);
        when(historyService.getHistoryNamesByUserId(userId)).thenReturn(mockHistoryList);

        ResponseEntity<List<String>> response = historyController.getListHistory(session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals(mockHistoryList, response.getBody());
    }

    @Test
    void testGetListHistory_Unauthorized() {
        when(session.getAttribute("userId")).thenReturn(null);

        ResponseEntity<List<String>> response = historyController.getListHistory(session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals(null, response.getBody());
    }

    @Test
    void testRenameHistory_Success() {
        Long userId = 1L;
        HistoryDto request = new HistoryDto();
        request.setHistoryName("oldHistory");
        request.setNewName("newHistory");

        when(session.getAttribute("userId")).thenReturn(userId);

        ResponseEntity<String> response = historyController.renameHistory(request, session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("History name updated successfully", response.getBody());
        verify(historyService, times(1)).updateHistoryName("oldHistory", "newHistory", userId);
    }

    @Test
    void testRenameHistory_Unauthorized() {
        HistoryDto request = new HistoryDto();
        request.setHistoryName("oldHistory");
        request.setNewName("newHistory");

        when(session.getAttribute("userId")).thenReturn(null);

        ResponseEntity<String> response = historyController.renameHistory(request, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("User not authenticated", response.getBody());
    }

    @Test
    void testDeleteHistory_Success() {
        Long userId = 1L;
        String historyName = "testHistory";
        HistoryDto request = new HistoryDto();
        request.setHistoryName(historyName);

        History mockHistory = new History();
        mockHistory.setId(1L);

        when(session.getAttribute("userId")).thenReturn(userId);
        when(historyService.getHistory(historyName, userId)).thenReturn(mockHistory);

        ResponseEntity<String> response = historyController.deleteHistory(request, session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("History deleted successfully", response.getBody());
        verify(jsonDataService, times(1)).deleteJsonData(mockHistory);
        verify(historyService, times(1)).deleteHistory(historyName, userId);
    }

    @Test
    void testDeleteHistory_Unauthorized() {
        HistoryDto request = new HistoryDto();
        request.setHistoryName("testHistory");

        when(session.getAttribute("userId")).thenReturn(null);

        ResponseEntity<String> response = historyController.deleteHistory(request, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("User not authenticated", response.getBody());
    }

    @Test
    void testGetHistory_Success() {
        Long userId = 1L;
        String historyName = "testHistory";
        HistoryDto request = new HistoryDto();
        request.setHistoryName(historyName);

        History mockHistory = new History();
        JsonData mockJsonData = new JsonData();
        mockJsonData.setAlignment("alignmentData");
        mockJsonData.setLinearDesign("linearDesignData");
        mockJsonData.setPdb("pdbData");

        when(session.getAttribute("userId")).thenReturn(userId);
        when(historyService.getHistory(historyName, userId)).thenReturn(mockHistory);
        when(jsonDataService.getJsonData(mockHistory)).thenReturn(mockJsonData);

        ResponseEntity<String> response = historyController.getHistory(request, session);

        Map<String, String> expectedJson = new HashMap<>();
        expectedJson.put("alignment", "alignmentData");
        expectedJson.put("linearDesign", "linearDesignData");
        expectedJson.put("pdb", "pdbData");
        String expectedResponse = new Gson().toJson(expectedJson);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals(expectedResponse, response.getBody());
    }

    @Test
    void testGetHistory_Unauthorized() {
        HistoryDto request = new HistoryDto();
        request.setHistoryName("testHistory");

        when(session.getAttribute("userId")).thenReturn(null);

        ResponseEntity<String> response = historyController.getHistory(request, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("User not authenticated", response.getBody());
    }

    @Test
    void testGetHistory_NoHistoryFound() {
        Long userId = 1L;
        String historyName = "testHistory";
        HistoryDto request = new HistoryDto();
        request.setHistoryName(historyName);

        History mockHistory = new History();

        when(session.getAttribute("userId")).thenReturn(userId);
        when(historyService.getHistory(historyName, userId)).thenReturn(mockHistory);
        when(jsonDataService.getJsonData(mockHistory)).thenReturn(null);

        ResponseEntity<String> response = historyController.getHistory(request, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("There is no history", response.getBody());
    }
}
