package VirusDecode.backend.service;

import VirusDecode.backend.dto.initialData.ReferenceDto;
import VirusDecode.backend.dto.initialData.fasta.VarientDto;
import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.entity.User;
import com.google.gson.Gson;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.*;

public class InitialDataServiceTest {

    @Mock
    private PythonScriptService pythonScriptService;

    @Mock
    private FastaFileService fastaFileService;

    @Mock
    private JsonDataService jsonDataService;

    @Mock
    private UserService userService;

    @Mock
    private HistoryService historyService;

    @InjectMocks
    private InitialDataService initialDataService;

    @BeforeEach
    public void setup() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    public void testProcessMetadata_Success() {
        // Given
        ReferenceDto request = new ReferenceDto();
        request.setSequenceId("seq123");

        when(pythonScriptService.executePythonScript("1", "seq123"))
                .thenReturn(ResponseEntity.ok("mockMetadataResponse"));

        // When
        ResponseEntity<String> response = initialDataService.processMetadata(request);

        // Then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("mockMetadataResponse", response.getBody());
        verify(pythonScriptService, times(1)).executePythonScript("1", "seq123");
    }

    @Test
    public void testProcessAlignment_Success() throws IOException {
        // Given
        VarientDto request = new VarientDto();
        request.setHistoryName("history1");
        request.setReferenceId("ref123");

        Long userId = 1L;
        User mockUser = new User();
        mockUser.setId(userId);

        History mockHistory = new History();
        JsonData mockJsonData = new JsonData();
        mockJsonData.setAlignment("{ \"alignment\": \"mockAlignment\" }");

        when(fastaFileService.saveFastaContent(request)).thenReturn("mockFastaContent");
        when(pythonScriptService.executePythonScript(eq("2"), eq("ref123"), eq("mockFastaContent")))
                .thenReturn(ResponseEntity.ok("{ \"alignment\": \"mockAlignment\" }"));
        when(userService.getUserById(userId)).thenReturn(Optional.of(mockUser));
        when(historyService.validateHistoryName("history1", userId)).thenReturn("validatedHistory1");
        when(historyService.createHistory(any(History.class))).thenReturn(mockHistory);
        when(jsonDataService.saveJsonData(any(JsonData.class))).thenReturn(mockJsonData);

        // When
        ResponseEntity<String> response = initialDataService.processAlignment(request, userId);

        // Then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        Map<String, String> expectedResponseMap = new HashMap<>();
        expectedResponseMap.put("historyName", "validatedHistory1");
        expectedResponseMap.put("alignment", "{ \"alignment\": \"mockAlignment\" }");
        String expectedResponse = new Gson().toJson(expectedResponseMap);
        assertEquals(expectedResponse, response.getBody());
        verify(fastaFileService, times(1)).saveFastaContent(request);
        verify(pythonScriptService, times(1)).executePythonScript(eq("2"), eq("ref123"), eq("mockFastaContent"));
        verify(userService, times(1)).getUserById(userId);
        verify(historyService, times(1)).validateHistoryName("history1", userId);
        verify(historyService, times(1)).createHistory(any(History.class));
        verify(jsonDataService, times(1)).saveJsonData(any(JsonData.class));
    }

    @Test
    public void testProcessAlignment_UserNotFound() throws IOException {
        // Given
        VarientDto request = new VarientDto();
        request.setHistoryName("history1");
        request.setReferenceId("ref123");

        Long userId = 1L;

        when(fastaFileService.saveFastaContent(request)).thenReturn("mockFastaContent");
        when(pythonScriptService.executePythonScript(eq("2"), eq("ref123"), eq("mockFastaContent")))
                .thenReturn(ResponseEntity.ok("{ \"alignment\": \"mockAlignment\" }"));
        when(userService.getUserById(userId)).thenReturn(Optional.empty());

        // When
        ResponseEntity<String> response = initialDataService.processAlignment(request, userId);

        // Then
        assertEquals(HttpStatus.NOT_FOUND, response.getStatusCode());
        assertEquals("User not found", response.getBody());
        verify(userService, times(1)).getUserById(userId);
        verify(historyService, never()).createHistory(any());
        verify(jsonDataService, never()).saveJsonData(any());
    }
    @Test
    public void testProcessAlignment_PythonScriptExecutionFailed() throws IOException {
        // Given
        VarientDto request = new VarientDto();
        request.setHistoryName("history1");
        request.setReferenceId("ref123");

        Long userId = 1L;

        when(fastaFileService.saveFastaContent(request)).thenReturn("mockFastaContent");
        when(pythonScriptService.executePythonScript(eq("2"), eq("ref123"), eq("mockFastaContent")))
                .thenReturn(ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Python script execution failed"));

        // When
        ResponseEntity<String> response = initialDataService.processAlignment(request, userId);

        // Then
        assertEquals(HttpStatus.INTERNAL_SERVER_ERROR, response.getStatusCode());
        assertEquals("Python script execution failed", response.getBody());
        verify(pythonScriptService, times(1)).executePythonScript(eq("2"), eq("ref123"), eq("mockFastaContent"));
        verify(userService, never()).getUserById(userId);
        verify(historyService, never()).createHistory(any());
        verify(jsonDataService, never()).saveJsonData(any());
    }


    @Test
    public void testProcessAlignment_FastaIOException() throws IOException {
        // Given
        VarientDto request = new VarientDto();
        request.setHistoryName("history1");
        request.setReferenceId("ref123");

        Long userId = 1L;

        when(fastaFileService.saveFastaContent(request)).thenThrow(new IOException("Mock IO Exception"));

        // When
        ResponseEntity<String> response = initialDataService.processAlignment(request, userId);

        // Then
        assertEquals(HttpStatus.INTERNAL_SERVER_ERROR, response.getStatusCode());
        assertEquals("Fasta 파일 저장에 문제 발생하였습니다.", response.getBody());
        verify(fastaFileService, times(1)).saveFastaContent(request);
        verify(pythonScriptService, never()).executePythonScript(anyString(), anyString(), anyString());
    }
}
