package VirusDecode.backend.service;

import VirusDecode.backend.dto.analysis.LinearDesignDto;
import VirusDecode.backend.dto.analysis.PdbDto;
import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.*;

public class AnalysisServiceTest {

    @Mock
    private JsonDataService jsonDataService;

    @Mock
    private PythonScriptService pythonScriptService;

    @Mock
    private HistoryService historyService;

    @InjectMocks
    private AnalysisService analysisService;

    @BeforeEach
    public void setup() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    public void testProcessLinearDesign_Success() {
        // Given
        Long userId = 1L;
        LinearDesignDto request = new LinearDesignDto();
        request.setGene("gene1");
        request.setVarientName("variant1");
        request.setStart(1);
        request.setEnd(10);
        request.setHistoryName("history1");

        History mockHistory = new History();
        JsonData mockJsonData = new JsonData();
        mockJsonData.setAlignment("{ \"alignment_index\": { \"gene1\": [0, 10] }, \"aligned_sequences\": { \"variant1\": \"ATCGATCGAA\" } }");

        when(historyService.getHistory("history1", userId)).thenReturn(mockHistory);
        when(jsonDataService.getJsonData(mockHistory)).thenReturn(mockJsonData);
        when(pythonScriptService.executePythonScript(eq("3"), anyString()))
                .thenReturn(ResponseEntity.ok("{\"linearDesign\": \"mockData\"}"));

        // When
        ResponseEntity<String> response = analysisService.processLinearDesign(request, userId);

        // Then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        verify(jsonDataService, times(1)).saveJsonData(mockJsonData);
        verify(pythonScriptService, times(1)).executePythonScript(eq("3"), anyString());
    }
    @Test
    public void testProcessLinearDesign_InvalidSequence() {
        // Given
        Long userId = 1L;
        LinearDesignDto request = new LinearDesignDto();
        request.setGene("gene1");
        request.setVarientName("variant1");
        request.setStart(1);
        request.setEnd(10);
        request.setHistoryName("history1");

        History mockHistory = new History();
        JsonData mockJsonData = new JsonData();
        mockJsonData.setAlignment("{ \"alignment_index\": { \"gene1\": [0, 10] }, \"aligned_sequences\": { \"variant1\": \"-----------\" } }"); // 빈 서열

        when(historyService.getHistory("history1", userId)).thenReturn(mockHistory);
        when(jsonDataService.getJsonData(mockHistory)).thenReturn(mockJsonData);

        // When
        ResponseEntity<String> response = analysisService.processLinearDesign(request, userId);

        // Then
        assertEquals(HttpStatus.BAD_REQUEST, response.getStatusCode());
        assertEquals("선택된 구간에 유효한 서열이 없습니다.", response.getBody());
        verify(pythonScriptService, times(0)).executePythonScript(anyString(), anyString());
        verify(jsonDataService, times(0)).saveJsonData(any(JsonData.class));
    }

    @Test
    public void testProcessLinearDesign_NoHistory() {
        // Given
        Long userId = 1L;
        LinearDesignDto request = new LinearDesignDto();
        request.setHistoryName("nonExistentHistory");

        when(historyService.getHistory("nonExistentHistory", userId)).thenReturn(null);

        // When
        ResponseEntity<String> response = analysisService.processLinearDesign(request, userId);

        // Then
        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("There is no history", response.getBody());
        verify(jsonDataService, never()).saveJsonData(any());
    }

    @Test
    public void testProcessPdb_Success() {
        // Given
        Long userId = 1L;
        PdbDto request = new PdbDto();
        request.setGene("gene1");
        request.setHistoryName("history1");

        History mockHistory = new History();
        JsonData mockJsonData = new JsonData();
        mockJsonData.setAlignment("{ \"alignment_index\": { \"gene1\": [0, 10] }, \"aligned_sequences\": { \"reference1\": \"ATCGATCGAA\" } }");
        mockJsonData.setReferenceId("reference1");

        when(historyService.getHistory("history1", userId)).thenReturn(mockHistory);
        when(jsonDataService.getJsonData(mockHistory)).thenReturn(mockJsonData);
        when(pythonScriptService.executePythonScript(eq("4"), anyString()))
                .thenReturn(ResponseEntity.ok("{\"pdb\": \"mockData\"}"));

        // When
        ResponseEntity<String> response = analysisService.processPdb(request, userId);

        // Then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        verify(jsonDataService, times(1)).saveJsonData(mockJsonData);
        verify(pythonScriptService, times(1)).executePythonScript(eq("4"), anyString());
    }

    @Test
    public void testProcessPdb_NoHistory() {
        // Given
        Long userId = 1L;
        PdbDto request = new PdbDto();
        request.setHistoryName("nonExistentHistory");

        when(historyService.getHistory("nonExistentHistory", userId)).thenReturn(null);

        // When
        ResponseEntity<String> response = analysisService.processPdb(request, userId);

        // Then
        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("There is no history", response.getBody());
        verify(jsonDataService, never()).saveJsonData(any());
    }

    @Test
    void testProcessPdb_UnsuccessfulScriptResponse() {
        // Arrange
        PdbDto request = new PdbDto();
        request.setGene("gene1");
        request.setHistoryName("history1");

        Long userId = 1L;
        History history = new History();
        history.setId(1L);

        JsonData jsonData = new JsonData();
        jsonData.setReferenceId("ref123");
        jsonData.setAlignment("{ \"alignment_index\": {\"gene1\": [1, 5]}, \"aligned_sequences\": {\"ref123\": \"ACTG-\"} }");

        when(historyService.getHistory("history1", userId)).thenReturn(history);
        when(jsonDataService.getJsonData(history)).thenReturn(jsonData);
        when(pythonScriptService.executePythonScript(eq("4"), anyString()))
                .thenReturn(ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Script error"));

        // Act
        ResponseEntity<String> response = analysisService.processPdb(request, userId);

        // Assert
        assertEquals(HttpStatus.INTERNAL_SERVER_ERROR, response.getStatusCode());
        assertEquals("Script error", response.getBody());
        verify(jsonDataService, never()).saveJsonData(any());
    }
    @Test
    void testProcessLinearDesign_UnsuccessfulScriptResponse() {
        // Arrange
        LinearDesignDto request = new LinearDesignDto();
        request.setGene("gene1");
        request.setVarientName("variant1");
        request.setStart(1);
        request.setEnd(5);
        request.setHistoryName("history1");

        Long userId = 1L;
        History history = new History();
        history.setId(1L);

        JsonData jsonData = new JsonData();
        jsonData.setAlignment("{ \"alignment_index\": {\"gene1\": [0, 5]}, \"aligned_sequences\": {\"variant1\": \"ACTGGGACTGGG--\"} }");

        when(historyService.getHistory("history1", userId)).thenReturn(history);
        when(jsonDataService.getJsonData(history)).thenReturn(jsonData);
        when(pythonScriptService.executePythonScript(eq("3"), anyString()))
                .thenReturn(ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Script error"));

        // Act
        ResponseEntity<String> response = analysisService.processLinearDesign(request, userId);

        // Assert
        assertEquals(HttpStatus.INTERNAL_SERVER_ERROR, response.getStatusCode());
        assertEquals("Script error", response.getBody());
        verify(jsonDataService, never()).saveJsonData(any());
    }

}
