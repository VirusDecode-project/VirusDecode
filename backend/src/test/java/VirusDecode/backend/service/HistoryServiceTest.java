package VirusDecode.backend.service;

import VirusDecode.backend.entity.History;
import VirusDecode.backend.repository.HistoryRepository;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.*;

class HistoryServiceTest {

    @Mock
    private HistoryRepository historyRepository;

    @InjectMocks
    private HistoryService historyService;

    @BeforeEach
    void setUp() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    void testCreateHistory() {
        // Given
        History history = new History();
        history.setId(1L);
        when(historyRepository.save(any(History.class))).thenReturn(history);

        // When
        History savedHistory = historyService.createHistory(history);

        // Then
        assertNotNull(savedHistory);
        assertEquals(1L, savedHistory.getId());
        verify(historyRepository, times(1)).save(history);
    }

    @Test
    void testGetHistory() {
        // Given
        String historyName = "history1";
        Long userId = 1L;
        History history = new History();
        history.setId(1L);
        history.setHistoryName(historyName);
        when(historyRepository.findByHistoryNameAndUserId(historyName, userId)).thenReturn(history);

        // When
        History retrievedHistory = historyService.getHistory(historyName, userId);

        // Then
        assertNotNull(retrievedHistory);
        assertEquals(1L, retrievedHistory.getId());
        assertEquals(historyName, retrievedHistory.getHistoryName());
        verify(historyRepository, times(1)).findByHistoryNameAndUserId(historyName, userId);
    }

    @Test
    void testGetHistoryNamesByUserId() {
        // Given
        Long userId = 1L;
        List<String> historyNames = Arrays.asList("history1", "history2", "history3");
        when(historyRepository.findHistoryNamesByUserId(userId)).thenReturn(historyNames);

        // When
        List<String> retrievedHistoryNames = historyService.getHistoryNamesByUserId(userId);

        // Then
        assertNotNull(retrievedHistoryNames);
        assertEquals(3, retrievedHistoryNames.size());
        assertEquals("history3", retrievedHistoryNames.get(0));  // Reverse order
        verify(historyRepository, times(1)).findHistoryNamesByUserId(userId);
    }

    @Test
    void testUpdateHistoryName() {
        // Given
        String historyName = "oldHistoryName";
        String newName = "newHistoryName";
        Long userId = 1L;

        // When
        historyService.updateHistoryName(historyName, newName, userId);

        // Then
        verify(historyRepository, times(1)).updateHistoryName(eq(historyName), eq(newName), eq(userId));
    }

    @Test
    void testDeleteHistory() {
        // Given
        String historyName = "historyToDelete";
        Long userId = 1L;

        // When
        historyService.deleteHistory(historyName, userId);

        // Then
        verify(historyRepository, times(1)).deleteByHistoryNameAndUserId(historyName, userId);
    }

    @Test
    void testValidateHistoryName() {
        // Given
        String historyName = "history1";
        Long userId = 1L;

        History existingHistory = new History();
        existingHistory.setHistoryName(historyName);
        when(historyRepository.findByHistoryNameAndUserId(eq("history1"), eq(userId))).thenReturn(existingHistory);
        when(historyRepository.findByHistoryNameAndUserId(eq("history1_1"), eq(userId))).thenReturn(null);

        // When
        String validatedName = historyService.validateHistoryName(historyName, userId);

        // Then
        assertEquals("history1_1", validatedName);
        verify(historyRepository, times(2)).findByHistoryNameAndUserId(anyString(), eq(userId));
    }
}
