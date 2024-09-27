package VirusDecode.backend.service;

import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.repository.JsonDataRepository;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

class JsonDataServiceTest {

    @Mock
    private JsonDataRepository jsonDataRepository;

    @InjectMocks
    private JsonDataService jsonDataService;

    @BeforeEach
    void setUp() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    void testSaveJsonData() {
        JsonData jsonData = new JsonData();
        when(jsonDataRepository.save(any(JsonData.class))).thenReturn(jsonData);

        JsonData result = jsonDataService.saveJsonData(jsonData);

        assertNotNull(result);
        verify(jsonDataRepository, times(1)).save(jsonData);
    }

    @Test
    void testGetJsonData() {
        String historyName = "testHistory";
        Long userId = 1L;
        JsonData jsonData = new JsonData();
        when(jsonDataRepository.findByHistoryNameAndUserId(historyName, userId)).thenReturn(jsonData);

        JsonData result = jsonDataService.getJsonData(historyName, userId);

        assertNotNull(result);
        assertEquals(jsonData, result);
        verify(jsonDataRepository, times(1)).findByHistoryNameAndUserId(historyName, userId);
    }

    @Test
    void testGetHistoryNamesByUserId() {
        Long userId = 1L;
        List<String> historyNames = Arrays.asList("History1", "History2", "History3");
        when(jsonDataRepository.findHistoryNamesByUserId(userId)).thenReturn(historyNames);

        List<String> result = jsonDataService.getHistoryNamesByUserId(userId);

        assertNotNull(result);
        assertEquals(3, result.size());
        assertEquals("History3", result.get(0)); // Collections.reverse가 작동하는지 확인
        verify(jsonDataRepository, times(1)).findHistoryNamesByUserId(userId);
    }

    @Test
    void testUpdateHistoryName() {
        String oldName = "oldName";
        String newName = "newName";
        Long userId = 1L;
        doNothing().when(jsonDataRepository).updateHistoryName(oldName, newName, userId);

        jsonDataService.updateHistoryName(oldName, newName, userId);

        verify(jsonDataRepository, times(1)).updateHistoryName(oldName, newName, userId);
    }

    @Test
    void testDeleteHistory() {
        String historyName = "testHistory";
        Long userId = 1L;
        doNothing().when(jsonDataRepository).deleteByHistoryNameAndUserId(historyName, userId);

        jsonDataService.deleteHistory(historyName, userId);

        verify(jsonDataRepository, times(1)).deleteByHistoryNameAndUserId(historyName, userId);
    }
}
