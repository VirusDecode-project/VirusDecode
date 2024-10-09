package VirusDecode.backend.service;

import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.repository.JsonDataRepository;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.ArgumentMatchers.any;
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
        // Given
        JsonData jsonData = new JsonData();
        jsonData.setId(1L);
        when(jsonDataRepository.save(any(JsonData.class))).thenReturn(jsonData);

        // When
        JsonData savedJsonData = jsonDataService.saveJsonData(jsonData);

        // Then
        assertNotNull(savedJsonData);
        assertEquals(1L, savedJsonData.getId());
        verify(jsonDataRepository, times(1)).save(jsonData);
    }

    @Test
    void testGetJsonData() {
        // Given
        History history = new History();
        history.setId(1L);
        JsonData jsonData = new JsonData();
        jsonData.setId(1L);
        jsonData.setHistory(history);
        when(jsonDataRepository.findByHistoryId(1L)).thenReturn(jsonData);

        // When
        JsonData retrievedJsonData = jsonDataService.getJsonData(history);

        // Then
        assertNotNull(retrievedJsonData);
        assertEquals(1L, retrievedJsonData.getId());
        verify(jsonDataRepository, times(1)).findByHistoryId(1L);
    }

    @Test
    void testDeleteJsonData() {
        // Given
        History history = new History();
        history.setId(1L);

        // When
        jsonDataService.deleteJsonData(history);

        // Then
        verify(jsonDataRepository, times(1)).deleteByHistoryId(1L);
    }
}
