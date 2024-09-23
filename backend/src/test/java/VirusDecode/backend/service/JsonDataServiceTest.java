package VirusDecode.backend.service;

import VirusDecode.backend.entity.JsonDataEntity;
import VirusDecode.backend.repository.JsonDataRepository;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import java.util.Optional;

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
        // Arrange
        String key = "testKey";
        String jsonData = "{\"name\":\"test\"}";
        JsonDataEntity jsonDataEntity = new JsonDataEntity(key, jsonData);

        // Act
        jsonDataService.saveJsonData(key, jsonData);

        // Assert
        verify(jsonDataRepository, times(1)).save(jsonDataEntity);
    }

    @Test
    void testGetJsonDataExists() {
        // Arrange
        String key = "testKey";
        String jsonData = "{\"name\":\"test\"}";
        JsonDataEntity jsonDataEntity = new JsonDataEntity(key, jsonData);
        when(jsonDataRepository.findById(key)).thenReturn(Optional.of(jsonDataEntity));

        // Act
        String result = jsonDataService.getJsonData(key);

        // Assert
        assertEquals(jsonData, result);
    }

    @Test
    void testGetJsonDataNotExists() {
        // Arrange
        String key = "nonExistentKey";
        when(jsonDataRepository.findById(key)).thenReturn(Optional.empty());

        // Act
        String result = jsonDataService.getJsonData(key);

        // Assert
        assertNull(result);
    }

    @Test
    void testGetAllJsonData() {
        // Act
        jsonDataService.getAllJsonData();

        // Assert
        verify(jsonDataRepository, times(1)).findAll();
    }

    @Test
    void testDeleteAllJsonData() {
        // Act
        jsonDataService.deleteAllJsonData();

        // Assert
        verify(jsonDataRepository, times(1)).deleteAll();
    }

    @Test
    void testDeleteJsonData() {
        // Arrange
        String key = "testKey";

        // Act
        jsonDataService.deleteJsonData(key);

        // Assert
        verify(jsonDataRepository, times(1)).deleteById(key);
    }

    @Test
    void testFindByIdExists() {
        // Arrange
        String key = "testKey";
        JsonDataEntity jsonDataEntity = new JsonDataEntity(key, "{\"name\":\"test\"}");
        when(jsonDataRepository.findById(key)).thenReturn(Optional.of(jsonDataEntity));

        // Act
        Optional<JsonDataEntity> result = jsonDataService.findById(key);

        // Assert
        assertTrue(result.isPresent());
        assertEquals(jsonDataEntity, result.get());
    }

    @Test
    void testFindByIdNotExists() {
        // Arrange
        String key = "nonExistentKey";
        when(jsonDataRepository.findById(key)).thenReturn(Optional.empty());

        // Act
        Optional<JsonDataEntity> result = jsonDataService.findById(key);

        // Assert
        assertFalse(result.isPresent());
    }
}
