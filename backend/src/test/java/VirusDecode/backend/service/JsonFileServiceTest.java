package VirusDecode.backend.service;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.MockedStatic;
import org.springframework.http.ResponseEntity;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mockStatic;

class JsonFileServiceTest {

    private static final Path testDir = Paths.get("").toAbsolutePath().resolve("data");
    private JsonFileService jsonFileService;

    @BeforeEach
    void setUp() {
        jsonFileService = new JsonFileService();
    }

    @Test
    void testReadJsonFileSuccess() throws IOException {
        // Arrange: Mocking Files.readString and Files.exists to simulate a successful file read
        try (MockedStatic<Files> filesMock = mockStatic(Files.class)) {
            Path jsonFilePath = testDir.resolve("test.json");
            String expectedJsonContent = "{\"key\":\"value\"}";

            filesMock.when(() -> Files.exists(jsonFilePath)).thenReturn(true);
            filesMock.when(() -> Files.readString(jsonFilePath)).thenReturn(expectedJsonContent);

            // Act
            ResponseEntity<String> response = JsonFileService.readJsonFile("test.json");

            // Assert
            assertEquals(200, response.getStatusCodeValue());
            assertEquals(expectedJsonContent, response.getBody());
        }
    }

    @Test
    void testReadJsonFileNotFound() {
        // Arrange: Mocking Files.exists to simulate file not found
        try (MockedStatic<Files> filesMock = mockStatic(Files.class)) {
            Path jsonFilePath = testDir.resolve("test.json");

            filesMock.when(() -> Files.exists(jsonFilePath)).thenReturn(false);

            // Act
            ResponseEntity<String> response = JsonFileService.readJsonFile("test.json");

            // Assert
            assertEquals(404, response.getStatusCodeValue());
            assertEquals("Result JSON file not found", response.getBody());
        }
    }

    @Test
    void testReadJsonFileIOException() throws IOException {
        // Arrange: Mocking Files.readString to throw an IOException
        try (MockedStatic<Files> filesMock = mockStatic(Files.class)) {
            Path jsonFilePath = testDir.resolve("test.json");

            filesMock.when(() -> Files.exists(jsonFilePath)).thenReturn(true);
            filesMock.when(() -> Files.readString(jsonFilePath)).thenThrow(new IOException("Test IO Exception"));

            // Act
            ResponseEntity<String> response = JsonFileService.readJsonFile("test.json");

            // Assert
            assertEquals(500, response.getStatusCodeValue());
            assertEquals("Error reading JSON file", response.getBody());
        }
    }
}
