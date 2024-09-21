package VirusDecode.backend.service;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.http.ResponseEntity;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.*;

class PythonScriptExecutorTest {

    private ProcessBuilder processBuilderMock;
    private Process processMock;
    private BufferedReader stdInputMock;
    private BufferedReader stdErrorMock;
    private PythonScriptExecutor pythonScriptExecutor;

    @BeforeEach
    void setUp() {
        // Mocking the ProcessBuilder, Process, and BufferedReader
        processBuilderMock = mock(ProcessBuilder.class);
        processMock = mock(Process.class);
        stdInputMock = mock(BufferedReader.class);
        stdErrorMock = mock(BufferedReader.class);
        pythonScriptExecutor = new PythonScriptExecutor();
    }

    @Test
    void testExecutePythonScriptSuccess() throws Exception {
        // Arrange: Create a mock InputStream for process output
        InputStream inputStreamMock = mock(InputStream.class);
        InputStream errorStreamMock = mock(InputStream.class);

        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(0);

        // Act
        ResponseEntity<String> response = pythonScriptExecutor.executePythonScript("1", "NC_045512");

        // Assert
        assertEquals(200, response.getStatusCodeValue());
    }

    @Test
    void testNotAvailableNucleotide() throws Exception {
        // Arrange: Create a mock InputStream for process output
        InputStream inputStreamMock = mock(InputStream.class);
        InputStream errorStreamMock = mock(InputStream.class);

        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(11);

        // Act
        ResponseEntity<String> response = pythonScriptExecutor.executePythonScript("1", "noExist");

        // Assert: Verify the error message for exit code 11
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("NCBI에 요청한 nucleotide ID가 존재하지 않습니다.", response.getBody());
    }

    @Test
    void testEmptyArgument() throws Exception {
        // Arrange: Create a mock InputStream for process output
        InputStream inputStreamMock = mock(InputStream.class);
        InputStream errorStreamMock = mock(InputStream.class);

        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(11);

        // Act
        ResponseEntity<String> response = pythonScriptExecutor.executePythonScript("1");

        // Assert: Verify the error message for exit code 11
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("전달된 인자가 부족합니다.", response.getBody());
    }

}
