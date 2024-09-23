package VirusDecode.backend.service;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.http.ResponseEntity;

import java.io.InputStream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.*;

class PythonScriptExecutorTest {

    private ProcessBuilder processBuilderMock;
    private Process processMock;
    private PythonScriptExecutor pythonScriptExecutor;

    @BeforeEach
    void setUp() {
        // Mocking the ProcessBuilder, Process, and BufferedReader
        processBuilderMock = mock(ProcessBuilder.class);
        processMock = mock(Process.class);
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
    void testExecutePythonScriptSuccess2() throws Exception {
        // Arrange: Create a mock InputStream for process output
        InputStream inputStreamMock = mock(InputStream.class);
        InputStream errorStreamMock = mock(InputStream.class);

        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(0);

        String fastaContent = "";
        String metadataJson = "{\n" +
                "    \"Sequence ID\": \"NC_001803.1\",\n" +
                "    \"Name\": \"NC_001803\",\n" +
                "    \"Description\": \"Respiratory syncytial virus, complete genome\",\n" +
                "    \"Length\": 15191\n" +
                "}";

        // Act
        ResponseEntity<String> response = pythonScriptExecutor.executePythonScript("2", metadataJson, fastaContent);

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
    void testForCase1() throws Exception {
        // Arrange: Create a mock InputStream for process output
        InputStream inputStreamMock = mock(InputStream.class);
        InputStream errorStreamMock = mock(InputStream.class);

        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(1);

        // Act
        ResponseEntity<String> response = pythonScriptExecutor.executePythonScript("1", anyString());

        // Assert: Verify the error message for exit code 1
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("필요한 파이썬 환경이 제대로 설치되지 않았습니다.\nVirusDecode Github를 참고하세요.", response.getBody());
    }


    @Test
    void testEmptyArgument() throws Exception {
        // Arrange: Create a mock InputStream for process output
        InputStream inputStreamMock = mock(InputStream.class);
        InputStream errorStreamMock = mock(InputStream.class);

        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(2);

        // Act
        ResponseEntity<String> response = pythonScriptExecutor.executePythonScript();

        // Assert: Verify the error message for exit code 11
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("전달된 인자가 부족합니다.", response.getBody());
    }

    @Test
    void testLinearDesignCase31() throws Exception {
        // Arrange: Create a mock InputStream for process output
        InputStream inputStreamMock = mock(InputStream.class);
        InputStream errorStreamMock = mock(InputStream.class);

        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(2);
        String metadataJson = "{\n" +
                "    \"Sequence ID\": \"NC_001803.1\",\n" +
                "    \"Name\": \"NC_001803\",\n" +
                "    \"Description\": \"Respiratory syncytial virus, complete genome\",\n" +
                "    \"Length\": 15191\n" +
                "}";

        String alignmentJson="{\"alignment_index\": {\"ORF1ab\": [0, 13]}, \"aligned_sequences\": {\"NC_001803.1\": \"MESLVPGFNEKTH\",\"MT576556.1\": \"-------------\"}}";

        // Act
        ResponseEntity<String> response = pythonScriptExecutor.executePythonScript("3", metadataJson, alignmentJson, "ORF1ab", "MT576556.1", "1", "10");

        // Assert: Verify the error message for exit code 11
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("선택된 구간에 유효한 서열이 없습니다.", response.getBody());
    }
}
