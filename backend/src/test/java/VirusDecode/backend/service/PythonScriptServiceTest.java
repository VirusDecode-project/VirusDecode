package VirusDecode.backend.service;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;
import org.springframework.http.ResponseEntity;
import java.io.ByteArrayInputStream;
import java.io.InputStream;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.*;

class PythonScriptServiceTest {
    private ProcessBuilder processBuilderMock;
    private Process processMock;
    private PythonScriptService pythonScriptService;
    private InputStream inputStreamMock;
    private InputStream errorStreamMock;
    private PythonScriptService pythonScriptServiceMock;
    @BeforeEach
    void setUp() {
        processBuilderMock = mock(ProcessBuilder.class);
        processMock = mock(Process.class);
        pythonScriptService = new PythonScriptService();
        pythonScriptServiceMock = Mockito.spy(new PythonScriptService());
        inputStreamMock = new ByteArrayInputStream("".getBytes());
        errorStreamMock = new ByteArrayInputStream("".getBytes());
    }

    @Test
    void 메타데이터_성공1() throws Exception {
        // Act
        ResponseEntity<String> response = pythonScriptService.executePythonScript("1", "NC_045512");

        // Assert
        assertEquals(200, response.getStatusCodeValue());
    }

    @Test
    void testExecutePythonScriptSuccess2() throws Exception {
        String fastaContent = "";
        String referenceId = "NC_001803.1";

        // Act
        ResponseEntity<String> response = pythonScriptService.executePythonScript("2", referenceId, fastaContent);

        // Assert
        assertEquals(200, response.getStatusCodeValue());
    }

    @Test
    void args누락() throws Exception {
        // Act
        ResponseEntity<String> response = pythonScriptService.executePythonScript();

        // Assert: Verify the error message for exit code 1
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("Error: No arguments provided.", response.getBody());
    }

    @Test
    void testNotAvailableNucleotide() throws Exception {
        // Act
        ResponseEntity<String> response = pythonScriptService.executePythonScript("1", "noExist");

        // Assert: Verify the error message for exit code 1
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("NCBI에 요청한 nucleotide ID가 존재하지 않습니다.", response.getBody());
    }


    @Test
    void testEmptyArgument() throws Exception {
        // Act
        ResponseEntity<String> response = pythonScriptService.executePythonScript("1");

        // Assert: Verify the error message for exit code 11
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("전달된 인자가 부족합니다.", response.getBody());
    }

    @Test
    void testExitCode1() throws Exception {
        // Mock 설정
        doReturn(processBuilderMock).when(pythonScriptServiceMock).createProcessBuilder(anyList());
        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(1);

        // Act: Python 스크립트 실행
        ResponseEntity<String> response = pythonScriptServiceMock.executePythonScript("2");

        // Assert: 결과 검증
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("필요한 파이썬 환경이 제대로 설치되지 않았습니다.\nVirusDecode Github를 참고하세요.", response.getBody());
    }
    @Test
    void testExitCode21() throws Exception {
        // Mock 설정
        doReturn(processBuilderMock).when(pythonScriptServiceMock).createProcessBuilder(anyList());
        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(21);

        String fastaContent = "";
        String referenceId = "NC_001803.1";

        // Act: Python 스크립트 실행
        ResponseEntity<String> response = pythonScriptServiceMock.executePythonScript("2", referenceId, fastaContent);

        // Assert: 결과 검증
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("MUSCLE 다중 서열 정리에 문제가 발생하였습니다.", response.getBody());
    }
    @Test
    void testExitCode32() throws Exception {
        // Mock 설정
        doReturn(processBuilderMock).when(pythonScriptServiceMock).createProcessBuilder(anyList());
        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(32);

        String fastaContent = "";
        String referenceId = "NC_001803.1";

        // Act: Python 스크립트 실행
        ResponseEntity<String> response = pythonScriptServiceMock.executePythonScript("2", referenceId, fastaContent);

        // Assert: 결과 검증
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("LinearDesign 실행파일이 정상적으로 만들어지지 않았습니다.\nLinux 또는 Max 사용자가 맞으신가요?", response.getBody());
    }
    @Test
    void testExitCode33() throws Exception {
        // Mock 설정
        doReturn(processBuilderMock).when(pythonScriptServiceMock).createProcessBuilder(anyList());
        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(33);

        String referenceId = "NC_001803.1";
        String alignmentJson="{\"alignment_index\": {\"ORF1ab\": [0, 13]}, \"aligned_sequences\": {\"NC_001803.1\": \"MESLVPGFNEKTH\",\"MT576556.1\": \"-------------\"}}";

        // Act
        ResponseEntity<String> response = pythonScriptServiceMock.executePythonScript("3", referenceId, alignmentJson, "ORF1ab", "NC_001803.1", "1", "10");

        // Assert: 결과 검증
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("LinearDesign 디렉토리가 원하는 위치가 존재하지 않습니다.", response.getBody());
    }
    @Test
    void testExitCode41() throws Exception {
        // Mock 설정
        doReturn(processBuilderMock).when(pythonScriptServiceMock).createProcessBuilder(anyList());
        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(41);


        String referenceId = "NC_001803.1";
        String alignmentJson="{\"alignment_index\": {\"ORF1ab\": [0, 13]}, \"aligned_sequences\": {\"NC_001803.1\": \"MESLVPGFNEKTH\",\"MT576556.1\": \"-------------\"}}";
        String gene = "ORF1ab";

        // Act: Python 스크립트 실행
        ResponseEntity<String> response = pythonScriptServiceMock.executePythonScript("4", referenceId, alignmentJson, gene);

        // Assert: 결과 검증
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("RCSB PDB 서버로부터 PDB ID 검색에 실패하였습니다.", response.getBody());
    }
    @Test
    void testExitCode42() throws Exception {
        // Mock 설정
        doReturn(processBuilderMock).when(pythonScriptServiceMock).createProcessBuilder(anyList());
        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(42);

        String referenceId = "NC_001803.1";
        String alignmentJson="{\"alignment_index\": {\"ORF1ab\": [0, 13]}, \"aligned_sequences\": {\"NC_001803.1\": \"MESLVPGFNEKTH\",\"MT576556.1\": \"-------------\"}}";
        String gene = "ORF1ab";

        // Act: Python 스크립트 실행
        ResponseEntity<String> response = pythonScriptServiceMock.executePythonScript("4", referenceId, alignmentJson, gene);

        // Assert: 결과 검증
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("RCSB PDB 서버의 문제로 3D viewer 데이터 로드에 실패하였습니다.", response.getBody());
    }
    @Test
    void testExitCodeDefault() throws Exception {
        // Mock 설정
        doReturn(processBuilderMock).when(pythonScriptServiceMock).createProcessBuilder(anyList());
        when(processBuilderMock.start()).thenReturn(processMock);
        when(processMock.getInputStream()).thenReturn(inputStreamMock);
        when(processMock.getErrorStream()).thenReturn(errorStreamMock);
        when(processMock.waitFor()).thenReturn(5);

        // Act: Python 스크립트 실행
        ResponseEntity<String> response = pythonScriptServiceMock.executePythonScript("1");

        // Assert: 결과 검증
        assertEquals(500, response.getStatusCodeValue());
        assertEquals("Error executing Python script: ", response.getBody());
    }
}
