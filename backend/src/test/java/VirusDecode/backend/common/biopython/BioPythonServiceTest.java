package VirusDecode.backend.common.biopython;

import VirusDecode.backend.common.exception.BioPythonException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.*;
import org.mockito.junit.jupiter.MockitoExtension;

import java.io.*;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

@ExtendWith(MockitoExtension.class)
class BioPythonServiceTest {

    @Spy
    @InjectMocks
    private BioPythonService bioPythonService;

    @BeforeEach
    void setUp() {
//        MockitoAnnotations.openMocks(this);
    }

    // Helper to mock Process with given stdout, stderr, and exit code
    private Process createMockProcess(String stdout, String stderr, int exitCode) throws Exception {
        Process mockProcess = mock(Process.class);

        InputStream stdOutStream = new ByteArrayInputStream(stdout.getBytes());
        InputStream errOutStream = new ByteArrayInputStream(stderr.getBytes());

        when(mockProcess.getInputStream()).thenReturn(stdOutStream);
        when(mockProcess.getErrorStream()).thenReturn(errOutStream);
        when(mockProcess.waitFor()).thenReturn(exitCode);

        return mockProcess;
    }

    @Test
    @DisplayName("정상 실행 시 BioPythonDto 반환")
    void executePythonScript_ShouldReturnDto_WhenExitCodeIsZero() throws Exception {
        Process mockProcess = createMockProcess("output success", "", 0);
        ProcessBuilder mockBuilder = mock(ProcessBuilder.class);
        when(mockBuilder.start()).thenReturn(mockProcess);

        doReturn(mockBuilder).when(bioPythonService).createProcessBuilder(anyList());

        BioPythonDto result = bioPythonService.executePythonScript("1", "Nc_045512");

        assertEquals(0, result.getExitCode());
        assertTrue(result.getOutput().contains("output success"));
    }

    @Test
    @DisplayName("비정상 종료 시 예외 발생 - ExitCode 11 (ID 없음)")
    void executePythonScript_ShouldThrow_WhenExitCodeIs11() throws Exception {
        Process mockProcess = createMockProcess("", "ID not found", 11);
        ProcessBuilder mockBuilder = mock(ProcessBuilder.class);
        when(mockBuilder.start()).thenReturn(mockProcess);

        doReturn(mockBuilder).when(bioPythonService).createProcessBuilder(anyList());

        BioPythonException ex = assertThrows(BioPythonException.class, () ->
                bioPythonService.executePythonScript("arg1")
        );

        assertTrue(ex.getMessage().contains("NCBI에 요청한 nucleotide ID가 존재하지 않습니다."));
    }

    @Test
    @DisplayName("예외 발생 시 BioPythonException으로 감싸서 반환")
    void executePythonScript_ShouldThrow_WhenIOExceptionOccurs() {
        ProcessBuilder mockBuilder = mock(ProcessBuilder.class);
        try {
            when(mockBuilder.start()).thenThrow(new IOException("파일 실행 실패"));
        } catch (IOException e) {
            // 무시
        }

        doReturn(mockBuilder).when(bioPythonService).createProcessBuilder(anyList());

        BioPythonException ex = assertThrows(BioPythonException.class, () ->
                bioPythonService.executePythonScript("arg1")
        );

        assertTrue(ex.getMessage().contains("파일 실행 실패"));
    }
}
