package VirusDecode.backend.controller;

import VirusDecode.backend.dto.ReferenceDTO;
import VirusDecode.backend.dto.VarientDTO;
import VirusDecode.backend.service.FastaFileService;
import VirusDecode.backend.service.PythonScriptExecutor;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;

import java.io.IOException;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

class InputSeqControllerTest {

    @Mock
    private PythonScriptExecutor pythonScriptExecutor;

    @Mock
    private FastaFileService fastaFileService;

    @InjectMocks
    private InputSeqController inputSeqController;

    @BeforeEach
    void setUp() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    void testGetMetadata_Success() {
        // Given
        ReferenceDTO referenceDTO = new ReferenceDTO();
        referenceDTO.setSequenceId("ABC123");

        when(pythonScriptExecutor.executePythonScript("1", "ABC123"))
                .thenReturn(ResponseEntity.ok("Python Script Executed Successfully"));

        // When
        ResponseEntity<String> response = inputSeqController.getMetadata(referenceDTO);

        // Then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("Python Script Executed Successfully", response.getBody());
    }

    @Test
    void testGetMetadata_Failure() {
        // Given
        ReferenceDTO referenceDTO = new ReferenceDTO();
        referenceDTO.setSequenceId("ABC123");

        when(pythonScriptExecutor.executePythonScript("1", "ABC123"))
                .thenReturn(ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Python Script Execution Failed"));

        // When
        ResponseEntity<String> response = inputSeqController.getMetadata(referenceDTO);

        // Then
        assertEquals(HttpStatus.INTERNAL_SERVER_ERROR, response.getStatusCode());
        assertEquals("Python Script Execution Failed", response.getBody());
    }

    @Test
    void testGetAlignment_Success() throws IOException {
        // Given
        VarientDTO varientDTO = new VarientDTO();
        when(fastaFileService.saveFastaContent(any())).thenReturn("FastaContent");
        when(pythonScriptExecutor.executePythonScript("2", "FastaContent"))
                .thenReturn(ResponseEntity.ok("Python Script Executed Successfully"));

        // When
        ResponseEntity<String> response = inputSeqController.getAlignment(varientDTO);

        // Then
        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("Python Script Executed Successfully", response.getBody());
    }

    @Test
    void testGetAlignment_Failure() throws IOException {
        // Given
        VarientDTO varientDTO = new VarientDTO();
        when(fastaFileService.saveFastaContent(any())).thenThrow(new IOException("Fasta file save error"));

        // When
        ResponseEntity<String> response = inputSeqController.getAlignment(varientDTO);

        // Then
        assertEquals(HttpStatus.INTERNAL_SERVER_ERROR, response.getStatusCode());
        assertEquals("Fasta 파일 저장에 문제 발생하였습니다.", response.getBody());
    }
}
