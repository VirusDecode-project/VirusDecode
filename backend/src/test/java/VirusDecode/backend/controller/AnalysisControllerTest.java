package VirusDecode.backend.controller;

import VirusDecode.backend.dto.analysis.LinearDesignDto;
import VirusDecode.backend.dto.analysis.PdbDto;
import VirusDecode.backend.service.AnalysisService;
import jakarta.servlet.http.HttpSession;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;

class AnalysisControllerTest {

    @Mock
    private AnalysisService analysisService;

    @Mock
    private HttpSession session;

    @InjectMocks
    private AnalysisController analysisController;

    @BeforeEach
    void setup() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    void testGetLinearDesign_Success() {
        Long userId = 1L;
        LinearDesignDto request = new LinearDesignDto();
        request.setGene("testGene");
        request.setVarientName("variant1");
        request.setStart(1);
        request.setEnd(10);
        request.setHistoryName("testHistory");

        when(session.getAttribute("userId")).thenReturn(userId);
        when(analysisService.processLinearDesign(any(LinearDesignDto.class), eq(userId)))
                .thenReturn(ResponseEntity.ok("Linear design processed successfully"));

        ResponseEntity<String> response = analysisController.getLinearDesign(request, session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("Linear design processed successfully", response.getBody());
    }

    @Test
    void testGetLinearDesign_Unauthorized() {
        when(session.getAttribute("userId")).thenReturn(null);

        LinearDesignDto request = new LinearDesignDto();
        request.setGene("testGene");
        request.setVarientName("variant1");
        request.setStart(1);
        request.setEnd(10);
        request.setHistoryName("testHistory");

        ResponseEntity<String> response = analysisController.getLinearDesign(request, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("User not authenticated", response.getBody());
    }

    @Test
    void testGetPdb_Success() {
        Long userId = 1L;
        PdbDto request = new PdbDto();
        request.setGene("testGene");
        request.setHistoryName("testHistory");

        when(session.getAttribute("userId")).thenReturn(userId);
        when(analysisService.processPdb(any(PdbDto.class), eq(userId)))
                .thenReturn(ResponseEntity.ok("PDB processed successfully"));

        ResponseEntity<String> response = analysisController.getPdb(request, session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals("PDB processed successfully", response.getBody());
    }

    @Test
    void testGetPdb_Unauthorized() {
        when(session.getAttribute("userId")).thenReturn(null);

        PdbDto request = new PdbDto();
        request.setGene("testGene");
        request.setHistoryName("testHistory");

        ResponseEntity<String> response = analysisController.getPdb(request, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("User not authenticated", response.getBody());
    }
}
