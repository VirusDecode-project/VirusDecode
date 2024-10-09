package VirusDecode.backend.controller;

import VirusDecode.backend.dto.initialData.ReferenceDto;
import VirusDecode.backend.dto.initialData.fasta.VarientDto;
import VirusDecode.backend.service.InitialDataService;
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

class InitialDataControllerTest {

    @Mock
    private InitialDataService initialDataService;

    @Mock
    private HttpSession session;

    @InjectMocks
    private InitialDataController initialDataController;

    @BeforeEach
    void setup() {
        MockitoAnnotations.openMocks(this);
    }

    @Test
    void testGetMetadata_Success() {
        ReferenceDto request = new ReferenceDto();
        request.setSequenceId("sequence123");
        String expectedResponse = "metadata response";

        when(initialDataService.processMetadata(any(ReferenceDto.class)))
                .thenReturn(ResponseEntity.ok(expectedResponse));

        ResponseEntity<String> response = initialDataController.getMetadata(request);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals(expectedResponse, response.getBody());
        verify(initialDataService, times(1)).processMetadata(any(ReferenceDto.class));
    }

    @Test
    void testGetAlignment_Success() {
        Long userId = 1L;
        VarientDto request = new VarientDto();
        String expectedResponse = "alignment response";

        when(session.getAttribute("userId")).thenReturn(userId);
        when(initialDataService.processAlignment(any(VarientDto.class), eq(userId)))
                .thenReturn(ResponseEntity.ok(expectedResponse));

        ResponseEntity<String> response = initialDataController.getAlignment(request, session);

        assertEquals(HttpStatus.OK, response.getStatusCode());
        assertEquals(expectedResponse, response.getBody());
        verify(initialDataService, times(1)).processAlignment(any(VarientDto.class), eq(userId));
    }

    @Test
    void testGetAlignment_Unauthorized() {
        when(session.getAttribute("userId")).thenReturn(null);

        VarientDto request = new VarientDto();
        ResponseEntity<String> response = initialDataController.getAlignment(request, session);

        assertEquals(HttpStatus.UNAUTHORIZED, response.getStatusCode());
        assertEquals("User not authenticated", response.getBody());
        verify(initialDataService, never()).processAlignment(any(VarientDto.class), anyLong());
    }
}
