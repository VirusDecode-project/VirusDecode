package VirusDecode.backend.analysis.controller;

import VirusDecode.backend.analysis.dto.LinearDesignDto;
import VirusDecode.backend.analysis.dto.PdbDto;
import VirusDecode.backend.analysis.service.AnalysisService;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.http.MediaType;
import org.springframework.test.web.servlet.MockMvc;

import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.*;

@WebMvcTest(AnalysisController.class)
class AnalysisControllerTest {

    @Autowired
    private MockMvc mockMvc;

    @MockBean
    private AnalysisService analysisService;

    private final ObjectMapper objectMapper = new ObjectMapper();

    @Test
    @DisplayName("LinearDesign 처리 - 성공")
    void getLinearDesign_ShouldReturnJsonString() throws Exception {
        LinearDesignDto dto = new LinearDesignDto();
        dto.setGene("testGene");
        dto.setVarientName("variant1");
        dto.setStart(1);
        dto.setEnd(10);
        dto.setHistoryName("testHistory");
//        dto.setSequence("AUGCGUAGC"); // 예시

        when(analysisService.processLinearDesign(any(LinearDesignDto.class), any(Long.class)))
                .thenReturn("{\"structure\": \"((...))\"}");

        mockMvc.perform(post("/api/analysis/linearDesign")
                        .sessionAttr("userId", 1L)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(objectMapper.writeValueAsString(dto)))
                .andExpect(status().isOk())
                .andExpect(content().json("{\"structure\": \"((...))\"}"));
    }

    @Test
    @DisplayName("PDB 처리 - 성공")
    void getPdb_ShouldReturnJsonString() throws Exception {
        PdbDto dto = new PdbDto();
        dto.setGene("testGene");
        dto.setHistoryName("testHistory");

        when(analysisService.processPdb(any(PdbDto.class), any(Long.class)))
                .thenReturn("{\"pdb\": \"ATOM ...\"}");

        mockMvc.perform(post("/api/analysis/pdb")
                        .sessionAttr("userId", 1L)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(objectMapper.writeValueAsString(dto)))
                .andExpect(status().isOk())
                .andExpect(content().json("{\"pdb\": \"ATOM ...\"}"));
    }
}
