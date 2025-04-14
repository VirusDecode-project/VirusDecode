package virusdecode.backend.bioinput.controller;

import virusdecode.backend.analysis.dto.AlignmentDto;
import virusdecode.backend.bioinput.dto.ReferenceDto;
import virusdecode.backend.bioinput.dto.VarientSequenceDto;
import virusdecode.backend.bioinput.entity.MetaData;
import virusdecode.backend.bioinput.service.BioInputService;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.http.MediaType;
import org.springframework.test.web.servlet.MockMvc;

import java.util.LinkedHashMap;
import java.util.Map;

import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.*;

@WebMvcTest(BioInputController.class)
class BioInputControllerTest {

    @Autowired
    private MockMvc mockMvc;

    @MockBean
    private BioInputService bioInputService;

    private final ObjectMapper objectMapper = new ObjectMapper();

    @Test
    @DisplayName("메타데이터 요청 - 성공")
    void getMetadata_ShouldReturnMetaDataString() throws Exception {
        ReferenceDto dto = new ReferenceDto("REF123");
        MetaData mockMetaData = new MetaData("{\"meta\":true}", "REF123");

        Mockito.when(bioInputService.getMetadata(any(ReferenceDto.class))).thenReturn(mockMetaData);

        mockMvc.perform(post("/api/inputSeq/metadata")
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(objectMapper.writeValueAsString(dto)))
                .andExpect(status().isOk())
                .andExpect(content().json("{\"meta\":true}"));
    }

    @Test
    @DisplayName("alignment 요청 - 성공")
    void getAlignment_ShouldReturnAlignmentDto() throws Exception {
        VarientSequenceDto dto = new VarientSequenceDto();
        dto.setHistoryName("testHistory");
        dto.setReferenceId("REF123");

        Map<String, String> sequences = new LinkedHashMap<>();
        sequences.put("var1", "ATCG");
        dto.setSequences(sequences);

        AlignmentDto responseDto = new AlignmentDto("{alignment:true}", "validatedHistory");

        Mockito.when(bioInputService.processAlignment(any(VarientSequenceDto.class), eq(1L)))
                .thenReturn(responseDto);

        mockMvc.perform(post("/api/inputSeq/alignment")
                        .sessionAttr("userId", 1L)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(objectMapper.writeValueAsString(dto)))
                .andExpect(status().isOk())
                .andExpect(jsonPath("$.alignment").value("{alignment:true}"))
                .andExpect(jsonPath("$.historyName").value("validatedHistory"));
    }
}
