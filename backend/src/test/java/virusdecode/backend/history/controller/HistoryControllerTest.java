package virusdecode.backend.history.controller;

import virusdecode.backend.analysis.entity.Analysis;
import virusdecode.backend.history.dto.HistoryDto;
import virusdecode.backend.history.entity.History;
import virusdecode.backend.history.service.HistoryService;
import virusdecode.backend.analysis.service.AnalysisService;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.http.MediaType;
import org.springframework.test.web.servlet.MockMvc;

import java.util.*;

import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.*;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.*;

@WebMvcTest(HistoryController.class)
class HistoryControllerTest {

    @Autowired
    private MockMvc mockMvc;

    @MockBean
    private HistoryService historyService;

    @MockBean
    private AnalysisService analysisService;

    private final ObjectMapper objectMapper = new ObjectMapper();

    @Test
    @DisplayName("히스토리 목록 조회 - 성공")
    void getListHistory_ShouldReturnHistoryList() throws Exception {
        List<String> mockList = Arrays.asList("History1", "History2");
        Mockito.when(historyService.getHistoryNamesByUserId(1L)).thenReturn(mockList);

        mockMvc.perform(get("/api/history/list")
                        .sessionAttr("userId", 1L))
                .andExpect(status().isOk())
                .andExpect(jsonPath("$[0]").value("History1"))
                .andExpect(jsonPath("$[1]").value("History2"));
    }

    @Test
    @DisplayName("히스토리 이름 변경 - 성공")
    void renameHistory_ShouldUpdateName() throws Exception {
        HistoryDto dto = new HistoryDto();
        dto.setHistoryName("OldName");
        dto.setNewName("NewName");

        mockMvc.perform(put("/api/history/rename")
                        .sessionAttr("userId", 1L)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(objectMapper.writeValueAsString(dto)))
                .andExpect(status().isOk())
                .andExpect(content().string("History name updated successfully"));

        Mockito.verify(historyService).updateHistoryName("OldName", "NewName", 1L);
    }

    @Test
    @DisplayName("히스토리 삭제 - 성공")
    void deleteHistory_ShouldDeleteSuccessfully() throws Exception {
        HistoryDto dto = new HistoryDto();
        dto.setHistoryName("HistoryToDelete");

        History mockHistory = new History();
        Mockito.when(historyService.getHistory("HistoryToDelete", 1L)).thenReturn(mockHistory);

        mockMvc.perform(delete("/api/history/delete")
                        .sessionAttr("userId", 1L)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(objectMapper.writeValueAsString(dto)))
                .andExpect(status().isOk())
                .andExpect(content().string("History deleted successfully"));

        Mockito.verify(analysisService).deleteAnalysisData(mockHistory);
        Mockito.verify(historyService).deleteHistory("HistoryToDelete", 1L);
    }

    @Test
    @DisplayName("히스토리 분석 데이터 조회 - 성공")
    void getHistory_ShouldReturnAnalysisData() throws Exception {
        HistoryDto dto = new HistoryDto();
        dto.setHistoryName("History1");

        History mockHistory = new History();
        Analysis mockAnalysis = new Analysis();
        mockAnalysis.setAlignment("alignmentData");
        mockAnalysis.setLinearDesign("linearData");
        mockAnalysis.setPdb("pdbData");

        Mockito.when(historyService.getHistory("History1", 1L)).thenReturn(mockHistory);
        Mockito.when(analysisService.getAnalysisData(mockHistory)).thenReturn(mockAnalysis);

        mockMvc.perform(post("/api/history/get")
                        .sessionAttr("userId", 1L)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(objectMapper.writeValueAsString(dto)))
                .andExpect(status().isOk())
                .andExpect(jsonPath("$.alignment").value("alignmentData"))
                .andExpect(jsonPath("$.linearDesign").value("linearData"))
                .andExpect(jsonPath("$.pdb").value("pdbData"));
    }
}
