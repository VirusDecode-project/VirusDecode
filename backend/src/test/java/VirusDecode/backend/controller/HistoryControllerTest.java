package VirusDecode.backend.controller;

import VirusDecode.backend.dto.HistoryDto;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.service.JsonDataService;
import com.google.gson.Gson;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.mock.web.MockHttpSession;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.setup.MockMvcBuilders;

import java.util.*;

import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.*;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.*;

class HistoryControllerTest {

    private MockMvc mockMvc;

    @Mock
    private JsonDataService jsonDataService;

    @InjectMocks
    private HistoryController historyController;

    private MockHttpSession mockSession;

    @BeforeEach
    void setUp() {
        MockitoAnnotations.openMocks(this);
        mockMvc = MockMvcBuilders.standaloneSetup(historyController).build();

        // Mock된 세션 생성
        mockSession = new MockHttpSession();
        mockSession.setAttribute("userId", 1L);
    }

    @Test
    void testGetListHistory_Unauthorized() throws Exception {
        mockMvc.perform(get("/history/list"))
                .andExpect(status().isUnauthorized());
    }

    @Test
    void testGetListHistory_Success() throws Exception {
        List<String> historyList = Arrays.asList("History1", "History2");
        when(jsonDataService.getHistoryNamesByUserId(anyLong())).thenReturn(historyList);

        mockMvc.perform(get("/history/list").session(mockSession))
                .andExpect(status().isOk())
                .andExpect(jsonPath("$[0]").value("History1"))
                .andExpect(jsonPath("$[1]").value("History2"));
    }

    @Test
    void testRenameHistory_Unauthorized() throws Exception {
        HistoryDto historyDto = new HistoryDto();
        historyDto.setHistoryName("OldHistory");
        historyDto.setNewName("NewHistory");

        mockMvc.perform(put("/history/rename")
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(new Gson().toJson(historyDto)))
                .andExpect(status().isUnauthorized());
    }

    @Test
    void testRenameHistory_Success() throws Exception {
        HistoryDto historyDto = new HistoryDto();
        historyDto.setHistoryName("OldHistory");
        historyDto.setNewName("NewHistory");

        mockMvc.perform(put("/history/rename")
                        .session(mockSession)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(new Gson().toJson(historyDto)))
                .andExpect(status().isOk())
                .andExpect(content().string("History name updated successfully"));

        verify(jsonDataService, times(1)).updateHistoryName("OldHistory", "NewHistory", 1L);
    }

    @Test
    void testDeleteHistory_Unauthorized() throws Exception {
        HistoryDto historyDto = new HistoryDto();
        historyDto.setHistoryName("HistoryToDelete");

        mockMvc.perform(delete("/history/delete")
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(new Gson().toJson(historyDto)))
                .andExpect(status().isUnauthorized());
    }

    @Test
    void testDeleteHistory_Success() throws Exception {
        HistoryDto historyDto = new HistoryDto();
        historyDto.setHistoryName("HistoryToDelete");

        mockMvc.perform(delete("/history/delete")
                        .session(mockSession)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(new Gson().toJson(historyDto)))
                .andExpect(status().isOk())
                .andExpect(content().string("History deleted successfully"));

        verify(jsonDataService, times(1)).deleteHistory("HistoryToDelete", 1L);
    }

    @Test
    void testGetHistory_Unauthorized() throws Exception {
        HistoryDto historyDto = new HistoryDto();
        historyDto.setHistoryName("SomeHistory");

        mockMvc.perform(post("/history/get")
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(new Gson().toJson(historyDto)))
                .andExpect(status().isUnauthorized());
    }

    @Test
    void testGetHistory_Success() throws Exception {
        HistoryDto historyDto = new HistoryDto();
        historyDto.setHistoryName("SomeHistory");

        JsonData jsonData = new JsonData();
        jsonData.setAlignment("AlignmentData");
        jsonData.setLinearDesign("LinearDesignData");
        jsonData.setPdb("PdbData");

        when(jsonDataService.getJsonData("SomeHistory", 1L)).thenReturn(jsonData);

        mockMvc.perform(post("/history/get")
                        .session(mockSession)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(new Gson().toJson(historyDto)))
                .andExpect(status().isOk())
                .andExpect(jsonPath("$.alignment").value("AlignmentData"))
                .andExpect(jsonPath("$.linearDesign").value("LinearDesignData"))
                .andExpect(jsonPath("$.pdb").value("PdbData"));
    }

    @Test
    void testGetHistory_NoHistory() throws Exception {
        HistoryDto historyDto = new HistoryDto();
        historyDto.setHistoryName("NonExistentHistory");

        when(jsonDataService.getJsonData("NonExistentHistory", 1L)).thenReturn(null);

        mockMvc.perform(post("/history/get")
                        .session(mockSession)
                        .contentType(MediaType.APPLICATION_JSON)
                        .content(new Gson().toJson(historyDto)))
                .andExpect(status().isUnauthorized())
                .andExpect(content().string("There is no history"));
    }
}
