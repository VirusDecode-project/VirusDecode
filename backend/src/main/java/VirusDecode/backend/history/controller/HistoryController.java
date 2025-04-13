package VirusDecode.backend.history.controller;

import VirusDecode.backend.history.dto.HistoryDto;
import VirusDecode.backend.history.entity.History;
import VirusDecode.backend.analysis.entity.Analysis;
import VirusDecode.backend.history.service.HistoryService;
import VirusDecode.backend.analysis.service.AnalysisService;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import java.util.*;

import static VirusDecode.backend.util.UserSessionUtil.getAuthenticatedUserId;

@RestController
@RequestMapping("/api/history")
public class HistoryController {
    private final AnalysisService analysisService;
    private final HistoryService historyService;

    @Autowired
    public HistoryController(AnalysisService analysisService, HistoryService historyService) {
        this.analysisService = analysisService;
        this.historyService = historyService;
    }

    @GetMapping("/list")
    public ResponseEntity<List<String>> getListHistory(HttpSession session) {
        Long userId = getAuthenticatedUserId(session);

        List<String> historyList = historyService.getHistoryNamesByUserId(userId);
        return ResponseEntity.ok(historyList);
    }

    @PutMapping("/rename")
    public ResponseEntity<String> renameHistory(@RequestBody HistoryDto request, HttpSession session) {
        Long userId = getAuthenticatedUserId(session);

        historyService.updateHistoryName(request.getHistoryName(), request.getNewName(), userId);
        return ResponseEntity.ok("History name updated successfully");
    }

    @DeleteMapping("/delete")
    public ResponseEntity<String> deleteHistory(@RequestBody HistoryDto request, HttpSession session) {
        Long userId = getAuthenticatedUserId(session);

        String historyName = request.getHistoryName();
        History history = historyService.getHistory(historyName, userId);
        analysisService.deleteAnalysisData(history);
        historyService.deleteHistory(historyName, userId);

        return ResponseEntity.ok("History deleted successfully");
    }

    @PostMapping("/get")
    public ResponseEntity<?> getHistory(@RequestBody HistoryDto request, HttpSession session) {
        Long userId = getAuthenticatedUserId(session);

        String historyName = request.getHistoryName();
        History history = historyService.getHistory(historyName, userId);
        Analysis analysis = analysisService.getAnalysisData(history);

        // JSON 데이터 통합
        Map<String, String> historyJson = new HashMap<>();
        historyJson.put("alignment", analysis.getAlignment());
        historyJson.put("linearDesign", analysis.getLinearDesign());
        historyJson.put("pdb", analysis.getPdb());


        return ResponseEntity.ok(historyJson);
    }
}