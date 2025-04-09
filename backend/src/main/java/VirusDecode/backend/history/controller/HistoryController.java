package VirusDecode.backend.history.controller;

import VirusDecode.backend.history.dto.HistoryDto;
import VirusDecode.backend.history.entity.History;
import VirusDecode.backend.analysis.entity.Analysis;
import VirusDecode.backend.history.service.HistoryService;
import VirusDecode.backend.analysis.service.AnalysisService;
import com.google.gson.Gson;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import java.util.*;

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
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body(null);
        }

        List<String> historyList = historyService.getHistoryNamesByUserId(userId);
        return ResponseEntity.ok(historyList);
    }

    @PutMapping("/rename")
    public ResponseEntity<String> renameHistory(@RequestBody HistoryDto request, HttpSession session) {
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }
        historyService.updateHistoryName(request.getHistoryName(), request.getNewName(), userId);
        return ResponseEntity.ok("History name updated successfully");
    }

    @DeleteMapping("/delete")
    public ResponseEntity<String> deleteHistory(@RequestBody HistoryDto request, HttpSession session) {
        String historyName = request.getHistoryName();
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }

        History history = historyService.getHistory(historyName, userId);
        if(history == null){
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history");
        }
        analysisService.deleteAnalysisData(history);
        historyService.deleteHistory(historyName, userId);

        return ResponseEntity.ok("History deleted successfully");
    }

    @PostMapping("/get")
    public ResponseEntity<String> getHistory(@RequestBody HistoryDto request, HttpSession session) {
        String historyName = request.getHistoryName();
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }

        History history = historyService.getHistory(historyName, userId);
        if(history == null){
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history");
        }

        Analysis analysis = analysisService.getAnalysisData(history);
        if (analysis == null){
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history");
        }

        // JSON 데이터 통합
        Map<String, String> combinedJson = new HashMap<>();
        combinedJson.put("alignment", analysis.getAlignment());
        combinedJson.put("linearDesign", analysis.getLinearDesign());
        combinedJson.put("pdb", analysis.getPdb());

        // JSON 문자열로 변환
        String jsonResponse = new Gson().toJson(combinedJson);

        return ResponseEntity.ok(jsonResponse);
    }
}