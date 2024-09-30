package VirusDecode.backend.controller;

import VirusDecode.backend.dto.HistoryDto;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.service.JsonDataService;
import com.google.gson.Gson;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import java.util.*;

@RestController
@RequestMapping("/history")
public class HistoryController {
    private final JsonDataService jsonDataService;

    @Autowired
    public HistoryController(JsonDataService jsonDataService) {
        this.jsonDataService = jsonDataService;
    }

    @GetMapping("/list")
    public ResponseEntity<List<String>> getListHistory() {
        List<String> historyList = jsonDataService.getHistoryNames();
        return ResponseEntity.ok(historyList);
    }

    @PutMapping("/rename")
    public ResponseEntity<String> renameHistory(@RequestBody HistoryDto request) {
        jsonDataService.updateHistoryName(request.getHistoryName(), request.getNewName());
        return ResponseEntity.ok("History name updated successfully");
    }

    @DeleteMapping("/delete")
    public ResponseEntity<String> deleteHistory(@RequestBody HistoryDto request) {
        String historyName = request.getHistoryName();

        jsonDataService.deleteHistory(historyName);
        return ResponseEntity.ok("History deleted successfully");
    }

    @PostMapping("/get")
    public ResponseEntity<String> getHistory(@RequestBody HistoryDto request) {
        String historyName = request.getHistoryName();
        JsonData jsonData = jsonDataService.getJsonData(historyName);
        if (jsonData == null){
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history");
        }

        // JSON 데이터 통합
        Map<String, String> combinedJson = new HashMap<>();
        combinedJson.put("alignment", jsonData.getAlignment());
        combinedJson.put("linearDesign", jsonData.getLinearDesign());
        combinedJson.put("pdb", jsonData.getPdb());

        // JSON 문자열로 변환
        String jsonResponse = new Gson().toJson(combinedJson);

        return ResponseEntity.ok(jsonResponse);
    }
}