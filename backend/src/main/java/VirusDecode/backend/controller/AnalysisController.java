package VirusDecode.backend.controller;

import VirusDecode.backend.dto.HistoryDto;
import VirusDecode.backend.dto.analysis.LinearDesignDto;
import VirusDecode.backend.dto.analysis.PdbDto;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.service.JsonDataService;
import VirusDecode.backend.service.PythonScriptService;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

@RestController
@RequestMapping("/api/analysis")
public class AnalysisController {
    private final PythonScriptService pythonScriptService;
    private final JsonDataService jsonDataService;

    @Autowired
    public AnalysisController(PythonScriptService pythonScriptService, JsonDataService jsonDataService) {
        this.pythonScriptService = pythonScriptService;
        this.jsonDataService = jsonDataService;
    }


    @PostMapping("/linearDesign")
    public ResponseEntity<String> getLinearDesign(@RequestBody LinearDesignDto request, HttpSession session) {
        String region = request.getRegion();
        String varientName = request.getVarientName();
        String start = String.valueOf(request.getStart());
        String end = String.valueOf(request.getEnd());
        String historyName = request.getHistoryName();
        String referenceId;
        String alignmentJson;

        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }
        JsonData jsonData = jsonDataService.getJsonData(historyName, userId);
        if (jsonData == null){
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history");
        }

        referenceId = jsonData.getReferenceId();
        alignmentJson = jsonData.getAlignment();

        ResponseEntity<String> scriptResponse = pythonScriptService.executePythonScript("3", referenceId, alignmentJson, region, varientName, start, end);
        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            String linearDesignJson = scriptResponse.getBody();
            jsonData.setLinearDesign(linearDesignJson);
            jsonDataService.saveJsonData(jsonData);
        }
        return scriptResponse;
    }

    @PostMapping("/pdb")
    public ResponseEntity<String> getPdb(@RequestBody PdbDto request, HttpSession session) {
        String gene = request.getGene();
        String historyName = request.getHistoryName();
        String referenceId;
        String alignmentJson;

        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }
        JsonData jsonData = jsonDataService.getJsonData(historyName, userId);
        if (jsonData == null){
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history");
        }

        referenceId = jsonData.getReferenceId();
        alignmentJson = jsonData.getAlignment();

        ResponseEntity<String> scriptResponse = pythonScriptService.executePythonScript("4", referenceId, alignmentJson, gene);
        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            String pdbJson = scriptResponse.getBody();
            jsonData.setPdb(pdbJson);
            jsonDataService.saveJsonData(jsonData);
        }
        return scriptResponse;
    }
}
