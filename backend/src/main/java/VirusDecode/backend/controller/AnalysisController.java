package VirusDecode.backend.controller;

import VirusDecode.backend.dto.HistoryDto;
import VirusDecode.backend.dto.analysis.LinearDesignDto;
import VirusDecode.backend.dto.analysis.PdbDto;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.service.JsonDataService;
import VirusDecode.backend.service.PythonScriptService;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import java.nio.file.Path;
import java.nio.file.Paths;

@RestController
@RequestMapping("/analysis")
public class AnalysisController {
    private final PythonScriptService pythonScriptService;
    private final JsonDataService jsonDataService;
    private static final Path currentDir = Paths.get("").toAbsolutePath();
    private static final Path HISTORY_DIR = currentDir.resolve("history");

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

        Long userId = (Long) session.getAttribute("userId");

        String referenceId;
        String alignmentJson;

        JsonData jsonData = jsonDataService.getJsonData(historyName, userId);
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

        Long userId = (Long) session.getAttribute("userId");
        String referenceId;
        String alignmentJson;

        JsonData jsonData = jsonDataService.getJsonData(historyName, userId);
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
