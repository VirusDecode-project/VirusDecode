package VirusDecode.backend.service;

import VirusDecode.backend.entity.History;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonArray;
import VirusDecode.backend.dto.analysis.LinearDesignDto;
import VirusDecode.backend.dto.analysis.PdbDto;
import VirusDecode.backend.entity.JsonData;

@Service
public class AnalysisService {
    private final JsonDataService jsonDataService;
    private final PythonScriptService pythonScriptService;
    private final HistoryService historyService;

    @Autowired
    public AnalysisService(JsonDataService jsonDataService, PythonScriptService pythonScriptService, HistoryService historyService) {
        this.jsonDataService = jsonDataService;
        this.pythonScriptService = pythonScriptService;
        this.historyService = historyService;
    }

    public ResponseEntity<String> processLinearDesign(LinearDesignDto request, Long userId) {
        String gene = request.getGene();
        String varientName = request.getVarientName();
        int start = request.getStart();
        int end = request.getEnd();
        String historyName = request.getHistoryName();

        History history = historyService.getHistory(historyName, userId);
        JsonData jsonData = jsonDataService.getJsonData(history);
        if (jsonData == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history");
        }

        String alignmentJson = jsonData.getAlignment();
        String aminoAcidSequence = extractAminoAcidSequence(alignmentJson, gene, varientName, start, end);

        if (aminoAcidSequence.length() == 0) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body("선택된 구간에 유효한 서열이 없습니다.");
        }

        ResponseEntity<String> scriptResponse = pythonScriptService.executePythonScript("3", aminoAcidSequence);
        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            String linearDesignJson = scriptResponse.getBody();
            jsonData.setLinearDesign(linearDesignJson);
            jsonDataService.saveJsonData(jsonData);
        }

        return scriptResponse;
    }


    public ResponseEntity<String> processPdb(PdbDto request, Long userId) {
        String gene = request.getGene();
        String historyName = request.getHistoryName();

        History history = historyService.getHistory(historyName, userId);
        JsonData jsonData = jsonDataService.getJsonData(history);
        if (jsonData == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("There is no history");
        }

        String referenceId = jsonData.getReferenceId();
        String alignmentJson = jsonData.getAlignment();

        String sequence = extractSequence(referenceId, alignmentJson, gene);

        ResponseEntity<String> scriptResponse = pythonScriptService.executePythonScript("4", sequence);
        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            String pdbJson = scriptResponse.getBody();
            jsonData.setPdb(pdbJson);
            jsonDataService.saveJsonData(jsonData);
        }

        return scriptResponse;
    }

    
    public String extractAminoAcidSequence(String alignmentJson, String gene, String varientName, int start, int end) {
        JsonObject jsonObject = JsonParser.parseString(alignmentJson).getAsJsonObject();
        JsonObject alignmentIndex = jsonObject.getAsJsonObject("alignment_index");

        JsonArray gRange = alignmentIndex.getAsJsonArray(gene);
        int startIdx = gRange.get(0).getAsInt();
        int endIdx = gRange.get(1).getAsInt();

        JsonObject alignedSequences = jsonObject.getAsJsonObject("aligned_sequences");
        String sequence = alignedSequences.get(varientName).getAsString();

        return sequence.substring(startIdx, endIdx).substring(start - 1, end).replace("-", "");
    }

    public String extractSequence(String referenceId, String alignmentJson, String gene) {
        JsonObject jsonObject = JsonParser.parseString(alignmentJson).getAsJsonObject();
        JsonObject alignmentIndex = jsonObject.getAsJsonObject("alignment_index");

        JsonArray gRange = alignmentIndex.getAsJsonArray(gene);
        int startIdx = gRange.get(0).getAsInt();
        int endIdx = gRange.get(1).getAsInt();

        JsonObject alignedSequences = jsonObject.getAsJsonObject("aligned_sequences");
        String sequence = alignedSequences.get(referenceId).getAsString();

        return sequence.substring(startIdx, endIdx).replace("-", "");
    }
}
