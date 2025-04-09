package VirusDecode.backend.analysis.service;

import VirusDecode.backend.common.biopython.BioPythonDto;
import VirusDecode.backend.analysis.repository.AnalysisRepository;
import VirusDecode.backend.common.biopython.BioPythonService;
import VirusDecode.backend.history.entity.History;
import VirusDecode.backend.history.service.HistoryService;
import jakarta.transaction.Transactional;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonArray;
import VirusDecode.backend.analysis.dto.LinearDesignDto;
import VirusDecode.backend.analysis.dto.PdbDto;
import VirusDecode.backend.analysis.entity.Analysis;

@Log4j2
@Service
public class AnalysisService {
    private final BioPythonService bioPythonService;
    private final HistoryService historyService;
    private final AnalysisRepository analysisRepository;

    @Autowired
    public AnalysisService(BioPythonService bioPythonService, HistoryService historyService, AnalysisRepository analysisRepository) {
        this.bioPythonService = bioPythonService;
        this.historyService = historyService;
        this.analysisRepository = analysisRepository;
    }

    public String processLinearDesign(LinearDesignDto request, Long userId) {
        String gene = request.getGene();
        String varientName = request.getVarientName();
        int start = request.getStart();
        int end = request.getEnd();
        String historyName = request.getHistoryName();

        History history = historyService.getHistory(historyName, userId);
        Analysis analysis = getAnalysisData(history);
        if (analysis == null) {
            log.error("There is no history");
            return null;
        }

        String alignmentJson = analysis.getAlignment();
        String aminoAcidSequence = extractAminoAcidSequence(alignmentJson, gene, varientName, start, end);

        if (aminoAcidSequence.length() == 0) {
            log.error("선택된 구간에 유효한 서열이 없습니다.");
            return null;
        }

        BioPythonDto response = bioPythonService.executePythonScript("3", aminoAcidSequence);
        if(response.isSuccess()){
            String linearDesignJson = response.getOutput();
            analysis.setLinearDesign(linearDesignJson);
            saveAnalysisData(analysis);
            return linearDesignJson;
        }

        return null;
    }


    public String processPdb(PdbDto request, Long userId) {
        String gene = request.getGene();
        String historyName = request.getHistoryName();

        History history = historyService.getHistory(historyName, userId);
        Analysis analysis = getAnalysisData(history);
        if (analysis == null) {
            log.error("There is no history");
            return null;
        }

        String referenceId = analysis.getReferenceId();
        String alignmentJson = analysis.getAlignment();

        String sequence = extractSequence(referenceId, alignmentJson, gene);

        BioPythonDto scriptResponse = bioPythonService.executePythonScript("4", sequence);
        if(scriptResponse.isSuccess()){
            String pdbJson = scriptResponse.getOutput();
            analysis.setPdb(pdbJson);
            saveAnalysisData(analysis);
            return pdbJson;
        }
        return null;
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



    @Transactional
    public Analysis saveAnalysisData(Analysis analysis) {
        return analysisRepository.save(analysis);
    }

    public Analysis getAnalysisData(History history) {
        Long historyId = history.getId();
        return analysisRepository.findByHistoryId(historyId);
    }
    @Transactional
    public void deleteAnalysisData(History history){
        Long historyId = history.getId();
        analysisRepository.deleteByHistoryId(historyId);
    }
}
