package virusdecode.backend.analysis.service;

import virusdecode.backend.analysis.exception.*;
import virusdecode.backend.analysis.exception.AnalysisNotFoundException;
import virusdecode.backend.analysis.exception.InvalidAnalysisDataException;
import virusdecode.backend.analysis.exception.InvalidSequenceRangeException;
import virusdecode.backend.common.biopython.BioPythonDto;
import virusdecode.backend.analysis.repository.AnalysisRepository;
import virusdecode.backend.common.biopython.BioPythonService;
import virusdecode.backend.history.entity.History;
import virusdecode.backend.history.service.HistoryService;
import jakarta.transaction.Transactional;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.stereotype.Service;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonArray;
import virusdecode.backend.analysis.dto.LinearDesignDto;
import virusdecode.backend.analysis.dto.PdbDto;
import virusdecode.backend.analysis.entity.Analysis;

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

        String alignmentJson = analysis.getAlignment();
        String aminoAcidSequence = extractAminoAcidSequence(alignmentJson, gene, varientName, start, end);

        BioPythonDto response = bioPythonService.executePythonScript("3", aminoAcidSequence);
        String linearDesignJson = response.getOutput();
        analysis.setLinearDesign(linearDesignJson);
        saveAnalysisData(analysis);
        return linearDesignJson;
    }


    public String processPdb(PdbDto request, Long userId) {
        String gene = request.getGene();
        String historyName = request.getHistoryName();

        History history = historyService.getHistory(historyName, userId);
        Analysis analysis = getAnalysisData(history);

        String referenceId = analysis.getReferenceId();
        String alignmentJson = analysis.getAlignment();

        String sequence = extractSequence(referenceId, alignmentJson, gene);

        BioPythonDto scriptResponse = bioPythonService.executePythonScript("4", sequence);
        String pdbJson = scriptResponse.getOutput();
        analysis.setPdb(pdbJson);
        saveAnalysisData(analysis);
        return pdbJson;
    }

    
    public String extractAminoAcidSequence(String alignmentJson, String gene, String varientName, int start, int end) {
        JsonObject jsonObject = JsonParser.parseString(alignmentJson).getAsJsonObject();
        JsonObject alignmentIndex = jsonObject.getAsJsonObject("alignment_index");

        JsonArray gRange = alignmentIndex.getAsJsonArray(gene);
        int startIdx = gRange.get(0).getAsInt();
        int endIdx = gRange.get(1).getAsInt();

        JsonObject alignedSequences = jsonObject.getAsJsonObject("aligned_sequences");
        String sequence = alignedSequences.get(varientName).getAsString();
        String aminoAcidSequence = sequence.substring(startIdx, endIdx).substring(start - 1, end).replace("-", "");
        if(aminoAcidSequence.isEmpty()){
            throw new InvalidSequenceRangeException("선택된 구간에 유효한 서열이 없습니다.");
        }
        return aminoAcidSequence;
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
    public void saveAnalysisData(Analysis analysis) {
        try {
            analysisRepository.save(analysis);
        } catch (DataIntegrityViolationException e) {
            throw new InvalidAnalysisDataException("분석 데이터 저장에 실패했습니다: " + e.getMessage());
        }
    }

    public Analysis getAnalysisData(History history) {
        Long historyId = history.getId();
        Analysis analysis = analysisRepository.findByHistoryId(historyId);
        if(analysis == null){
            throw new AnalysisNotFoundException("분석 정보를 찾을 수 없습니다.");
        }
        return analysis;
    }

    @Transactional
    public void deleteAnalysisData(History history){
        Long historyId = history.getId();
        analysisRepository.deleteByHistoryId(historyId);
    }
}
