package VirusDecode.backend.analysis.service;

import VirusDecode.backend.analysis.dto.LinearDesignDto;
import VirusDecode.backend.analysis.dto.PdbDto;
import VirusDecode.backend.analysis.entity.Analysis;
import VirusDecode.backend.analysis.exception.AnalysisNotFoundException;
import VirusDecode.backend.analysis.exception.InvalidSequenceRangeException;
import VirusDecode.backend.analysis.repository.AnalysisRepository;
import VirusDecode.backend.common.biopython.BioPythonDto;
import VirusDecode.backend.common.biopython.BioPythonService;
import VirusDecode.backend.history.entity.History;
import VirusDecode.backend.history.service.HistoryService;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.*;
import org.mockito.junit.jupiter.MockitoExtension;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

@ExtendWith(MockitoExtension.class)
class AnalysisServiceTest {

    @InjectMocks
    private AnalysisService analysisService;

    @Mock
    private BioPythonService bioPythonService;

    @Mock
    private HistoryService historyService;

    @Mock
    private AnalysisRepository analysisRepository;

    @BeforeEach
    void setUp() {
//        MockitoAnnotations.openMocks(this);
    }

    @Test
    @DisplayName("extractAminoAcidSequence 정상 동작")
    void extractAminoAcidSequence_ShouldReturnValidSequence() {
        String alignmentJson = """
            {
              "alignment_index": {
                "P": [0, 2],
                "geneX": [2, 12]
              },
              "aligned_sequences": {
                "ref1": "ATGTTTGTTCTT",
                "variant1": "M-AT-GC---ET"
              }
            }
            """;

        String result = analysisService.extractAminoAcidSequence(alignmentJson, "geneX", "variant1", 1, 6);
        assertEquals("ATGC", result);
    }

    @Test
    @DisplayName("extractAminoAcidSequence - 유효하지 않은 범위 예외")
    void extractAminoAcidSequence_ShouldThrowException_WhenEmptyResult() {
        String alignmentJson = """
            {
              "alignment_index": {
                "geneX": [0, 4],
                "S": [4, 12]
                },
              "aligned_sequences": {
                "ref1": "ATGTTTGTTCTT",
                "variant1": "------------"
                }
            }
            """;

        assertThrows(InvalidSequenceRangeException.class, () ->
                analysisService.extractAminoAcidSequence(alignmentJson, "geneX", "variant1", 1, 3)
        );
    }

    @Test
    @DisplayName("extractSequence 정상 동작")
    void extractSequence_ShouldReturnCleanedSequence() {
        String alignmentJson = """
            {
              "alignment_index": {
                "geneX": [0, 4],
                "S": [4, 12]
                },
              "aligned_sequences": {
                "ref1": "ATGTTTGTTCTT",
                "variant1": "ATG---------"
                }
            }
            """;

        String result = analysisService.extractSequence("ref1", alignmentJson, "geneX");
        assertEquals("ATGT", result);
    }

    @Test
    @DisplayName("getAnalysisData - 존재하지 않을 때 예외 발생")
    void getAnalysisData_ShouldThrowException_WhenNotFound() {
        History history = new History();
        history.setId(99L);

        when(analysisRepository.findByHistoryId(99L)).thenReturn(null);

        assertThrows(AnalysisNotFoundException.class, () -> {
            analysisService.getAnalysisData(history);
        });
    }

    @Test
    @DisplayName("processLinearDesign - 정상 처리")
    void processLinearDesign_ShouldReturnJson() {
        LinearDesignDto dto = new LinearDesignDto();
        dto.setGene("geneX");
        dto.setVarientName("variant1");
        dto.setStart(1);
        dto.setEnd(4);
        dto.setHistoryName("hist1");

        History history = new History();
        Analysis analysis = new Analysis();
        analysis.setAlignment("""
            {
              "alignment_index": {
                "geneX": [0, 4],
                "S": [4, 12]
                },
              "aligned_sequences": {
                "ref1": "ATGTTTGTTCTT",
                "variant1": "ATG---------"
                }
            }
            """);

        BioPythonDto mockResponse = new BioPythonDto();
        mockResponse.setOutput("{\"structure\": \"...\"}");

        when(historyService.getHistory("hist1", 1L)).thenReturn(history);
        when(analysisRepository.findByHistoryId(any())).thenReturn(analysis);
        when(bioPythonService.executePythonScript(eq("3"), anyString())).thenReturn(mockResponse);

        String result = analysisService.processLinearDesign(dto, 1L);

        assertEquals("{\"structure\": \"...\"}", result);
        verify(analysisRepository).save(any());
    }

    @Test
    @DisplayName("processPdb - 정상 처리")
    void processPdb_ShouldReturnJson() {
        PdbDto dto = new PdbDto();
        dto.setGene("geneX");
        dto.setHistoryName("hist1");

        History history = new History();
        Analysis analysis = new Analysis();
        analysis.setReferenceId("ref1");
        analysis.setAlignment("""
            {
              "alignment_index": {"geneX": [0, 4]},
              "aligned_sequences": {"ref1": "A-T-G"}
            }
            """);

        BioPythonDto mockResponse = new BioPythonDto();
        mockResponse.setOutput("{\"pdb\": \"ATOM ...\"}");

        when(historyService.getHistory("hist1", 1L)).thenReturn(history);
        when(analysisRepository.findByHistoryId(any())).thenReturn(analysis);
        when(bioPythonService.executePythonScript(eq("4"), anyString())).thenReturn(mockResponse);

        String result = analysisService.processPdb(dto, 1L);

        assertEquals("{\"pdb\": \"ATOM ...\"}", result);
        verify(analysisRepository).save(any());
    }

    @Test
    @DisplayName("deleteAnalysisData - 정상 호출")
    void deleteAnalysisData_ShouldCallDeleteByHistoryId() {
        History history = new History();
        history.setId(123L);

        analysisService.deleteAnalysisData(history);

        verify(analysisRepository).deleteByHistoryId(123L);
    }
}
