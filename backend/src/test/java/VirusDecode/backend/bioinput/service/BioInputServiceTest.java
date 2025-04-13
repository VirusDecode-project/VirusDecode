package VirusDecode.backend.bioinput.service;

import VirusDecode.backend.analysis.dto.AlignmentDto;
import VirusDecode.backend.analysis.entity.Analysis;
import VirusDecode.backend.analysis.service.AnalysisService;
import VirusDecode.backend.bioinput.dto.ReferenceDto;
import VirusDecode.backend.bioinput.dto.VarientSequenceDto;
import VirusDecode.backend.bioinput.entity.FastaFile;
import VirusDecode.backend.bioinput.entity.MetaData;
import VirusDecode.backend.bioinput.exception.FastaFileSaveFailException;
import VirusDecode.backend.bioinput.repository.BioInputRepository;
import VirusDecode.backend.common.biopython.BioPythonDto;
import VirusDecode.backend.common.biopython.BioPythonService;
import VirusDecode.backend.history.entity.History;
import VirusDecode.backend.history.service.HistoryService;
import VirusDecode.backend.user.entity.User;
import VirusDecode.backend.user.service.UserService;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.*;
import org.mockito.junit.jupiter.MockitoExtension;

import java.util.*;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

@ExtendWith(MockitoExtension.class)
class BioInputServiceTest {

    @InjectMocks
    private BioInputService bioInputService;

    @Mock
    private BioPythonService bioPythonService;

    @Mock
    private UserService userService;

    @Mock
    private HistoryService historyService;

    @Mock
    private BioInputRepository bioInputRepository;

    @Mock
    private AnalysisService analysisService;

    @BeforeEach
    void setUp() {
//        MockitoAnnotations.openMocks(this);
    }

    @Test
    @DisplayName("getMetadata - 기존 메타데이터가 있으면 반환")
    void getMetadata_ShouldReturnExistingMetaData() {
        ReferenceDto dto = new ReferenceDto("REF123");
        MetaData mockMeta = new MetaData("json-data", "REF123");

        when(bioInputRepository.findById("REF123")).thenReturn(Optional.of(mockMeta));

        MetaData result = bioInputService.getMetadata(dto);

        assertEquals(mockMeta, result);
    }

    @Test
    @DisplayName("getMetadata - 없으면 바이오파이썬 실행 후 저장")
    void getMetadata_ShouldCallPythonAndSave_WhenNotExists() {
        ReferenceDto dto = new ReferenceDto("REF123");

        when(bioInputRepository.findById("REF123")).thenReturn(Optional.empty());
        when(bioPythonService.executePythonScript("1", "REF123"))
                .thenReturn(new BioPythonDto("{metadata: true}"));

        MetaData result = bioInputService.getMetadata(dto);

        assertEquals("{metadata: true}", result.getMetadata());
        verify(bioInputRepository).save(any(MetaData.class));
    }

    @Test
    @DisplayName("processAlignment - 정상 처리")
    void processAlignment_ShouldReturnAlignmentDto() {
        VarientSequenceDto dto = new VarientSequenceDto();
        dto.setHistoryName("myHistory");
        dto.setReferenceId("REF123");

        Map<String, String> seqs = new HashMap<>();
        seqs.put("var1", "ATCG");
        dto.setSequences(seqs);

        User mockUser = new User();
        History mockHistory = new History(mockUser, "validatedHistory");

        when(userService.getUserById(1L)).thenReturn(mockUser);
        when(historyService.validateHistoryName("myHistory", 1L)).thenReturn("validatedHistory");
        when(bioPythonService.executePythonScript(eq("2"), eq("REF123"), anyString()))
                .thenReturn(new BioPythonDto("{alignment: true}"));

        AlignmentDto result = bioInputService.processAlignment(dto, 1L);

        assertEquals("{alignment: true}", result.getAlignment());
        assertEquals("validatedHistory", result.getHistoryName());
        verify(analysisService).saveAnalysisData(any(Analysis.class));
    }

    @Test
    @DisplayName("makeFastaContent - 시퀀스 및 파일 포함")
    void makeFastaContent_ShouldReturnFastaFormattedText() {
        VarientSequenceDto dto = new VarientSequenceDto();

        Map<String, String> seqMap = new LinkedHashMap<>();
        seqMap.put("var 1", "ATCG");
        seqMap.put("var 2", "GGGG");
        dto.setSequences(seqMap);

        FastaFile file = new FastaFile();
        file.setContent(">file_seq\nAACCGG");
        dto.setFiles(List.of(file));

        String result = bioInputService.makeFastaContent(dto);

        String expected = """
                >var1
                ATCG
                >var2
                GGGG
                >file_seq
                AACCGG
                """.trim();

        assertEquals(expected.replaceAll("\\s+", ""), result.replaceAll("\\s+", ""));
    }

    @Test
    @DisplayName("makeFastaContent - 예외 발생 시 FastaFileSaveFailException")
    void makeFastaContent_ShouldThrow_WhenUnexpectedException() {
        VarientSequenceDto dto = mock(VarientSequenceDto.class);
        when(dto.getSequences()).thenThrow(new RuntimeException("Unexpected"));

        assertThrows(FastaFileSaveFailException.class, () -> {
            bioInputService.makeFastaContent(dto);
        });
    }
}
