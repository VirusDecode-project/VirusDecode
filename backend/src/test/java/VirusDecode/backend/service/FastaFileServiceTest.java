package VirusDecode.backend.service;

import VirusDecode.backend.dto.FastaFileDTO;
import VirusDecode.backend.dto.VarientDTO;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

class FastaFileServiceTest {

    private FastaFileService fastaFileService;

    @BeforeEach
    void setUp() {
        fastaFileService = new FastaFileService();
    }


    @Test
    void testSaveFastaContentWithNullSequencesAndFiles() throws IOException {
        // Arrange: sequences와 files가 null인 경우
        VarientDTO dto = new VarientDTO();
        dto.setSequences(null);
        dto.setFiles(null);

        // Act
        String result = fastaFileService.saveFastaContent(dto);

        // Assert: 빈 문자열을 기대
        assertEquals("", result);
    }

    @Test
    void testSaveFastaContentWithEmptySequencesAndFiles() throws IOException {
        // Arrange: sequences와 files가 빈 경우
        VarientDTO dto = new VarientDTO();
        dto.setSequences(new HashMap<>());
        dto.setFiles(new ArrayList<>());

        // Act
        String result = fastaFileService.saveFastaContent(dto);

        // Assert: 빈 문자열을 기대
        assertEquals("", result);
    }

    @Test
    void testSaveFastaContentWithSequencesContainingEmptySequenceData() throws IOException {
        // Arrange: 빈 시퀀스 데이터 포함
        VarientDTO dto = new VarientDTO();
        Map<String, String> sequences = new LinkedHashMap<>();
        sequences.put("sequence1", " ");
        dto.setSequences(sequences);

        // Act
        String result = fastaFileService.saveFastaContent(dto);

        // Assert: 빈 문자열을 기대 (공백 시퀀스 무시)
        assertEquals("", result);
    }

    @Test
    void testSaveFastaContentWithSequencesContainingValidData() throws IOException {
        // Arrange: 올바른 시퀀스 데이터
        VarientDTO dto = new VarientDTO();
        Map<String, String> sequences = new LinkedHashMap<>();
        sequences.put("sequence1", "ATCG");
        sequences.put("sequence 2", "GCTA");  // 시퀀스 이름에 공백 포함
        dto.setSequences(sequences);

        // Act
        String result = fastaFileService.saveFastaContent(dto);

        // Assert
        String expected = ">sequence1\nATCG\n>sequence2\nGCTA\n";  // sequenceName에서 공백 제거됨
        assertEquals(expected, result);
    }

    @Test
    void testSaveFastaContentWithFiles() throws IOException {
        // Arrange: 파일 데이터 포함
        VarientDTO dto = new VarientDTO();
        FastaFileDTO fileDTO = new FastaFileDTO();
        fileDTO.setName("file1");
        fileDTO.setContent(">FileSequence\nACTGACTG");

        List<FastaFileDTO> files = new ArrayList<>();
        files.add(fileDTO);
        dto.setFiles(files);

        // Act
        String result = fastaFileService.saveFastaContent(dto);

        // Assert
        String expected = ">FileSequence\nACTGACTG\n";
        assertEquals(expected, result);
    }

    @Test
    void testSaveFastaContentWithSequencesAndFiles() throws IOException {
        // Arrange: 시퀀스와 파일 데이터가 모두 있는 경우
        VarientDTO dto = new VarientDTO();

        Map<String, String> sequences = new LinkedHashMap<>();
        sequences.put("sequence1", "ATCG");
        dto.setSequences(sequences);

        FastaFileDTO fileDTO = new FastaFileDTO();
        fileDTO.setName("file1");
        fileDTO.setContent(">FileSequence\nACTGACTG");

        List<FastaFileDTO> files = new ArrayList<>();
        files.add(fileDTO);
        dto.setFiles(files);

        // Act
        String result = fastaFileService.saveFastaContent(dto);

        // Assert
        String expected = ">sequence1\nATCG\n>FileSequence\nACTGACTG\n";
        assertEquals(expected, result);
    }

}
