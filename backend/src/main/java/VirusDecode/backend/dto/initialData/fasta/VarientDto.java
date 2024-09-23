package VirusDecode.backend.dto.initialData.fasta;

import lombok.Getter;
import lombok.Setter;

import java.util.List;
import java.util.Map;

@Getter
@Setter
public class VarientDto {
    private String referenceSequenceId;
    private Map<String, String> sequences;
    private List<FastaFileDto> files;
}
