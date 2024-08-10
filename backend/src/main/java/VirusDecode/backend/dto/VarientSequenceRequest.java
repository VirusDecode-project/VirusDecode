package VirusDecode.backend.dto;

import lombok.Getter;
import lombok.Setter;

import java.util.List;
import java.util.Map;

@Getter
@Setter
public class VarientSequenceRequest {
    private String referenceSequenceId;
    private Map<String, String> sequences;
    private List<FastaFileDTO> files;
}
