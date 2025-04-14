package virusdecode.backend.bioinput.dto;

import virusdecode.backend.bioinput.entity.FastaFile;
import lombok.Getter;
import lombok.Setter;

import java.util.List;
import java.util.Map;

@Getter
@Setter
public class VarientSequenceDto {
    private String referenceSequenceId;
    private Map<String, String> sequences;
    private List<FastaFile> files;
    private String historyName;
    private String referenceId;
}
