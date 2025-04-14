package virusdecode.backend.analysis.dto;

import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class AlignmentDto {
    String alignment;
    String historyName;

    public AlignmentDto(String alignment, String historyName) {
        this.alignment = alignment;
        this.historyName = historyName;
    }
}
