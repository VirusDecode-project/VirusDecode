package virusdecode.backend.analysis.dto;

import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class LinearDesignDto {
    private String gene;
    private String varientName;
    private int start;
    private int end;
    private String historyName;
}



