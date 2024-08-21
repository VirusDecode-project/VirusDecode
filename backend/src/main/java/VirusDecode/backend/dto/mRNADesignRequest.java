package VirusDecode.backend.dto;

import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class mRNADesignRequest {
    private String region;
    private String varientName;
    private int start;
    private int end;
}



