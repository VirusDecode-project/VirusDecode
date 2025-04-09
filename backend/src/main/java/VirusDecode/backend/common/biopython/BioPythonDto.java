package VirusDecode.backend.common.biopython;

import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class BioPythonDto {
    private final int exitCode;
    private final String output;
    private final String errorOutput;

    public BioPythonDto(int exitCode, String output, String errorOutput) {
        this.exitCode = exitCode;
        this.output = output;
        this.errorOutput = errorOutput;
    }


    public boolean isSuccess() {
        return exitCode == 0;
    }
}
