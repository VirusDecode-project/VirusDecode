package virusdecode.backend.common.biopython;

import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class BioPythonDto {
    private int exitCode;
    private String output;
    private String errorOutput;

    public BioPythonDto() {
    }

    public BioPythonDto(String output){
        this.output = output;
    }
    public BioPythonDto(int exitCode, String output, String errorOutput) {
        this.exitCode = exitCode;
        this.output = output;
        this.errorOutput = errorOutput;
    }

    public boolean isSuccess() {
        return exitCode == 0;
    }
}
