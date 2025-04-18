package virusdecode.backend.bioinput.dto;

import lombok.Getter;
import lombok.Setter;

// Request Body를 받을 DTO 클래스 정의
@Getter
@Setter
public class ReferenceDto {
    private String sequenceId;

    public ReferenceDto() {
    }

    public ReferenceDto(String sequenceId) {
        this.sequenceId = sequenceId;
    }
}