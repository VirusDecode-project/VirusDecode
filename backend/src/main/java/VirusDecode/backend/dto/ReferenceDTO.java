package VirusDecode.backend.dto;

import lombok.Getter;
import lombok.Setter;

// Request Body를 받을 DTO 클래스 정의
@Getter
@Setter
public class ReferenceDTO {
    private String sequenceId;
}