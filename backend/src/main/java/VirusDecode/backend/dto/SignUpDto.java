package VirusDecode.backend.dto;

import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class SignUpDto {
    private String firstName;
    private String lastName;
    private String loginId;
    private String password;
}
