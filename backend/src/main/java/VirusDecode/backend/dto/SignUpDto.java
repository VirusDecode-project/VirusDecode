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

    public SignUpDto() {
    }

    public SignUpDto(String firstName, String lastName, String loginId, String password) {
        this.firstName = firstName;
        this.lastName = lastName;
        this.loginId = loginId;
        this.password = password;
    }
}
