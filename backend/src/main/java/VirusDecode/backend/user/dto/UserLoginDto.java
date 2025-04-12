package VirusDecode.backend.user.dto;

import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
public class UserLoginDto {
    private String loginId;
    private String password;

    public UserLoginDto() {
    }

    public UserLoginDto(String loginId, String password) {
        this.loginId = loginId;
        this.password = password;
    }
}
