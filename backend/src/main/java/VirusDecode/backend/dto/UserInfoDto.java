package VirusDecode.backend.dto;

import lombok.Getter;

@Getter
public class UserInfoDto {
    private String loginId;
    private String userName;

    public UserInfoDto() {
    }

    public UserInfoDto(String userName) {
        this.userName = userName;
    }

    public UserInfoDto(String loginId, String userName) {
        this.loginId = loginId;
        this.userName = userName;
    }
}
