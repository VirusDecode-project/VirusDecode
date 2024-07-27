package VirusDecode.backend.dto;

import org.springframework.security.core.GrantedAuthority;
import org.springframework.security.oauth2.core.user.OAuth2User;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

public class CustomOAuth2User implements OAuth2User {

    private final OAuth2Response oAuth2Response;
    private final String role;

    public CustomOAuth2User(OAuth2Response oAuth2Response, String role){
        this.oAuth2Response=oAuth2Response;
        this.role = role;
    }

    @Override
    public Map<String, Object> getAttributes() {  // attribute 는 리소스 서버로부터 넘어온 모든 유저정보
//        return Map.of();
        return null;
    }

    @Override
    public Collection<? extends GrantedAuthority> getAuthorities() {
        Collection<GrantedAuthority> collection = new ArrayList<>();
        collection.add(new GrantedAuthority() {
            @Override
            public String getAuthority() {
                return role;
            }
        });
        return collection;
    }

    @Override
    public String getName() {
        return oAuth2Response.getName();
    }

    // id로 사용할 값을 만들어줌 ( response 에 있는 Name 을 사용하지 않는 것.. )
    public String getUsername(){
        return oAuth2Response.getProvider()+" "+oAuth2Response.getProviderId();
    }
}
