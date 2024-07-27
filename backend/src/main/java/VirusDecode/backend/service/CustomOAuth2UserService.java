package VirusDecode.backend.service;

import VirusDecode.backend.dto.CustomOAuth2User;
import VirusDecode.backend.dto.GoogleResponse;
import VirusDecode.backend.dto.OAuth2Response;
import VirusDecode.backend.dto.Role;
import VirusDecode.backend.entity.UserEntity;
import VirusDecode.backend.repository.UserRepository;
import org.springframework.security.oauth2.client.userinfo.DefaultOAuth2UserService;
import org.springframework.security.oauth2.client.userinfo.OAuth2UserRequest;
import org.springframework.security.oauth2.core.OAuth2AuthenticationException;
import org.springframework.security.oauth2.core.user.OAuth2User;
import org.springframework.stereotype.Service;

@Service
public class CustomOAuth2UserService extends DefaultOAuth2UserService {

    private final UserRepository userRepository;

    public CustomOAuth2UserService(UserRepository userRepository){
        this.userRepository = userRepository;
    }

    @Override
    public OAuth2User loadUser(OAuth2UserRequest userRequest) throws OAuth2AuthenticationException{
        // OAuth2User
        OAuth2User oAuth2User = super.loadUser(userRequest);
        System.out.println(oAuth2User.getAttributes());


        String registrationId = userRequest.getClientRegistration().getRegistrationId();

        // OAuth2Response : OAuth2User 의 각 특성들을 적절한 형태로 ( google, naver,   ... ) 분리하고 이를 가져오는 메서드를 구현한 객체
        // ( OAuth2User 를 OAuth2Response 객체에 담는다... )
        OAuth2Response oAuth2Response = null;
        if(registrationId.equals("google")){
            oAuth2Response = new GoogleResponse(oAuth2User.getAttributes());

        }else{
            return null;
        }


        // Repository 유저 정보 저장 & 확인
        String username = oAuth2Response.getProvider()+" "+oAuth2Response.getProviderId();
        UserEntity existData = userRepository.findByUsername(username);
        String role = "ROLE_USER";
        if (existData == null) {

            UserEntity userEntity = new UserEntity();
            userEntity.setUsername(username);
            userEntity.setEmail(oAuth2Response.getEmail());
            userEntity.setRole(role);

            userRepository.save(userEntity);
        }
        else {

            existData.setUsername(username);
            existData.setEmail(oAuth2Response.getEmail());

            role = existData.getRole();

            userRepository.save(existData);
        }

        return new CustomOAuth2User(oAuth2Response, role);
    }



}
