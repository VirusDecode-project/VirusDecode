package VirusDecode.backend.config;

import VirusDecode.backend.service.CustomOAuth2UserService;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.security.config.Customizer;
import org.springframework.security.config.annotation.web.builders.HttpSecurity;
import org.springframework.security.config.annotation.web.configuration.EnableWebSecurity;
import org.springframework.security.config.annotation.web.configurers.AbstractHttpConfigurer;
import org.springframework.security.web.SecurityFilterChain;


@EnableWebSecurity
@Configuration
public class SecurityConfig{

    private final CustomOAuth2UserService customOAuth2UserService;

    public SecurityConfig(CustomOAuth2UserService customOAuth2UserService){
        this.customOAuth2UserService = customOAuth2UserService;
    }

    @Bean
    public SecurityFilterChain filterChain(HttpSecurity http) throws Exception {
        http
                .csrf(AbstractHttpConfigurer::disable);
        http
                .formLogin(AbstractHttpConfigurer::disable);
        http
                .httpBasic(AbstractHttpConfigurer::disable);
        http
                .oauth2Login((oauth2)->oauth2
                        .defaultSuccessUrl("http://localhost:3000/login-success", true)  // 로그인 성공시 ( response 로 전달 )
                        .failureUrl("http://localhost:3000/login-failure")  // 로그인 실패시
                        .userInfoEndpoint((userInfoEndpointConfig ->
                                userInfoEndpointConfig.userService(customOAuth2UserService))));
        http
                .authorizeHttpRequests((auth)->auth
                        .requestMatchers("/", "/oauth2/**","/login/**").permitAll()
                        .anyRequest().authenticated());   // request 에 대해 모두 인증(로그인)이 필요함


        return http.build();
    }
}
