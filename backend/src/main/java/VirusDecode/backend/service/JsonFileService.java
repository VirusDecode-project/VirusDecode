package VirusDecode.backend.service;

import org.springframework.core.io.ClassPathResource;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;

import java.io.IOException;
import java.nio.file.Files;

@Service
public class JsonFileService {

    // 주어진 파일 경로의 JSON 파일을 읽어들여 그 내용을 반환하는 메서드
    public ResponseEntity<String> readJsonFile(String filePath) {
        try {
            // ClassPathResource를 사용하여 클래스패스 내의 JSON 파일을 로드
            ClassPathResource resource = new ClassPathResource(filePath);
            String jsonContent = Files.readString(resource.getFile().toPath());  // 파일 내용을 문자열로 읽기

            // JSON 데이터를 HTTP 응답으로 반환
            return ResponseEntity.ok(jsonContent);
        } catch (IOException e) {
            e.printStackTrace();
            // 파일 읽기 중 오류가 발생한 경우, 상태 코드 500과 오류 메시지를 반환
            return ResponseEntity.status(500).body("Error reading JSON file");
        }
    }
}
