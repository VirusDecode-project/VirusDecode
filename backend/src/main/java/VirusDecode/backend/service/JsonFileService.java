package VirusDecode.backend.service;

import org.springframework.core.io.ClassPathResource;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;

@Service
public class JsonFileService {

    private static final Logger logger = LoggerFactory.getLogger(JsonFileService.class);

    // 주어진 파일 경로의 JSON 파일을 읽어들여 그 내용을 반환하는 메서드
    public static ResponseEntity<String> readJsonFile(String jsonFilePath) {
        try {
            // ClassPathResource를 사용하여 클래스패스 내의 JSON 파일을 로드
            ClassPathResource resource = new ClassPathResource(jsonFilePath);
            if (!resource.exists()) {
                logger.error("결과 JSON 파일을 찾을 수 없음: {}", jsonFilePath);
                return ResponseEntity.status(404).body("Result JSON file not found");
            }

            String jsonContent = Files.readString(resource.getFile().toPath());  // 파일 내용을 문자열로 읽기

            logger.info("Read JSON content: {}", jsonContent);

            // JSON 데이터를 HTTP 응답으로 반환
            return ResponseEntity.ok(jsonContent);
        }catch (IOException e) {
            logger.error("IO 오류가 발생했습니다: {}", e.getMessage());
            return ResponseEntity.status(500).body("Error reading JSON file");
        }
    }
}
