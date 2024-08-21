package VirusDecode.backend.controller;
import org.springframework.core.io.ClassPathResource;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

@RestController
@RequestMapping("/analysis")
public class analysisController {

    @PostMapping("/re-alignment")
    public ResponseEntity<String> sendJsonFile() {
        try {
            // ClassPathResource를 사용하여 클래스패스 내에서 JSON 파일을 읽음
            ClassPathResource resource = new ClassPathResource("bioinformatics/test_without_bio/analyze.json");
            Path filePath = resource.getFile().toPath();
            String jsonContent = Files.readString(filePath);

            // JSON 데이터를 응답으로 반환
            return ResponseEntity.ok(jsonContent);

        } catch (IOException e) {
            e.printStackTrace();
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Error reading JSON file");
        }
    }

}
