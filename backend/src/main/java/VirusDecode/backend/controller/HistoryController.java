package VirusDecode.backend.controller;

import VirusDecode.backend.dto.HistoryDTO;
import VirusDecode.backend.service.JsonFileService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.util.Arrays;
import java.util.List;
import java.util.Comparator;
import java.util.stream.Collectors;

import com.fasterxml.jackson.databind.ObjectMapper;

@RestController
@RequestMapping("/history")
public class HistoryController {
    private final JsonFileService jsonFileService;

    @Autowired
    public HistoryController(JsonFileService jsonFileService) {
        this.jsonFileService = jsonFileService;
    }
    private static final Path currentDir = Paths.get("").toAbsolutePath();

    @PostMapping("/create")
    public ResponseEntity<String> createHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        try {
            createHistory(historyName);
            return ResponseEntity.ok("History created successfully.");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to create history: " + e.getMessage());
        }
    }

    @PostMapping("/delete")
    public ResponseEntity<String> deleteHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        try {
            deleteHistory(historyName);
            return ResponseEntity.ok("History deleted successfully.");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to delete history: " + e.getMessage());
        }
    }

    @PostMapping("/rename")
    public ResponseEntity<String> renameHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        String newName = request.getNewName();
        try {
            renameHistory(historyName, newName);
            return ResponseEntity.ok("History renamed successfully.");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to rename history: " + e.getMessage());
        }
    }

    @PostMapping("/get")
    public ResponseEntity<String> getHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        try {
            getHistory(historyName);
            return jsonFileService.readJsonFile("alignment_data.json");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to get history: " + e.getMessage());
        }
    }

    @GetMapping("/list")
    public ResponseEntity<String> listHistory() {
        try {
            String jsonHistoryList = getHistoryList();
            return ResponseEntity.ok(jsonHistoryList);
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to list history: " + e.getMessage());
        }
    }

    @GetMapping("/deleteData")
    public ResponseEntity<String> deleteData() {
        try {
            getDeleteData();
            return ResponseEntity.ok("Data delete successfully.");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to delete data: " + e.getMessage());
        }
    }


    private void createHistory(String historyName) throws IOException {
        Path newDir = currentDir.resolve("history").resolve(historyName);
        Files.createDirectories(newDir);

        Path dataDir = currentDir.resolve("data");

        if (Files.exists(dataDir)) {
            try (DirectoryStream<Path> stream = Files.newDirectoryStream(dataDir)) {
                for (Path filePath : stream) {
                    if (Files.isRegularFile(filePath)) {
                        Files.copy(filePath, newDir.resolve(filePath.getFileName()), StandardCopyOption.REPLACE_EXISTING);
                    }
                }
            }
        }
    }

    private void deleteHistory(String historyName) throws IOException {
        Path newDir = currentDir.resolve("history").resolve(historyName);

        if (Files.exists(newDir)) {
            Files.walk(newDir)
                    .sorted(Comparator.reverseOrder())
                    .map(Path::toFile)
                    .forEach(File::delete);
        }
    }

    private void renameHistory(String historyName, String newHistoryName) throws IOException {
        Path oldDir = currentDir.resolve("history").resolve(historyName);
        Path newDir = currentDir.resolve("history").resolve(newHistoryName);

        if (Files.exists(oldDir)) {
            Files.move(oldDir, newDir, StandardCopyOption.REPLACE_EXISTING);
        }
    }

    private void getHistory(String historyName) throws IOException {
        Path sourceDir = currentDir.resolve("history").resolve(historyName);
        Path dataDir = currentDir.resolve("data");

        if (Files.exists(sourceDir)) {
            if (Files.exists(dataDir)) {
                try (DirectoryStream<Path> stream = Files.newDirectoryStream(dataDir)) {
                    for (Path filePath : stream) {
                        if (Files.isRegularFile(filePath)) {
                            Files.delete(filePath);
                        }
                    }
                }
            }
            Files.createDirectories(dataDir);

            try (DirectoryStream<Path> stream = Files.newDirectoryStream(sourceDir)) {
                for (Path filePath : stream) {
                    if (Files.isRegularFile(filePath)) {
                        Files.copy(filePath, dataDir.resolve(filePath.getFileName()), StandardCopyOption.REPLACE_EXISTING);
                    }
                }
            }
        }
    }

    private String getHistoryList() throws IOException {
        Path historyDir = currentDir.resolve("history");

        if (Files.notExists(historyDir)) {
            Files.createDirectories(historyDir);
        }

        List<String> historyList = Arrays.stream(new File(historyDir.toString()).listFiles())
                .sorted(Comparator.comparingLong(File::lastModified))  // 마지막 수정 시간을 기준으로 정렬
                .map(File::getName)  // 파일 이름만 가져오기
                .collect(Collectors.toList());

        // 리스트를 JSON 문자열로 변환
        ObjectMapper objectMapper = new ObjectMapper();
        return objectMapper.writerWithDefaultPrettyPrinter().writeValueAsString(historyList);
    }

    private void getDeleteData() throws IOException {
        Path newDir = currentDir.resolve("data");


        if (Files.exists(newDir)) {
            Files.walk(newDir)
                    .sorted(Comparator.reverseOrder())
                    .map(Path::toFile)
                    .forEach(File::delete);
        }

    }
}