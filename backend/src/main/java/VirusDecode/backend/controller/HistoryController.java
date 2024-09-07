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
    private static final Path HISTORY_DIR = currentDir.resolve("history");
    private static final Path DATA_DIR = currentDir.resolve("data");

    @PostMapping("/create")
    public ResponseEntity<String> createHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        try {
                performCreateHistory(historyName);
            return ResponseEntity.ok("History created successfully.");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to create history: " + e.getMessage());
        }
    }

    @DeleteMapping("/delete")
    public ResponseEntity<String> deleteHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        try {
            performDeleteHistory(historyName);
            return ResponseEntity.ok("History deleted successfully.");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to delete history: " + e.getMessage());
        }
    }

    @PutMapping("/rename")
    public ResponseEntity<String> renameHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        String newName = request.getNewName();
        try {
            performRenameHistory(historyName, newName);
            return ResponseEntity.ok("History renamed successfully.");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to rename history: " + e.getMessage());
        }
    }

    @PostMapping("/get")
    public ResponseEntity<String> getHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        try {
            performGetHistory(historyName);
            return jsonFileService.readJsonFile("alignment_data.json");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to get history: " + e.getMessage());
        }
    }

    @GetMapping("/list")
    public ResponseEntity<String> listHistory() {
        try {
            String jsonHistoryList = performHistoryList();
            return ResponseEntity.ok(jsonHistoryList);
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to list history: " + e.getMessage());
        }
    }

    @GetMapping("/deleteData")
    public ResponseEntity<String> deleteData() {
        try {
            performDeleteData();
            return ResponseEntity.ok("Data delete successfully.");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to delete data: " + e.getMessage());
        }
    }


    private void performCreateHistory(String historyName) throws IOException {
        Path newDir = HISTORY_DIR.resolve(historyName);
        Files.createDirectories(newDir);

        if (Files.exists(DATA_DIR)) {
            try (DirectoryStream<Path> stream = Files.newDirectoryStream(DATA_DIR)) {
                for (Path filePath : stream) {
                    if (Files.isRegularFile(filePath)) {
                        Files.copy(filePath, newDir.resolve(filePath.getFileName()), StandardCopyOption.REPLACE_EXISTING);
                    }
                }
            }
        }
    }

    private void performDeleteHistory(String historyName) throws IOException {
        Path newDir = HISTORY_DIR.resolve(historyName);

        if (Files.exists(newDir)) {
            Files.walk(newDir)
                    .sorted(Comparator.reverseOrder())
                    .map(Path::toFile)
                    .forEach(File::delete);
        }
    }

    private void performRenameHistory(String historyName, String newHistoryName) throws IOException {
        Path oldDir = HISTORY_DIR.resolve(historyName);
        Path newDir = HISTORY_DIR.resolve(newHistoryName);

        if (Files.exists(oldDir)) {
            Files.move(oldDir, newDir, StandardCopyOption.REPLACE_EXISTING);
        }
    }

    private void performGetHistory(String historyName) throws IOException {
        Path sourceDir = HISTORY_DIR.resolve(historyName);

        if (Files.exists(sourceDir)) {
            if (Files.exists(DATA_DIR)) {
                try (DirectoryStream<Path> stream = Files.newDirectoryStream(DATA_DIR)) {
                    for (Path filePath : stream) {
                        if (Files.isRegularFile(filePath)) {
                            Files.delete(filePath);
                        }
                    }
                }
            }
            Files.createDirectories(DATA_DIR);

            try (DirectoryStream<Path> stream = Files.newDirectoryStream(sourceDir)) {
                for (Path filePath : stream) {
                    if (Files.isRegularFile(filePath)) {
                        Files.copy(filePath, DATA_DIR.resolve(filePath.getFileName()), StandardCopyOption.REPLACE_EXISTING);
                    }
                }
            }
        }
    }

    private String performHistoryList() throws IOException {
        if (Files.notExists(HISTORY_DIR)) {
            Files.createDirectories(HISTORY_DIR);
        }

        List<String> historyList = Arrays.stream(new File(HISTORY_DIR.toString()).listFiles())
                .sorted(Comparator.comparingLong(File::lastModified))  // 마지막 수정 시간을 기준으로 정렬
                .map(File::getName)  // 파일 이름만 가져오기
                .collect(Collectors.toList());

        // 리스트를 JSON 문자열로 변환
        ObjectMapper objectMapper = new ObjectMapper();
        return objectMapper.writerWithDefaultPrettyPrinter().writeValueAsString(historyList);
    }

    private void performDeleteData() throws IOException {
        if (Files.exists(DATA_DIR)) {
            Files.walk(DATA_DIR)
                    .sorted(Comparator.reverseOrder())
                    .map(Path::toFile)
                    .forEach(File::delete);
        }

    }
}