package VirusDecode.backend.controller;

import VirusDecode.backend.dto.HistoryDto;
import VirusDecode.backend.entity.JsonDataEntity;
import VirusDecode.backend.service.JsonDataService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;

import com.fasterxml.jackson.databind.ObjectMapper;

@RestController
@RequestMapping("/history")
public class HistoryController {
    private final JsonDataService jsonDataService;
    @Autowired
    public HistoryController(JsonDataService jsonDataService) {
        this.jsonDataService = jsonDataService;
    }

    private static final Path currentDir = Paths.get("").toAbsolutePath();
    private static final Path HISTORY_DIR = currentDir.resolve("history");

    @PostMapping("/create")
    public ResponseEntity<String> createHistory(@RequestBody HistoryDto request) {
        String historyName = request.getHistoryName();
        try {
            String createdHistoryName = performCreateHistory(historyName);
            return ResponseEntity.ok(createdHistoryName);  // Return the created directory name
        } catch (FileAlreadyExistsException e) {
            return ResponseEntity.status(409).body("History already exists: " + historyName);
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to create history: " + e.getMessage());
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Unexpected error: " + e.getMessage());
        }
    }

    @DeleteMapping("/delete")
    public ResponseEntity<String> deleteHistory(@RequestBody HistoryDto request) {
        String historyName = request.getHistoryName();
        try {
            performDeleteHistory(historyName);
            return ResponseEntity.ok("History deleted successfully.");
        } catch (NoSuchFileException e) {
            return ResponseEntity.status(404).body("History not found: " + historyName);
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to delete history: " + e.getMessage());
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Unexpected error: " + e.getMessage());
        }
    }

    @PutMapping("/rename")
    public ResponseEntity<String> renameHistory(@RequestBody HistoryDto request) {
        String historyName = request.getHistoryName();
        String newName = request.getNewName();
        try {
            performRenameHistory(historyName, newName);
            return ResponseEntity.ok("History renamed successfully.");
        } catch (NoSuchFileException e) {
            return ResponseEntity.status(404).body("History not found: " + historyName);
        } catch (FileAlreadyExistsException e) {
            return ResponseEntity.status(409).body("History with new name already exists: " + newName);
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to rename history: " + e.getMessage());
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Unexpected error: " + e.getMessage());
        }
    }

    @PostMapping("/get")
    public ResponseEntity<String> getHistory(@RequestBody HistoryDto request) {
        String historyName = request.getHistoryName();
        try {
            Path sourceDir = HISTORY_DIR.resolve(historyName);
            if (Files.notExists(sourceDir)) {
                return ResponseEntity.status(404).body("History not found: " + historyName);
            }
            performGetHistory(sourceDir);
            return ResponseEntity.ok(jsonDataService.getJsonData("alignment"));
        } catch (NoSuchFileException e) {
            return ResponseEntity.status(404).body("History file not found: " + historyName);
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to get history: " + e.getMessage());
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Unexpected error: " + e.getMessage());
        }
    }

    @GetMapping("/list")
    public ResponseEntity<String> listHistory() {
        try {
            String jsonHistoryList = performHistoryList();
            return ResponseEntity.ok(jsonHistoryList);
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to list history: " + e.getMessage());
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Unexpected error: " + e.getMessage());
        }
    }

    @GetMapping("/deleteData")
    public ResponseEntity<String> deleteData() {
        try {
            performDeleteData();
            return ResponseEntity.ok("Data delete successfully.");
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to delete data: " + e.getMessage());
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Unexpected error: " + e.getMessage());
        }
    }

    @PostMapping("/save")
    public ResponseEntity<String> saveHistory(@RequestBody HistoryDto request) {
        String historyName = request.getHistoryName();
        try {
            performSaveHistory(historyName);
            return ResponseEntity.ok("History save successfully.");  // Return the created directory name
        } catch (FileAlreadyExistsException e) {
            return ResponseEntity.status(409).body("History already exists: " + historyName);
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Failed to create history: " + e.getMessage());
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Unexpected error: " + e.getMessage());
        }
    }

    @GetMapping("/checkFiles")
    public ResponseEntity<Map<String, Boolean>> checkFiles() {
        Map<String, Boolean> fileStatus = new HashMap<>();

        // Repository에서 JSON 데이터 확인
        boolean pdbExists = jsonDataService.findByName("pdb").isPresent();
        boolean linearDesignExists = jsonDataService.findByName("linearDesign").isPresent();

        // 상태를 Map에 저장
        fileStatus.put("pdb.json", pdbExists);
        fileStatus.put("linearDesign.json", linearDesignExists);

        return ResponseEntity.ok(fileStatus);
    }

    private String performCreateHistory(String historyName) throws IOException {
        Path newDir = HISTORY_DIR.resolve(historyName);

        // 동일한 이름이 존재하면 숫자를 붙임
        int counter = 1;
        while (Files.exists(newDir)) {
            String newHistoryName = historyName + "_" + counter;
            newDir = HISTORY_DIR.resolve(newHistoryName);
            counter++;
        }

        Files.createDirectories(newDir);

        // repository에서 JSON 데이터를 가져와 newDir에 파일로 저장
        Iterable<JsonDataEntity> jsonDataEntities = jsonDataService.getAllJsonData();  // 모든 JSON 데이터 가져오기
        for (JsonDataEntity jsonDataEntity : jsonDataEntities) {
            // 각 JSON 데이터를 파일로 저장
            Path jsonFilePath = newDir.resolve(jsonDataEntity.getName() + ".json");  // ID를 파일명으로 사용
            Files.write(jsonFilePath, jsonDataEntity.getJsonData().getBytes());  // JSON 데이터를 파일로 저장
        }




        return newDir.getFileName().toString();  // Return the directory name
    }

    private void performDeleteHistory(String historyName) throws IOException {
        Path newDir = HISTORY_DIR.resolve(historyName);

        if (Files.exists(newDir)) {
            Files.walk(newDir)
                    .sorted(Comparator.reverseOrder())
                    .map(Path::toFile)
                    .forEach(File::delete);
        }else {
            throw new NoSuchFileException("History directory not found: " + historyName);
        }
    }

    private void performRenameHistory(String historyName, String newHistoryName) throws IOException {
        Path oldDir = HISTORY_DIR.resolve(historyName);
        Path newDir = HISTORY_DIR.resolve(newHistoryName);

        if (Files.exists(oldDir)) {
            Files.move(oldDir, newDir, StandardCopyOption.REPLACE_EXISTING);
        }else {
            throw new NoSuchFileException("History directory not found: " + historyName);
        }
    }

    private void performGetHistory(Path sourceDir) throws IOException {
        // sourceDir이 존재하는지 확인
        if (Files.exists(sourceDir)) {

            // 먼저 repository의 기존 데이터를 삭제
            jsonDataService.deleteAllJsonData();  // 기존 데이터 모두 삭제

            // sourceDir에서 JSON 파일들을 읽어와 repository에 저장
            try (DirectoryStream<Path> stream = Files.newDirectoryStream(sourceDir)) {
                for (Path filePath : stream) {
                    if (Files.isRegularFile(filePath) && filePath.toString().endsWith(".json")) {  // JSON 파일만 처리
                        // 파일 내용을 읽기
                        String jsonContent = new String(Files.readAllBytes(filePath));

                        // 파일 이름을 ID로 사용 (예: fileName.json -> fileName)
                        String fileName = filePath.getFileName().toString();
                        String name = fileName.substring(0, fileName.lastIndexOf('.'));

                        // repository에 저장
                        jsonDataService.saveJsonData(name, jsonContent);
                    }
                }
            }
        } else {
            throw new NoSuchFileException("Source directory not found: " + sourceDir.getFileName().toString());
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
        jsonDataService.deleteAllJsonData();
    }

    private void performSaveHistory(String historyName) throws IOException {
        // 새로운 history 디렉토리 경로 설정
        Path newDir = HISTORY_DIR.resolve(historyName);

        // newDir 디렉토리가 없으면 생성
        if (!Files.exists(newDir)) {
            Files.createDirectories(newDir);
        }

        // repository에서 JSON 데이터를 가져와 newDir에 파일로 저장
        Iterable<JsonDataEntity> jsonDataEntities = jsonDataService.getAllJsonData();  // 모든 JSON 데이터 가져오기
        for (JsonDataEntity jsonDataEntity : jsonDataEntities) {
            // 각 JSON 데이터를 파일로 저장
            Path jsonFilePath = newDir.resolve(jsonDataEntity.getName() + ".json");  // ID를 파일명으로 사용
            Files.write(jsonFilePath, jsonDataEntity.getJsonData().getBytes());  // JSON 데이터를 파일로 저장
        }
    }
}