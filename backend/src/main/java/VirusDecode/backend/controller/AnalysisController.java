package VirusDecode.backend.controller;

import VirusDecode.backend.dto.LinearDesignDTO;
import VirusDecode.backend.dto.PdbDTO;
import VirusDecode.backend.entity.JsonDataEntity;
import VirusDecode.backend.service.JsonDataService;
import VirusDecode.backend.service.PythonScriptExecutor;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.io.IOException;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

@RestController
@RequestMapping("/analysis")
public class AnalysisController {
    private final PythonScriptExecutor pythonScriptExecutor;
    private final JsonDataService jsonDataService;
    private static final Path currentDir = Paths.get("").toAbsolutePath();
    private static final Path HISTORY_DIR = currentDir.resolve("history");

    @Autowired
    public AnalysisController(PythonScriptExecutor pythonScriptExecutor, JsonDataService jsonDataService) {
        this.pythonScriptExecutor = pythonScriptExecutor;
        this.jsonDataService = jsonDataService;
    }


    // /analysis/mrnadesign 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/linearDesign")
    public ResponseEntity<String> getLinearDesign(@RequestBody LinearDesignDTO request) {
        // mRNA 디자인 요청에서 필요한 데이터를 추출
        String region = request.getRegion();
        String varientName = request.getVarientName();
        String start = String.valueOf(request.getStart());
        String end = String.valueOf(request.getEnd());
        String historyName = request.getHistoryName();

        String metadataJson = jsonDataService.getJsonData("metadata");
        String alignmentJson = jsonDataService.getJsonData("alignment");

        ResponseEntity<String> scriptResponse = pythonScriptExecutor.executePythonScript("3", metadataJson, alignmentJson, region, varientName, start, end);
        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            jsonDataService.saveJsonData("linearDesign", scriptResponse.getBody());
            try {
                performSaveHistory(historyName);
            } catch (Exception e) {
                return ResponseEntity.status(500).body("Unexpected error during save history");
            }
        }
        return scriptResponse;
    }

    // /analysis/pdb 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/pdb")
    public ResponseEntity<String> getPdb(@RequestBody PdbDTO request) {
        String gene = request.getGene();
        String historyName = request.getHistoryName();
        String metadataJson = jsonDataService.getJsonData("metadata");
        String alignmentJson = jsonDataService.getJsonData("alignment");

        ResponseEntity<String> scriptResponse = pythonScriptExecutor.executePythonScript("4", metadataJson, alignmentJson, gene);
        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            jsonDataService.saveJsonData("pdb", scriptResponse.getBody());
            try {
                performSaveHistory(historyName);
            } catch (Exception e) {
                return ResponseEntity.status(500).body("Unexpected error during save history");
            }
        }
        return scriptResponse;
    }

    // /analysis/re-alignment 엔드포인트에 대한 POST 요청 처리
    @GetMapping("/re-alignment")
    public ResponseEntity<String> reGetAlignment() {
        return ResponseEntity.ok(jsonDataService.getJsonData("alignment"));
    }

    // /analysis/re-mrnadesign 엔드포인트에 대한 POST 요청 처리
    @GetMapping("/re-linearDesign")
    public ResponseEntity<String> reGetLinearDesign() {
        return ResponseEntity.ok(jsonDataService.getJsonData("linearDesign"));
    }

    // /analysis/render3d 엔드포인트에 대한 GET 요청 처리
    @GetMapping("/re-pdb")
    public ResponseEntity<String> reGetPdb() {
        return ResponseEntity.ok(jsonDataService.getJsonData("pdb"));
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
