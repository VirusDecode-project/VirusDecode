package VirusDecode.backend.service;

import VirusDecode.backend.dto.initialData.ReferenceDto;
import VirusDecode.backend.dto.initialData.fasta.VarientDto;
import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.entity.User;
import com.google.gson.Gson;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

@Service
public class InitialDataService {
    private final PythonScriptService pythonScriptService;
    private final FastaFileService fastaFileService;
    private final JsonDataService jsonDataService;
    private final UserService userService;
    private final HistoryService historyService;

    @Autowired
    public InitialDataService(PythonScriptService pythonScriptService, FastaFileService fastaFileService, JsonDataService jsonDataService, UserService userService, HistoryService historyService) {
        this.pythonScriptService = pythonScriptService;
        this.fastaFileService = fastaFileService;
        this.jsonDataService = jsonDataService;
        this.userService = userService;
        this.historyService = historyService;
    }

    public ResponseEntity<String> processMetadata(ReferenceDto request) {
        String sequenceId = request.getSequenceId();
        return pythonScriptService.executePythonScript("1", sequenceId);
    }

    public ResponseEntity<String> processAlignment(VarientDto request, Long userId) {
        String historyName = request.getHistoryName();
        String referenceId = request.getReferenceId();

        try {
            String fastaContent = fastaFileService.saveFastaContent(request);
            ResponseEntity<String> scriptResponse = pythonScriptService.executePythonScript("2", referenceId, fastaContent);
            if (!scriptResponse.getStatusCode().is2xxSuccessful()) {
                return scriptResponse;
            }

            String alignmentJson = scriptResponse.getBody();

            Optional<User> userOptional = userService.getUserById(userId);
            if (!userOptional.isPresent()) {
                return ResponseEntity.status(HttpStatus.NOT_FOUND).body("User not found");
            }

            User user = userOptional.get();

            // historyName 중복 체크
            String validatedHistoryName = historyService.validateHistoryName(historyName, userId);

            History history = new History();
            history.setUser(user);
            history.setHistoryName(validatedHistoryName);
            historyService.createHistory(history);

            // JsonData 생성 및 저장
            JsonData jsonData = new JsonData();
            jsonData.setReferenceId(referenceId);
            jsonData.setAlignment(alignmentJson);
            jsonData.setHistory(history);
            jsonDataService.saveJsonData(jsonData);

            // JSON 응답 생성
            Map<String, String> combinedJson = new HashMap<>();
            combinedJson.put("alignment", jsonData.getAlignment());
            combinedJson.put("historyName", validatedHistoryName);

            String jsonResponse = new Gson().toJson(combinedJson);
            return ResponseEntity.ok(jsonResponse);
        } catch (IOException e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Fasta 파일 저장에 문제 발생하였습니다.");
        }
    }
}
