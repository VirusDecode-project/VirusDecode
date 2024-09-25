package VirusDecode.backend.controller;

import VirusDecode.backend.dto.initialData.ReferenceDto;
import VirusDecode.backend.dto.initialData.fasta.VarientDto;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.entity.User;
import VirusDecode.backend.service.*;
import com.google.gson.Gson;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

@RestController
@RequestMapping("/inputSeq")
public class InitialDataController {

    private final PythonScriptService pythonScriptService;
    private final FastaFileService fastaFileService;  // Fasta 파일 처리를 위한 서비스 주입
    private final JsonDataService jsonDataService;
    private final UserService userService;

    @Autowired
    public InitialDataController(PythonScriptService pythonScriptService, FastaFileService fastaFileService, JsonDataService jsonDataService, UserService userService) {
        this.pythonScriptService = pythonScriptService;
        this.fastaFileService = fastaFileService;
        this.jsonDataService = jsonDataService;
        this.userService = userService;
    }

    @PostMapping("/metadata")
    public ResponseEntity<String> getMetadata(@RequestBody ReferenceDto request) {
        String sequenceId = request.getSequenceId();
        ResponseEntity<String> scriptResponse = pythonScriptService.executePythonScript("1", sequenceId);
        return scriptResponse;
    }

    @PostMapping("/alignment")
    public ResponseEntity<String> getAlignment(@RequestBody(required = false) VarientDto request, HttpSession session) {
        String historyName = request.getHistoryName();
        String referenceId = request.getReferenceId();
        try {
            String fastaContent = fastaFileService.saveFastaContent(request);

            ResponseEntity<String> scriptResponse = pythonScriptService.executePythonScript("2", referenceId, fastaContent);
            String alignmentJson = scriptResponse.getBody();
            if (scriptResponse.getStatusCode().is2xxSuccessful()) {
                Long userId = (Long) session.getAttribute("userId");
                Optional<User> userOptional = userService.getUserById(userId);
                if (!userOptional.isPresent()) {
                    return ResponseEntity.status(HttpStatus.NOT_FOUND).body("User not found");
                }
                User user = userOptional.get();

                // historyName 중복 체크
                String originalHistoryName = historyName;
                int counter = 1;
                while (jsonDataService.getJsonData(historyName, userId) != null) {
                    historyName = originalHistoryName + "_" + counter;
                    counter++;
                }

                JsonData jsonData = new JsonData();
                jsonData.setReferenceId(referenceId);
                jsonData.setAlignment(alignmentJson);
                jsonData.setHistoryName(historyName);
                jsonData.setUser(user);
                jsonDataService.saveJsonData(jsonData);


                Map<String, String> combinedJson = new HashMap<>();
                combinedJson.put("alignment", jsonData.getAlignment());
                combinedJson.put("historyName", historyName);

                String jsonResponse = new Gson().toJson(combinedJson);
                return ResponseEntity.ok(jsonResponse);
            }
            else {
                return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Python script execution failed");
            }

        } catch (IOException e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Fasta 파일 저장에 문제 발생하였습니다.");
        }
    }



}