package VirusDecode.backend.controller;

import VirusDecode.backend.dto.initialData.ReferenceDto;
import VirusDecode.backend.dto.initialData.fasta.VarientDto;
import VirusDecode.backend.service.FastaFileService;
import VirusDecode.backend.service.JsonDataService;
import VirusDecode.backend.service.PythonScriptService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.io.*;
@RestController
@RequestMapping("/inputSeq")
public class InitialDataController {

    private final PythonScriptService pythonScriptService;
    private final FastaFileService fastaFileService;  // Fasta 파일 처리를 위한 서비스 주입
    private final JsonDataService jsonDataService;

    @Autowired
    public InitialDataController(PythonScriptService pythonScriptService, FastaFileService fastaFileService, JsonDataService jsonDataService) {
        this.pythonScriptService = pythonScriptService;
        this.fastaFileService = fastaFileService;
        this.jsonDataService = jsonDataService;
    }

    // /inputSeq/reference 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/metadata")
    public ResponseEntity<String> getMetadata(@RequestBody ReferenceDto request) {
        String sequenceId = request.getSequenceId();
        ResponseEntity<String> scriptResponse = pythonScriptService.executePythonScript("1", sequenceId);

        // Memory Repository
        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            jsonDataService.saveJsonData("metadata", scriptResponse.getBody());
        }
        return scriptResponse;

    }

    // /inputSeq/alignment 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/alignment")
    public ResponseEntity<String> getAlignment(@RequestBody(required = false) VarientDto request) {
        try {
            String fastaContent = fastaFileService.saveFastaContent(request);
            String metadataJson = jsonDataService.getJsonData("metadata");

            ResponseEntity<String> scriptResponse = pythonScriptService.executePythonScript("2", metadataJson, fastaContent);

            // Memory Repository
            if (scriptResponse.getStatusCode().is2xxSuccessful()) {
                jsonDataService.saveJsonData("alignment", scriptResponse.getBody());
            }
            return scriptResponse;

        } catch (IOException e) {
            return ResponseEntity.status(500).body("Fasta 파일 저장에 문제 발생하였습니다.");
        }
    }



}