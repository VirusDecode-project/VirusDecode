package VirusDecode.backend.controller;

import VirusDecode.backend.dto.ReferenceDTO;
import VirusDecode.backend.dto.VarientDTO;
import VirusDecode.backend.service.FastaFileService;
import VirusDecode.backend.service.JsonFileService;
import VirusDecode.backend.service.PythonScriptExecutor;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.io.*;
@RestController
@RequestMapping("/inputSeq")
public class InputSeqController {

    private final PythonScriptExecutor pythonScriptExecutor;
    private FastaFileService fastaFileService;  // Fasta 파일 처리를 위한 서비스 주입

    @Autowired
    public InputSeqController(PythonScriptExecutor pythonScriptExecutor, FastaFileService fastaFileService) {
        this.pythonScriptExecutor = pythonScriptExecutor;
        this.fastaFileService = fastaFileService;
    }

    // /inputSeq/reference 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/reference")
    public ResponseEntity<String> getMetadata(@RequestBody ReferenceDTO request) {
        String sequenceId = request.getSequenceId();
        ResponseEntity<String> scriptResponse = pythonScriptExecutor.executePythonScript("1", sequenceId);

        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            return JsonFileService.readJsonFile("metadata.json");
        } else {
            return scriptResponse;
        }
    }

    // /inputSeq/alignment 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/alignment")
    public ResponseEntity<String> getAlignment(@RequestBody(required = false) VarientDTO request) {
        try {
            String fastaContent = fastaFileService.saveFastaContent(request);
            ResponseEntity<String> scriptResponse = pythonScriptExecutor.executePythonScript("2", fastaContent);
            if (scriptResponse.getStatusCode().is2xxSuccessful()) {
                return JsonFileService.readJsonFile("alignment.json");
            } else {
                return scriptResponse;
            }
        } catch (IOException e) {
            return ResponseEntity.status(500).body("Fasta 파일 저장에 문제 발생하였습니다.");
        }
    }
}