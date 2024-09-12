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
        String sequenceId = request.getSequenceId();  // 요청에서 시퀀스 ID 추출
        pythonScriptExecutor.executePythonScript( "1", sequenceId);
        return JsonFileService.readJsonFile("metadata.json");
    }

    // /inputSeq/alignment 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/alignment")
    public ResponseEntity<String> getAlignment(@RequestBody(required = false) VarientDTO request) {
        try {
            String savedFilePath = fastaFileService.saveFastaContent(request);
            pythonScriptExecutor.executePythonScript("2");
            return JsonFileService.readJsonFile("alignment_data.json");
        } catch (IOException e) {
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("파일 저장 중 오류 발생");
        }
    }
}