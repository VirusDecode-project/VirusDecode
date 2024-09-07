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
        return pythonScriptExecutor.executePythonScript( "metadata.json","1", sequenceId);
    }

    // /inputSeq/alignment 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/alignment")
    public ResponseEntity<String> getAlignment(@RequestBody(required = false) VarientDTO request) {
        try {
            // 서비스 호출을 통해 사용자 입력을 파일로 저장
            String savedFilePath = fastaFileService.saveFastaContent(request);

            // 결과를 상태 코드 200 OK와 함께 반환
            return pythonScriptExecutor.executePythonScript("alignment_data.json", "2");
        } catch (IOException e) {
            // 파일 저장 중 오류가 발생한 경우 상태 코드 500과 함께 오류 메시지 반환
//            Map<String, Object> errorResponse = new HashMap<>();
//            errorResponse.put("status", "error");
//            errorResponse.put("message", "파일 저장 중 오류 발생: " + e.getMessage());

            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("파일 저장 중 오류 발생");
        }
    }
}