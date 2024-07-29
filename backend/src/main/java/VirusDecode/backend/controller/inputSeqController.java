package VirusDecode.backend.controller;


import lombok.Getter;
import lombok.Setter;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.*;

@RestController
@RequestMapping("/inputSeq")
public class inputSeqController {

    @PostMapping("/reference")
    public ResponseEntity<String> processDone(@RequestBody SequenceRequest request) {
        String sequenceId = request.getSequenceId();
        // 여기서 sequenceId를 사용하여 필요한 처리를 수행합니다.
        System.out.println("Processing DONE for sequence ID: " + sequenceId);

        return ResponseEntity.ok("DONE processing completed for sequence ID: " + sequenceId);
    }

}

// Request Body를 받을 DTO 클래스 정의
@Getter
@Setter
class SequenceRequest {
    private String sequenceId;
}

