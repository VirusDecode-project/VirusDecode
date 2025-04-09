package VirusDecode.backend.bioinput.controller;

import VirusDecode.backend.analysis.dto.AlignmentDto;
import VirusDecode.backend.bioinput.entity.MetaData;
import VirusDecode.backend.bioinput.service.BioInputService;
import VirusDecode.backend.bioinput.dto.ReferenceDto;
import VirusDecode.backend.bioinput.dto.VarientSequenceDto;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

@RestController
@RequestMapping("/api/inputSeq")
public class BioInputController {

    private final BioInputService bioInputService;
    @Autowired
    public BioInputController(BioInputService bioInputService) {
        this.bioInputService = bioInputService;
    }

    @PostMapping("/metadata")
    public ResponseEntity<String> getMetadata(@RequestBody ReferenceDto referenceDto) {
        MetaData metaData = bioInputService.getMetadata(referenceDto);

        if(metaData ==null){
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("메타데이터를 불러오는데 실패하였습니다.");
        }else{
            return ResponseEntity.ok(metaData.getMetadata());
        }
    }

    @PostMapping("/alignment")
    public ResponseEntity<?> getAlignment(@RequestBody(required = false) VarientSequenceDto varientSequenceDto, HttpSession session) {
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }

        AlignmentDto alignmentDto = bioInputService.processAlignment(varientSequenceDto, userId);
        if(alignmentDto==null){
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Fasta 파일 저장에 문제 발생하였습니다.");
        }else{
            return ResponseEntity.ok(alignmentDto);
        }
    }


}