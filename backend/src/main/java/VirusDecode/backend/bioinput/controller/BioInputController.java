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

import static VirusDecode.backend.util.UserSessionUtil.getAuthenticatedUserId;

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
        return ResponseEntity.ok(metaData.getMetadata());
    }

    @PostMapping("/alignment")
    public ResponseEntity<?> getAlignment(@RequestBody(required = false) VarientSequenceDto varientSequenceDto, HttpSession session) {
        Long userId = getAuthenticatedUserId(session);

        AlignmentDto alignmentDto = bioInputService.processAlignment(varientSequenceDto, userId);
        return ResponseEntity.ok(alignmentDto);
    }
}