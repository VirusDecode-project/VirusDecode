package virusdecode.backend.bioinput.controller;

import virusdecode.backend.analysis.dto.AlignmentDto;
import virusdecode.backend.bioinput.entity.MetaData;
import virusdecode.backend.bioinput.service.BioInputService;
import virusdecode.backend.bioinput.dto.ReferenceDto;
import virusdecode.backend.bioinput.dto.VarientSequenceDto;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import static virusdecode.backend.util.UserSessionUtil.getAuthenticatedUserId;

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