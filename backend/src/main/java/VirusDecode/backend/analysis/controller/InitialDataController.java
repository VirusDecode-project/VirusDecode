package VirusDecode.backend.analysis.controller;

import VirusDecode.backend.analysis.service.InitialDataService;
import VirusDecode.backend.analysis.dto.initialData.ReferenceDto;
import VirusDecode.backend.analysis.dto.initialData.fasta.VarientDto;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

@RestController
@RequestMapping("/api/inputSeq")
public class InitialDataController {

    private final InitialDataService initialDataService;
    @Autowired
    public InitialDataController(InitialDataService initialDataService) {
        this.initialDataService = initialDataService;
    }

    @PostMapping("/metadata")
    public ResponseEntity<String> getMetadata(@RequestBody ReferenceDto request) {
        return initialDataService.processMetadata(request);
    }

    @PostMapping("/alignment")
    public ResponseEntity<String> getAlignment(@RequestBody(required = false) VarientDto request, HttpSession session) {
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }

        return initialDataService.processAlignment(request, userId);
    }



}