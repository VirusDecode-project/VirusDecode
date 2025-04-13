package VirusDecode.backend.analysis.controller;

import VirusDecode.backend.analysis.dto.LinearDesignDto;
import VirusDecode.backend.analysis.dto.PdbDto;
import VirusDecode.backend.analysis.service.AnalysisService;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import static VirusDecode.backend.util.UserSessionUtil.getAuthenticatedUserId;

@RestController
@RequestMapping("/api/analysis")
public class AnalysisController {
    private final AnalysisService analysisService;

    @Autowired
    public AnalysisController(AnalysisService analysisService) {
        this.analysisService = analysisService;
    }

    @PostMapping("/linearDesign")
    public ResponseEntity<String> getLinearDesign(@RequestBody LinearDesignDto request, HttpSession session) {
        Long userId = getAuthenticatedUserId(session);

        String linearDesignJson = analysisService.processLinearDesign(request, userId);
        return ResponseEntity.ok(linearDesignJson);
    }

    @PostMapping("/pdb")
    public ResponseEntity<String> getPdb(@RequestBody PdbDto request, HttpSession session) {
        Long userId = getAuthenticatedUserId(session);

        String pdbJson = analysisService.processPdb(request, userId);
        return ResponseEntity.ok(pdbJson);
    }
}
