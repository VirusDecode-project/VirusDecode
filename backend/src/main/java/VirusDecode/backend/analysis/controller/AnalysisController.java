package VirusDecode.backend.analysis.controller;

import VirusDecode.backend.analysis.dto.LinearDesignDto;
import VirusDecode.backend.analysis.dto.PdbDto;
import VirusDecode.backend.analysis.service.AnalysisService;
import jakarta.servlet.http.HttpSession;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

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
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }
        String linearDesignJson = analysisService.processLinearDesign(request, userId);
        if(linearDesignJson==null){
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("데이터 분석 실패");
        }
        return ResponseEntity.ok(linearDesignJson);
    }

    @PostMapping("/pdb")
    public ResponseEntity<String> getPdb(@RequestBody PdbDto request, HttpSession session) {
        Long userId = (Long) session.getAttribute("userId");
        if (userId == null) {
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("User not authenticated");
        }

        String pdbJson = analysisService.processPdb(request, userId);
        if(pdbJson==null){
            return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body("데이터 분석 실패");
        }
        return ResponseEntity.ok(pdbJson);
    }
}
