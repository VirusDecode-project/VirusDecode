package VirusDecode.backend.controller;

import VirusDecode.backend.dto.LinearDesignDTO;
import VirusDecode.backend.dto.PdbDTO;
import VirusDecode.backend.service.JsonDataService;
import VirusDecode.backend.service.PythonScriptExecutor;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import java.util.Map;

@RestController
@RequestMapping("/analysis")
public class AnalysisController {
    private final PythonScriptExecutor pythonScriptExecutor;
    private JsonDataService jsonDataService;

    @Autowired
    public AnalysisController(PythonScriptExecutor pythonScriptExecutor, JsonDataService jsonDataService) {
        this.pythonScriptExecutor = pythonScriptExecutor;
        this.jsonDataService = jsonDataService;
    }


    // /analysis/mrnadesign 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/mrnadesign")
    public ResponseEntity<String> getMRNAdesign(@RequestBody LinearDesignDTO request) {
        // mRNA 디자인 요청에서 필요한 데이터를 추출
        String region = request.getRegion();
        String varientName = request.getVarientName();
        String start = String.valueOf(request.getStart());
        String end = String.valueOf(request.getEnd());

        String referenceId = jsonDataService.parseSequenceIdFromMetadata(jsonDataService.getJsonData("metadata"));
        String alignmentIndex = jsonDataService.jsonParser(jsonDataService.getJsonData("alignment"), "alignment_index");
        String alignedSequences = jsonDataService.jsonParser(jsonDataService.getJsonData("alignment"), "aligned_sequences");
        System.out.println(referenceId);
        System.out.println(alignmentIndex);
        System.out.println(alignedSequences);

        ResponseEntity<String> scriptResponse = pythonScriptExecutor.executePythonScript("3", referenceId, alignmentIndex, alignedSequences, region, varientName, start, end);
        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            jsonDataService.saveJsonData("linearDesign", scriptResponse.getBody());
        }
        return scriptResponse;
    }

    // /analysis/pdb 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/pdb")
    public ResponseEntity<String> getPdb(@RequestBody PdbDTO request) {
        String gene = request.getGene();
        String referenceId = jsonDataService.parseSequenceIdFromMetadata(jsonDataService.getJsonData("metadata"));
        String alignmentIndex = jsonDataService.jsonParser(jsonDataService.getJsonData("alignment"), "alignment_index");
        String alignedSequences = jsonDataService.jsonParser(jsonDataService.getJsonData("alignment"), "aligned_sequences");
        System.out.println(referenceId);
        System.out.println(alignmentIndex);
        System.out.println(alignedSequences);
        ResponseEntity<String> scriptResponse = pythonScriptExecutor.executePythonScript("4", referenceId, alignmentIndex, alignedSequences, gene);
        if (scriptResponse.getStatusCode().is2xxSuccessful()) {
            jsonDataService.saveJsonData("pdb", scriptResponse.getBody());
        }
        return scriptResponse;
    }

    // /analysis/re-alignment 엔드포인트에 대한 POST 요청 처리
    @GetMapping("/re-alignment")
    public ResponseEntity<String> sendJsonFile() {
        return ResponseEntity.ok(jsonDataService.getJsonData("alignment"));
    }

    // /analysis/re-mrnadesign 엔드포인트에 대한 POST 요청 처리
    @GetMapping("/re-mrnadesign")
    public ResponseEntity<String> re_mrna_design_json() {
        return ResponseEntity.ok(jsonDataService.getJsonData("linearDesign"));
    }

    // /analysis/render3d 엔드포인트에 대한 GET 요청 처리
    @GetMapping("/re-render3d")
    public ResponseEntity<String> return_pdb_list() {
        return ResponseEntity.ok(jsonDataService.getJsonData("pdb"));
    }
}
