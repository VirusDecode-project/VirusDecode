package VirusDecode.backend.controller;

import VirusDecode.backend.dto.LinearDesignDTO;
import VirusDecode.backend.dto.PdbDTO;
import VirusDecode.backend.service.PythonScriptExecutor;
import VirusDecode.backend.service.JsonFileService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import java.util.Map;

@RestController
@RequestMapping("/analysis")
public class AnalysisController {
    private final PythonScriptExecutor pythonScriptExecutor;
    private final JsonFileService jsonFileService;

    @Autowired
    public AnalysisController(PythonScriptExecutor pythonScriptExecutor, JsonFileService jsonFileService) {
        this.pythonScriptExecutor = pythonScriptExecutor;
        this.jsonFileService = jsonFileService;
    }


    // /analysis/mrnadesign 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/mrnadesign")
    public ResponseEntity<String> getMRNAdesign(@RequestBody LinearDesignDTO request) {
        // mRNA 디자인 요청에서 필요한 데이터를 추출
        String region = request.getRegion();
        String varientName = request.getVarientName();
        String start = String.valueOf(request.getStart());
        String end = String.valueOf(request.getEnd());

        return pythonScriptExecutor.executePythonScript("3", region, varientName, start, end);
    }

    // /analysis/pdb 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/pdb")
    public ResponseEntity<String> getPdb(@RequestBody PdbDTO request) {
        String gene = request.getGene();
        return pythonScriptExecutor.executePythonScript("4", gene);
    }

    // /analysis/re-alignment 엔드포인트에 대한 POST 요청 처리
    @GetMapping("/re-alignment")
    public ResponseEntity<String> sendJsonFile() {
        return jsonFileService.readJsonFile("alignment_data.json");
    }

    // /analysis/re-mrnadesign 엔드포인트에 대한 POST 요청 처리
    @GetMapping("/re-mrnadesign")
    public ResponseEntity<String> re_mrna_design_json() {
        return jsonFileService.readJsonFile("linearDesign_data.json");
    }

    // /analysis/render3d 엔드포인트에 대한 GET 요청 처리
    @GetMapping("/re-render3d")
    public ResponseEntity<String> return_pdb_list() {
        return jsonFileService.readJsonFile("pdb_data.json");
    }
}
