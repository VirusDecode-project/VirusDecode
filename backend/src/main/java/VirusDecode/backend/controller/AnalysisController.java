package VirusDecode.backend.controller;

import VirusDecode.backend.dto.LinearDesignDTO;
import VirusDecode.backend.service.PythonScriptExecutor;
import VirusDecode.backend.service.JsonFileService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import java.util.Map;

@RestController
@RequestMapping("/analysis")
public class AnalysisController {

    // /analysis/mrnadesign 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/mrnadesign")
    public ResponseEntity<String> getMRNAdesign(@RequestBody LinearDesignDTO request) {
        // mRNA 디자인 요청에서 필요한 데이터를 추출
        String region = request.getRegion();
        String varientName = request.getVarientName();
        String start = String.valueOf(request.getStart());
        String end = String.valueOf(request.getEnd());

        // Python 스크립트를 실행하여 mRNA 디자인 결과를 얻음
//        Map<String, Object> mRNAResult = PythonScriptExecutor.executePythonScriptAsObjectMap("bioinformatics/data/linearDesign_data.json","3", region, varientName, start, end);
        return PythonScriptExecutor.executePythonScript("bioinformatics/data/linearDesign_data.json","3", region, varientName, start, end);
    }

    // /analysis/re-alignment 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/re-alignment")
    public ResponseEntity<String> sendJsonFile() {
        // JSON 파일을 읽어서 반환
        return JsonFileService.readJsonFile("bioinformatics/data/alignment_data.json");
    }

    // /analysis/re-mrnadesign 엔드포인트에 대한 POST 요청 처리
    @PostMapping("/re-mrnadesign")
    public ResponseEntity<String> re_mrna_design_json() {
        // JSON 파일을 읽어서 반환
        return JsonFileService.readJsonFile("bioinformatics/data/linearDesign_data.json");
    }

    // /analysis/render3d 엔드포인트에 대한 GET 요청 처리
    @GetMapping("/render3d")
    public ResponseEntity<String> return_pdb_list() {
        // JSON 파일을 읽어서 반환
        return JsonFileService.readJsonFile("bioinformatics/data/pdb_data.json");
    }
}
