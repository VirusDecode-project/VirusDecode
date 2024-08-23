package VirusDecode.backend.controller;

import VirusDecode.backend.dto.mRNADesignRequest;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.springframework.core.io.ClassPathResource;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


@RestController
@RequestMapping("/analysis")
public class analysisController {

    @PostMapping("/re-alignment")
    public ResponseEntity<String> sendJsonFile() {
        try {
            // ClassPathResource를 사용하여 클래스패스 내에서 JSON 파일을 읽음

            // GK - test path 주석
//            ClassPathResource resource = new ClassPathResource("bioinformatics/test_without_bio/analyze.json");
            ClassPathResource resource = new ClassPathResource("bioinformatics/data/alignment_data.json");

            Path filePath = resource.getFile().toPath();
            String jsonContent = Files.readString(filePath);

            // JSON 데이터를 응답으로 반환
            return ResponseEntity.ok(jsonContent);

        } catch (IOException e) {
            e.printStackTrace();
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Error reading JSON file");
        }
    }

    @PostMapping("/mrnadesign")
    public ResponseEntity<Map> getMRNAdesign(@RequestBody mRNADesignRequest request) {
        // 0. request 데이터 처리 //
        String region = request.getRegion();
        String varientName = request.getVarientName();
        int start = request.getStart();
        int end = request.getEnd();

        // 1. 데이터가 잘 추출되었는지 확인하는 디버깅 메시지 추가
        System.out.println("Region: " + region);
        System.out.println("Varient Name: " + varientName);
        System.out.println("Start Index: " + start);
        System.out.println("End Index: " + end);

        // 2. 파이썬 스크립트 실행 & 결과값(json) client로 전달 //
        Map<String, Object> mRNAResult = new HashMap<>();
        mRNAResult = mRNADesignPy(region, varientName, start, end);

        return ResponseEntity.ok(mRNAResult);
    }

    @PostMapping("/re-mrnadesign")
    public ResponseEntity<String> re_mrna_design_json() {
        try {
            // ClassPathResource를 사용하여 클래스패스 내에서 JSON 파일을 읽음
            // GK - test path 주석
//            ClassPathResource resource = new ClassPathResource("bioinformatics/test_without_bio/mRNA.json");
            ClassPathResource resource = new ClassPathResource("bioinformatics/data/linearDesign_data.json");


            Path filePath = resource.getFile().toPath();
            String jsonContent = Files.readString(filePath);

            // JSON 데이터를 응답으로 반환
            return ResponseEntity.ok(jsonContent);

        } catch (IOException e) {
            e.printStackTrace();
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Error reading JSON file");
        }
    }

    @GetMapping("/render3d")
    public ResponseEntity<Map<String, String>> return_pdb_list() {

        // 2. 파이썬 스크립트 실행 & 결과값(json) client로 전달 //
        Map<String, String> pdbList = new HashMap<>();
        pdbList = pdbListPy();

        return ResponseEntity.ok(pdbList);
    }



    public static Map<String, Object> mRNADesignPy(String region, String varientName, int start, int end ) {
        Map<String, Object> mRNADesignResult = new HashMap<>();
        try {
            // 파이썬 스크립트 경로를 ClassPathResource를 사용하여 얻기

            // GK - test path 주석
//            ClassPathResource resource = new ClassPathResource("bioinformatics/test_without_bio/test.py");
            ClassPathResource resource = new ClassPathResource("bioinformatics/virusdecode.py");

            String scriptPath = resource.getFile().getAbsolutePath();

            // 파이썬 스크립트와 인자를 설정
            String[] command = new String[]{"python3", scriptPath, "3", region, varientName, String.valueOf(start), String.valueOf(end)};

            // ProcessBuilder를 사용하여 프로세스를 시작
            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream(true);
            Process process = pb.start();


            // GK
            process.waitFor();

            ClassPathResource resultResource = new ClassPathResource("bioinformatics/data/linearDesign_data.json");

            File jsonFile = resultResource.getFile(); // 파일 객체로 변환
            Path filePath = jsonFile.toPath();  // 파일 경로로 변환
            // JSON 파일의 내용을 문자열로 읽음
            String jsonContent = new String(Files.readAllBytes(filePath));

            // GK - Debug
            System.out.println(jsonContent);

            // JSON 문자열을 Map 객체로 변환
            ObjectMapper objectMapper = new ObjectMapper();
            mRNADesignResult = objectMapper.readValue(jsonContent, HashMap.class);
            System.out.println(mRNADesignResult);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return mRNADesignResult;
    }

    public static Map<String, String> pdbListPy() {
        Map<String, String> pdbList = null;
        try {
            // 파이썬 스크립트 경로 설정

            // GK - test path 주석
//            ClassPathResource resource = new ClassPathResource("bioinformatics/test_without_bio/test.py");
            ClassPathResource resource = new ClassPathResource("bioinformatics/virusdecode.py");

            String scriptPath = resource.getFile().getAbsolutePath();

            // 파이썬 스크립트 실행 명령어 설정
            String[] command = new String[]{"python3", scriptPath,"4"};

            // ProcessBuilder를 사용해 파이썬 스크립트 실행
            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream(true);
            Process process = pb.start();

            // GK
            process.waitFor();

            ClassPathResource resultResource = new ClassPathResource("bioinformatics/data/pdb_data.json");

            File jsonFile = resultResource.getFile(); // 파일 객체로 변환
            Path filePath = jsonFile.toPath();  // 파일 경로로 변환
            // JSON 파일의 내용을 문자열로 읽음
            String jsonContent = new String(Files.readAllBytes(filePath));


            // JSON 문자열을 Map 객체로 변환
            ObjectMapper objectMapper = new ObjectMapper();
            pdbList = objectMapper.readValue(jsonContent, new TypeReference<Map<String, String>>() {});

            // GK - Debug
            System.out.println(pdbList);

        } catch (Exception e) {
            e.printStackTrace();
        }
        return pdbList;
    }

}
