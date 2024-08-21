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
    public ResponseEntity<Map<String, Object>> getMRNAdesign(@RequestBody(required = false) mRNADesignRequest request) {
        StringBuilder varientContent = new StringBuilder();

        // 0. request 데이터 처리 //
        String region = request.getRegion();
        String varientName = request.getVarientName();
        int start = request.getStart();
        int end = request.getEnd();

        // 1. input ( 사용자 입력 ) 저장 //
        varientContent.append("{\n");
        varientContent.append("\"region\": \"").append(region).append("\",\n");
        varientContent.append("\"varientName\": \"").append(varientName).append("\",\n");
        varientContent.append("\"start\": ").append(start).append(",\n");
        varientContent.append("\"end\": ").append(end).append("\n");
        varientContent.append("}");
        // 파일 저장 경로 설정 (예: 현재 작업 디렉토리 내의 저장 위치)
        String currentDir = System.getProperty("user.dir");  // 현재 작업 디렉토리 경로
//        String jsonFilePath = Paths.get(currentDir, "backend/src/main/resources/bioinformatics/User_input_data/mrnadesign_request.json").toString();

        // GK - 경로 수정: ClassPathResource build 디렉토리 내에서 경로 검색
        String jsonFilePath = Paths.get(currentDir, "build/resources/main/bioinformatics/User_input_data/mrnadesign_request.json").toString();

        // 파일로 저장
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(jsonFilePath))) {
            writer.write(varientContent.toString());
        } catch (IOException e) {
            // 오류 발생 시, Map에 오류 메시지를 담아 반환
            Map<String, Object> errorResponse = new HashMap<>();
            errorResponse.put("status", "error");
            errorResponse.put("message", "파일 저장 중 오류 발생: " + e.getMessage());

            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body(errorResponse);
        }

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

            // 프로세스의 출력을 읽기 위한 BufferedReader
            BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
            StringBuilder output = new StringBuilder();

            String line;
            while ((line = in.readLine()) != null) {
                output.append(line);
            }
            in.close();

            // 파이썬 스크립트 출력 결과 확인
            // System.out.println("Output for python script: " + output);

            // JSON 문자열을 Map 객체로 변환
            ObjectMapper objectMapper = new ObjectMapper();
            mRNADesignResult = objectMapper.readValue(output.toString(), HashMap.class);

            // 프로세스가 완료될 때까지 대기
            process.waitFor();
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

            // 파이썬 스크립트 출력을 읽기
            BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
            StringBuilder output = new StringBuilder();
            String line;
            while ((line = in.readLine()) != null) {
                output.append(line);
            }
            in.close();

            // 프로세스가 완료될 때까지 대기
            process.waitFor();

            // 출력된 JSON 문자열을 Map<String, String>로 변환
            ObjectMapper objectMapper = new ObjectMapper();
            pdbList = objectMapper.readValue(output.toString(), new TypeReference<Map<String, String>>() {});

        } catch (Exception e) {
            e.printStackTrace();
        }
        return pdbList;
    }

}
