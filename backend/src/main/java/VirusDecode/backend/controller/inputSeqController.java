package VirusDecode.backend.controller;


import VirusDecode.backend.dto.FastaFileDTO;
import VirusDecode.backend.dto.ReferenceSequenceRequest;
import VirusDecode.backend.dto.VarientSequenceRequest;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.springframework.core.io.ClassPathResource;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

@RestController
@RequestMapping("/inputSeq")
public class inputSeqController {

    @PostMapping("/reference")
    public ResponseEntity<Map> getMetadata(@RequestBody ReferenceSequenceRequest request) {
        // 1. 사용자 입력(Sequence id) 추출 - (저장X)
        String sequenceId = request.getSequenceId();

        // 2. 파이썬 스크립트 실행 & 결과값(json) client로 전달 //
        Map<String, Object> metadata = new HashMap<>();
        metadata = referenceIdPy(sequenceId);

        // GK - 비정상 nucleotide 값 처리
        if (metadata == null || metadata.isEmpty()) {
            // 상태 코드 204 No Content 반환
            return ResponseEntity.status(HttpStatus.NO_CONTENT).body(null);
        }

        // 메타데이터가 정상적인 경우 상태 코드 200 OK와 함께 반환
        return ResponseEntity.ok(metadata);
    }

    @PostMapping("/alignment")
    public ResponseEntity<Map<String, Object>> getAlignment(@RequestBody(required = false) VarientSequenceRequest request) {
        StringBuilder fastaContent = new StringBuilder();

        // 0. fasta 파일들 + sequence 입력들 --> fasta 파일 1개 //
        // 시퀀스 데이터가 있는 경우에만 FASTA 형식으로 변환
        if (request.getSequences() != null && !request.getSequences().isEmpty()) {
            for (Map.Entry<String, String> entry : request.getSequences().entrySet()) {
                String sequenceName = entry.getKey();
                String sequenceData = entry.getValue();

                // GK - 시퀀스 데이터가 빈 문자열이 아닌 경우에만 처리
                sequenceName = sequenceName.replace(" ", "");   // 공백 제거
                if (sequenceData != null && !sequenceData.trim().isEmpty()) {
                    fastaContent.append(">").append(sequenceName).append("\n");
                    fastaContent.append(sequenceData).append("\n");
                }
            }
        }
        // 파일 데이터가 있는 경우에만 파일 내용 파싱 후 추가
        if (request.getFiles() != null && !request.getFiles().isEmpty()) {
            for (FastaFileDTO file : request.getFiles()) {
                fastaContent.append(file.getContent()).append("\n");
            }
        }

        // 1. input ( 사용자 입력 ) 저장 //
        // 파일 저장 경로 설정
        String currentDir = System.getProperty("user.dir");  // 현재 작업 디렉토리 경로
//        String varient_fasta_Path = Paths.get(currentDir, "backend/src/main/resources/bioinformatics/User_input_data/varient_for_alignment.fasta").toString();

        // GK - 경로 수정: ClassPathResource build 디렉토리 내에서 경로 검색
        String varient_fasta_Path = Paths.get(currentDir, "build/resources/main/bioinformatics/data/combined.fasta").toString();
        // 파일로 저장
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(varient_fasta_Path))) {
            writer.write(fastaContent.toString());
        } catch (IOException e) {
            // 오류 발생 시, Map에 오류 메시지를 담아 반환
            Map<String, Object> errorResponse = new HashMap<>();
            errorResponse.put("status", "error");
            errorResponse.put("message", "파일 저장 중 오류 발생: " + e.getMessage());

            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body(errorResponse);
        }


        // 2. 파이썬 스크립트 실행 & 결과값(json) client로 전달 //
        Map<String, Object> alignmentResult = new HashMap<>();
        alignmentResult = alignmentPy(varient_fasta_Path);

        return ResponseEntity.ok(alignmentResult);
    }


    public static Map<String, Object> referenceIdPy(String referenceId) {
        Map<String, Object> metadata = new HashMap<>();
        try {
            // 1. 파이썬 스크립트 실행 //
            // 파이썬 스크립트 경로를 ClassPathResource를 사용하여 얻기

            // GK - test path 주석
//            ClassPathResource resource = new ClassPathResource("bioinformatics/test_without_bio/test.py");
            ClassPathResource resource = new ClassPathResource("bioinformatics/virusdecode.py");

            String scriptPath = resource.getFile().getAbsolutePath();
            // 파이썬 스크립트와 인자를 설정
            String[] command = new String[]{"python3", scriptPath, "1", referenceId};
            // ProcessBuilder를 사용하여 프로세스를 시작
            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream(true);
            Process process = pb.start();


            // 2. 파이썬 스크립트 출력값 읽어오기 & 반환 //
            // 프로세스의 출력을 읽기 위한 BufferedReader
            BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
            StringBuilder output = new StringBuilder();
            String line;
            while ((line = in.readLine()) != null) {
                output.append(line);
            }
            in.close();

            // GK - Debug
            System.out.println(output);

            // GK - 비정상 nucleotide 값 처리
            if(output == null || output.isEmpty()){
                // GK - Debug
                System.out.println("Empty metadata");
                return null;
            }

            // JSON 문자열을 Map 객체로 변환
            ObjectMapper objectMapper = new ObjectMapper();
            metadata = objectMapper.readValue(output.toString(), HashMap.class);
            // 프로세스가 완료될 때까지 대기
            process.waitFor();

        } catch (Exception e) {
            e.printStackTrace();
        }
        return metadata;
    }


    public static Map<String, Object> alignmentPy(String fastaString) {
        Map<String, Object> analyzeResult = new HashMap<>();
        try {

            // 1. 파이썬 스크립트 실행 //
            // 파이썬 스크립트 경로를 ClassPathResource를 사용하여 얻기

            // GK - test path 주석
//            ClassPathResource resource = new ClassPathResource("bioinformatics/test_without_bio/test.py");
            ClassPathResource resource = new ClassPathResource("bioinformatics/virusdecode.py");

            String scriptPath = resource.getFile().getAbsolutePath();
            // 파이썬 스크립트와 인자를 설정
            String[] command = new String[]{"python3", scriptPath, "2", fastaString};
            // ProcessBuilder를 사용하여 프로세스를 시작
            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream(true);
            Process process = pb.start();

            // GK - 프로세스가 완료될 때까지 대기
            process.waitFor();

            // 2. 파이썬 결과 파일( 저장된 ) 반환 //
            // 파이썬 스크립트 실행 결과 파일 불러오기 & 보내기

            // GK - test path 주석
//            ClassPathResource result_resource = new ClassPathResource("bioinformatics/test_without_bio/analyze.json");
            ClassPathResource result_resource = new ClassPathResource("bioinformatics/data/alignment_data.json");

            File jsonFile = result_resource.getFile(); // 파일 객체로 변환
            Path filePath = jsonFile.toPath();  // 파일 경로로 변환
            // JSON 파일의 내용을 문자열로 읽음
            String jsonContent = new String(Files.readAllBytes(filePath));

            // GK - Debug
            System.out.println(jsonContent);

            // JSON 문자열을 Map 객체로 변환
            ObjectMapper objectMapper = new ObjectMapper();
            analyzeResult = objectMapper.readValue(jsonContent, HashMap.class);

        } catch (Exception e) {
            e.printStackTrace();
        }
        return analyzeResult;
    }



}


