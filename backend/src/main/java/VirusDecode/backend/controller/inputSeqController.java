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
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

@RestController
@RequestMapping("/inputSeq")
public class inputSeqController {

    @PostMapping("/reference")
    public ResponseEntity<Map> processDone(@RequestBody ReferenceSequenceRequest request) {
        String sequenceId = request.getSequenceId();
        // 여기서 sequenceId를 사용하여 필요한 처리를 수행합니다.
        System.out.println("Processing DONE for sequence ID: " + sequenceId);

        Map<String, Object> metadata = new HashMap<>();
        metadata = referenceIdPy(sequenceId);

        // Frontend로 전달 값
        return ResponseEntity.ok(metadata);
    }

    @PostMapping("/analyze")
    public ResponseEntity<Map<String, Object>> createFasta(@RequestBody(required = false) VarientSequenceRequest request) {
        StringBuilder fastaContent = new StringBuilder();

        // 시퀀스 데이터가 있는 경우에만 FASTA 형식으로 변환
        if (request.getSequences() != null && !request.getSequences().isEmpty()) {
            for (Map.Entry<String, String> entry : request.getSequences().entrySet()) {
                String sequenceName = entry.getKey();
                String sequenceData = entry.getValue();
                fastaContent.append(">").append(sequenceName).append("\n");
                fastaContent.append(sequenceData).append("\n");
            }
        }

        // 파일 데이터가 있는 경우에만 파일 내용 파싱 후 추가
        if (request.getFiles() != null && !request.getFiles().isEmpty()) {
            for (FastaFileDTO file : request.getFiles()) {
                fastaContent.append(file.getContent()).append("\n");
            }
        }

        // 고유한 파일 이름 생성 (UUID 또는 타임스탬프 사용)
        String uniqueFileName = "sequences_" + UUID.randomUUID().toString() + ".fasta";
        String outputPath = "./backend/src/main/resources/bioinformatics/test_temp_data_repository/" + uniqueFileName;

        // 파일로 저장
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(Paths.get(outputPath).toFile()))) {
            writer.write(fastaContent.toString());
        } catch (IOException e) {
            // 오류 발생 시, Map에 오류 메시지를 담아 반환
            Map<String, Object> errorResponse = new HashMap<>();
            errorResponse.put("status", "error");
            errorResponse.put("message", "파일 저장 중 오류 발생: " + e.getMessage());

            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body(errorResponse);
        }

        // 파이썬 스크립트 실행 (데이터가 없는 경우에도 파이썬 스크립트가 필요하다면 여전히 호출)
        Map<String, Object> analyzeResult = new HashMap<>();
        analyzeResult = analyzePy(outputPath);

        return ResponseEntity.ok(analyzeResult);
    }


    // reference id --> metadata ( MAP )
    public static Map<String, Object> referenceIdPy(String referenceId) {
        Map<String, Object> metadata = new HashMap<>();
        try {
            // 파이썬 스크립트 경로를 ClassPathResource를 사용하여 얻기
            ClassPathResource resource = new ClassPathResource("bioinformatics/test_without_bio/test.py");
            String scriptPath = resource.getFile().getAbsolutePath();

            // 파이썬 스크립트와 인자를 설정
            String[] command = new String[]{"python3", scriptPath, "1", referenceId};

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
            System.out.println("Output for python script: " + output);

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

    //
    public static Map<String, Object> analyzePy(String fastaFilePath) {
        Map<String, Object> analyzeResult = new HashMap<>();
        try {
            // 파이썬 스크립트 경로를 ClassPathResource를 사용하여 얻기
            ClassPathResource resource = new ClassPathResource("bioinformatics/test_without_bio/test.py");
            String scriptPath = resource.getFile().getAbsolutePath();

            // 파이썬 스크립트와 인자를 설정
            String[] command = new String[]{"python3", scriptPath, "2", fastaFilePath};

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
            analyzeResult = objectMapper.readValue(output.toString(), HashMap.class);

            // 프로세스가 완료될 때까지 대기
            process.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return analyzeResult;
    }


}


