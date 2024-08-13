package VirusDecode.backend.controller;


import VirusDecode.backend.dto.ReferenceSequenceRequest;
import VirusDecode.backend.dto.VarientSequenceRequest;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.springframework.core.io.ClassPathResource;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

@RestController
@RequestMapping("/inputSeq")
public class inputSeqController {

    @PostMapping("/reference")
    public ResponseEntity<String> processDone(@RequestBody ReferenceSequenceRequest request) {
        String sequenceId = request.getSequenceId();
        // 여기서 sequenceId를 사용하여 필요한 처리를 수행합니다.
        System.out.println("Processing DONE for sequence ID: " + sequenceId);

        Map<String, Object> metadata = new HashMap<>();
        metadata = referenceIdPy(sequenceId);

        // Frontend로 전달 값
        return ResponseEntity.ok(metadata.toString());
    }

    @PostMapping("/analyze")
    public ResponseEntity<String> handleSequences(@RequestBody VarientSequenceRequest request) {
        // sequence의 value들을 처리
        StringBuilder processedMessage = new StringBuilder("Processed sequences: ");

        Map<String, String> sequences = request.getSequences();
        sequences.forEach((key, value) -> {
            processedMessage.append(key).append(": ").append(value).append("; ");
        });

        // 처리된 결과를 JSON 형식으로 클라이언트에게 전송
        return ResponseEntity.ok(processedMessage.toString());
    }



    // reference id --> metadata ( MAP )
    public static Map<String, Object> referenceIdPy(String referenceId) {
        Map<String, Object> metadata = new HashMap<>();
        try {
            // 파이썬 스크립트 경로를 ClassPathResource를 사용하여 얻기
            ClassPathResource resource = new ClassPathResource("bioinformatics/alignment_and_translation.py");
            String scriptPath = resource.getFile().getAbsolutePath();

            // 파이썬 스크립트와 인자를 설정
            String[] command = new String[]{"python3", scriptPath, referenceId};

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


}


