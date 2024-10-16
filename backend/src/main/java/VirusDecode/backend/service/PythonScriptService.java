package VirusDecode.backend.service;

import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Service
public class PythonScriptService {
    private static final String pythonScriptPath = Paths.get("").toAbsolutePath().resolve("virusdecode.py").normalize().toString();
    private static final Logger logger = LoggerFactory.getLogger(PythonScriptService.class);

    // ProcessBuilder를 생성하는 메서드를 따로 분리하여 테스트에서 모의 가능하도록 설계
    protected ProcessBuilder createProcessBuilder(List<String> command) {
        return new ProcessBuilder(command);
    }

    public ResponseEntity<String> executePythonScript(String... args) {
        try {
            List<String> command = new ArrayList<>();
            command.add("python3");
            command.add(pythonScriptPath);
            command.addAll(Arrays.asList(args));

            ProcessBuilder pb = createProcessBuilder(command);
            Process process = pb.start();

            StringBuilder output = new StringBuilder();
            StringBuilder errorOutput = new StringBuilder();

            try (BufferedReader stdInput = new BufferedReader(new InputStreamReader(process.getInputStream()));
                 BufferedReader stdError = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {

                String line;
                while ((line = stdInput.readLine()) != null) {
                    output.append(line).append("\n");
                }

                while ((line = stdError.readLine()) != null) {
                    errorOutput.append(line).append("\n");
                }
            }

            int exitCode = process.waitFor();
            // 파이썬 오류 코드 처리
            if (exitCode != 0) {
                logger.error("Python script 종료 코드: {}", exitCode);
                if (!errorOutput.isEmpty()) {
                    logger.error("Python script error output: \n{}", errorOutput);
                }
                return switch (exitCode) {
                    case 1 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("필요한 파이썬 환경이 제대로 설치되지 않았습니다.\nVirusDecode Github를 참고하세요.");
                    case 2 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("전달된 인자가 부족합니다.");
                    case 11 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("NCBI에 요청한 nucleotide ID가 존재하지 않습니다.");
                    case 21 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("MUSCLE 다중 서열 정리에 문제가 발생하였습니다.");
                    case 22 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("입력하신 서열 정보가 올바르지 않습니다. A, T, C, 그리고 G만 허용됩니다.");
                    case 31 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("서버의 메모리 부족 문제로 LinearDesign 실행 중 문제가 발생하였습니다. 더 짧은 구간을 선택해 주세요.");
                    case 32 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("LinearDesign 실행파일이 정상적으로 만들어지지 않았습니다.\n서버는 Linux 또는 Mac에서 실행해야 하며, 도커를 사용하신 경우 아키텍쳐 문제일 수 있습니다.");
                    case 33 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("LinearDesign 디렉토리가 원하는 위치가 존재하지 않습니다.");
                    case 41 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("RCSB PDB 서버로부터 PDB ID 검색에 실패하였습니다.");
                    case 42 -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("RCSB PDB 서버의 문제로 3D viewer 데이터 로드에 실패하였습니다.");
                    default -> ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("Error executing Python script: " + errorOutput);
                };
            }
//            if (!output.isEmpty()) {
//                logger.info("Python script output: \n{}", output);
//            }
            return ResponseEntity.ok(output.toString());
        } catch (Exception e) {
            logger.error("Python Script 실행 과정에서 알 수 없는 오류가 발생했습니다: {}", e.getMessage());
            return ResponseEntity.status(HttpStatus.INTERNAL_SERVER_ERROR).body("An unknown error occurred during Python Script execution.");
        }
    }
}
