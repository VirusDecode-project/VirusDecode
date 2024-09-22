package VirusDecode.backend.service;

import org.springframework.core.io.ClassPathResource;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

@Service
public class PythonScriptExecutor {
    private static final Path currentDir = Paths.get("").toAbsolutePath();
    private static final Path pythonScriptPath = currentDir.resolve("../virusdecode.py").normalize();
    private static final Logger logger = LoggerFactory.getLogger(PythonScriptExecutor.class);

    public static ResponseEntity<String> executePythonScript(String... args) {
        try {
            List<String> command = new ArrayList<>();
            // 운영 체제에 따라 python3 또는 python 사용
            String os = System.getProperty("os.name").toLowerCase();
            if (os.contains("win")) {
                command.add("python");  // Windows에서는 python 사용
            } else {
                command.add("python3"); // Linux, MacOS에서는 python3 사용
            }
//            command.add("python3");
            command.add(pythonScriptPath.toString());
            command.add(currentDir.toString());

            if (args != null) {
                for (String arg : args) {
                    command.add(arg);
                }
            }

            ProcessBuilder pb = new ProcessBuilder(command);
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
                if (!output.isEmpty()) {
                    logger.info("Python script output: \n{}", output);
                }
                if (!errorOutput.isEmpty()) {
                    logger.error("Python script error output: \n{}", errorOutput);
                }
                return switch (exitCode) {
                    case 1 -> ResponseEntity.status(500).body("필요한 파이썬 환경이 제대로 설치되지 않았습니다.\nVirusDecode Github를 참고하세요.");
                    case 2 -> ResponseEntity.status(500).body("전달된 인자가 부족합니다.");
                    case 11 -> ResponseEntity.status(500).body("NCBI에 요청한 nucleotide ID가 존재하지 않습니다.");
                    case 21 -> ResponseEntity.status(500).body("MUSCLE 다중 서열 정리에 문제가 발생하였습니다.");
                    case 31 -> ResponseEntity.status(500).body("선택된 구간에 유효한 서열이 없습니다.");
                    case 32 -> ResponseEntity.status(500).body("LinearDesign 실행파일이 정상적으로 만들어지지 않았습니다.\nLinux 또는 Max 사용자가 맞으신가요?");
                    case 33 -> ResponseEntity.status(500).body("LinearDesign 디렉토리가 원하는 위치가 존재하지 않습니다.");
                    case 41 -> ResponseEntity.status(500).body("RCSB PDB 서버로부터 PDB ID 검색에 실패하였습니다.");
                    case 42 -> ResponseEntity.status(500).body("RCSB PDB 서버의 문제로 3D viewer 데이터 로드에 실패하였습니다.");
                    default -> ResponseEntity.status(500).body("Error executing Python script: " + errorOutput);
                };
            }
//            test
            if (!output.isEmpty()) {
                logger.info("Python script output: \n{}", output);
            }
            return ResponseEntity.ok(output.toString());
        } catch (IOException e) {
            logger.error("IO 오류가 발생했습니다: {}", e.getMessage());
            return ResponseEntity.status(500).body("Error executing Python script");
        } catch (InterruptedException e) {
            logger.error("프로세스가 중단되었습니다: {}", e.getMessage());
            Thread.currentThread().interrupt();
            return ResponseEntity.status(500).body("Process was interrupted");
        } catch (Exception e) {
            logger.error("알 수 없는 오류가 발생했습니다: {}", e.getMessage());
            return ResponseEntity.status(500).body("An unknown error occurred");
        }
    }
}
