package VirusDecode.backend.service;

import org.springframework.core.io.ClassPathResource;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

@Service
public class PythonScriptExecutor {

    private static final Logger logger = LoggerFactory.getLogger(PythonScriptExecutor.class);

    public static ResponseEntity<String> executePythonScript(String jsonFilePath, String... args) {
        try {
            // Python 스크립트 리소스를 로드
            ClassPathResource pythonResource = new ClassPathResource("bioinformatics/virusdecode.py");
            if (!pythonResource.exists()) {
                logger.error("Python 스크립트를 찾을 수 없음: {}", pythonResource.getPath());
                return ResponseEntity.status(404).body("Python script not found");
            }


            // Python 스크립트의 절대 경로 가져오기
            String pythonScriptPath = pythonResource.getFile().getAbsolutePath();

            // 명령어 리스트에 스크립트 경로 및 인자 추가
            List<String> command = new ArrayList<>();
            command.add("python3");
            command.add(pythonScriptPath);

            // 제공된 인자들을 명령어 리스트에 추가
            if (args != null) {
                for (String arg : args) {
                    command.add(arg);
                }
            }

            // 프로세스 시작
            ProcessBuilder pb = new ProcessBuilder(command);
//            pb.redirectErrorStream(true);  // 오류 스트림을 표준 출력 스트림으로 리다이렉트
            Process process = pb.start();

            // 프로세스의 표준 출력과 오류 출력 읽기
            try (BufferedReader stdInput = new BufferedReader(new InputStreamReader(process.getInputStream()));
                 BufferedReader stdError = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {

                String line;
                StringBuilder output = new StringBuilder();
                while ((line = stdInput.readLine()) != null) {
                    output.append(line).append("\n");
                }

                StringBuilder errorOutput = new StringBuilder();
                while ((line = stdError.readLine()) != null) {
                    errorOutput.append(line).append("\n");
                }

                if (!output.isEmpty()) {
                    logger.info("Python script output: \n{}", output);
                }
                if (!errorOutput.isEmpty()) {
                    logger.error("Python script error output: \n{}", errorOutput);
                }
            }

            // 프로세스 완료 대기
            int exitCode = process.waitFor();

            if (exitCode != 0) {
                logger.error("Python script 종료 코드: {}", exitCode);
                return ResponseEntity.status(500).body("Error executing Python script");
            }

            return JsonFileService.readJsonFile(jsonFilePath);
        } catch (IOException e) {
            logger.error("IO 오류가 발생했습니다: {}", e.getMessage());
            return ResponseEntity.status(500).body("Error executing Python script");
        } catch (InterruptedException e) {
            logger.error("프로세스가 중단되었습니다: {}", e.getMessage());
            Thread.currentThread().interrupt();  // 인터럽트 상태 복원
            return ResponseEntity.status(500).body("Process was interrupted");
        } catch (Exception e) {
            logger.error("알 수 없는 오류가 발생했습니다: {}", e.getMessage());
            return ResponseEntity.status(500).body("An unknown error occurred");
        }
    }
}
