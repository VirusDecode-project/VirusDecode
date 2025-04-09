package VirusDecode.backend.common.biopython;

import lombok.extern.log4j.Log4j2;
import org.springframework.stereotype.Service;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Log4j2
@Service
public class BioPythonService {
    private static final String pythonScriptPath = Paths.get("").toAbsolutePath().resolve("../bioinformatics/main.py").normalize().toString();

    // ProcessBuilder를 생성하는 메서드를 따로 분리하여 테스트에서 모의 가능하도록 설계
    public ProcessBuilder createProcessBuilder(List<String> command) {
        return new ProcessBuilder(command);
    }

    public BioPythonDto executePythonScript(String... args) {
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
            return new BioPythonDto(exitCode, output.toString(), errorOutput.toString());


        } catch (Exception e) {
            log.error("알 수 없는 오류가 발생했습니다: {}", e.getMessage());
            return new BioPythonDto(-1, "", e.getMessage());
        }
    }
}
