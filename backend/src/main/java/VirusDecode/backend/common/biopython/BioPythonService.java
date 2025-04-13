package VirusDecode.backend.common.biopython;

import VirusDecode.backend.common.exception.BioPythonException;
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
            // 정상 실행
            if (exitCode == 0) {
                return new BioPythonDto(exitCode, output.toString(), errorOutput.toString());
            } else {
                String errorMessage = switch (exitCode) {
                    case 1 -> "필요한 파이썬 환경이 제대로 설치되지 않았습니다.";
                    case 11 -> "NCBI에 요청한 nucleotide ID가 존재하지 않습니다.";
                    case 21 -> "MUSCLE 다중 서열 정리에 문제가 발생하였습니다.";
                    case 22 -> "입력하신 서열 정보가 올바르지 않습니다. A, T, C, 그리고 G만 허용됩니다.";
                    case 31 -> "서버 메모리 부족으로 LinearDesign 실행 중 문제가 발생하였습니다.";
                    case 32 -> "LinearDesign 실행파일이 정상적으로 만들어지지 않았습니다.";
                    case 33 -> "LinearDesign 디렉토리가 원하는 위치에 없습니다.";
                    case 41 -> "PDB ID 검색 실패.";
                    case 42 -> "3D viewer 데이터 로드 실패.";
                    default -> "Unknown error (Exit Code: " + exitCode + "): " + errorOutput.toString();
                };
                throw new BioPythonException("바이오 파이썬 실행 중 오류가 발생했습니다.\n" + errorMessage);
            }

        } catch (Exception e) {
            throw new BioPythonException("바이오 파이썬 실행 중 오류가 발생했습니다. " + e.getMessage());
        }
    }
}
