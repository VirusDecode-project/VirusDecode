package VirusDecode.backend.service;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.springframework.core.io.ClassPathResource;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

@Service
public class PythonScriptExecutor {

    // Python 스크립트를 실행하고 결과를 Map<String, Object>로 반환하는 메서드
    public static Map<String, Object> executePythonScriptAsObjectMap(String resultJsonPath, String... args) {
        return executePythonScript(resultJsonPath, new TypeReference<Map<String, Object>>() {}, args);
    }

    // Python 스크립트를 실행하고 결과를 Map<String, String>으로 반환하는 메서드
    public static Map<String, String> executePythonScriptAsStringMap(String resultJsonPath, String... args) {
        return executePythonScript(resultJsonPath, new TypeReference<Map<String, String>>() {}, args);
    }

    // 제네릭 메서드로, Python 스크립트를 실행하고 결과를 특정 타입의 Map으로 반환
    private static <T> T executePythonScript(String resultJsonPath, TypeReference<T> typeRef, String... args) {
        T result = null;  // 결과를 저장할 변수
        try {
            // Python 스크립트 리소스를 로드
            ClassPathResource resource = new ClassPathResource("bioinformatics/virusdecode.py");

            if (!resource.exists()) {
                throw new FileNotFoundException("Python 스크립트를 찾을 수 없음: " + resource.getPath());
            }

            // Python 스크립트의 절대 경로 가져오기
            String scriptPath = resource.getFile().getAbsolutePath();

            // 명령어 리스트에 스크립트 경로 및 인자 추가
            List<String> command = new ArrayList<>();
            command.add("python3");
            command.add(scriptPath);

            // 제공된 인자들을 명령어 리스트에 추가
            if (args != null) {
                for (String arg : args) {
                    command.add(arg);
                }
            }

            // 프로세스 시작
            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream(true);  // 오류 스트림을 표준 출력 스트림으로 리다이렉트
            Process process = pb.start();

            // 프로세스 완료 대기
            process.waitFor();

            // 결과 JSON 파일 로드
            ClassPathResource resultResource = new ClassPathResource(resultJsonPath);

            if (!resultResource.exists()) {
                throw new FileNotFoundException("결과 JSON 파일을 찾을 수 없음: " + resultResource.getPath());
            }

            File jsonFile = resultResource.getFile();
            Path filePath = jsonFile.toPath();
            String jsonContent = new String(Files.readAllBytes(filePath));  // JSON 파일 내용을 문자열로 읽기

            // JSON 문자열을 제공된 타입 참조에 따라 Map으로 변환
            ObjectMapper objectMapper = new ObjectMapper();
            result = objectMapper.readValue(jsonContent, typeRef);

            // 디버깅 용도로 결과 출력
            System.out.println(result);
        } catch (Exception e) {
            e.printStackTrace();  // 예외 발생 시 스택 트레이스 출력
        }
        return result;  // 결과 반환
    }
}
