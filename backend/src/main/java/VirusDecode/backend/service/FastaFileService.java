package VirusDecode.backend.service;

import VirusDecode.backend.dto.FastaFileDTO;
import VirusDecode.backend.dto.VarientDTO;
import org.springframework.stereotype.Service;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

@Service
public class FastaFileService {
    private static final Path currentDir = Paths.get("").toAbsolutePath();
    // FASTA 형식의 콘텐츠를 저장하는 메서드
    public String saveFastaContent(VarientDTO request) throws IOException {
        StringBuilder fastaContent = new StringBuilder();  // FASTA 콘텐츠를 저장할 StringBuilder 객체

        // 시퀀스 데이터가 있는 경우에만 FASTA 형식으로 변환하여 추가
        if (request.getSequences() != null && !request.getSequences().isEmpty()) {
            for (Map.Entry<String, String> entry : request.getSequences().entrySet()) {
                String sequenceName = entry.getKey();
                String sequenceData = entry.getValue();

                // 시퀀스 이름에서 공백 제거 및 시퀀스 데이터가 빈 문자열이 아닌 경우에만 처리
                sequenceName = sequenceName.replace(" ", "");
                if (sequenceData != null && !sequenceData.trim().isEmpty()) {
                    fastaContent.append(">").append(sequenceName).append("\n");
                    fastaContent.append(sequenceData).append("\n");
                }
            }
        }

        // 파일 데이터가 있는 경우에만 파일 내용을 파싱하여 추가
        if (request.getFiles() != null && !request.getFiles().isEmpty()) {
            for (FastaFileDTO file : request.getFiles()) {
                fastaContent.append(file.getContent()).append("\n");
            }
        }

        // 파일 저장 경로 설정
//        String currentDir = System.getProperty("user.dir");  // 현재 작업 디렉토리 경로
//        String varient_fasta_Path = Paths.get(currentDir, "build/resources/main/bioinformatics/data/combined.fasta").toString();
        String varient_fasta_Path = currentDir.resolve("data/combined.fasta").toString();
        // FASTA 콘텐츠를 파일로 저장
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(varient_fasta_Path))) {
            writer.write(fastaContent.toString());  // 파일에 콘텐츠 쓰기
        }

        return varient_fasta_Path; // 저장된 파일 경로를 반환
    }
}
