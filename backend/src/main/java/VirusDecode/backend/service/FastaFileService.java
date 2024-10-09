package VirusDecode.backend.service;

import VirusDecode.backend.dto.initialData.fasta.FastaFileDto;
import VirusDecode.backend.dto.initialData.fasta.VarientDto;
import org.springframework.stereotype.Service;
import java.io.*;
import java.util.Map;

@Service
public class FastaFileService {
    // FASTA 형식의 콘텐츠를 저장하는 메서드
    public String saveFastaContent(VarientDto request) throws IOException {
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
            for (FastaFileDto file : request.getFiles()) {
                fastaContent.append(file.getContent()).append("\n");
            }
        }

        return fastaContent.toString();
    }
}
