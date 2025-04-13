package VirusDecode.backend.bioinput.service;

import VirusDecode.backend.analysis.service.AnalysisService;
import VirusDecode.backend.bioinput.entity.FastaFile;
import VirusDecode.backend.bioinput.exception.FastaFileSaveFailException;
import VirusDecode.backend.common.biopython.BioPythonService;
import VirusDecode.backend.user.service.UserService;
import VirusDecode.backend.common.biopython.BioPythonDto;
import VirusDecode.backend.analysis.dto.AlignmentDto;
import VirusDecode.backend.bioinput.dto.ReferenceDto;
import VirusDecode.backend.bioinput.dto.VarientSequenceDto;
import VirusDecode.backend.history.entity.History;
import VirusDecode.backend.analysis.entity.Analysis;
import VirusDecode.backend.bioinput.entity.MetaData;
import VirusDecode.backend.user.entity.User;
import VirusDecode.backend.history.service.HistoryService;
import VirusDecode.backend.bioinput.repository.BioInputRepository;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import java.io.IOException;
import java.util.Map;
import java.util.Optional;

@Log4j2
@Service
public class BioInputService {
    private final BioPythonService bioPythonService;
    private final UserService userService;
    private final HistoryService historyService;
    private final BioInputRepository bioInputRepository;
    private final AnalysisService analysisService;

    @Autowired
    public BioInputService(BioPythonService bioPythonService, UserService userService, HistoryService historyService, BioInputRepository bioInputRepository, AnalysisService analysisService) {
        this.bioPythonService = bioPythonService;
        this.userService = userService;
        this.historyService = historyService;
        this.bioInputRepository = bioInputRepository;
        this.analysisService = analysisService;
    }

    public MetaData getMetadata(ReferenceDto referenceDto) {
        String sequenceId = referenceDto.getSequenceId();

        Optional<MetaData> existingData = bioInputRepository.findById(sequenceId);
        if (existingData.isPresent()) {
            return existingData.get();
        }

        BioPythonDto response = bioPythonService.executePythonScript("1", sequenceId);
        MetaData metaData = new MetaData(response.getOutput(), sequenceId);
        bioInputRepository.save(metaData);
        return metaData;
    }

    public AlignmentDto processAlignment(VarientSequenceDto varientSequenceDto, Long userId) {
        String historyName = varientSequenceDto.getHistoryName();
        String referenceId = varientSequenceDto.getReferenceId();

        String fastaContent = makeFastaContent(varientSequenceDto);
        BioPythonDto response = bioPythonService.executePythonScript("2", referenceId, fastaContent);
        String alignmentJson = response.getOutput();

        User user = userService.getUserById(userId);

        // historyName 중복 체크
        String validatedHistoryName = historyService.validateHistoryName(historyName, userId);

        History history = new History(user, validatedHistoryName);
        historyService.createHistory(history);

        // JsonData 생성 및 저장
        Analysis analysis = new Analysis(referenceId, alignmentJson, history);
        analysisService.saveAnalysisData(analysis);

        // JSON 응답 생성
        return new AlignmentDto(analysis.getAlignment(), validatedHistoryName);
    }


    // FASTA 형식의 콘텐츠를 저장하는 메서드
    public String makeFastaContent(VarientSequenceDto request) {
        try{

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
                for (FastaFile file : request.getFiles()) {
                    fastaContent.append(file.getContent()).append("\n");
                }
            }

            return fastaContent.toString();
        }catch (Exception e){
            throw new FastaFileSaveFailException("Fasta 파일 저장에 문제 발생하였습니다.");
        }
    }
}
