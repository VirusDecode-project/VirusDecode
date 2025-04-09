package VirusDecode.backend.bioinput.service;

import VirusDecode.backend.analysis.service.AnalysisService;
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
import java.util.Optional;

@Log4j2
@Service
public class BioInputService {
    private final BioPythonService bioPythonService;
    private final FastaFileService fastaFileService;
    private final UserService userService;
    private final HistoryService historyService;
    private final BioInputRepository bioInputRepository;
    private final AnalysisService analysisService;

    @Autowired
    public BioInputService(BioPythonService bioPythonService, FastaFileService fastaFileService, UserService userService, HistoryService historyService, BioInputRepository bioInputRepository, AnalysisService analysisService) {
        this.bioPythonService = bioPythonService;
        this.fastaFileService = fastaFileService;
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

        if(response.isSuccess()){
            MetaData metaData = new MetaData(response.getOutput(), sequenceId);
            bioInputRepository.save(metaData);
            return metaData;
        }else{
            log.error(errorMessageByCode(response.getExitCode(), response.getErrorOutput()));
            return null;
        }
    }

    public AlignmentDto processAlignment(VarientSequenceDto varientSequenceDto, Long userId) {
        String historyName = varientSequenceDto.getHistoryName();
        String referenceId = varientSequenceDto.getReferenceId();

        try {
            String fastaContent = fastaFileService.saveFastaContent(varientSequenceDto);
            BioPythonDto response = bioPythonService.executePythonScript("2", referenceId, fastaContent);

            if(!response.isSuccess()){
                log.error(errorMessageByCode(response.getExitCode(), response.getErrorOutput()));
                System.out.println("response is not success");
                return null;
            }

            String alignmentJson = response.getOutput();

            Optional<User> userOptional = userService.getUserById(userId);
            if (!userOptional.isPresent()) {
                log.error("User not found");
                return null;
            }

            User user = userOptional.get();

            // historyName 중복 체크
            String validatedHistoryName = historyService.validateHistoryName(historyName, userId);

            History history = new History(user, validatedHistoryName);
            historyService.createHistory(history);

            // JsonData 생성 및 저장
            Analysis analysis = new Analysis(referenceId, alignmentJson, history);
            analysisService.saveAnalysisData(analysis);

            // JSON 응답 생성
            return new AlignmentDto(analysis.getAlignment(), validatedHistoryName);

        } catch (IOException e) {
            log.error("Fasta 파일 저장에 문제 발생하였습니다. " + e.getMessage());
            return null;
        }
    }



    private String errorMessageByCode(int exitCode, String errorOutput) {
        return switch (exitCode) {
            case 1 -> "필요한 파이썬 환경이 제대로 설치되지 않았습니다.";
            case 11 -> "NCBI에 요청한 nucleotide ID가 존재하지 않습니다.";
            case 21 -> "MUSCLE 다중 서열 정리에 문제가 발생하였습니다.";
            case 22 -> "입력하신 서열 정보가 올바르지 않습니다. A, T, C, 그리고 G만 허용됩니다.";
            case 31 -> "서버 메모리 부족으로 LinearDesign 실행 중 문제가 발생하였습니다.";
            case 32 -> "LinearDesign 실행파일이 정상적으로 만들어지지 않았습니다.";
            case 33 -> "LinearDesign 디렉토리가 원하는 위치에 없습니다.";
            case 41 -> "PDB ID 검색 실패.";
            case 42 -> "3D viewer 데이터 로드 실패.";
            default -> "Unknown error: " + errorOutput;
        };
    }


}
