package VirusDecode.backend.history.service;

import VirusDecode.backend.history.entity.History;
import VirusDecode.backend.history.exception.HistoryNotFoundException;
import VirusDecode.backend.history.repository.HistoryRepository;
import VirusDecode.backend.user.entity.User;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.junit.jupiter.MockitoExtension;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

@ExtendWith(MockitoExtension.class)
public class HistoryServiceTest {
    // @Mock으로 만들어진 historyRepository를 자동으로 주입
    @InjectMocks
    private HistoryService historyService;

    @Mock
    private HistoryRepository historyRepository;

    private User user;
    private History history;

    @BeforeEach
    void setup(){
        user = new User("John", "Doe", "testuser", "encodedpassword", "USER");
        user.setId(1L);
        history = new History(user, "testHistoryName");
    }

    @Test
    @DisplayName("히스토리 생성 - 성공")
    void createHistory_Success(){
        when(historyRepository.save(history)).thenReturn(history);

        History result = historyService.createHistory(history);

        assertEquals(history, result);
        verify(historyRepository).save(history);
    }

    @Test
    @DisplayName("히스토리 조회 - 성공")
    void getHistory_Success(){
        when(historyRepository.findByHistoryNameAndUserId("sample", 1L)).thenReturn(history);

        History result = historyService.getHistory("sample", 1L);

        assertEquals(history, result);
    }

    @Test
    @DisplayName("히스토리 조회 - 실패")
    void getHistory_ShouldThrowException_WhenNotFound(){
        when(historyRepository.findByHistoryNameAndUserId("not_exist", 1L)).thenReturn(null);

        assertThrows(HistoryNotFoundException.class, ()->
            historyService.getHistory("not_exist", 1L)
        );
    }

    @Test
    @DisplayName("히스토리 이름 목록 가져오기 - 최신 순 정렬 확인")
    void getHistoryNamesByUserId_ShouldReturnReversedList() {
        List<String> historyNames = Arrays.asList("h1", "h2", "h3");
        when(historyRepository.findHistoryNamesByUserId(1L)).thenReturn(historyNames);

        List<String> result = historyService.getHistoryNamesByUserId(1L);

        assertEquals(Arrays.asList("h3", "h2", "h1"), result);
    }

    @Test
    @DisplayName("히스토리 이름 업데이트 호출 확인")
    void updateHistoryName_ShouldCallRepositoryMethod() {
        historyService.updateHistoryName("old", "new", 1L);

        verify(historyRepository).updateHistoryName("old", "new", 1L);
    }

    @Test
    @DisplayName("히스토리 삭제 호출 확인")
    void deleteHistory_ShouldCallRepositoryMethod() {
        historyService.deleteHistory("deleteMe", 1L);

        verify(historyRepository).deleteByHistoryNameAndUserId("deleteMe", 1L);
    }

    @Test
    @DisplayName("히스토리 이름 중복 확인 및 새 이름 반환")
    void validateHistoryName_ShouldReturnUniqueName() {
        when(historyRepository.findByHistoryNameAndUserId("history", 1L)).thenReturn(new History());
        when(historyRepository.findByHistoryNameAndUserId("history_1", 1L)).thenReturn(new History());
        when(historyRepository.findByHistoryNameAndUserId("history_2", 1L)).thenReturn(null);

        String uniqueName = historyService.validateHistoryName("history", 1L);

        assertEquals("history_2", uniqueName);
    }

    @Test
    @DisplayName("히스토리 이름이 중복되지 않으면 그대로 반환")
    void validateHistoryName_ShouldReturnOriginalName_IfNotExists() {
        when(historyRepository.findByHistoryNameAndUserId("uniqueName", 1L)).thenReturn(null);

        String result = historyService.validateHistoryName("uniqueName", 1L);

        assertEquals("uniqueName", result);
    }
}
