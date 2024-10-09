package VirusDecode.backend.service;

import VirusDecode.backend.entity.History;
import VirusDecode.backend.repository.HistoryRepository;
import jakarta.transaction.Transactional;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.Collections;
import java.util.List;

@Service
public class HistoryService {

    @Autowired
    private HistoryRepository historyRepository;

    @Transactional
    public History createHistory(History history){
        return historyRepository.save(history);
    }

    public History getHistory(String historyName, Long userId){
        return historyRepository.findByHistoryNameAndUserId(historyName, userId);
    }


    public List<String> getHistoryNamesByUserId(Long userId) {
        List<String> historyNames = historyRepository.findHistoryNamesByUserId(userId);
        Collections.reverse(historyNames);
        return historyNames;
    }

    @Transactional
    public void updateHistoryName(String historyName, String newName, Long userId) {
        historyRepository.updateHistoryName(historyName, newName, userId);
    }

    @Transactional
    public void deleteHistory(String historyName, Long userId){
        historyRepository.deleteByHistoryNameAndUserId(historyName, userId);
    }


    public String validateHistoryName(String historyName, Long userId) {
        String originalHistoryName = historyName;
        int counter = 1;
        while (getHistory(historyName, userId) != null) {
            historyName = originalHistoryName + "_" + counter;
            counter++;
        }
        return historyName;
    }
}
