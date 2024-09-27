package VirusDecode.backend.service;

import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.repository.JsonDataRepository;
import jakarta.transaction.Transactional;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

@Service
public class JsonDataService {

    @Autowired
    private JsonDataRepository jsonDataRepository;

    @Transactional
    public JsonData saveJsonData(JsonData jsonData) {
        return jsonDataRepository.save(jsonData);
    }

    public JsonData getJsonData(String historyName, Long userId) {
        return jsonDataRepository.findByHistoryNameAndUserId(historyName, userId);
    }

    public List<String> getHistoryNamesByUserId(Long userId) {
        List<String> historyNames = jsonDataRepository.findHistoryNamesByUserId(userId);
        Collections.reverse(historyNames);
        return historyNames;
    }

    @Transactional
    public void updateHistoryName(String historyName, String newName, Long userId) {
        jsonDataRepository.updateHistoryName(historyName, newName, userId);
    }

    @Transactional
    public void deleteHistory(String historyName, Long userId){
        jsonDataRepository.deleteByHistoryNameAndUserId(historyName, userId);
    }

}
