package VirusDecode.backend.service;

import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.repository.JsonDataRepository;
import jakarta.transaction.Transactional;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

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

    public JsonData getJsonData(String historyName) {
        return jsonDataRepository.findByHistoryName(historyName);
    }

    public List<String> getHistoryNames() {
        List<String> historyNames = jsonDataRepository.findHistoryNames();
        Collections.reverse(historyNames);
        return historyNames;
    }

    @Transactional
    public void updateHistoryName(String historyName, String newName) {
        jsonDataRepository.updateHistoryName(historyName, newName);
    }

    @Transactional
    public void deleteHistory(String historyName){
        jsonDataRepository.deleteByHistoryName(historyName);
    }

}
