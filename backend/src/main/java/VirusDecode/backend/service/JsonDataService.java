package VirusDecode.backend.service;

import VirusDecode.backend.entity.History;
import VirusDecode.backend.entity.JsonData;
import VirusDecode.backend.repository.JsonDataRepository;
import jakarta.transaction.Transactional;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

@Service
public class JsonDataService {

    @Autowired
    private JsonDataRepository jsonDataRepository;

    @Transactional
    public JsonData saveJsonData(JsonData jsonData) {
        return jsonDataRepository.save(jsonData);
    }

    public JsonData getJsonData(History history) {
        Long historyId = history.getId();
        return jsonDataRepository.findByHistoryId(historyId);
    }
    @Transactional
    public void deleteJsonData(History history){
        Long historyId = history.getId();
        jsonDataRepository.deleteByHistoryId(historyId);
    }

}
