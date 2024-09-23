package VirusDecode.backend.service;

import VirusDecode.backend.entity.JsonDataEntity;
import VirusDecode.backend.repository.JsonDataRepository;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.Optional;

@Service
public class JsonDataService {

    private final JsonDataRepository jsonDataRepository;

    @Autowired
    public JsonDataService(JsonDataRepository jsonDataRepository) {
        this.jsonDataRepository = jsonDataRepository;
    }

    // Store JSON data
    public void saveJsonData(String name, String jsonData) {
        JsonDataEntity jsonDataEntity = new JsonDataEntity();
        jsonDataEntity.setName(name);
        jsonDataEntity.setJsonData(jsonData);
        jsonDataRepository.save(jsonDataEntity);
    }

    // Retrieve JSON data
    public String getJsonData(String name) {
        Optional<JsonDataEntity> jsonDataEntity = jsonDataRepository.findByName(name);
        return jsonDataEntity.map(JsonDataEntity::getJsonData).orElse(null);
    }

    // Retrieve all JSON data
    public Iterable<JsonDataEntity> getAllJsonData() {
        return jsonDataRepository.findAll();
    }

    public void deleteAllJsonData() {
        jsonDataRepository.deleteAll();  // 기존 데이터 삭제
    }

    // Delete JSON data by key
    public void deleteJsonData(String name) {
        jsonDataRepository.deleteByName(name);
    }

    public Optional<JsonDataEntity> findByName(String name) {
        return jsonDataRepository.findByName(name);
    }
}
