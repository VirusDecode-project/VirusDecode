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
    public void saveJsonData(String key, String jsonData) {
        JsonDataEntity jsonDataEntity = new JsonDataEntity(key, jsonData);
        jsonDataRepository.save(jsonDataEntity);
    }

    // Retrieve JSON data
    public String getJsonData(String key) {
        Optional<JsonDataEntity> jsonDataEntity = jsonDataRepository.findById(key);
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
    public void deleteJsonData(String key) {
        jsonDataRepository.deleteById(key);
    }

    public Optional<JsonDataEntity> findById(String id) {
        return jsonDataRepository.findById(id);
    }
}
