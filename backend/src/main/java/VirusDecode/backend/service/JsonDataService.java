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

    // metadata JSON에서 Sequence ID를 추출하는 간단한 메서드
    public String parseSequenceIdFromMetadata(String jsonString) {
        try {
            ObjectMapper objectMapper = new ObjectMapper();
            JsonNode jsonNode = objectMapper.readTree(jsonString);
            return jsonNode.get("Sequence ID").asText();  // Sequence ID 값 추출
        } catch (Exception e) {
            throw new RuntimeException("Failed to parse Sequence ID from metadata", e);
        }
    }

    // metadata JSON에서 Sequence ID를 추출하는 간단한 메서드
    public String jsonParser(String jsonString, String jsonKey) {
        try {
            ObjectMapper objectMapper = new ObjectMapper();
            JsonNode rootNode = objectMapper.readTree(jsonString);
            JsonNode jsonNode = rootNode.get(jsonKey);
            return jsonNode.toString();
        } catch (Exception e) {
            throw new RuntimeException("Failed to parse json data", e);
        }
    }
}
