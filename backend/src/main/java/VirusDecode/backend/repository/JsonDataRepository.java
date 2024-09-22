package VirusDecode.backend.repository;

import VirusDecode.backend.entity.JsonDataEntity;

import java.util.Optional;

public interface JsonDataRepository {

    JsonDataEntity save(JsonDataEntity jsonDataEntity); // Save or update JSON data

    Optional<JsonDataEntity> findById(String id); // Find JSON data by ID

    Iterable<JsonDataEntity> findAll(); // Get all stored JSON data

    void deleteById(String id); // Delete JSON data by ID

    void deleteAll(); // Delete all JSON data (명시적으로 추가)
}
