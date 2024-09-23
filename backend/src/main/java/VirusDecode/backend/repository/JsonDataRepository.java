package VirusDecode.backend.repository;

import VirusDecode.backend.entity.JsonDataEntity;

import java.util.Optional;

public interface JsonDataRepository {

    JsonDataEntity save(JsonDataEntity jsonDataEntity); // Save or update JSON data

//    Optional<JsonDataEntity> findById(Long id); // Find JSON data by ID
    Optional<JsonDataEntity> findByName(String name); // Find JSON data by Name

    Iterable<JsonDataEntity> findAll(); // Get all stored JSON data

//    void deleteById(Long id); // Delete JSON data by ID

    void deleteByName(String name); // Delete JSON data by Name

    void deleteAll(); // Delete all JSON data (명시적으로 추가)
}
