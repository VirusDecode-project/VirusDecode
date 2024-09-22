package VirusDecode.backend.repository;

import VirusDecode.backend.entity.JsonDataEntity;
import org.springframework.stereotype.Repository;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

@Repository
public class InMemoryJsonDataRepository implements JsonDataRepository {

    private final Map<String, JsonDataEntity> jsonDataStore = new HashMap<>();

    @Override
    public JsonDataEntity save(JsonDataEntity jsonDataEntity) {
        jsonDataStore.put(jsonDataEntity.getId(), jsonDataEntity);
        return jsonDataEntity;
    }

    @Override
    public Optional<JsonDataEntity> findById(String id) {
        return Optional.ofNullable(jsonDataStore.get(id));
    }

    @Override
    public Iterable<JsonDataEntity> findAll() {
        return jsonDataStore.values();
    }

    @Override
    public void deleteById(String id) {
        jsonDataStore.remove(id);
    }

    // deleteAll() 메서드 구현
    @Override
    public void deleteAll() {
        jsonDataStore.clear();  // 저장된 모든 데이터를 삭제
    }
}
