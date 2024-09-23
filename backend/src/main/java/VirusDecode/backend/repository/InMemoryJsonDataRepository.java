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
        jsonDataStore.put(jsonDataEntity.getName(), jsonDataEntity);
        return jsonDataEntity;
    }

//    @Override
//    public Optional<JsonDataEntity> findById(Long id) {
//        return Optional.ofNullable(jsonDataStore.get(id));
//    }

    @Override
    public Optional<JsonDataEntity> findByName(String name) {
        return Optional.ofNullable(jsonDataStore.get(name));
    }

    @Override
    public Iterable<JsonDataEntity> findAll() {
        return jsonDataStore.values();
    }

//    @Override
//    public void deleteById(String id) {
//        jsonDataStore.remove(id);
//    }
    @Override
    public void deleteByName(String name) {
        jsonDataStore.remove(name);
    }

    // deleteAll() 메서드 구현
    @Override
    public void deleteAll() {
        jsonDataStore.clear();  // 저장된 모든 데이터를 삭제
    }
}
