//package VirusDecode.backend.repository;
//
//import VirusDecode.backend.entity.JsonDataEntity;
//import org.springframework.beans.factory.annotation.Autowired;
//import org.springframework.stereotype.Repository;
//
//import java.util.Optional;
//
//@Repository
//public class JpaJsonDataRepository implements JsonDataRepository {
//
//    private final org.springframework.data.jpa.repository.JpaRepository<JsonDataEntity, String> jpaRepository;
//
//    @Autowired
//    public JpaJsonDataRepository(org.springframework.data.jpa.repository.JpaRepository<JsonDataEntity, String> jpaRepository) {
//        this.jpaRepository = jpaRepository;
//    }
//
//    @Override
//    public JsonDataEntity save(JsonDataEntity jsonDataEntity) {
//        return jpaRepository.save(jsonDataEntity);
//    }
//
//    @Override
//    public Optional<JsonDataEntity> findById(String id) {
//        return jpaRepository.findById(id);
//    }
//
//    @Override
//    public Iterable<JsonDataEntity> findAll() {
//        return jpaRepository.findAll();
//    }
//
//    @Override
//    public void deleteById(String id) {
//        jpaRepository.deleteById(id);
//    }
//}
