package VirusDecode.backend.repository;

import VirusDecode.backend.entity.JsonDataEntity;
import jakarta.persistence.TypedQuery;
import org.springframework.stereotype.Repository;
import jakarta.persistence.EntityManager;

import java.util.List;
import java.util.Optional;

@Repository
public class JpaJsonDataRepository implements JsonDataRepository {

    private final EntityManager em;

    public JpaJsonDataRepository(EntityManager em) {
        this.em = em;
    }
    @Override
    public JsonDataEntity save(JsonDataEntity jsonDataEntity) {
        TypedQuery<JsonDataEntity> query = em.createQuery(
                "SELECT j FROM JsonDataEntity j WHERE j.name = :name", JsonDataEntity.class);
        query.setParameter("name", jsonDataEntity.getName());
        List<JsonDataEntity> existingEntities = query.getResultList();

        if (!existingEntities.isEmpty()) {
            JsonDataEntity existingEntity = existingEntities.get(0);
            existingEntity.setJsonData(jsonDataEntity.getJsonData());
            return em.merge(existingEntity);
        } else {
            em.persist(jsonDataEntity);
            return jsonDataEntity;
        }
    }

    @Override
    public Optional<JsonDataEntity> findByName(String name) {
        TypedQuery<JsonDataEntity> query = em.createQuery(
                "SELECT j FROM JsonDataEntity j WHERE j.name = :name", JsonDataEntity.class);
        query.setParameter("name", name);

        List<JsonDataEntity> results = query.getResultList();

        // 결과가 있으면 Optional에 값을 담아 반환, 없으면 Optional.empty() 반환
        return results.isEmpty() ? Optional.empty() : Optional.of(results.get(0));
    }

    @Override
    public Iterable<JsonDataEntity> findAll() {
        return em.createQuery("SELECT j FROM JsonDataEntity j", JsonDataEntity.class).getResultList();
    }
    @Override
    public void deleteByName(String name) {
        em.createQuery("DELETE FROM JsonDataEntity j WHERE j.name = :name")
                .setParameter("name", name)
                .executeUpdate();
    }
    @Override
    public void deleteAll() {
        em.createQuery("DELETE FROM JsonDataEntity").executeUpdate();
    }
}
