package VirusDecode.backend.repository;

import VirusDecode.backend.entity.JsonData;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.data.jpa.repository.Modifying;
import org.springframework.data.jpa.repository.Query;
import org.springframework.data.repository.query.Param;

import java.util.List;
import java.util.Optional;

public interface JsonDataRepository extends JpaRepository<JsonData, Long> {
    JsonData findByHistoryId(Long historyId);

    @Modifying
    @Query("DELETE FROM JsonData j WHERE j.history.id = :historyId")
    void deleteByHistoryId(@Param("historyId") Long historyId);

}
