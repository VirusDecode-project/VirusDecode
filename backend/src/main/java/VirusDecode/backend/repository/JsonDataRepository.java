package VirusDecode.backend.repository;

import VirusDecode.backend.entity.JsonData;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.data.jpa.repository.Modifying;
import org.springframework.data.jpa.repository.Query;
import org.springframework.data.repository.query.Param;

import java.util.List;
import java.util.Optional;

public interface JsonDataRepository extends JpaRepository<JsonData, Long> {
//    @Query("SELECT j FROM JsonData j WHERE j.historyName = :historyName AND j.user.id = :userId")
//    Optional<JsonData> findByHistoryNameAndUserId(@Param("historyName") String historyName, @Param("userId") Long userId);

    @Query("SELECT j FROM JsonData j WHERE j.historyName = :historyName AND j.user.id = :userId")
    JsonData findByHistoryNameAndUserId(@Param("historyName") String historyName, @Param("userId") Long userId);

    @Query("SELECT j.historyName FROM JsonData j WHERE j.user.id = :userId")
    List<String> findHistoryNamesByUserId(@Param("userId") Long userId);

    @Modifying
    @Query("UPDATE JsonData j SET j.historyName = :newName WHERE j.historyName = :historyName AND j.user.id = :userId")
    void updateHistoryName(@Param("historyName") String historyName, @Param("newName") String newName, @Param("userId") Long userId);

    @Modifying
    @Query("DELETE FROM JsonData j WHERE j.historyName = :historyName AND j.user.id = :userId")
    void deleteByHistoryNameAndUserId(@Param("historyName") String historyName, @Param("userId") Long userId);

}
