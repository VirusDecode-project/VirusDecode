package VirusDecode.backend.repository;

import VirusDecode.backend.entity.History;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.data.jpa.repository.Modifying;
import org.springframework.data.jpa.repository.Query;
import org.springframework.data.repository.query.Param;

import java.util.List;

public interface HistoryRepository extends JpaRepository<History, Long> {
    @Query("SELECT h FROM History h WHERE h.historyName = :historyName AND h.user.id = :userId")
    History findByHistoryNameAndUserId(@Param("historyName") String historyName, @Param("userId") Long userId);

    @Query("SELECT j.historyName FROM History j WHERE j.user.id = :userId")
    List<String> findHistoryNamesByUserId(@Param("userId") Long userId);

    @Modifying
    @Query("UPDATE History j SET j.historyName = :newName WHERE j.historyName = :historyName AND j.user.id = :userId")
    void updateHistoryName(@Param("historyName") String historyName, @Param("newName") String newName, @Param("userId") Long userId);

    @Modifying
    @Query("DELETE FROM History j WHERE j.historyName = :historyName AND j.user.id = :userId")
    void deleteByHistoryNameAndUserId(@Param("historyName") String historyName, @Param("userId") Long userId);
}
