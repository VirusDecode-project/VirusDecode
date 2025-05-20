package virusdecode.backend.analysis.repository;

import virusdecode.backend.analysis.entity.Analysis;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.data.jpa.repository.Modifying;
import org.springframework.data.jpa.repository.Query;
import org.springframework.data.repository.query.Param;

public interface AnalysisRepository extends JpaRepository<Analysis, Long> {
    Analysis findByHistoryId(Long historyId);

    @Modifying
    @Query("DELETE FROM Analysis j WHERE j.history.id = :historyId")
    void deleteByHistoryId(@Param("historyId") Long historyId);

}
