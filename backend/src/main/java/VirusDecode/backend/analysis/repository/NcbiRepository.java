package VirusDecode.backend.analysis.repository;

import VirusDecode.backend.analysis.entity.NcbiData;
import org.springframework.data.jpa.repository.JpaRepository;

public interface NcbiRepository extends JpaRepository<NcbiData, String> {

}
