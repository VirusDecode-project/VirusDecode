package VirusDecode.backend.repository;

import VirusDecode.backend.entity.NcbiData;
import org.springframework.data.jpa.repository.JpaRepository;

public interface NcbiRepository extends JpaRepository<NcbiData, String> {

}
