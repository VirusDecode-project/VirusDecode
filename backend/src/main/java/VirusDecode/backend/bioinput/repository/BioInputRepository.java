package VirusDecode.backend.bioinput.repository;

import VirusDecode.backend.bioinput.entity.MetaData;
import org.springframework.data.jpa.repository.JpaRepository;

public interface BioInputRepository extends JpaRepository<MetaData, String> {

}
