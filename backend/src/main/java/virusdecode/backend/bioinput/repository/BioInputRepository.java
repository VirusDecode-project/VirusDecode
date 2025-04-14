package virusdecode.backend.bioinput.repository;

import virusdecode.backend.bioinput.entity.MetaData;
import org.springframework.data.jpa.repository.JpaRepository;

public interface BioInputRepository extends JpaRepository<MetaData, String> {

}
