package VirusDecode.backend.entity;

import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
@Entity
@Table(name="ncbi_data")
public class NcbiData {

    @Id
    @Column(nullable = false)
    private String referenceId;

    @Column(columnDefinition = "TEXT")
    private String metadata;
}
