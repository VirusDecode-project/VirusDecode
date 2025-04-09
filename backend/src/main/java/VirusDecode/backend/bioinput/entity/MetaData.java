package VirusDecode.backend.bioinput.entity;

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
public class MetaData {

    @Id
    @Column(nullable = false)
    private String referenceId;

    @Column(columnDefinition = "TEXT")
    private String metadata;

    public MetaData() {
    }

    public MetaData(String metadata, String referenceId) {
        this.metadata = metadata;
        this.referenceId = referenceId;
    }
}
