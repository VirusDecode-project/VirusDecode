package VirusDecode.backend.analysis.entity;

import VirusDecode.backend.history.entity.History;
import jakarta.persistence.*;
import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
@Entity
@Table(name="json_data")
public class Analysis {
    @Id
    @GeneratedValue(strategy = GenerationType.IDENTITY)
    private Long id;

    @Column(nullable = false)
    private String referenceId;

    @Column(columnDefinition = "MEDIUMTEXT")
    private String alignment;

    @Column(columnDefinition = "TEXT")
    private String linearDesign;

    @Column(columnDefinition = "TEXT")
    private String pdb;

    @ManyToOne
    @JoinColumn(nullable = false, name = "history_id")
    private History history;

    public Analysis() {
    }

    public Analysis(String referenceId, String alignment, History history) {
        this.referenceId = referenceId;
        this.alignment = alignment;
        this.history = history;
    }
}
