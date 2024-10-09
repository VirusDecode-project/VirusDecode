package VirusDecode.backend.entity;

import jakarta.persistence.*;
import lombok.Getter;
import lombok.Setter;

@Getter
@Setter
@Entity
@Table(name="json_data")
public class JsonData {
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
}
