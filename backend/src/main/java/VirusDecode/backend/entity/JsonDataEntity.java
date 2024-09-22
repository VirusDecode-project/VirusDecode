package VirusDecode.backend.entity;

import jakarta.persistence.Entity;
import jakarta.persistence.Id;
import jakarta.persistence.Column;
import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;

@Entity
@Data
@NoArgsConstructor
@AllArgsConstructor
public class JsonDataEntity {

    @Id
    private String id;  // This will be the unique key to identify each JSON object (e.g., "metadata", "alignment")

    @Column(columnDefinition = "TEXT")  // Store JSON data as a large text field in the database
    private String jsonData;  // This will store the actual JSON string

}
