package VirusDecode.backend.entity;

import jakarta.persistence.*;
//import lombok.AllArgsConstructor;
import lombok.Data;
//import lombok.NoArgsConstructor;

@Entity
//@Data
//@NoArgsConstructor
//@AllArgsConstructor
public class JsonDataEntity {

    @Id @GeneratedValue(strategy = GenerationType.IDENTITY)
    private Long id;

    private String name;  // This will be the unique key to identify each JSON object (e.g., "metadata", "alignment")

    @Column(columnDefinition = "TEXT")  // Store JSON data as a large text field in the database
    private String jsonData;  // This will store the actual JSON string

    public Long getId() {
        return id;
    }

    public void setId(Long id) {
        this.id = id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setJsonData(String jsonData){
        this.jsonData=jsonData;
    }
    public String getJsonData(){
        return jsonData;
    }
}
