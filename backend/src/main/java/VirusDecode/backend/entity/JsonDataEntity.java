package VirusDecode.backend.entity;

import jakarta.persistence.*;

@Entity
@Table(name = "history_data")
public class HistoryData {

    @Id @GeneratedValue(strategy = GenerationType.IDENTITY)
    private Long id;

    @Id
    @Column(name = "data_type", nullable = false)
    private String dataType;  // This will be the unique key to identify each JSON object (e.g., "metadata", "alignment")

//    @Column(name = "data_value", columnDefinition = "TEXT")
    @Column(name = "data_value", columnDefinition = "json")
    private String dataValue;  // This will store the actual JSON string

    @ManyToOne
    @JoinColumn(name = "history_id")
    private History history;  // 해당 데이터가 속하는 히스토리


    public Long getId() {
        return id;
    }

    public void setId(Long id) {
        this.id = id;
    }

    public String getDataType() {
        return dataType;
    }

    public void setDataType(String dataType) {
        this.dataType = dataType;
    }

    public void setDataValue(String dataValue){
        this.dataValue = dataValue;
    }
    public String getDataValue(){
        return dataValue;
    }
    public History getHistory() {
        return history;
    }

    public void setHistory(History history) {
        this.history = history;
    }
}
