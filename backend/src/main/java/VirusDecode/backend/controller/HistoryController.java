package VirusDecode.backend.controller;

import VirusDecode.backend.dto.HistoryDTO;
import VirusDecode.backend.dto.LinearDesignDTO;
import VirusDecode.backend.service.JsonFileService;
import VirusDecode.backend.service.PythonScriptExecutor;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;

@RestController
@RequestMapping("/history")
public class HistoryController {

    @PostMapping("/create")
    public ResponseEntity<String> createHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        return PythonScriptExecutor.executePythonScriptWithoutJson("history.py",  "1", historyName);
    }

    @PostMapping("/delete")
    public ResponseEntity<String> deleteHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        return PythonScriptExecutor.executePythonScriptWithoutJson("history.py", "2", historyName);
    }

    @PostMapping("/rename")
    public ResponseEntity<String> renameHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        String newName = request.getNewName();
        return PythonScriptExecutor.executePythonScriptWithoutJson("history.py",  "3", historyName, newName);
    }

    @PostMapping("/get")
    public ResponseEntity<String> getHistory(@RequestBody HistoryDTO request) {
        String historyName = request.getHistoryName();
        return PythonScriptExecutor.executePythonScript("history.py", "bioinformatics/data/alignment_data.json", "4", historyName);
    }

    @GetMapping("/list")
    public ResponseEntity<String> listHistory() {
        return PythonScriptExecutor.executePythonScript("history.py","bioinformatics/history_list.json", "5");
    }
}
