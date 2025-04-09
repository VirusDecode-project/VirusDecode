package VirusDecode.backend.common.config;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import java.util.HashMap;
import java.util.Map;

@RestController
@RequestMapping("/api/config")
public class ConfigController {

    @Value("${max.interval.length}")
    private int maxIntervalLength;

    @GetMapping("/max-interval")
    public Map<String, Object> getMaxIntervalLength() {
        Map<String, Object> response = new HashMap<>();
        response.put("MAX_INTERVAL_LENGTH", maxIntervalLength);
        return response;
    }
}
