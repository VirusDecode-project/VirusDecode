package VirusDecode.backend.controller;

import VirusDecode.backend.config.ConfigController;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.util.ReflectionTestUtils;

import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.result.MockMvcResultHandlers.print;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.jsonPath;
import static org.hamcrest.Matchers.is;

@WebMvcTest(ConfigController.class)
public class ConfigControllerTest {

    @Autowired
    private MockMvc mockMvc;

    @Autowired
    private ConfigController configController;

    @BeforeEach
    public void setup() {
        // ReflectionTestUtils를 사용해 private 필드에 값을 주입
        ReflectionTestUtils.setField(configController, "maxIntervalLength", 25);
    }

    @Test
    public void testGetMaxIntervalLength() throws Exception {
        mockMvc.perform(get("/api/config/max-interval"))
                .andExpect(status().isOk())
                .andDo(print())
                .andExpect(jsonPath("$.MAX_INTERVAL_LENGTH", is(25)));  // 25는 설정된 테스트 값
    }
}
