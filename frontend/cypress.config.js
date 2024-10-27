import { defineConfig } from "cypress";

export default defineConfig({
  e2e: {
    setupNodeEvents(on, config) {
      // implement node event listeners here
    },
    video: false,
    defaultCommandTimeout: 20000,  // 20초 타임아웃 설정
  },
  projectId: "virusdecode",
});
