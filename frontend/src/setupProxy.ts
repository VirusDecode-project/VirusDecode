// src/setupProxy.ts
import { createProxyMiddleware } from "http-proxy-middleware";
import { Express } from "express";

const setupProxy = (app: Express) => {
    app.use(
        "/api",
        createProxyMiddleware({
            target: "http://localhost:8080", // 비즈니스 서버 URL 설정
            changeOrigin: true,
        })
    );
};

export default setupProxy;
