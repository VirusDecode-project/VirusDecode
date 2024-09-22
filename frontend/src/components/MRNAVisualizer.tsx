import React, { useEffect, useRef } from 'react';
import { FornaContainer } from 'fornac';

interface RNAVisualizerProps {
    sequence: string;
    structure: string;
}

const RNAVisualizer: React.FC<RNAVisualizerProps> = ({ sequence, structure }) => {
    const containerRef = useRef<HTMLDivElement | null>(null);

    useEffect(() => {
        if (containerRef.current) {
            const options = {
                applyForce: false,
                allowPanningAndZooming: true,
                circularizeExternal: true,
                labelInterval: 10,
                initialSize: [600, 500],
            };

            const container = new FornaContainer(containerRef.current, options);
            container.addRNA(structure, { sequence: sequence });

            // 마우스 휠 이벤트 핸들러 추가
            const handleWheel = (event: WheelEvent) => {
                event.preventDefault(); // 페이지 스크롤 방지
            };

            const currentContainer = containerRef.current;
            currentContainer?.addEventListener('wheel', handleWheel);

            return () => {
                container.clearNodes();
                currentContainer?.removeEventListener('wheel', handleWheel); // 이벤트 리스너 제거
            };
        }
    }, [sequence, structure]);

    return <div ref={containerRef} style={{ width: '600px', height: '500px', overflow: 'hidden' }} />;
};

export default RNAVisualizer;
