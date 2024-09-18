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
                'applyForce': false,
                'allowPanningAndZooming': true,
                'circularizeExternal': true,
                'labelInterval': 10,
                'initialSize':[600,500]
            };
    
            const container = new FornaContainer(containerRef.current, options);
    
            container.addRNA(structure, { sequence: sequence });
    
            return () => {
                container.clearNodes();
            };
        }
        
    }, [sequence, structure]);

    return <div ref={containerRef} style={{ width: '600px', height: '500px' }} />;
};

export default RNAVisualizer;

