import React, { useEffect, useRef } from 'react';
import { FornaContainer } from 'fornac';

const RNAVisualizer = ({ sequence, structure }) => {
    const containerRef = useRef(null);

    useEffect(() => {
        const options = {
            'applyForce': false,
            'allowPanningAndZooming': true,
            'circularizeExternal': true,
            'labelInterval': 10,
            'initialSize':[500,500]
        };

        const container = new FornaContainer(containerRef.current, options);

        container.addRNA(structure, { sequence: sequence });

        return () => {
            container.clearNodes();
        };
    }, [sequence, structure]);

    return <div ref={containerRef} style={{ width: '500px', height: '500px' }} />;
};

export default RNAVisualizer;

