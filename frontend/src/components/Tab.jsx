import React, { useState } from 'react';
import Alignment from '../pages/Alignment';
import MRNAdesign from '../pages/MRNAdesign';
import Render3D from '../pages/Render3D';

function Tab({ tab, setTab, responseData }) {
    const [modalRegion, setModalRegion] = useState('');

    const handleModalRegion = (region) => {
        setModalRegion(region);
    };

    return (
        <div>
            {/* 현재 탭 상태에 따라 다른 컴포넌트 렌더링 */}
            {tab === 0 && (
                <Alignment
                    responseData={responseData}
                    setTab={setTab}
                    onRegionUpdate={handleModalRegion}
                />
            )}
            {tab === 1 && <MRNAdesign />}
            {tab === 2 && <Render3D region={modalRegion} />}
        </div>
    );
}

export default Tab;
