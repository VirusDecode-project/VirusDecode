import React, { useState } from 'react';
import { useLocation } from 'react-router-dom';
import { Nav } from 'react-bootstrap';
import '../styles/Analysis.css';
import Alignment from './Alignment';
import MRNAdesign from './MRNAdesign';
import Render3D from './Render3D';

function Analysis({ tab, setTab, mRNAReceived, setMRNAReceived, PDBReceived, setPDBReceived, workingHistory }) {
  const location = useLocation();
  const responseData = location.state?.responseData || null;
  const [modalRegion, setModalRegion] = useState('');
  const handleModalRegion = (region) => {
      setModalRegion(region);
  };

  return (
    <div>
      <div className={"analysis-container"}>
        <>
          <Nav variant="tabs" defaultActiveKey="link0" className="justify-content-start">
            <Nav.Item>
              <Nav.Link eventKey="link0" active={tab === 0} onClick={() => setTab(0)}>Alignment</Nav.Link>
            </Nav.Item>
            <Nav.Item>
              <Nav.Link eventKey="link1" active={tab === 1} onClick={() => setTab(1)} disabled={!mRNAReceived}>mRNA design</Nav.Link>
            </Nav.Item>
            <Nav.Item>
              <Nav.Link eventKey="link2" active={tab === 2} onClick={() => setTab(2)} disabled={!PDBReceived}>3D viewer</Nav.Link>
            </Nav.Item>
          </Nav>
          <div>
            {tab === 0 && (
                <Alignment
                    responseData={responseData}
                    setMRNAReceived={setMRNAReceived}
                    setPDBReceived={setPDBReceived}
                    setTab={setTab}
                    onRegionUpdate={handleModalRegion}
                    workingHistory={workingHistory}
                />
            )}
            {tab === 1 && <MRNAdesign />}
            {tab === 2 && <Render3D region={modalRegion} />}
        </div>
        </>
      </div>
    </div>
  );
}

export default Analysis;
