import React, { useState } from 'react';
import { useLocation } from 'react-router-dom';
import { Nav } from 'react-bootstrap';
import './analysis.css';
import Alignment from '../pages/Alignment';
import MRNAdesign from '../pages/MRNAdesign';
import Render3D from '../pages/Render3D';

function Analysis({ tab, setTab, mRNAReceived, setMRNAReceived, PDBReceived, setPDBReceived }) {
  
  // const [mRNAReceived, setMRNAReceived] = useState(false);
  // const [PDBReceived, setPDBReceived] = useState(false);
  const [show, setShow] = useState(false);
  const location = useLocation();
  const responseData = location.state?.responseData || null;
  const [modalRegion, setModalRegion] = useState('');
  const handleModalRegion = (region) => {
      setModalRegion(region);
  };

  return (
    <div>
      <div className={`analysis-container ${show ? 'shrink' : ''}`}>
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
