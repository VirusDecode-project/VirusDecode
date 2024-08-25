import React, { useState, useEffect } from 'react';
import { useLocation } from 'react-router-dom';
import { Nav } from 'react-bootstrap';
import './analysis.css';
import Alignment from './Alignment';
import MRNAdesign from './MRNAdesign';
import Render3D from './Render3D';

function Analysis() {
  let [tab, setTab] = useState(0);
  const [show, setShow] = useState(false);
  const location = useLocation();
  const responseData = location.state?.responseData || null;

  return (
    <div>
      <div className={`analysis-container ${show ? 'shrink' : ''}`}>
        <>
          <Nav variant="tabs" defaultActiveKey="link0" className="justify-content-start">
            <Nav.Item>
              <Nav.Link eventKey="link0" active={tab === 0} onClick={() => setTab(0)}>Alignment</Nav.Link>
            </Nav.Item>
            <Nav.Item>
              <Nav.Link eventKey="link1" active={tab === 1} onClick={() => setTab(1)}>mRNA design</Nav.Link>
            </Nav.Item>
            <Nav.Item>
              <Nav.Link eventKey="link2" active={tab === 2} onClick={() => setTab(2)}>3D viewer</Nav.Link>
            </Nav.Item>

          </Nav>
          <Tab tab={tab} setTab={setTab} responseData={responseData} />
        </>
      </div>
    </div>
  );
}

function Tab({ tab, setTab, responseData }) {
  const [modalRegion, setModalRegion] = useState('');

  const handleModalRegion = (region) => {
    setModalRegion(region);
  }
  return (
    <div>
      {/* 현재 탭 상태에 따라 다른 컴포넌트 렌더링 */}
      {tab === 0 && <Alignment responseData={responseData} setTab={setTab} onRegionUpdate={handleModalRegion} />}
      {tab === 1 && <MRNAdesign />}
      {tab === 2 && <Render3D region={modalRegion} />}
    </div>
  );
}

export default Analysis;