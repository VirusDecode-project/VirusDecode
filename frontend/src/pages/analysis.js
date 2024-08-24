import React, { useState, useEffect } from 'react';
import { useLocation } from 'react-router-dom';
import { Nav } from 'react-bootstrap';
import { useNavigate } from 'react-router-dom';
import loadingImage from '../image/loading.png';
import './analysis.css';
import Alignment from './Alignment';
import MRNAdesign from './MRNAdesign';
import Render3D from './Render3D';

function Analysis() {
  let [tab, setTab] = useState(0);
  let [isLoading, setIsLoading] = useState(true);
  let [loadingText, setLoadingText] = useState('Analyzing');
  let navigate = useNavigate();
  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);
  const [show, setShow] = useState(false);
  /* parkki */
  const [responseData, setResponseData] = useState(null); // 서버에서 받은 데이터를 상태로 관리
  /* parkki */

  useEffect(() => {
    const interval = setInterval(() => {
      setLoadingText((prev) => {
        if (prev === 'Analyzing...') return 'Analyzing';
        return prev + '.';
      });
    }, 500);

    setTimeout(() => {
      setIsLoading(false);
      clearInterval(interval);
    }, 3000);

    return () => clearInterval(interval);
  }, []);

  useEffect(() => {
    // 서버로 POST 요청을 보내고 데이터를 받아옴
    const fetchData = async () => {
      try {
        const response = await fetch('http://localhost:8080/analysis/re-alignment', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
        });

        if (!response.ok) {
          throw new Error('Network response was not ok');
        }

        const data = await response.json(); // JSON 데이터로 변환
        setResponseData(data); // 데이터를 상태로 저장
        setIsLoading(false); // 로딩 상태 해제
      } catch (error) {
        console.error('Error fetching data:', error);
        setIsLoading(false); // 오류 발생 시에도 로딩 상태 해제
      }
    };

    fetchData(); // 데이터를 가져오는 함수 호출
  }, []); // 컴포넌트가 마운트될 때 한 번만 실행

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