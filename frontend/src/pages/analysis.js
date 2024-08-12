import React, { useState, useEffect } from 'react';
import { Nav } from 'react-bootstrap';
import { useNavigate } from 'react-router-dom';
import loadingImage from '../image/loading.png';
import './analysis.css';
import Alignment from './Alignment';

function Analysis() {
  let [tab, setTab] = useState(0); // 현재 선택된 탭 상태
  let [isLoading, setIsLoading] = useState(true); // 로딩 상태 관리
  let [loadingText, setLoadingText] = useState('Analyzing'); // 로딩 텍스트 상태
  let navigate = useNavigate(); // 페이지 이동을 위한 useNavigate 훅
  const handleClose = () => setShow(false); // 모달 닫기 
  const handleShow = () => setShow(true); // 모달 열기 
  const [show, setShow] = useState(false); // 모달의 표시 상태 

  useEffect(() => {
    // 로딩 창 점(...) 나타내는 부분 
    const interval = setInterval(() => {
      setLoadingText((prev) => {
        if (prev === 'Analyzing...') return 'Analyzing';
        return prev + '.';
      });
    }, 500);

    // 3초 로딩 지연
    setTimeout(() => {
      setIsLoading(false);
      clearInterval(interval); // 점 나타내기?
    }, 3000);

    
    return () => clearInterval(interval);
  }, []);

  return (
    <div>
      {isLoading ? (
        // 로딩 화면
        <div className="loading-container">
          <img src={loadingImage} alt="Loading..." className="loading-image" />
          <div className="loading-text">{loadingText}</div>
        </div>
      ) : (
        <div className={`analysis-container ${show ? 'shrink' : ''}`}>
          <>
            <Nav variant="tabs" defaultActiveKey="link0" className="justify-content-center">
              <Nav.Item>
                <Nav.Link eventKey="link0" onClick={() => setTab(0)}>Alignment</Nav.Link>
              </Nav.Item>
              <Nav.Item>
                <Nav.Link eventKey="link1" onClick={() => setTab(1)}>mRNA design</Nav.Link>
              </Nav.Item>
              <Nav.Item>
                <Nav.Link eventKey="link2" onClick={() => setTab(2)}>3D viewer</Nav.Link>
              </Nav.Item>
            </Nav>
            <Tab tab={tab} /> {/* 선택된 탭에 따라 다른 내용을 표시 */}
          </>
        </div>
      )}
    </div>
  );
}

function Tab(props) {
  // alignment 화면 스택바 예시 데이터
  const [data, setData] = useState([
    { label: 'ORF1ab', value: 50, color: 'rgba(144, 238, 144, 0.6)' },
    { label: 'S', value: 20, color: 'rgba(255, 99, 132, 0.6)' },
    { label: 'ORF3a', value: 10, color: 'rgba(255, 182, 193, 0.6)' },
    { label: 'E', value: 5, color: 'rgba(255, 255, 102, 0.6)' },
    { label: 'M', value: 10, color: 'rgba(135, 206, 235, 0.6)' },
    { label: 'ORF7b', value: 5, color: 'rgba(255, 160, 122, 0.6)' },
  ]);

  return (
    <div>
      {/* 탭 별 화면 연결하는 부분*/}
      {props.tab === 0 && <Alignment data={data} />}
      {props.tab === 1 && <div>mRNA design에 관한 내용을 나타냅니다.<br />(예정)</div>}
      {props.tab === 2 && <div>3D viewer가 실행되는 페이지 입니다.<br />(베타서비스)</div>}
    </div>
  );
}

export default Analysis;
