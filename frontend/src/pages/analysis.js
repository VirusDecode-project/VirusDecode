import React, { useState, useEffect } from 'react';
import { Nav } from 'react-bootstrap';
import { useNavigate } from 'react-router-dom';
import loadingImage from './loading.png';
import './analysis.css';
import Alignment from './Alignment';

function Analysis() {
  let [tab, setTab] = useState(0);
  let [isLoading, setIsLoading] = useState(true);
  let [loadingText, setLoadingText] = useState('Analyzing');
  let navigate = useNavigate();
  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);
  const [show, setShow] = useState(false);

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

  return (
    <div>
      {isLoading ? (
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
            <Tab tab={tab} />
          </>
        </div>
      )}
    </div>
  );
}

function Tab(props) {
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
      {props.tab === 0 && <Alignment data={data} />}
      {props.tab === 1 && <div>mRNA design에 관한 내용을 나타냅니다.<br />(예정)</div>}
      {props.tab === 2 && <div>3D viewer가 실행되는 페이지 입니다.<br />(베타서비스)</div>}
      
    </div>
  );
}

export default Analysis;
