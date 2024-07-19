import { Nav } from 'react-bootstrap';
import React, { useState, useEffect } from 'react';
import loadingImage from './loading.png';
import './analysis.css';

function Analysis() {
    let [tab, setTab] = useState(0)
    let [isLoading, setIsLoading] = useState(true);
    let [loadingText, setLoadingText] = useState('Analyzing');

    useEffect(() => {
        const interval = setInterval(() => {
            setLoadingText((prev) => {
                if (prev === 'Analyzing...') return 'Analyzing';
                return prev + '.';
            });
        }, 500); // 500ms 간격으로 텍스트 업데이트

        setTimeout(() => {
            setIsLoading(false);
            clearInterval(interval); // 로딩 완료 시 인터벌 정리
        }, 5000); // 페이지가 처음 로딩될 때 3초 동안 로딩 이미지 표시

        return () => clearInterval(interval); // 컴포넌트 언마운트 시 인터벌 정리
    }, []);

    return (
        <div>
            {isLoading ? (
                <div className="loading-container">
                    <img src={loadingImage} alt="Loading..." className="loading-image" />
                    <div className="loading-text">{loadingText}</div>
                </div>
            ) : (
                <>
                    <Nav variant="tabs" defaultActiveKey="link0" className="justify-content-center">
                        <Nav.Item>
                            <Nav.Link eventKey="link0" onClick={() => setTab(0)}>Alignment</Nav.Link>
                        </Nav.Item>
                        <Nav.Item>
                            <Nav.Link eventKey="link1" onClick={() => setTab(1)}>Gene</Nav.Link>
                        </Nav.Item>
                        <Nav.Item>
                            <Nav.Link eventKey="link2" onClick={() => setTab(2)}>Mutation</Nav.Link>
                        </Nav.Item>
                        <Nav.Item>
                            <Nav.Link eventKey="link3" onClick={() => setTab(3)}>Pathogenecity</Nav.Link>
                        </Nav.Item>
                        <Nav.Item>
                            <Nav.Link eventKey="link4" onClick={() => setTab(4)}>Analyze</Nav.Link>
                        </Nav.Item>
                    </Nav>
                    <Tab tab={tab} />
                </>
            )}
        </div>

    );
}

function Tab(props) {
    if (props.tab == 0) {
        return <div>Alignment을 분석한 결과에 대한 내용입니다.<br />(오픈예정)</div>
    } else if (props.tab == 1) {
        return <div>Gene에 관한 내용을 나타냅니다.<br />(예정)</div>
    } else if (props.tab == 2) {
        return <div>Mutation에 관한 정보<br />(베타서비스중)</div>
    } else if (props.tab == 3) {
        return <div>Pathogenecity를 분석한 내용입니다.<br />(추후오픈예정)</div>
    } else if (props.tab == 4) {
        return <div>Analyze 결과 입니다.<br />(유료서비스)</div>
    }
}

export default Analysis;