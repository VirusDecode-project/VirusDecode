import React, { useState, useEffect } from 'react';
import { NavLink, Route, Routes, Navigate, useLocation } from 'react-router-dom';
import './Tabs.css';
import Sequence from './tab_pages/Sequence';
import Gene from './tab_pages/Gene';
import Mutation from './tab_pages/Mutation';
import Pathogenicity from './tab_pages/Pathogenicity';
import Analyze from './tab_pages/Analyze';

// 탭 데이터 정의
const tabData = [
  { path: '/Tabs/Sequence', label: 'Sequence', Component: Sequence },
  { path: '/Tabs/Gene', label: 'Gene', Component: Gene },
  { path: '/Tabs/Mutation', label: 'Mutation', Component: Mutation },
  { path: '/Tabs/Pathogenicity', label: 'Pathogenicity', Component: Pathogenicity },
  { path: '/Tabs/Analyze', label: 'Analyze', Component: Analyze }
];

function Tabs() {
  const location = useLocation(); // 현재 위치 저장
  //기본 경로 '/Tabs/Sequence'로 설정
  const initialActiveTab = location.pathname === "/" ? "/Tabs/Sequence" : location.pathname;
  const [activeTab, setActiveTab] = useState(initialActiveTab); // 활성 탭 상태 기록

  // location.pathname이 변경될 때마다 activeTab 업데이트
  useEffect(() => {
    setActiveTab(location.pathname);
  }, [location.pathname]);

  return (
    <div>
      <nav className="tabs">
        {tabData.map(({ path, label }) => (
          <NavLink
            key={path} // 각 NavLink에 고유한 키 부여
            to={path} // 경로 설정
            // 활성화여부에 따라 'tab active' 또는 'tab'  적용
            className={({ isActive }) => (isActive ? 'tab active' : 'tab')}
            onClick={() => setActiveTab(path)} // 탭 클릭 시 activeTab 업데이트
          >
            {label}
          </NavLink>
        ))}
      </nav>
      <Routes>
        {/* 기본 경로로 이동 시 '/Tabs/Sequence'로 돌아옴 */}
        <Route path="/" element={<Navigate to="/Tabs/Sequence" />} />
        {tabData.map(({ path, Component }) => (
          // 각 경로 페이지 렌더링
          <Route key={path} path={path} element={<Component />} />
        ))}
      </Routes>
    </div>
  );
};

export default Tabs;
