import React, { useState } from 'react';
import { Routes, Route } from 'react-router-dom';
import Sidebar from './components/Sidebar';
import NavigationBar from './components/NavBar';
import './App.css';
import SeqInput from'./SeqInput';
import Tabs from './Tabs';
import LoadingPage from './LoadingPage';


function App() {
  const [isSidebarOpen, setSidebarOpen] = useState(false);
  const toggleSidebar = () => {
    setSidebarOpen(!isSidebarOpen);
  };

  return (
      <div className="App">
        {/*사이드바가 여닫힘 state 갱신 */}
        <Sidebar isOpen={isSidebarOpen} toggleSidebar={toggleSidebar} />
        {/*최상단 네비게이션 바 출력 */}
        <div className={`content-container ${isSidebarOpen ? 'shifted' : ''}`}>
          <NavigationBar isSidebarOpen={isSidebarOpen} toggleSidebar={toggleSidebar} />
          {/* 페이지 내부 컨텐츠 출력 */}
          <div className="main-content">
            <Routes>
            <Route path="/" element={<SeqInput />} />
            <Route path="/Tabs/*" element={<Tabs />} />
            <Route path="/LoadingPage" element={<LoadingPage />} />
              {/*루트경로(/)에서 시퀀스 입력페이지 렌더링 */}  
              <Route path="/" element={<SeqInput />}/>
              {/*사이드 바 내부 요소 클릭 시 해당 레퍼런스 결과 창으로 이동 */}            
              <Route path="/page0" element={<><h1>Reference1</h1><p>Reference1 분석 결과</p></>} />
              <Route path="/page1" element={<><h1>Reference2</h1><p>Reference2 분석 결과</p></>} />
              <Route path="/page2" element={<><h1>Reference3</h1><p>Reference3 분석 결과</p></>} />
              { /*--Todo--받아오는 데이터의 로그 수만큼 route 추가하는 함수 만들기--데이터 파일..??*/ }
            </Routes>
          </div>
        </div>
      </div>
  );
}

export default App;

