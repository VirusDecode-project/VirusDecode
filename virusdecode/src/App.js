import React, { useState } from 'react';
import { Routes, Route } from 'react-router-dom';
import Sidebar from './components/Sidebar';
import NavigationBar from './components/NavBar';
import './App.css';
import SeqInput from './SeqInput';
import Tabs from './Tabs';
import LoadingPage from './LoadingPage';

function App() {
  const [isSidebarOpen, setSidebarOpen] = useState(false);
  const toggleSidebar = () => {
    setSidebarOpen(!isSidebarOpen);
  };

  return (
    <div className="App">
      <Sidebar isOpen={isSidebarOpen} toggleSidebar={toggleSidebar} />
      <div className={`content-container ${isSidebarOpen ? 'shifted' : ''}`}>
        <NavigationBar isSidebarOpen={isSidebarOpen} toggleSidebar={toggleSidebar} />
        <div className="main-content">
          <Routes>
            <Route path="/" element={<SeqInput />} />
            <Route path="Tabs/*" element={<Tabs />} />
            <Route path="LoadingPage" element={<LoadingPage />} />
            <Route path="page0" element={<><h1>Reference1</h1><p>Reference1 분석 결과</p></>} />
            <Route path="page1" element={<><h1>Reference2</h1><p>Reference2 분석 결과</p></>} />
            <Route path="page2" element={<><h1>Reference3</h1><p>Reference3 분석 결과</p></>} />
          </Routes>
        </div>
      </div>
    </div>
  );
}

export default App;
