// SeqInput 페이지 위에 띄우는 방식으로 수정예정

import React from 'react';
import './LoadingPage.css'; 


const LoadingPage = () => {
  return (
    <div className="loading-container">
      <img className="loading-icon" src="/loadingimage.png" alt="Loading icon" />
      <div className="loading-text">Analyzing...</div>
    </div>
  );
};

export default LoadingPage;