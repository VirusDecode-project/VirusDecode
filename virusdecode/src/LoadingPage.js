import React from 'react';
import './LoadingPage.css'; 
import loadingimage from './loadingimage.png';

const LoadingPage = () => {
  return (
    <div className="loading-container">
      <img className="loading-icon" src={loadingimage} alt="Loading icon" />
      <div className="loading-text">Analyzing...</div>
    </div>
  );
};

export default LoadingPage;
