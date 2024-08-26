// Loading.js
import React, { useState, useEffect } from 'react';
import loadingImage from '../image/loading.png'; // Adjust the path as necessary
import './Loading.css'; // Optional: You can add styles for the loading component

const Loading = ({text}) => {
  const [loadingText, setLoadingText] = useState(text);

  useEffect(() => {
    let intervalId = setInterval(() => {
      setLoadingText((prev) => {
        if (prev.endsWith('...')) return loadingText;
        return prev + '.';
      });
    }, 500);

    return () => {
      clearInterval(intervalId);
    };
  }, [loadingText]);

  return (
    <div className="loading-container">
      <img src={loadingImage} alt="Loading..." className="loading-image" />
      <div className="loading-text">{loadingText}</div>
    </div>
  );
};

export default Loading;
