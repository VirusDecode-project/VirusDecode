import React from 'react';
import { useNavigate } from 'react-router-dom';
import '../styles/Home.css';

const Home: React.FC = () => {
  let navigate = useNavigate();
  return (
    <div>
      <div className="main-bg"></div>
      <div className="text-box">
        <p>
          Decode the virus’s genetic code, analyze its mutations, and determine the vaccine sequence.
        </p>
      </div>
      <button className="decode-button" onClick={() => navigate("/inputSeq")}>
        Try Decoding
      </button>
    </div>
  );
};

export default Home;
