// Home.jsx
import React from 'react';
import { useNavigate } from 'react-router-dom';
import '../styles/Home.css';

const Home = ({ setHistory, setShow, setMRNAReceived, setPDBReceived }) => {
  let navigate = useNavigate();

  return (
    <div>
      <div className="main-bg"></div>
      <div className="text-box">
        <p>
          Decode the virus’s genetic code, analyze its mutations, and determine the vaccine sequence.
        </p>
      </div>
      <button
        className="decode-button"
        onClick={() => {
          const fetchHistory = async () => {
            try {
              const serverResponse = await fetch("http://localhost:8080/history/list");
              if (!serverResponse.ok) {
                throw new Error("Failed to fetch history list");
              }
              const responseData = await serverResponse.json();
              setHistory(responseData);
            } catch (error) {
              console.error("Error fetching history:", error);
            }
          };

          fetchHistory();
          navigate("inputSeq");
          setShow(true); // Make sure the sidebar shows after navigation
          setMRNAReceived(false);
          setPDBReceived(false);
        }}
      >
        Try Decoding
      </button>
    </div>
  );
};

export default Home;
