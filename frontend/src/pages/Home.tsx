import React, { useEffect, Dispatch, SetStateAction } from 'react';
import { useNavigate } from 'react-router-dom';
import '../styles/Home.css';

interface HomeProps {
  history: string[],
  setHistory: Dispatch<SetStateAction<string[]>>;
  setShow: Dispatch<SetStateAction<boolean>>;
  setMRNAReceived: Dispatch<SetStateAction<boolean>>;
  setPDBReceived: Dispatch<SetStateAction<boolean>>;
}

const Home: React.FC<HomeProps> = ({ history, setHistory, setShow, setMRNAReceived, setPDBReceived }) => {
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
              ///////////////test signUp
              const signUpData = { firstName:"John", lastName:"Doe", loginId: "johndoe123", password: "1234" };
              const signUpResponse = await fetch('http://localhost:8080/auth/signup', {
                method: 'POST',
                credentials: 'include',
                headers: {
                  'Content-Type': 'application/json',
                },
                body: JSON.stringify(signUpData),
              });

              if (!signUpResponse.ok) {
                const errorMessage = await signUpResponse.text();
                throw new Error(errorMessage);
              }
              await signUpResponse.text();

              ////////////////////////test
              ///////////////test Login
              const loginData = { loginId: "johndoe123", password: "1234" };
              const loginResponse = await fetch('http://localhost:8080/auth/login', {
                method: 'POST',
                credentials: 'include',
                headers: {
                  'Content-Type': 'application/json',
                },
                body: JSON.stringify(loginData),
              });

              if (!loginResponse.ok) {
                const errorMessage = await loginResponse.text();
                throw new Error(errorMessage);
              }
              await loginResponse.text();
              ////////////////////////test

              const serverResponse = await fetch("http://localhost:8080/history/list", {
                method: 'GET',
                credentials: 'include',
              });
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
          setShow(true); // Make sure the sidebar shows after navigation
          setMRNAReceived(false);
          setPDBReceived(false);
          navigate("inputSeq");
        }}
      >
        Try Decoding
      </button>
    </div>
  );
};

export default Home;
