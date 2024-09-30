import React, { useState, Dispatch, SetStateAction } from "react";
import { useNavigate } from "react-router-dom";
import "../styles/Login.css"

interface LoginProps {
  history: string[],
  setHistory: Dispatch<SetStateAction<string[]>>;
  setShow: Dispatch<SetStateAction<boolean>>;
  setMRNAReceived: Dispatch<SetStateAction<boolean>>;
  setPDBReceived: Dispatch<SetStateAction<boolean>>;
}

const Login: React.FC<LoginProps> = ({history, setHistory, setShow, setMRNAReceived, setPDBReceived}) => {
  let navigate = useNavigate();
  const [loginId, setLoginId] = useState<string | null>(null);
  const [password, setPassword] = useState<string | null>(null);

  const handleLoginSubmit = async (event: React.FormEvent) => {
    event.preventDefault();
    try {
      const loginResponse = await fetch(`${process.env.REACT_APP_BACKEND_URL}/auth/login`, {
        method: "POST",
        credentials: 'include',
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          loginId,
          password,
        }),
      });

      if (loginResponse.ok) {
        const fetchHistory = async () => {
          try {
            const historyResponse = await fetch(`${process.env.REACT_APP_BACKEND_URL}/history/list`, {
              method: 'GET',
              credentials: 'include',
            });
            if (!historyResponse.ok) {
              throw new Error("Failed to fetch history list");
            }
            const responseData = await historyResponse.json();
            console.log(responseData);
            setHistory(responseData);
          } catch (error) {
            console.error("Error fetching history:", error);
          }
        };
        fetchHistory();
        setShow(true); // Make sure the sidebar shows after navigation
        setMRNAReceived(false);
        setPDBReceived(false);
        navigate("/inputSeq");
      } else {
        alert("잘못된 ID나 비밀번호입니다.");
      }
    } catch (error) {
      console.error("로그인 요청 중 에러 발생:", error);
      alert("로그인 중 문제가 발생했습니다. 다시 시도해주세요.");
    }
  };
  
  return (
    <div className="login-wrapper">
      <p className="loginTitle">Login</p>
      <form className="loginForm" onSubmit={handleLoginSubmit}>
        <div className="input-container">
          <input 
            className="loginInput"
            type="text" 
            name="loginId" 
            placeholder=" "
            onChange={(e) => setLoginId(e.target.value)}
          />
          <label className="loginLabel" htmlFor="id">ID</label>
        </div>
        <div className="input-container">
          <input 
            className="loginInput"
            type="password" 
            name="password" 
            placeholder=" "
            onChange={(e) => setPassword(e.target.value)}
          />
          <label className="loginLabel" htmlFor="password">Password</label>
        </div>
        <button className="loginBtn" type="submit">
          Login
        </button>
      </form>
      <div className="linkBtns">
        <button className="gotoSignupBtn" onClick={() => navigate("/signup")}>
          Signup
        </button>
        <button className="forgotPasswordBtn" onClick={() => navigate("/signup")}>
          Forgot password?
        </button>
      </div>
    </div>
  );
}

export default Login;