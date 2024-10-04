import React, { useCallback, useEffect, useRef, Dispatch, SetStateAction } from 'react';
import { useNavigate } from "react-router-dom";
import '../styles/Modal.css';
import { useRecoilState } from "recoil";
import { authState } from '../state/authState';

interface LoginModalProps {
  isOpen: boolean;
  onClose: () => void;
  setUserName:Dispatch<SetStateAction<string | null>>;
}

const LoginModal: React.FC<LoginModalProps> = ({ isOpen, onClose, setUserName }) => {
  const navigate = useNavigate();
  const [isLoggedIn, setIsLoggedIn] = useRecoilState(authState);
  const loginId = 'virusdecode';
  const password = 'virusdecode';
  const modalRef = useRef<HTMLDivElement>(null); // 모달을 참조하기 위한 useRef
  
  const login = useCallback(async () => {
    try {
      const loginResponse = await fetch(`${process.env.REACT_APP_BACKEND_URL}/api/auth/login`, {
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

      if (loginResponse.ok) {const fetchName = async () => {
        try {
          const nameResponse = await fetch(`${process.env.REACT_APP_BACKEND_URL}/api/auth/userinfo`, { 
            method: "POST",
            credentials: 'include',
            headers: {
              "Content-Type": "application/json",
            },
          });
      
          if (nameResponse.ok) {
            const responseData = await nameResponse.text();
            setUserName(responseData);
          }else{
            const errorMessage = await nameResponse.text();
            throw new Error(errorMessage);
          }
        } catch (error) {
          // if(error instanceof Error){
              // window.alert(error.message);
            // }
            setUserName(null);
          }
        };
        fetchName();
        setIsLoggedIn(true);
        navigate("/inputSeq");
      } else {
        const errorMessage = await loginResponse.text();
        throw new Error(errorMessage);
      }
    } catch (error) {
      if(error instanceof Error){
        window.alert(error.message);
      }
    }
  },[setIsLoggedIn]);

  const handleStayLoggedOut = () => {
    onClose();
    login();
  };

  useEffect(() => {
    const handleEsc = (event: KeyboardEvent) => {
        if (event.key === 'Escape') {
            onClose();
        }
    };

    const handleClickOutside = (event: MouseEvent) => {
        if (modalRef.current && !modalRef.current.contains(event.target as Node)) {
            onClose();
        }
    };

    if (isOpen) {
        document.addEventListener('keydown', handleEsc);
        document.addEventListener('mousedown', handleClickOutside); 
    }

    return () => {
        document.removeEventListener('keydown', handleEsc);
        document.removeEventListener('mousedown', handleClickOutside);
    };
}, [isOpen, onClose]);

  if (!isOpen) {
    return null; 
  }

  return (
      <div className="login-modal-overlay">
          <div className="login-modal" ref={modalRef}>
              <button className="close-button" onClick={onClose}>×</button>
              <p className='welcomeTitle'>Welcome to VirusDecode</p>
              <p className='welcomeText'>Log in to get your <br/> virus analysis records.</p>
              <button className='loginBtn_modal' onClick={() => navigate("/login")}>Log in</button>
              <button className='signupBtn_modal' onClick={() => navigate("/signup")}>Sign up</button>
              <button className='stayLoggedOutBtn' onClick={(handleStayLoggedOut)}>stay logged out</button>
          </div>
      </div>
  );
}

export default LoginModal;