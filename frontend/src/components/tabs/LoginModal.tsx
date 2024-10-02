import React, { useEffect, useRef } from 'react';
import { useNavigate } from "react-router-dom";
import '../../styles/Modal.css'

interface LoginModalProps {
    isOpen: boolean;
    onClose: () => void;
}

const LoginModal: React.FC<LoginModalProps> = ({ isOpen, onClose }) => {
  let navigate = useNavigate();
  const modalRef = useRef<HTMLDivElement>(null); // 모달을 참조하기 위한 useRef
  if (!isOpen) {
      return null; //여기에 test id 등록하는 코드 넣기
  }
  //useEffect(() => {
  //     const handleEsc = (event: KeyboardEvent) => {
  //         if (event.key === 'Escape') {
  //             onClose();
  //         }
  //     };

  //     const handleClickOutside = (event: MouseEvent) => {
  //         if (modalRef.current && !modalRef.current.contains(event.target as Node)) {
  //             onClose();
  //         }
  //     };

  //     if (isOpen) {
  //         document.addEventListener('keydown', handleEsc);
  //         document.addEventListener('mousedown', handleClickOutside); // 모달 바깥 클릭 감지
  //     }

  //     // Clean up 이벤트 리스너
  //     return () => {
  //         document.removeEventListener('keydown', handleEsc);
  //         document.removeEventListener('mousedown', handleClickOutside);
  //     };
  // }, [isOpen, onClose]);

  // if (!isOpen) {
  //     return null;
  // }

  return (
      <div className="login-modal-overlay">
          <div className="login-modal" ref={modalRef}>
              <p className='welcomeTitle'>Welcome to VirusDecode</p>
              <p className='welcomeText'>Log in to get your <br/> virus analysis records.</p>
              <button className='loginBtn_modal' onClick={() => navigate("/login")}>Log in</button>
              <button className='signupBtn_modal' onClick={() => navigate("/signup")}>Sign up</button>
              <button className='stayLoggedOutBtn' onClick={(e) => onClose()}>stay logged out</button>
          </div>
      </div>
  );
}

export default LoginModal;
