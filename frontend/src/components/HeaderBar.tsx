import React, { useEffect, useState, Dispatch, SetStateAction, useRef } from 'react';
import { Row, Col } from "react-bootstrap";
import historyIcon from '../assets/history.png';
import editIcon from '../assets/edit.png';
import logo from '../assets/logo.png';
import { NavigateFunction } from 'react-router-dom';
import userIcon from '../assets/user.png';

interface HeaderBarProps {
  show: boolean;
  isHome: boolean;
  handleShow: () => void;
  handleEditClick: () => void;
  navigate: NavigateFunction;
  userName: string | null;
  setUserName: Dispatch<SetStateAction<string | null>>;
  handleError: (message: string) => void;
}

const HeaderBar: React.FC<HeaderBarProps> = ({ show, isHome, handleShow, handleEditClick, navigate, userName, setUserName, handleError }) => {
  const [loginId, setLoginId] = useState<string | null>(null);
  const [isUserInfoOpen, setIsUserInfoOpen] = useState<boolean>(false);
  const userInfoRef = useRef<HTMLDivElement | null>(null);
  const iconRef = useRef<HTMLImageElement | null>(null);

  const handleRestart = () => {
    navigate('/');
  };

  const toggleUserInfoOpen = () => {
    if (userName != null) {
      setIsUserInfoOpen(!isUserInfoOpen);
    }
  }

  useEffect(() => {
    const fetchName = async () => {
      try {
        const nameResponse = await fetch(`/api/auth/userinfo`, {
          method: "POST",
          credentials: 'include',
          headers: {
            "Content-Type": "application/json",
          },
        });

        if (nameResponse.ok) {
          const responseData = await nameResponse.json();
          setUserName(responseData.userName);
          setLoginId(responseData.loginId);
        } else {
          const errorMessage = await nameResponse.text();
          throw new Error(errorMessage);
        }
      } catch (error) {
        setUserName(null);
        setLoginId(null);
      }
    };
    fetchName();
  }, [userName]);

  const handleLogout = async (event: React.MouseEvent<HTMLButtonElement>) => {
    event.preventDefault();
    try {
      const nameResponse = await fetch(`/api/auth/logout`, {
        method: "POST",
        credentials: 'include',
        headers: {
          "Content-Type": "application/json",
        },
      });

      if (nameResponse.ok) {
        setUserName(null)
        navigate("/");
        setIsUserInfoOpen(false);
      } else {
        const errorMessage = await nameResponse.text();
        throw new Error(errorMessage);
      }
    } catch (error) {
      if (error instanceof Error) {
        // window.alert(error.message);
        handleError(error.message);
      }
    }
  };

  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (userInfoRef.current && !userInfoRef.current.contains(event.target as Node) &&
          iconRef.current && !iconRef.current.contains(event.target as Node)) 
      {setIsUserInfoOpen(false);}
    };

    if (isUserInfoOpen) {
      document.addEventListener('mousedown', handleClickOutside);
    } else {
      document.removeEventListener('mousedown', handleClickOutside);
    }

    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, [isUserInfoOpen]);

  return (
    <div className="header-bar">
      <div className="header-left">
        {!show && !isHome && (
          <>
            <img
              onClick={handleShow}
              style={{ cursor: 'pointer' }}
              src={historyIcon}
              alt="History"
              className="history-icon"
            />
            <img
              src={editIcon}
              alt="Edit"
              className="edit-icon"
              onClick={handleEditClick}
              style={{ cursor: 'pointer' }}
            />
          </>
        )}
        <span
          className="logo-text"
          onClick={handleRestart}
          style={{ cursor: 'pointer' }}
        >
          {isHome && (
            <img
              alt="VirusDecode Logo"
              src={logo}
              width="30"
              height="30"
              className="d-inline-block align-top"
              style={{ marginRight: '10px' }}
            />
          )}
          VirusDecode
        </span>
      </div>
      <div className="userInfo-wrapper" ref={userInfoRef}>
        <img
          src={userIcon}
          className="user-icon"
          style={{ cursor: "pointer" }}
          alt="User Icon"
          onClick={toggleUserInfoOpen}
          ref={iconRef}
        />
        {isUserInfoOpen && (
          <div className="userInfo-menu">
            <div className="userInfo-display">
              <Row>
                <Col xs="auto" className="text-center">
                  <div>Name:</div>
                  <div>ID:</div>
                </Col>
                <Col xs="auto">
                  <div>{userName}</div>
                  <div>{loginId}</div>
                </Col>
              </Row>
            </div>
            <div className="divider"></div>
            <button className="logoutBtn text-center" onClick={handleLogout}>logout</button>
          </div>
        )}
      </div>
    </div>
  );
};

export default HeaderBar;
