import React, { useCallback, useEffect, useState, Dispatch, SetStateAction,  } from 'react';
import historyIcon from '../assets/history.png';
import editIcon from '../assets/edit.png';
import logo from '../assets/logo.png';
import { NavigateFunction } from 'react-router-dom';

interface HeaderBarProps {
  show: boolean;
  isHome: boolean;
  handleShow: () => void;
  handleEditClick: () => void;
  navigate: NavigateFunction;
  userName: string | null;
  setUserName:Dispatch<SetStateAction<string | null>>;
}

const HeaderBar: React.FC<HeaderBarProps> = ({ show, isHome, handleShow, handleEditClick, navigate, userName, setUserName })=> {
  // const [userName, setUserName] = useState<string | null>(null);
  const [logOutIsOpen, setLogOutIsOpen] = useState<boolean>(false);

  const handleRestart = () => {
    navigate('/');
  };

  const toggleSignOutIsOpen = () => {
    setLogOutIsOpen(!logOutIsOpen);
  }

  useEffect(() => {
    const fetchName = async () => {
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
  },[userName]);

  const handleLogout = async(event: React.MouseEvent<HTMLButtonElement>) => {
    event.preventDefault();
    try {
      const nameResponse = await fetch(`${process.env.REACT_APP_BACKEND_URL}/api/auth/logout`, {
        method: "POST",
        credentials: 'include',
        headers: {
          "Content-Type": "application/json",
        },
      });
  
      if (nameResponse.ok) {
        setUserName(null)
        navigate("/");
        toggleSignOutIsOpen();
      }else{
        const errorMessage = await nameResponse.text();
        throw new Error(errorMessage);
      }
    } catch (error) {
      if(error instanceof Error){
        window.alert(error.message);
      }
    }
  };

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
      <button className="username-display" onClick={toggleSignOutIsOpen}>{userName}</button>
      {logOutIsOpen &&  (
        <button className="logoutBtn" onClick={handleLogout}>logout</button>
      )}
    </div>
  );
};

export default HeaderBar;
