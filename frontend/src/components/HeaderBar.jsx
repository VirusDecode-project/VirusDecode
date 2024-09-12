// HeaderBar.jsx
import React from 'react';
import historyIcon from '../assets/history.png';
import editIcon from '../assets/edit.png';
import logo from '../assets/logo.png';

const HeaderBar = ({ show, isHome, handleShow, handleEditClick, navigate }) => {
  const handleRestart = async() => {
    try {
      const serverResponse = await fetch("http://localhost:8080/history/deleteData");

      if (!serverResponse.ok) {
        const errorMessage = await serverResponse.text();
        throw new Error(errorMessage);
      }

      await serverResponse.text();
    } catch (error) {
      console.error("An error occurred while fetching history details: ", error.message);
    }
    navigate('/');
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
    </div>
  );
};

export default HeaderBar;
