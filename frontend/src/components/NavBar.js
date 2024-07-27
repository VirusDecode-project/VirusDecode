import React from 'react';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faBars } from '@fortawesome/free-solid-svg-icons';
import { Link } from 'react-router-dom';
import './NavBar.css';

const NavigationBar = ({ isSidebarOpen, toggleSidebar }) => {
  return (
    <nav className={`navigation-bar ${isSidebarOpen ? 'shifted' : ''}`}>
      {!isSidebarOpen && (
        <button className="toggle-button" onClick={toggleSidebar}>
          <FontAwesomeIcon icon={faBars} />
        </button>
      )}
      <div>
        <Link to="/" className="nav-title">VirusDecode</Link>
      </div>
    </nav>
  );
};

export default NavigationBar;
