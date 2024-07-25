import React from 'react';
import { Link } from 'react-router-dom';
import './Sidebar.css';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faBars, faPen } from '@fortawesome/free-solid-svg-icons';

const Sidebar = ({ isOpen, toggleSidebar }) => {
  return (
    <div className={`sidebar ${isOpen ? 'open' : 'closed'}`}>
      <div className='side-top'>
        {isOpen && (
          <button className="toggle-button-opened" onClick={toggleSidebar}>
            <FontAwesomeIcon icon={faBars} />
          </button>
        )}
        <button className="new-input">
          <Link to="/app" className="new-input">
            <FontAwesomeIcon icon={faPen} />
          </Link>
        </button>
      </div>
      <div className="sidebar-content">
        <ul>
          <li className="day_after"><p/>yesterday</li>
          <li><Link to="/app/page0">Reference1</Link></li>
          <li><Link to="/app/page1">Reference2</Link></li>
          <li><Link to="/app/page2">Reference3</Link></li>
        </ul>
      </div>
    </div>
  );
};

export default Sidebar;
