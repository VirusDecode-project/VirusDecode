import React from 'react';
import { Link } from 'react-router-dom';
import './Sidebar.css'; 
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faBars, faPen } from '@fortawesome/free-solid-svg-icons';

const Sidebar = ({ isOpen, toggleSidebar }) => {
  return (
    <div className={`sidebar ${isOpen ? 'open' : 'closed'}`}> {/* 사이드바열림 닫힘 클래스 구분 */}
      <div className='side-top'>
        {isOpen && ( 
          <button className="toggle-button-opened" onClick={toggleSidebar}>
            <FontAwesomeIcon icon={faBars} /> 
          </button>
        )}
        <button className="new-input">
          <Link to="/" className="new-input">
            <FontAwesomeIcon icon={faPen} /> {/* 새 시퀀스 입력 버튼 */}
          </Link>
        </button>
      </div>
      <div className="sidebar-content"> 
        <ul>
          <li class="day_after"><p/>yesterday</li> {/* 날짜 */}
          <li><Link to="/page0">Reference1</Link></li> {/* 분석결과 링크 연결 */}
          <li><Link to="/page1">Reference2</Link></li> 
          <li><Link to="/page2">Reference3</Link></li> 
          <li><p/>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li>
          <li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li>
          <li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li><li>aa</li>
          {/* 데이터파일에서 정보 가져오면 기록 리스트 추가하는 로직 만들기  */}
        </ul>
      </div>
    </div>
  );
};

export default Sidebar;
