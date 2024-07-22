import './App.css';
import 'bootstrap/dist/css/bootstrap.min.css';
import { Navbar, Container } from 'react-bootstrap';
import logo from './logo.png';//이거 퍼블릭 폴더에 넣기
import { Routes, Route, useNavigate } from 'react-router-dom'

import InputSeq from './pages/inputSeq.js'
import Analysis from './pages/analysis.js'
import { React, useState } from 'react';
import historyIcon from './history.png';
import editIcon from './edit.png';

function App() {

  let navigate = useNavigate();

  const [show, setShow] = useState(false);
  const [isHome, setIsHome] = useState(true);



  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);

  const handleNavigate = (path) => {
    if (path === '/') {
      setIsHome(true);
      setShow(false); // 홈으로 갈 때는 사이드바를 닫습니다.
    } else {
      setIsHome(false);
      setShow(true); // 다른 페이지로 갈 때는 사이드바를 엽니다.
    }
    navigate(path);
  }

  return (
    <div className="App">

      <div className={`sidebar ${show ? 'show' : ''}`}>
        <div className="sidebar-header">

          <img className="history-icon" src={historyIcon} onClick={handleClose} style={{ cursor: 'pointer' }} />
          <img src={editIcon} alt="Edit" className="edit-icon" />
        </div>
        <div className="sidebar-body">
          <div>Today</div>
          <div>Reference1</div>
          <div>Reference2</div>
          <br />
          <div>Previous 7 days</div>
          <div>Reference1</div>
          <div>Reference2</div>
          <div>Reference3</div>

        </div>
      </div>


      <div className={`content-container ${show ? 'shrink' : ''}`}>
        <div className="header-bar">
          {!show && !isHome && (
            <>

              <img onClick={handleShow} style={{ cursor: 'pointer' }} src={historyIcon} alt="History" className="history-icon" />
              <img src={editIcon} alt="Edit" className="edit-icon" />
            </>
          )}
          <span className='logo-text' onClick={() => { handleNavigate('/') }} style={{ cursor: 'pointer' }}>
            {isHome &&
              <img
                alt="VirusDecode Logo"
                src={logo}
                width="30"
                height="30"
                className="d-inline-block align-top"
                style={{ marginRight: '10px' }}
              />}VirusDecode</span>
        </div>
        <Routes>
          <Route path='/' element={
            <div>
              <div className='main-bg'></div>
              <div className='text-box'>
                <p style={{ fontSize: '20px' }}>Decode the virus’s genetic code,<br />
                  analyze its mutations,<br />
                  and determine the vaccine sequence.</p>
              </div>
              <button className="image-button" onClick={() => { handleNavigate('inputSeq') }}></button>
            </div>} />

          <Route path="/inputSeq" element={<InputSeq />} />
          <Route path="/analysis" element={<Analysis />} />

        </Routes>
      </div>




    </div>
  );
}

export default App;
