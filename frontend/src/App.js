import './App.css';
import 'bootstrap/dist/css/bootstrap.min.css';
import { Navbar, Container } from 'react-bootstrap';
import logo from './logo.png';//이거 퍼블릭 폴더에 넣기
import { Routes, Route, useNavigate } from 'react-router-dom'

import InputSeq from './pages/inputSeq.js'
import Analysis from './pages/analysis.js'
import LoginButtonComponent from './loginButtonComponent';
import LoginFailurePage from './pages/loginFailurePage';
import LoginSuccessPage from './pages/loginSuccessPage';
import React from 'react';

// const InputSeq = React.lazy(() => import('./pages/inputSeq'));
// const Analysis = React.lazy(() => import('./pages/analysis'));

function App() {

  let navigate = useNavigate();

  return (
    <div className="App">



      <Routes>
        <Route path='/' element={
          <div>
            <Navbar bg="white" data-bs-theme="white">
              <Container className="custom-container">
                <Navbar.Brand style={{ cursor: 'pointer' }} onClick={() => {
                  navigate('/')
                }}>
                  <img
                    alt="VirusDecode Logo"
                    src={logo}
                    width="30"
                    height="30"
                    className="d-inline-block align-top"
                  />{' '}
                  <span className="ibm-plex-mono-regular">VirusDecode</span>
                </Navbar.Brand>
              </Container>
            </Navbar>
            <div className='main-bg'></div>
            <div className='text-box'>
              <p>Decode the virus’s genetic code,<br />
                analyze its mutations,<br />
                and determine the vaccine sequence.</p>
            </div>
            <LoginButtonComponent />      {/* login Button js파일로 따로 구현 */}
          </div>} />

        <Route path="/login-success" element={<LoginSuccessPage />} />   {/* server 로부터 path 값 전달받음 (성공) */}
        <Route path="/login-failure" element={<LoginFailurePage />} />  {/* server 로부터 path 값 전달받음 (실패) */}
        {/* <Route path="/inputSeq" element={<InputSeq />} />
        <Route path="/analysis" element={<Analysis />} /> */}

      </Routes>



    </div>
  );
}

export default App;
