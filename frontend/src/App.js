import './App.css';
import 'bootstrap/dist/css/bootstrap.min.css';
import { Navbar, Container } from 'react-bootstrap';
import logo from './logo.png';//이거 퍼블릭 폴더에 넣기
import { Routes, Route, Link, useNavigate, Outlet } from 'react-router-dom'

import InputSeq from './pages/inputSeq.js'
import Analysis from './pages/analysis.js'
import React from 'react';

// const InputSeq = React.lazy(() => import('./pages/inputSeq'));
// const Analysis = React.lazy(() => import('./pages/analysis'));

function App() {

  let navigate = useNavigate();

  return (
    <div className="App">

      <Navbar bg="white" data-bs-theme="white">
        <Container>
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

      <Routes>
        <Route path='/' element={<div>
          <div className='main-bg'></div>
          <div className='text-box'>
            Decode the virus’s genetic code,<br />
            analyze its mutations,<br />
            and determine the vaccine sequence.
          </div>
          <button className="image-button" onClick={() => { navigate('inputSeq') }}></button>
        </div>} />

        <Route path="/inputSeq" element={<InputSeq />}/>
        <Route path="/analysis" element={<Analysis />}/>

      </Routes>




    </div>
  );
}

export default App;
