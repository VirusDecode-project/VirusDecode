import "./App.css";
import "bootstrap/dist/css/bootstrap.min.css";
import logo from "./image/logo.png";
import { Routes, Route, useNavigate, useLocation } from "react-router-dom";
import InputSeq from "./pages/inputSeq.js";
import Analysis from "./pages/analysis.js";
import { React, useState, useEffect } from "react";
import historyIcon from "./image/history.png";
import editIcon from "./image/edit.png";

function App() {
  let navigate = useNavigate();
  let location = useLocation();

  const [show, setShow] = useState(false);
  const [isHome, setIsHome] = useState(true);
  const [username, setUsername] = useState(null); // 사용자 이름 상태 추가

  useEffect(() => {
    if (location.pathname === "/") {
      setIsHome(true);
      setShow(false);
    } else {
      setIsHome(false);
      setShow(true);
    }
  }, [location.pathname]);

  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);

  const handleNavigate = (path) => {
    navigate(path);
  };

  return (
    <div className="App">
      <div className={`sidebar ${show ? "show" : ""}`}>
        <div className="sidebar-header">
          <img
            className="history-icon"
            src={historyIcon}
            onClick={handleClose}
            style={{ cursor: "pointer" }}
          />
          <img src={editIcon} alt="Edit" className="edit-icon" />
        </div>
        <div className="sidebar-body">
          <div>Today</div>
          <div>Reference1</div>
          <div>Reference2</div>
          <div>Reference3</div>
          <br />
          <div>Yesterday</div>
          <div>Reference1</div>
          <div>Reference2</div>
          <div>Reference3</div>
          <div>Reference4</div>
          <div>Reference5</div>
          <div>Reference6</div>
          <div>Reference7</div>

          <br />
          <div>Previous 7days</div>
          <div>Reference1</div>
          <div>Reference2</div>
          <div>Reference3</div>
          <div>Reference4</div>
          <div>Reference5</div>
          <div>Reference6</div>
          <div>Reference7</div>
          <div>Reference8</div>
          <div>Reference9</div>
          <div>Referencea</div>
          <div>Referenceb</div>
          <div>Referencec</div>
        </div>
      </div>

      <div className={`content-container ${show ? "shrink" : ""}`}>
        <div className="header-bar">
        <div className="header-left">
          {!show && !isHome && (
            <>
              <img
                onClick={handleShow}
                style={{ cursor: "pointer" }}
                src={historyIcon}
                alt="History"
                className="history-icon"
              />
              <img src={editIcon} alt="Edit" className="edit-icon" />
            </>
          )}
          <span
            className="logo-text"
            onClick={() => {
              handleNavigate("/");
            }}
            style={{ cursor: "pointer" }}
          >
            {isHome && (
              <img
                alt="VirusDecode Logo"
                src={logo}
                width="30"
                height="30"
                className="d-inline-block align-top"
                style={{ marginRight: "10px" }}
              />
            )}
            VirusDecode
          </span>
          </div>
          {username && <span className="username-display">{username}</span>}
        </div>
        <Routes>
          <Route
            path="/"
            element={
              <div>
                <div className="main-bg"></div>
                <div className="text-box">
                  <p style={{ fontSize: "20px" }}>
                    Decode the virus’s genetic code,
                    <br />
                    analyze its mutations,
                    <br />
                    and determine the vaccine sequence.
                  </p>
                </div>
                <button
                  className="image-button"
                  onClick={() => {
                    handleNavigate("inputSeq");
                  }}
                ></button>
              </div>
            }
          />

          <Route
            path="/inputSeq"
            element={<InputSeq setUsername={setUsername} />}
          />
          <Route path="/analysis" element={<Analysis />} />
        </Routes>
      </div>
    </div>
  );
}

export default App;
