// App.js
import React, { useState, useEffect, useRef } from "react";
import { Routes, Route, useNavigate, useLocation } from "react-router-dom";
import "./styles/App.css";
import "bootstrap/dist/css/bootstrap.min.css";
import Home from "./pages/Home";
import InputSeq from "./pages/InputSeq";
import Analysis from "./pages/Analysis";
import Sidebar from "./components/Sidebar";
import HeaderBar from "./components/HeaderBar";

function App() {
  let navigate = useNavigate();
  let location = useLocation();
  const [show, setShow] = useState(false);
  const [isHome, setIsHome] = useState(true);
  const [history, setHistory] = useState([]);
  const [mRNAReceived, setMRNAReceived] = useState(false);
  const [PDBReceived, setPDBReceived] = useState(false);
  const [tab, setTab] = useState(0);
  const [workingHistory, setWorkingHistory] = useState(null);

  useEffect(() => {
    if (location.pathname === "/") {
      setIsHome(true);
      setShow(false);
    } else {
      setIsHome(false);
      setShow(true);
    }
  }, [location.pathname]);

  useEffect(() => {
    const fetchHistory = async () => {
      try {
        const serverResponse = await fetch("http://localhost:8080/history/list");
        if (!serverResponse.ok) {
          throw new Error("Failed to fetch history list");
        }
        const responseData = await serverResponse.json();
        setHistory(responseData.reverse());
      } catch (error) {
        console.error("Error fetching history:", error);
      }
    };

    fetchHistory();
  }, [navigate]);

  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);

  return (
    <div className="App">
      <Sidebar
        show={show}
        handleClose={handleClose}
        history={history}
        setHistory={setHistory}
        navigate={navigate}
        setMRNAReceived={setMRNAReceived}
        setPDBReceived={setPDBReceived}
        setTab={setTab}
        workingHistory={workingHistory}
        setWorkingHistory={setWorkingHistory}
      />

      <div className={`content-container ${show ? "shrink" : ""}`}>
        <HeaderBar
          show={show}
          isHome={isHome}
          handleShow={handleShow}
          navigate={navigate}
        />
        <Routes>
          <Route
            path="/"
            element={<Home setHistory={setHistory} setShow={setShow} setMRNAReceived={setMRNAReceived} setPDBReceived={setPDBReceived} />}
          />
          <Route
            path="/inputSeq"
            element={<InputSeq setTab={setTab} setWorkingHistory={setWorkingHistory}/>}
          />
          <Route
            path="/analysis"
            element={<Analysis
              tab={tab}
              setTab={setTab}
              mRNAReceived={mRNAReceived}
              setMRNAReceived={setMRNAReceived}
              PDBReceived={PDBReceived}
              setPDBReceived={setPDBReceived}
              workingHistory={workingHistory}
            />}
          />
        </Routes>
      </div>
    </div>
  );
}

export default App;
