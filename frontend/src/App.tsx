import React, { useState, useEffect } from "react";
import { Routes, Route, useNavigate, useLocation } from "react-router-dom";
import "./styles/App.css";
import "bootstrap/dist/css/bootstrap.min.css";
import Home from "./pages/Home";
import Login from "./pages/Login";
import Signup from "./pages/Signup";
import InputSeq from "./pages/InputSeq";
import Analysis from "./pages/Analysis";
import Sidebar from "./components/Sidebar";
import HeaderBar from "./components/HeaderBar";
import { MRNAData, AlignmentData, PDBResponse } from './components/types';

function App() {
  let navigate = useNavigate();
  let location = useLocation();
  const [show, setShow] = useState(false);
  const [isHome, setIsHome] = useState(true);
  const [history, setHistory] = useState<string[]>([]);
  const [mRNAReceived, setMRNAReceived] = useState(false);
  const [PDBReceived, setPDBReceived] = useState(false);
  const [tab, setTab] = useState(0);
  const [workingHistory, setWorkingHistory] = useState("");
  const [linearDesignData, setLinearDesignData] = useState<MRNAData | null>(null);
  const [PDBids, setPDBids] = useState<string[]>([]);
  const [PDBInfo, setPDBInfo] = useState<string[]>([]);
  const [selectedPDBid, setSelectedPDBid] = useState("");
  const [alignmentData, setAlignmentData] = useState<AlignmentData>({
    aligned_sequences: {},
    alignment_index: {},
  });
  

  useEffect(() => {
    if (location.pathname === "/" || location.pathname === "/login" || location.pathname === "/signup") {
      setIsHome(true);
      setShow(false);
    } else {
      setIsHome(false);
      setShow(true);
    }
  }, [location.pathname]);



  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);
  const handleEditClick = () => { return; };


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
        setLinearDesignData={setLinearDesignData}
        setPDBids={setPDBids}
        setPDBInfo={setPDBInfo}
        setSelectedPDBid={setSelectedPDBid}
        setAlignmentData={setAlignmentData}
      />

      <div className={`content-container ${show ? "shrink" : ""}`}>
        <HeaderBar
          show={show}
          isHome={isHome}
          handleShow={handleShow}
          handleEditClick={handleEditClick}
          navigate={navigate}
        />
        <Routes>
          <Route
            path="/"
            element={<Home/>}
          />
          <Route
            path="/login"
            element={<Login
              history={history}
              setHistory={setHistory} 
              setShow={setShow} 
              setMRNAReceived={setMRNAReceived} 
              setPDBReceived={setPDBReceived}
            />}
          />
          <Route
            path="/signup"
            element={<Signup/>}
          />
          <Route
            path="/inputSeq"
            element={<InputSeq
              setTab={setTab}
              setWorkingHistory={setWorkingHistory}
              workingHistory={workingHistory}
              setMRNAReceived={setMRNAReceived}
              setPDBReceived={setPDBReceived}
              setAlignmentData={setAlignmentData}
              setHistory={setHistory}
            />}
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
              setWorkingHistory={setWorkingHistory}
              linearDesignData={linearDesignData}
              setLinearDesignData={setLinearDesignData}
              PDBids={PDBids}
              setPDBids={setPDBids}
              PDBInfo={PDBInfo}
              setPDBInfo={setPDBInfo}
              selectedPDBid={selectedPDBid}
              setSelectedPDBid={setSelectedPDBid}
              alignmentData={alignmentData}
              setHistory={setHistory}
              setAlignmentData={setAlignmentData}
            />}
          />
        </Routes>
      </div>
    </div>
  );
}

export default App;
