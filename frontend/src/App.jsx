import "./App.css";
import "bootstrap/dist/css/bootstrap.min.css";
import logo from "./image/logo.png";
import { Routes, Route, useNavigate, useLocation } from "react-router-dom";
import InputSeq from "./pages/inputSeq";
import Analysis from "./pages/analysis";
import { React, useState, useEffect } from "react";
import historyIcon from "./image/history.png";
import editIcon from "./image/edit.png";

function App() {

  let navigate = useNavigate();
  let location = useLocation();

  const [show, setShow] = useState(false);
  const [isHome, setIsHome] = useState(true);
  const [username, setUsername] = useState(null); // User name state
  const [history, setHistory] = useState([]); // State for history items
  const [activeHistoryItem, setActiveHistoryItem] = useState(null); // State for active history item

  useEffect(() => {
    if (location.pathname === "/") {
      setIsHome(true);
      setShow(false);
    } else {
      setIsHome(false);
      setShow(true);
    }
  }, [location.pathname]);
  // 백엔드에서 history 리스트 가져오기
  useEffect(() => {
    const fetchHistory = async () => {
      try {
        const response = await fetch("http://localhost:8080/history/list");
        if (!response.ok) {
          throw new Error('Failed to fetch history list');
        }
        const data = await response.json();
        setHistory(data); // 가져온 데이터를 history 상태로 설정
      } catch (error) {
        console.error("Error fetching history:", error);
      }
    };

    fetchHistory();
  }, []);
  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);

  const handleNavigate = (path) => {
    navigate(path);
  };

  const handleEditClick = async () => {
    // Prompt the user for a name
    const name = prompt("Enter a name for the new history item:");
    if (name) {
      // Backend 호출
      const requestData = { historyName: name };
      try {
        // GK - Loading 위치 이동
        const serverResponse = await fetch(
          "http://localhost:8080/history/create",
          {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
            },
            body: JSON.stringify(requestData),
          }
        );

        if (!serverResponse.ok) {
          const errorMessage = await serverResponse.text();
          throw new Error(errorMessage);
        }

        await serverResponse.text()
        setHistory((prevHistory) => [...prevHistory, name]); // Add the new history item
        navigate("/inputSeq"); // Navigate to the inputSeq page

      } catch (error) {
        console.error("An error occurred during the request: ", error.message);
      }
    }
  };
  const handleHistoryClick = async (index) => {
    // 선택한 history의 세부 정보를 가져오기 위한 함수
    const historyName = history[index]; // index를 사용하여 선택된 history 항목의 이름을 가져옴
    const requestData = {
      historyName: historyName
    };
    try {
      const serverResponse = await fetch(
        `http://localhost:8080/history/get`,
        {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify(requestData),
        }
      );

      if (!serverResponse.ok) {
        const errorMessage = await serverResponse.text();
        throw new Error(errorMessage);
      }

      const responseData = await serverResponse.json();
      navigate("/analysis", { state: { responseData } });

    } catch (error) {
      console.error("An error occurred while fetching history details: ", error.message);
    }
  };

  const handleEllipsisClick = (index) => {
    // Toggle active history item to show/hide options
    setActiveHistoryItem(activeHistoryItem === index ? null : index);
  };

  const handleRename = async (index) => {
    const newName = prompt("Enter a new name for the history item:");
    if (newName) {

      const historyName = history[index]; // 선택한 history의 현재 이름 가져오기
      const requestData = {
        historyName: historyName, // 기존 이름을 포함
        newName: newName, // 새 이름
      };

      try {
        // GK - Loading 위치 이동
        const serverResponse = await fetch(
          "http://localhost:8080/history/rename",
          {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
            },
            body: JSON.stringify(requestData),
          }
        );

        if (!serverResponse.ok) {
          const errorMessage = await serverResponse.text();
          throw new Error(errorMessage);
        }

        await serverResponse.text()
        setHistory((prevHistory) =>
          prevHistory.map((item, i) => (i === index ? newName : item))
        );

      } catch (error) {
        console.error("An error occurred during the request: ", error.message);
      }
    }
    setActiveHistoryItem(null); // Hide options after renaming
  };

  const handleDelete = async (index) => {
    if (window.confirm("Are you sure you want to delete this history item?")) {
      const currentName = history[index]; // 선택한 history의 현재 이름 가져오기
      const requestData = { historyName: currentName };

      try {
        const serverResponse = await fetch(
          "http://localhost:8080/history/delete",
          {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
            },
            body: JSON.stringify(requestData), // 요청 데이터로 JSON 변환
          }
        );

        if (!serverResponse.ok) {
          const errorMessage = await serverResponse.text();
          throw new Error(errorMessage);
        }

        await serverResponse.text()
        setHistory((prevHistory) =>
          prevHistory.filter((_, i) => i !== index)
        ); // Remove the deleted history item from state
      } catch (error) {
        console.error("An error occurred during the request: ", error.message);
      }
    }
    setActiveHistoryItem(null); // Hide options after deleting
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
          <img
            src={editIcon}
            alt="Edit"
            className="edit-icon"
            onClick={handleEditClick} // Add onClick handler
            style={{ cursor: "pointer" }}
          />
        </div>
        <div className="sidebar-body">
          {/* Render history items */}
          <div className="history-list">
            {history.map((item, index) => (
              <div key={index} className="history-item">
                <span onClick={() => handleHistoryClick(index)} style={{ cursor: "pointer" }}>{item}</span> {/* index를 전달하도록 수정 */}
                <button
                  className="ellipsis-button"
                  onClick={() => handleEllipsisClick(index)}
                >
                  ...
                </button>
                {activeHistoryItem === index && (
                  <div className="options-menu">
                    <button
                      className="option-button"
                      onClick={() => handleRename(index)}
                    >
                      Rename
                    </button>
                    <button
                      className="option-button"
                      onClick={() => handleDelete(index)}
                    >
                      Delete
                    </button>
                  </div>
                )}
              </div>
            ))}
          </div>
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
                <img
                  src={editIcon}
                  alt="Edit"
                  className="edit-icon"
                  onClick={handleEditClick} // Add onClick handler here too
                  style={{ cursor: "pointer" }}
                />
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
