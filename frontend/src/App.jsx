import "./App.css";
import "bootstrap/dist/css/bootstrap.min.css";
import logo from "./image/logo.png";
import { Routes, Route, useNavigate, useLocation } from "react-router-dom";
import InputSeq from "./pages/inputSeq";
import Analysis from "./pages/analysis";
import React, { useState, useEffect, useRef } from "react";
import historyIcon from "./image/history.png";
import editIcon from "./image/edit.png";

function App() {
  let navigate = useNavigate();
  let location = useLocation();

  const [show, setShow] = useState(false);
  const [isHome, setIsHome] = useState(true);
  const [username, setUsername] = useState(null);
  const [history, setHistory] = useState([]);
  const [activeHistoryItem, setActiveHistoryItem] = useState(null);
  const [menuPosition, setMenuPosition] = useState({ top: 0, left: 0 });

  const optionsMenuRef = useRef(null);

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
        const response = await fetch("http://localhost:8080/history/list");
        if (!response.ok) {
          throw new Error("Failed to fetch history list");
        }
        const data = await response.json();
        setHistory(data);
      } catch (error) {
        console.error("Error fetching history:", error);
      }
    };

    fetchHistory();
  }, []);

  useEffect(() => {
    const handleClickOutside = (event) => {
      if (
        optionsMenuRef.current &&
        !optionsMenuRef.current.contains(event.target)
      ) {
        setActiveHistoryItem(null); // Deactivate the history item
      }
    };

    document.addEventListener("mousedown", handleClickOutside);

    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, []);

  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);

  const handleNavigate = (path) => {
    navigate(path);
  };

  const handleEditClick = async () => {
    const name = prompt("Enter a name for the new history item:");
    if (name) {
      const requestData = { historyName: name };
      try {
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

        await serverResponse.text();
        setHistory((prevHistory) => [...prevHistory, name]);
        navigate("/inputSeq");
      } catch (error) {
        console.error("An error occurred during the request: ", error.message);
      }
    }
  };

  const handleHistoryClick = async (index) => {
    const historyName = history[index];
    const requestData = {
      historyName: historyName,
    };
    try {
      const serverResponse = await fetch(`http://localhost:8080/history/get`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(requestData),
      });

      if (!serverResponse.ok) {
        const errorMessage = await serverResponse.text();
        throw new Error(errorMessage);
      }

      const responseData = await serverResponse.json();
      navigate("/analysis", { state: { responseData } });
    } catch (error) {
      console.error(
        "An error occurred while fetching history details: ",
        error.message
      );
    }
  };

  const handleEllipsisClick = (e, index) => {
    e.stopPropagation(); // Prevents handleHistoryClick from being triggered
    const { top, left } = e.currentTarget.getBoundingClientRect(); // Get the position of the clicked button
    setMenuPosition({
      top: top + window.scrollY,
      left: left + e.currentTarget.offsetWidth + 10,
    }); // Set menu position
    setActiveHistoryItem(index);
  };

  const handleRename = async (index) => {
    const newName = prompt("Enter a new name for the history item:");
    if (newName) {
      const historyName = history[index];
      const requestData = {
        historyName: historyName,
        newName: newName,
      };

      try {
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

        await serverResponse.text();
        setHistory((prevHistory) =>
          prevHistory.map((item, i) => (i === index ? newName : item))
        );
      } catch (error) {
        console.error("An error occurred during the request: ", error.message);
      }
    }
    setActiveHistoryItem(null);
  };

  const handleDelete = async (index) => {
    if (window.confirm("Are you sure you want to delete this history item?")) {
      const currentName = history[index];
      const requestData = { historyName: currentName };

      try {
        const serverResponse = await fetch(
          "http://localhost:8080/history/delete",
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

        await serverResponse.text();
        setHistory((prevHistory) => prevHistory.filter((_, i) => i !== index));
      } catch (error) {
        console.error("An error occurred during the request: ", error.message);
      }
    }
    setActiveHistoryItem(null);
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
            onClick={handleEditClick}
            style={{ cursor: "pointer" }}
          />
        </div>
        <div className="sidebar-body">
          <br />
          <div className="history-list">
            {history.map((item, index) => (
              <div
                key={index}
                className="history-item"
                onClick={() => handleHistoryClick(index)}
                style={{ cursor: "pointer" }}
              >
                <span>{item}</span>
                <div className="ellipsis-container">
                  <button
                    className="ellipsis-button"
                    onClick={(e) => handleEllipsisClick(e, index)} // Pass event and index
                  >
                    ...
                  </button>
                </div>
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
                  onClick={handleEditClick}
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
                    const fetchHistory = async () => {
                      try {
                        const response = await fetch("http://localhost:8080/history/list");
                        if (!response.ok) {
                          throw new Error("Failed to fetch history list");
                        }
                        const data = await response.json();
                        setHistory(data);
                      } catch (error) {
                        console.error("Error fetching history:", error);
                      }
                    };
                
                    fetchHistory();
                    
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

      {/* Modal for options */}
      {activeHistoryItem !== null && (
        <div
          className="options-modal"
          ref={optionsMenuRef}
          style={{ top: menuPosition.top, left: menuPosition.left }} // Use dynamic positioning
        >
          <div className="options-menu">
            <button
              className="option-button"
              onClick={(e) => {
                e.stopPropagation();
                handleRename(activeHistoryItem);
              }}
            >
              Rename
            </button>
            <button
              className="option-button"
              onClick={(e) => {
                e.stopPropagation();
                handleDelete(activeHistoryItem);
              }}
            >
              Delete
            </button>
          </div>
        </div>
      )}
    </div>
  );
}

export default App;
