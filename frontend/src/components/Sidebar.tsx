import React, { Dispatch, useState, useRef, useEffect, SetStateAction, MouseEvent } from "react";
import historyIcon from "../assets/history.png";
import editIcon from "../assets/edit.png";
import renameIcon from "../assets/rename.png";
import deleteIcon from "../assets/delete.png";
import HistoryList from "./HistoryList";
import CreateModal from "./CreateModal";
import RenameModal from "./RenameModal";
import DeleteModal from "./DeleteModal";
import { NavigateFunction } from 'react-router-dom';
import { MRNAData, AlignmentData, PDBResponse } from '../components/types';

interface SidebarProps {
  show: boolean;
  handleClose: () => void;
  history: string[],
  setHistory: Dispatch<SetStateAction<string[]>>,
  navigate: NavigateFunction;
  setMRNAReceived: Dispatch<SetStateAction<boolean>>;
  setPDBReceived: Dispatch<SetStateAction<boolean>>;
  setTab: Dispatch<SetStateAction<number>>;
  workingHistory: string;
  setWorkingHistory: Dispatch<SetStateAction<string>>;
  setLinearDesignData: Dispatch<SetStateAction<MRNAData | null>>;
  setPDBids: Dispatch<SetStateAction<string[]>>;
  setPDBInfo: Dispatch<SetStateAction<string[]>>;
  setSelectedPDBid: Dispatch<SetStateAction<string>>;
  setAlignmentData: Dispatch<SetStateAction<AlignmentData>>;
}

const Sidebar: React.FC<SidebarProps> = ({
  show,
  handleClose,
  history,
  setHistory,
  navigate,
  setMRNAReceived,
  setPDBReceived,
  setTab,
  workingHistory,
  setWorkingHistory,
  setLinearDesignData,
  setPDBids,
  setPDBInfo,
  setSelectedPDBid,
  setAlignmentData,
}) => {
  const [showEditModal, setShowEditModal] = useState(false);
  const [showRenameModal, setShowRenameModal] = useState(false);
  const [showDeleteModal, setShowDeleteModal] = useState(false);
  const [activeHistoryItem, setActiveHistoryItem] = useState<number | null>(null);
  const [menuPosition, setMenuPosition] = useState({ top: 0, left: 0 });
  const [showOptionsModal, setShowOptionsModal] = useState(false);
  const optionsMenuRef = useRef<HTMLDivElement>(null); // 옵션 메뉴를 감지할 참조 생성

  useEffect(() => {
    const handleClickOutside: EventListener = (event) => {
      if (
        optionsMenuRef.current &&
        !optionsMenuRef.current.contains(event.target as Node)
      ) {
        setShowOptionsModal(false);
      }
    };
    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, []);

  const handleEditClick = () => {
    setShowEditModal(true); // 편집 모달 열기
  };
  const checkPDBFileExists = async (url: string) => {
    try {
      const response = await fetch(url, { method: 'HEAD' });
      return response.ok;
    } catch (error) {
      console.error('Error checking PDB file existence:', error);
      return false;
    }
  };

  const handleHistoryClick = async (index: number) => {
    const historyName = history[index];
    const requestData = { historyName: historyName };

    try {
      // Fetch history details
      const serverResponse = await fetch(`/api/history/get`, {
        method: "POST",
        credentials: 'include',
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
      console.log("History details fetched successfully: ", responseData);

      if (responseData.alignment) {
        const alignmentJson = JSON.parse(responseData.alignment);
        setAlignmentData(alignmentJson);
      } else {
        console.warn('Alignment data is null or undefined');
      }

      if (responseData.linearDesign) {
        const linearDesignJson = JSON.parse(responseData.linearDesign);
        setLinearDesignData(linearDesignJson);
        setMRNAReceived(true);
      } else {
        setMRNAReceived(false);
        setTab(0)
      }

      if (responseData.pdb) {
        const pdbJson: PDBResponse = JSON.parse(responseData.pdb);
        const keys = Object.keys(pdbJson);
        setPDBids(keys);
        const values = Object.values(pdbJson);
        setPDBInfo(values);
        if (keys.length > 0) {
          setPDBReceived(true);
          for (let i = 0; i < keys.length; i++) {
            const exist = await (checkPDBFileExists(`https://files.rcsb.org/download/${keys[i]}.pdb`))
            if (exist) {
              setSelectedPDBid(keys[i]);
              break;
            }
          }
        }
      } else {
        setPDBReceived(false);
      }
      setWorkingHistory(historyName);
      navigate('/analysis');
    } catch (error) {
      if (error instanceof Error) {
        console.error("An error occurred: ", error.message);
      }
    }
  };

  const handleEllipsisClick = (e: MouseEvent<HTMLButtonElement>, index: number) => {
    e.stopPropagation(); // Prevents handleHistoryClick from being triggered
    const { top, left } = e.currentTarget.getBoundingClientRect(); // Get the position of the clicked button
    setMenuPosition({
      top: top + window.scrollY,
      left: left + e.currentTarget.offsetWidth + 10,
    }); // Set menu position
    setActiveHistoryItem(index);
    setShowOptionsModal(true);
  };

  const handleRename = async (newName: string) => {
    if (newName && activeHistoryItem !== null) {
      const historyName = history[activeHistoryItem];
      const requestData = {
        historyName: historyName,
        newName: newName,
      };

      try {
        const serverResponse = await fetch(
          `/api/history/rename`,
          {
            method: "PUT",
            credentials: 'include',
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
          prevHistory.map((item, i) =>
            i === activeHistoryItem ? newName : item
          )
        );

        if (historyName === workingHistory) {
          setWorkingHistory(newName);
        }
        setShowRenameModal(false); // 이름 변경 모달 닫기
      } catch (error) {
        if (error instanceof Error) {
          console.error("An error occurred during the request: ", error.message);
        }
      }
    } else {
      console.log("No name provided or no active item selected.");
    }
  };

  const handleDelete = async () => {
    if (activeHistoryItem !== null) {
      const currentName = history[activeHistoryItem];
      const requestData = { historyName: currentName };

      try {
        const serverResponse = await fetch(
          `/api/history/delete`,
          {
            method: "DELETE",
            credentials: 'include',
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
          prevHistory.filter((_, i) => i !== activeHistoryItem)
        );
      } catch (error) {
        if (error instanceof Error) {
          console.error("An error occurred during the request: ", error.message);
        }
      }
    }
    setShowDeleteModal(false); // 삭제 확인 모달 닫기
  };

  return (
    <div className={`sidebar ${show ? "show" : ""}`}>
      <div className="sidebar-header">
        <img
          className="history-icon"
          src={historyIcon}
          onClick={handleClose}
          style={{ cursor: "pointer" }}
          alt="History Icon"
        />
        <img
          src={editIcon}
          className="edit-icon"
          onClick={handleEditClick}
          style={{ cursor: "pointer" }}
          alt="Edit Icon"
        />
      </div>
      <div className="sidebar-body">
        <br />
        <HistoryList
          history={history}
          handleHistoryClick={handleHistoryClick}
          handleEllipsisClick={handleEllipsisClick}
        />
        {showOptionsModal && (
          <div
            className="options-modal"
            style={{ top: menuPosition.top, left: menuPosition.left }} // Use dynamic positioning
            ref={optionsMenuRef} // ref 추가
          >
            <div className="options-menu">
              <button
                className="option-button"
                onClick={(e) => {
                  e.stopPropagation();
                  setShowOptionsModal(false);
                  setShowRenameModal(true);
                }}
              >
                <img
                  src={renameIcon}
                  style={{ cursor: "pointer" }}
                  alt="rename Icon"
                />
                Rename
              </button>
              <button
                className="option-button"
                onClick={(e) => {
                  e.stopPropagation();
                  setShowOptionsModal(false);
                  setShowDeleteModal(true);
                }}
              >
                <div>
                  <img
                    src={deleteIcon}
                    style={{ cursor: "pointer" }}
                    alt="delete Icon"
                  />
                  Delete
                </div>
              </button>
            </div>
          </div>
        )}
      </div>
      {/* 커스텀 모달들을 추가 */}
      <CreateModal
        show={showEditModal}
        onClose={() => setShowEditModal(false)}
      />
      <RenameModal
        show={showRenameModal}
        onClose={() => setShowRenameModal(false)}
        onRename={handleRename}
      />
      <DeleteModal
        show={showDeleteModal}
        onClose={() => setShowDeleteModal(false)}
        onDelete={handleDelete}
      />
    </div>
  );
};

export default Sidebar;
