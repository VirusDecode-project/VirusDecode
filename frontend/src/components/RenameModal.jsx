import React, { useState, useRef, useEffect } from 'react';

const RenameModal = ({ show, onClose, onRename }) => {
  const [newName, setNewName] = useState('');
  const modalRef = useRef(null); // 모달 컨텐츠를 참조할 ref

  const handleRename = () => {
    if (newName.trim()) {  // 이름이 비어 있지 않은 경우에만 실행
      onRename(newName);  // 새로운 이름을 부모 컴포넌트로 전달
      setNewName('');  // 입력창 초기화
      onClose();  // 모달 닫기
    } else {
      alert("Please enter correct name.");  // 빈 입력 방지
    }
  };

  useEffect(() => {
    const handleClickOutside = (event) => {
      if (modalRef.current && !modalRef.current.contains(event.target)) {
        onClose(); // 모달 외부를 클릭했을 때 모달 닫기
        setNewName('');
      }
    };

    if (show) {
      document.addEventListener('mousedown', handleClickOutside);
    }

    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, [show, onClose]);

  const handleClose = () => {
    onClose();
    setNewName('');
  }

  if (!show) return null;  // 모달이 보이지 않으면 렌더링하지 않음

  return (
    <div className="history-modal-overlay">
      <div className="history-modal-content" ref={modalRef}>
        <h2>Please enter a new name.</h2>
        <input
          type="text"
          value={newName}
          onChange={(e) => setNewName(e.target.value)}
          placeholder="Enter new name"
        />
        <div className="history-modal-buttons">
          <button className="modal-next-button" onClick={handleRename}>Save</button>
          <button className="modal-close-button" onClick={handleClose}>Cancle</button>
        </div>
      </div>
    </div>
  );
};

export default RenameModal;
