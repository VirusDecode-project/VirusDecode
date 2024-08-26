// DeleteModal.jsx
import React, { useRef, useEffect } from 'react';

const DeleteModal = ({ show, onClose, onDelete }) => {
  const modalRef = useRef(null); // 모달 컨텐츠를 참조할 ref

  useEffect(() => {
    const handleClickOutside = (event) => {
      if (modalRef.current && !modalRef.current.contains(event.target)) {
        onClose(); // 모달 외부를 클릭했을 때 모달 닫기
      }
    };

    if (show) {
      document.addEventListener('mousedown', handleClickOutside);
    }

    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, [show, onClose]);

  if (!show) return null; // 모달이 보이지 않으면 렌더링하지 않음

  return (
    <div className="modal-overlay">
      <div className="modal-content" ref={modalRef}>
        <h2>Are you sure you want to delete this history?</h2>
        <div className="modal-buttons">
          <button onClick={onDelete}>Delete</button>
          <button onClick={onClose}>Cancle</button>
        </div>
      </div>
    </div>
  );
};

export default DeleteModal;
