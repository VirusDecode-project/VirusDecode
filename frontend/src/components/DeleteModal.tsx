// DeleteModal.tsx
import React, { useRef, useEffect } from 'react';

interface DeleteModalProps {
  show: boolean;
  onClose: () => void;
  onDelete: () => void;
}

const DeleteModal: React.FC<DeleteModalProps> = ({ show, onClose, onDelete }) => {
  const modalRef = useRef<HTMLDivElement>(null); // 모달 컨텐츠를 참조할 ref

  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (modalRef.current && !modalRef.current.contains(event.target as Node)) {
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
    <div className="history-modal-overlay">
      <div className="history-modal-content" ref={modalRef}>
        <h2>Are you sure you want to delete this history?</h2>
        <div className="history-modal-buttons">
          <button className="modal-next-button" onClick={onDelete}>Delete</button>
          <button className="modal-close-button" onClick={onClose}>Cancle</button>
        </div>
      </div>
    </div>
  );
};

export default DeleteModal;
