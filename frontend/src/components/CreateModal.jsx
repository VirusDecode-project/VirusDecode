import React, { useRef, useEffect } from 'react';
import { useNavigate } from 'react-router-dom'; // Import the useNavigate hook

const CreateModal = ({ show, onClose }) => {
  const modalRef = useRef(null);
  const navigate = useNavigate(); // Initialize the navigate function

  useEffect(() => {
    const handleClickOutside = (event) => {
      if (modalRef.current && !modalRef.current.contains(event.target)) {
        onClose();
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
  }

  const handleStartNew = () => {
    navigate('/InputSeq'); // Navigate to the InputSeq page
    onClose(); // Close the modal after navigation
  };

  if (!show) return null;

  return (
    <div className="history-modal-overlay">
      <div className="history-modal-content" ref={modalRef}>
        <h2>Your history was saved. Do you want to restart?</h2>

        <div className="history-modal-buttons">
          <button className="modal-next-button" onClick={handleStartNew}>restart</button> {/* New button */}
          <button className="modal-close-button"onClick={handleClose}>Cancel</button>
        </div>
      </div>
    </div>
  );
};

export default CreateModal;
