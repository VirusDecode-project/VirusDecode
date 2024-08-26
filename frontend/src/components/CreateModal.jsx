import React, { useState, useRef, useEffect } from 'react';
import { useNavigate } from 'react-router-dom'; // Import the useNavigate hook

const CreateModal = ({ show, onClose, onSave }) => {
  const [inputValue, setInputValue] = useState('');
  const modalRef = useRef(null);
  const navigate = useNavigate(); // Initialize the navigate function

  useEffect(() => {
    const handleClickOutside = (event) => {
      if (modalRef.current && !modalRef.current.contains(event.target)) {
        onClose();
        setInputValue('');
      }
    };

    if (show) {
      document.addEventListener('mousedown', handleClickOutside);
    }

    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, [show, onClose]);

  const handleSave = () => {
    onSave(inputValue);
    setInputValue('');
    onClose();
  };

  const handleClose = () => {
    onClose();
    setInputValue('');
  }

  const handleStartNew = () => {
    navigate('/InputSeq'); // Navigate to the InputSeq page
    setInputValue('');
    onClose(); // Close the modal after navigation
  };

  if (!show) return null;

  return (
    <div className="history-modal-overlay">
      <div className="history-modal-content" ref={modalRef}>
        <h2>Do you want to save this to your history?</h2>
        <input
          type="text"
          value={inputValue}
          onChange={(e) => setInputValue(e.target.value)}
          placeholder="Enter name to save"
        />
        <div className="history-modal-buttons">
          <button className="modal-next-button" onClick={handleStartNew}>New</button> {/* New button */}
          <button className="modal-next-button" onClick={handleSave}>Save</button>
          <button className="modal-close-button"onClick={handleClose}>Cancel</button>
        </div>
      </div>
    </div>
  );
};

export default CreateModal;
