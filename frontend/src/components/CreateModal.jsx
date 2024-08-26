import React, { useState, useRef, useEffect } from 'react';
import { useNavigate } from 'react-router-dom'; // Import the useNavigate hook

const CustomModal = ({ show, onClose, onSave }) => {
  const [inputValue, setInputValue] = useState('');
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

  const handleSave = () => {
    onSave(inputValue);
    onClose();
  };

  const handleStartNew = () => {
    navigate('/InputSeq'); // Navigate to the InputSeq page
    onClose(); // Close the modal after navigation
  };

  if (!show) return null;

  return (
    <div className="modal-overlay">
      <div className="modal-content" ref={modalRef}>
        <h2>Do you want to save this to your history for future reference?</h2>
        <input
          type="text"
          value={inputValue}
          onChange={(e) => setInputValue(e.target.value)}
          placeholder="Enter name to save"
        />
        <div className="modal-buttons">
          <button onClick={handleStartNew}>New</button> {/* New button */}
          <button onClick={handleSave}>Save</button>
          <button onClick={onClose}>Cancel</button>
        </div>
      </div>
    </div>
  );
};

export default CustomModal;
