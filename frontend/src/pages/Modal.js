import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import './Modal.css';

const Modal = ({ isOpen, onClose, selectedSequence, onIndicesChange, sequences }) => {
  const [startIndex, setStartIndex] = useState('');
  const [endIndex, setEndIndex] = useState('');
  const [selectedGenome, setSelectedGenome] = useState(selectedSequence ? selectedSequence.label : '');

  useEffect(() => {
    // Reset the indices when the selected genome changes
    setStartIndex('');
    setEndIndex('');
  }, [selectedGenome]);

  if (!isOpen) return null;

  const handleGenomeChange = (e) => {
    setSelectedGenome(e.target.value);
  };

  const handleStartIndexChange = (e) => {
    const value = parseInt(e.target.value, 10);
    if (value >= 1 && value <= getMaxEndIndex()) {
      setStartIndex(value);
      if (endIndex < value) {
        setEndIndex(value); // Ensure endIndex is always greater or equal to startIndex
      }
    }
  };

  const handleEndIndexChange = (e) => {
    const value = parseInt(e.target.value, 10);
    if (value >= startIndex && value <= getMaxEndIndex()) {
      setEndIndex(value);
    }
  };

  const getMaxEndIndex = () => {
    const selectedSeq = sequences.find(seq => seq.label === selectedGenome);
    return selectedSeq ? selectedSeq.sequence.length : 0;
  };

  const handleNext = () => {
    onIndicesChange(startIndex, endIndex);
    onClose();
  };

  return (
    <div className="modal-overlay">
      <div className="modal-content">
        <h2>Select Indices for {selectedGenome}</h2>
        <div className="modal-inputs">
          <label>
            Genome:
            <select value={selectedGenome} onChange={handleGenomeChange}>
              {sequences.map((seq) => (
                <option key={seq.label} value={seq.label}>
                  {seq.label}
                </option>
              ))}
            </select>
          </label>
          <label>
            Start Index:
            <input
              type="number"
              value={startIndex}
              onChange={handleStartIndexChange}
              min="1"
              max={getMaxEndIndex()}
            />
          </label>
          <label>
            End Index:
            <input
              type="number"
              value={endIndex}
              onChange={handleEndIndexChange}
              min={startIndex}
              max={getMaxEndIndex()}
            />
          </label>
        </div>
        <button className="modal-close-button" onClick={onClose}>
          Cancel
        </button>
        <button className="modal-next-button" onClick={handleNext} disabled={!startIndex || !endIndex}>
          <Link to={`/${selectedGenome.replace(/\s+/g, '-')}`}>
            Next
          </Link>
        </button>
      </div>
    </div>
  );
};

export default Modal;
