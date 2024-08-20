import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import './Modal.css';

const Modal = ({ isOpen, onClose, selectedSequence, onIndicesChange, sequences }) => {
  const [startIndex, setStartIndex] = useState('');
  const [endIndex, setEndIndex] = useState('');
  const [selectedGenome, setSelectedGenome] = useState(selectedSequence ? selectedSequence.label : '');
  const [error, setError] = useState('');

  useEffect(() => {
    setStartIndex('');
    setEndIndex('');
    setError('');
  }, [selectedGenome]);

  if (!isOpen) return null;

  const handleGenomeChange = (e) => {
    setSelectedGenome(e.target.value);
  };

  const handleStartIndexChange = (e) => {
    setStartIndex(e.target.value);
  };

  const handleEndIndexChange = (e) => {
    setEndIndex(e.target.value);
  };

  const getMaxEndIndex = () => {
    const selectedSeq = sequences.find(seq => seq.label === selectedGenome);
    return selectedSeq ? selectedSeq.sequence.length : 0;
  };

  const handleNext = () => {
    const start = parseInt(startIndex, 10);
    const end = parseInt(endIndex, 10);
    const maxEndIndex = getMaxEndIndex();

    if (isNaN(start) || isNaN(end)) {
      setError('Please enter valid indices.');
      return;
    }

    if (start < 1 || start > maxEndIndex) {
      setError(`Start index must be between 1 and ${maxEndIndex}.`);
      return;
    }

    if (end < start || end > maxEndIndex) {
      setError(`End index must be between ${start} and ${maxEndIndex}.`);
      return;
    }

    onIndicesChange(start, end);
    onClose();
  };

  return (
    <div className="modal-overlay">
      <div className="modal-content">
        <h2>Select Amino Acid Number for {selectedGenome}</h2>
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
            />
          </label>
          <label>
            End Index:
            <input
              type="number"
              value={endIndex}
              onChange={handleEndIndexChange}
              min={startIndex}
            />
          </label>
        </div>
        {error && <div className="error-message">{error}</div>}
        <div className="modal-button-group">
          <button className="modal-close-button" onClick={onClose}>
            Cancel
          </button>
          <button className="modal-next-button" onClick={handleNext}>
            <Link to={`/${selectedGenome.replace(/\s+/g, '-')}`}>
              Convert
            </Link>
          </button>
        </div>
      </div>
    </div>
  );
};

export default Modal;
