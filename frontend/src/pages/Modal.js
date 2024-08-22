import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import './Modal.css';

const Modal = ({ isOpen, onClose, selectedSequence, onIndicesChange, sequences }) => {
  const [startIndex, setStartIndex] = useState('');
  const [endIndex, setEndIndex] = useState('');
  const [selectedGenome, setSelectedGenome] = useState(selectedSequence ? selectedSequence.label : '');
  const [error, setError] = useState('');
  const navigate = useNavigate(); // useNavigate 훅 사용

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

    // 조건 검사 통과 여부 확인
    if (
      isNaN(start) ||
      isNaN(end) ||
      start < 1 ||
      start > maxEndIndex ||
      end < start ||
      end > maxEndIndex
    ) {
      setError('Please enter valid indices.');
      return;
    }

    onIndicesChange(start, end);
    onClose();

    // 라우팅을 프로그래밍적으로 수행
    navigate(`/${selectedGenome.replace(/\s+/g, '-')}`);
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
            Convert
          </button>
        </div>
      </div>
    </div>
  );
};

export default Modal;
