import React, { useState, useEffect } from 'react';
import './Modal.css';

const Modal = ({ isOpen, onClose, sequences, alignmentIndex, modalData, setTab }) => {
  const [startIndex, setStartIndex] = useState('');
  const [endIndex, setEndIndex] = useState('');
  const [selectedGenome, setSelectedGenome] = useState('');
  const [selectedRegion, setSelectedRegion] = useState('');
  const [error, setError] = useState('');

  useEffect(() => {
    if (modalData) {
      setSelectedGenome(modalData.genome || '');
      setSelectedRegion(modalData.protein || '');
    }
  }, [modalData]);

  if (!isOpen) return null;

  const handleGenomeChange = (e) => {
    setSelectedGenome(e.target.value);
  };

  const handleRegionChange = (e) => {
    setSelectedRegion(e.target.value);
  };

  const handleStartIndexChange = (e) => {
    setStartIndex(e.target.value);
  };

  const handleEndIndexChange = (e) => {
    setEndIndex(e.target.value);
  };

  const getMaxEndIndex = () => {
    if (selectedRegion) {
      return alignmentIndex[selectedRegion][1] - alignmentIndex[selectedRegion][0];
    }
    const selectedSeq = sequences.find(seq => seq.label === selectedGenome);
    return selectedSeq ? selectedSeq.sequence.length : 0;
  };


  // 백엔드 요청을 처리하는 함수
const handleConvertButton = async () => {
  // 데이터 객체 생성
  const data = {
    region: selectedRegion,   // region 필드에 해당하는 값
    varientName: selectedGenome, // varientName 필드에 해당하는 값
    start: parseInt(startIndex, 10),  // start 필드에 해당하는 값 (숫자형 변환)
    end: parseInt(endIndex, 10),      // end 필드에 해당하는 값 (숫자형 변환)
  };

  console.log("Data being sent to backend:", data); // 데이터 전송 전 로그 출력
  
  try {
    const response = await fetch('http://localhost:8080/analysis/mrnadesign', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(data), // JSON으로 변환된 데이터 객체
    });

    if (!response.ok) {
      throw new Error('Network response was not ok');
    }

    const responseData = await response.json();
    console.log('Response from server:', responseData);
  } catch (error) {
    console.error('Error sending data:', error);
    // 사용자에게 에러 메시지를 표시하는 로직을 추가할 수 있음
  }
};


const handleNext = async () => {
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

    // 백엔드로 데이터 전송
    await handleConvertButton();
    setTab(1);
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
            Protein Region:
            <select value={selectedRegion} onChange={handleRegionChange}>
              <option value="">Select a region</option>
              {Object.keys(alignmentIndex).map((region) => (
                <option key={region} value={region}>
                  {region}
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
