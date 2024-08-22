import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import './Modal.css';

const Modal = ({ isOpen, onClose, selectedSequence, sequences, selectedRegion}) => {
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

  // 백엔드 요청을 처리하는 함수
const handleConvertButton = async () => {
  // 데이터 객체 생성
  const data = {
    region: selectedRegion,   // region 필드에 해당하는 값
    varientName: selectedSequence.label, // varientName 필드에 해당하는 값
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

    onClose();

    // 백엔드로 데이터 전송
    await handleConvertButton(); // handleConvertButton 호출 추가

    // 라우팅을 프로그래밍적으로 수행
    // navigate(`/${selectedGenome.replace(/\s+/g, '-')}`);
    // GK
    // changeTab(1)
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
