import React, { useState, useEffect } from 'react';
import '../styles/Modal.css';

const Modal = ({ onRegionUpdate, isOpen, onClose, sequences, alignmentIndex, modalData, setTab, setIsLoading, setMRNAReceived, setPDBReceived }) => {
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

  const handleConvertButton = async () => {
    const mRnaData = {
      region: selectedRegion,
      varientName: selectedGenome,
      start: parseInt(startIndex, 10),
      end: parseInt(endIndex, 10),
    };
    onRegionUpdate(mRnaData.region);
    console.log("Data being sent to backend:", mRnaData);

    try {
      setIsLoading(true);
      const serverResponse = await fetch('http://localhost:8080/analysis/mrnadesign', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(mRnaData),
      });

      if (!serverResponse.ok) {
        const errorMessage = await serverResponse.text();
        throw new Error(errorMessage);
      }

      const responseData = await serverResponse.json();
      setIsLoading(false);
      setTab(1);
      setMRNAReceived(true);
      console.log('Response from server:', responseData);
    } catch (error) {
      console.error("An error occurred during the request: ", error.message);
    }


    const pdbData = { gene: selectedRegion }
    try {
      const serverResponse = await fetch('http://localhost:8080/analysis/pdb', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(pdbData),
      });

      if (!serverResponse.ok) {
        const errorMessage = await serverResponse.text();
        throw new Error(errorMessage);
      }

      await serverResponse.text();
      setPDBReceived(true);
    } catch (error) {
      console.error("An error occurred during the request: ", error.message);
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

    await handleConvertButton();
    onClose();
  };

  return (
    <div className="modal-overlay">
      <div className="modal-content">
        <h2>Select Amino Acid Number for {selectedGenome}</h2>
        <div className="modal-inputs">
          <label>
            Sublineage:
            <select value={selectedGenome} onChange={handleGenomeChange}>
              {sequences.map((seq) => (
                <option key={seq.label} value={seq.label}>
                  {seq.label}
                </option>
              ))}
            </select>
          </label>
          <label>
            Select Coding Sequence:
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
            Start Amino Acid Position:
            <input
              type="number"
              value={startIndex}
              onChange={handleStartIndexChange}
              min="1"
            />
          </label>
          <label>
            End Amino Acid Position:
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
