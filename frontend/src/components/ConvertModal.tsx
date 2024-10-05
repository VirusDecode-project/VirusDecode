import React, { Dispatch, SetStateAction, useState, useEffect, ChangeEvent } from 'react';
import '../styles/Modal.css';
import {Sequence, AlignmentIndex, ModalData, MRNAData, PDBResponse} from './types';

interface ConvertModalProps {
  onRegionUpdate: (region: string) => void;
  isOpen: boolean;
  onClose: () => void;
  sequences: Sequence[];
  alignmentIndex: AlignmentIndex;
  modalData: ModalData;
  setTab: Dispatch<SetStateAction<number>>;
  setIsLoading: Dispatch<SetStateAction<boolean>>;
  setMRNAReceived: Dispatch<SetStateAction<boolean>>;
  setPDBReceived: Dispatch<SetStateAction<boolean>>;
  workingHistory: string;
  setLinearDesignData:Dispatch<SetStateAction<MRNAData | null>>;
  setPDBids: Dispatch<SetStateAction<string[]>>;
  setPDBInfo: Dispatch<SetStateAction<string[]>>;
  setSelectedPDBid: Dispatch<SetStateAction<string>>;
}



const ConvertModal: React.FC<ConvertModalProps> = ({ onRegionUpdate, isOpen, onClose, sequences, alignmentIndex, modalData, setTab, setIsLoading, setMRNAReceived, setPDBReceived, workingHistory, setLinearDesignData, setPDBids, setPDBInfo, setSelectedPDBid }) => {
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

  const handleGenomeChange = (e: ChangeEvent<HTMLSelectElement>) => {
    setSelectedGenome(e.target.value);
  };

  const handleRegionChange = (e: ChangeEvent<HTMLSelectElement>) => {
    setSelectedRegion(e.target.value);
  };

  const handleStartIndexChange = (e: ChangeEvent<HTMLInputElement>) => {
    setStartIndex(e.target.value);
  };

  const handleEndIndexChange = (e: ChangeEvent<HTMLInputElement>) => {
    setEndIndex(e.target.value);
  };

  const getMaxEndIndex = () => {
    if (selectedRegion) {
      return alignmentIndex[selectedRegion][1] - alignmentIndex[selectedRegion][0];
    }
    const selectedSeq = sequences.find(seq => seq.label === selectedGenome);
    return selectedSeq ? selectedSeq.sequence.length : 0;
  };

  const checkPDBFileExists = async (url: string) => {
    try {
      const response = await fetch(url, { method: 'HEAD' });
      return response.ok;
    } catch (error) {
      console.error('Error checking PDB file existence:', error);
      return false;
    }
  };
  const handleConvertButton = async () => {
    setMRNAReceived(false);
    setPDBReceived(false);
    const mRnaData = {
      region: selectedRegion,
      varientName: selectedGenome,
      start: parseInt(startIndex, 10),
      end: parseInt(endIndex, 10),
      historyName: workingHistory,
    };
    onRegionUpdate(mRnaData.region);
    console.log("Data being sent to backend:", mRnaData);
  
    try {
      setIsLoading(true);
  
      // 1. mRNA 디자인 POST 요청
      const mRnaResponse = await fetch(`/api/analysis/linearDesign`, {
        method: 'POST',
        credentials: 'include',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(mRnaData),
      });
  
      if (!mRnaResponse.ok) {
        const errorMessage = await mRnaResponse.text();
        throw new Error(errorMessage);
      }
  
      const linearDesignResponse = await mRnaResponse.json();

      setMRNAReceived(true);
      setLinearDesignData(linearDesignResponse);
      setTab(1);
      setIsLoading(false);
  
      // 2. PDB 디자인 POST 요청
      const pdbData = { gene: selectedRegion, historyName: workingHistory};
      const pdbResponse = await fetch(`/api/analysis/pdb`, {
        method: 'POST',
        credentials: 'include',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(pdbData),
      });
  
      if (!pdbResponse.ok) {
        const errorMessage = await pdbResponse.text();
        throw new Error(errorMessage);
      }
  
      const responseData: PDBResponse = await pdbResponse.json();

        const keys = Object.keys(responseData);
        setPDBids(keys);
        const values = Object.values(responseData);
        setPDBInfo(values);
        
        if (keys.length > 0) {
          setPDBReceived(true);
          for (let i = 0; i < keys.length; i++){
            const exist = await (checkPDBFileExists(`https://files.rcsb.org/download/${keys[i]}.pdb`))
            if (exist) {
              setSelectedPDBid(keys[i]);
              break;
            }
          }
        }

    } catch (error) {
      if (error instanceof Error){
        console.error("An error occurred during the request: ", error.message);
        window.alert(error.message);
      }
    } finally {
      setIsLoading(false);
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

export default ConvertModal;
