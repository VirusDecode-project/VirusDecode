import React, { Dispatch, SetStateAction, useState, useEffect } from 'react';
import '../styles/Alignment.css';
import ConvertModal from './ConvertModal';
import helpIcon from '../assets/help.png';
import HelpModal from './HelpModal';
import SequenceDisplay from './SequenceDisplay'
import { Sequence, AlignmentData, MRNAData } from '../components/types';

interface ProteinSeqProps {
  onRegionUpdate: (region: string) => void;
  selectedRegion: string;
  setSelectedRegion: Dispatch<SetStateAction<string>>;
  alignmentData: AlignmentData;
  setTab: Dispatch<SetStateAction<number>>;
  setIsLoading: Dispatch<SetStateAction<boolean>>;
  setMRNAReceived: Dispatch<SetStateAction<boolean>>;
  setPDBReceived: Dispatch<SetStateAction<boolean>>;
  workingHistory: string;
  setLinearDesignData: Dispatch<SetStateAction<MRNAData | null>>;
  setPDBids: Dispatch<SetStateAction<string[]>>;
  setPDBInfo: Dispatch<SetStateAction<string[]>>;
  setSelectedPDBid: Dispatch<SetStateAction<string>>;
}

const ProteinSeq: React.FC<ProteinSeqProps> = ({ onRegionUpdate, selectedRegion, setSelectedRegion, alignmentData, setTab, setIsLoading, setMRNAReceived, setPDBReceived, workingHistory, setLinearDesignData, setPDBids, setPDBInfo, setSelectedPDBid }) => {
  const [sequences, setSequences] = useState<Sequence[]>([]);
  const [displayedSequences, setDisplayedSequences] = useState<Sequence[]>([]);
  const [selectedSequence, setSelectedSequence] = useState<Sequence | null>(null);
  const [isModalOpen, setModalOpen] = useState(false);
  const [modalData, setModalData] = useState({ genome: '', protein: '' });
  const [isHelpModalOpen, setHelpModalOpen] = useState(false);

  const SEQUENCES_TO_SHOW = 10; // Number of sequences to show initially and incrementally

  useEffect(() => {
    if (alignmentData) {
      const sequenceData: Sequence[] = Object.entries(alignmentData.aligned_sequences).map(([label, sequence]) => ({
        label,
        sequence,
      }));
      setSequences(sequenceData);
      setDisplayedSequences(sequenceData.slice(0, SEQUENCES_TO_SHOW)); // Show initial set of sequences
    }
  }, [alignmentData]);

  const handleShowMore = () => {
    setDisplayedSequences(prevDisplayed => [
      ...prevDisplayed,
      ...sequences.slice(prevDisplayed.length, prevDisplayed.length + SEQUENCES_TO_SHOW),
    ]);
  };

  if (!displayedSequences.length || !alignmentData.alignment_index[selectedRegion]) {
    return <p>Loading sequences...</p>;
  }

  const referenceSequence = sequences[0].sequence;
  const regionIndices = alignmentData.alignment_index[selectedRegion];
  const regionSequence = referenceSequence.slice(regionIndices[0], regionIndices[1]);

  const handleSequenceClick = (sequence: Sequence) => {
    setSelectedSequence(sequence);
    setModalData({
      genome: sequence.label,
      protein: selectedRegion,
    });
    setModalOpen(true);
  };

  const toggleHelpModal = () => {
    setHelpModalOpen(!isHelpModalOpen);
  };

  return (
    <div className="protein-sequence-container">
      <div className="region-selector">
        <div className='region-selector-wrapper'>
          <div className='region-selector-container'>
            <img
              className="help-icon"
              src={helpIcon}
              onClick={toggleHelpModal}
              style={{ cursor: "pointer" }}
              alt="Help"
            />
            <label htmlFor="region-select">Select Coding Sequence: </label>
            <select id="region-select" value={selectedRegion} onChange={(e) => setSelectedRegion(e.target.value)}>
              {Object.keys(alignmentData.alignment_index).map(region => (
                <option key={region} value={region}>
                  {region}
                </option>
              ))}
            </select>
          </div>
          <HelpModal
            isOpen={isHelpModalOpen}
            onClose={toggleHelpModal}
          />
        </div>
      </div>
      <div>
        <SequenceDisplay
          sequences={displayedSequences} // Display the sliced list of sequences
          referenceSequence={regionSequence}
          onSequenceClick={handleSequenceClick}
          selectedSequence={selectedSequence}
          regionIndices={regionIndices}
          selectedRegion={selectedRegion}
        />
        {displayedSequences.length < sequences.length && (
          <button className="show-more-button" onClick={handleShowMore}>
            Show More
          </button>
        )}
      </div>
      <ConvertModal
        onRegionUpdate={onRegionUpdate}
        isOpen={isModalOpen}
        onClose={() => setModalOpen(false)}
        sequences={sequences}
        alignmentIndex={alignmentData.alignment_index}
        modalData={modalData}
        setTab={setTab}
        setIsLoading={setIsLoading}
        setMRNAReceived={setMRNAReceived}
        setPDBReceived={setPDBReceived}
        workingHistory={workingHistory}
        setLinearDesignData={setLinearDesignData}
        setPDBids={setPDBids}
        setPDBInfo={setPDBInfo}
        setSelectedPDBid={setSelectedPDBid}
      />
    </div>
  );
};

export default ProteinSeq;