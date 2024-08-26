import React, { useState, useEffect } from 'react';
import '../styles/Alignment.css';
import ConvertModal from '../components/ConvertModal';
import helpIcon from '../assets/help.png';
import HelpModal from '../components/HelpModal';
import Loading from '../components/Loading';

let lastHue = 0;

const generatePastelColor = () => {
  const minDifference = 60;
  let hue;
  do {
    hue = Math.floor(Math.random() * 360);
  } while (Math.abs(hue - lastHue) < minDifference);

  lastHue = hue;
  const saturation = 70 + Math.floor(Math.random() * 10);
  const lightness = 85 + Math.floor(Math.random() * 10);
  return `hsl(${hue}, ${saturation}%, ${lightness}%)`;
};

function Alignment({ responseData, setTab, onRegionUpdate, setMRNAReceived, setPDBReceived, workingHistory }) {
  const [chartData, setChartData] = useState([]);
  const [selectedRegion, setSelectedRegion] = useState('');
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    if (responseData) {
      const firstRegion = Object.keys(responseData.alignment_index)[0];
      setSelectedRegion(firstRegion);
      try {
        const data = Object.entries(responseData.alignment_index).map(([label, [start, end]]) => {
          const value = end - start;
          return {
            label,
            value,
            color: generatePastelColor(),
            start,
            end,
          };
        });
        setChartData(data);
      } catch (error) {
        console.error('Error processing response data:', error);
      }
    }
  }, [responseData]);

  return (
    <div>
      {isLoading ? (
        <Loading text="Converting" />
      ) : (
        <div>
          <div className="stacked-bar-chart">
            <StackedBar data={chartData} onBarClick={setSelectedRegion} />
          </div>
          {responseData && (
            <ProteinSeq
              selectedRegion={selectedRegion}
              setSelectedRegion={setSelectedRegion}
              responseData={responseData}
              setTab={setTab}
              setIsLoading={setIsLoading}
              onRegionUpdate={onRegionUpdate}
              setMRNAReceived={setMRNAReceived}
              setPDBReceived={setPDBReceived}
              workingHistory={workingHistory}
            />
          )}
        </div>
      )}
    </div>
  );
}

export default Alignment;

const StackedBar = ({ data, onBarClick }) => {
  const [tooltip, setTooltip] = useState({ visible: false, text: '', x: 0, y: 0 });

  const handleMouseOver = (e, item) => {
    setTooltip({
      visible: true,
      text: `${item.label}: ${item.start + 1}~${item.end}`,
      x: e.clientX,
      y: e.clientY,
    });
  };

  const handleMouseOut = () => {
    setTooltip({ visible: false, text: '', x: 0, y: 0 });
  };

  const totalValue = data.reduce((acc, item) => acc + item.value, 0);

  return (
    <div className='stacked-bar-container'>
      <div className="stacked-bar">
        {data.map((item, index) => {
          const segmentWidthPercentage = (item.value / totalValue) * 100;
          const segmentWidthInPixels = (segmentWidthPercentage / 100) * window.innerWidth;

          return (
            <div
              key={index}
              className="stacked-bar-segment"
              style={{
                width: `${segmentWidthPercentage}%`,
                backgroundColor: item.color,
              }}
              onClick={() => onBarClick(item.label)}
              onMouseOver={(e) => handleMouseOver(e, item)}
              onMouseOut={handleMouseOut}
            >
              {segmentWidthInPixels > 50 && (
                <span className="stacked-bar-label">
                  {item.label}
                </span>
              )}
            </div>
          );
        })}
        {tooltip.visible && (
          <div
            className="custom-tooltip"
            style={{ top: tooltip.y + 10, left: tooltip.x + 10 }}
          >
            {tooltip.text}
          </div>
        )}
      </div>
    </div>
  );
};

const ProteinSeq = ({ onRegionUpdate, selectedRegion, setSelectedRegion, responseData, setTab, setIsLoading, setMRNAReceived, setPDBReceived, workingHistory }) => {
  const [sequences, setSequences] = useState([]);
  const [displayedSequences, setDisplayedSequences] = useState([]);
  const [selectedSequence, setSelectedSequence] = useState(null);
  const [isModalOpen, setModalOpen] = useState(false);
  const [modalData, setModalData] = useState({ genome: '', protein: '' });
  const [isHelpModalOpen, setHelpModalOpen] = useState(false);
  const [showMore, setShowMore] = useState(false); // New state for show more

  const SEQUENCES_TO_SHOW = 10; // Number of sequences to show initially and incrementally

  useEffect(() => {
    if (responseData) {
      const sequenceData = Object.entries(responseData.aligned_sequences).map(([label, sequence]) => ({
        label,
        sequence,
      }));
      setSequences(sequenceData);
      setDisplayedSequences(sequenceData.slice(0, SEQUENCES_TO_SHOW)); // Show initial set of sequences
    }
  }, [responseData]);

  const handleShowMore = () => {
    setDisplayedSequences(prevDisplayed => [
      ...prevDisplayed,
      ...sequences.slice(prevDisplayed.length, prevDisplayed.length + SEQUENCES_TO_SHOW),
    ]);
  };

  if (!displayedSequences.length || !responseData.alignment_index[selectedRegion]) {
    return <p>Loading sequences...</p>;
  }

  const referenceSequence = sequences[0].sequence;
  const regionIndices = responseData.alignment_index[selectedRegion];
  const regionSequence = referenceSequence.slice(regionIndices[0], regionIndices[1]);

  const handleSequenceClick = (sequence) => {
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
              {Object.keys(responseData.alignment_index).map(region => (
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
        alignmentIndex={responseData.alignment_index}
        modalData={modalData}
        setTab={setTab}
        setIsLoading={setIsLoading}
        setMRNAReceived={setMRNAReceived}
        setPDBReceived={setPDBReceived}
        workingHistory={workingHistory}
      />
    </div>
  );
};

const SequenceDisplay = ({ sequences, referenceSequence, onSequenceClick, selectedSequence, regionIndices, selectedRegion }) => {
  const maxLabelWidth = Math.max(
    ...sequences.map(seq => seq.label.length)
  ) + 1;

  return (
    <div className="sequence-container">
      {splitSequences(sequences, referenceSequence, regionIndices).map((chunk, chunkIndex) => (
        <div key={chunkIndex} className="sequence-chunk">
          {chunk.map((seq, index) => (
            <div
              key={index}
              className={`sequence ${selectedSequence && seq.label === selectedSequence.label ? 'selected' : ''}`}
              onClick={() => onSequenceClick(seq)}
              style={index === 0 ? { borderBottom: '2px solid #aaaaaa', paddingBottom: '6px', marginBottom: '6px' } : {}}
            >
              <div
                className="sequence-label"
                style={{ width: `${maxLabelWidth}ch` }}
              >
                {seq.label}
              </div>
              <div className="sequence-boxes">
                {seq.lines.map((line, lineIndex) => (
                  <div key={lineIndex} className="sequence-line">
                    {line.map((charObj, charIndex) => (
                      <div
                        key={charIndex}
                        className={`sequence-box ${charObj.char === '-' ? 'gap' : ''} ${charObj.different ? 'different' : ''} ${charObj.char === '' ? 'empty' : ''}`}
                        style={charObj.char === '' ? { border: 'none' } : {}}
                      >
                        {charObj.char}
                      </div>
                    ))}
                  </div>
                ))}
              </div>
            </div>
          ))}
          <div className="sequence-indexes" style={{ marginLeft: `${maxLabelWidth + 1}ch` }}>
            {Array.from({ length: 6 }, (_, i) => {
              const startPos = i === 0 ? chunkIndex * 50 + 1 : ((chunkIndex * 50) + (i * 10)).toString().padStart(2, '0');
              return (
                <div key={i} className="sequence-index">
                  {startPos || ' '}
                </div>
              );
            })}
          </div>
        </div>
      ))}
    </div>
  );
};

const splitSequences = (sequences, referenceSequence, regionIndices) => {
  const chunkSize = 50;
  const result = [];
  const numChunks = Math.ceil((regionIndices[1] - regionIndices[0]) / chunkSize);

  for (let chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
    const chunk = sequences.map((seq) => ({
      label: seq.label,
      lines: [
        seq.sequence
          .slice(regionIndices[0] + chunkIndex * chunkSize, regionIndices[0] + (chunkIndex + 1) * chunkSize)
          .split('')
          .map((char, index) => {
            const globalIndex = regionIndices[0] + chunkIndex * chunkSize + index;
            const isInRange = globalIndex >= regionIndices[0] && globalIndex < regionIndices[1];
            return {
              char: isInRange ? char : '',
              different: isInRange && char !== referenceSequence[chunkIndex * chunkSize + index]
            };
          })
      ]
    }));
    result.push(chunk);
  }

  return result;
};
