import Modal from './Modal';
import React, { useState, useEffect } from 'react';
import './ProteinSeq.css';

const generatePastelColor = () => {
  const hue = Math.floor(Math.random() * 360);
  const saturation = 70 + Math.floor(Math.random() * 10);
  const lightness = 85 + Math.floor(Math.random() * 10);
  return `hsl(${hue}, ${saturation}%, ${lightness}%)`;
};

function Alignment() {
    const [chartData, setChartData] = useState([]);
    const [selectedRegion, setSelectedRegion] = useState("ORF1ab");

    useEffect(() => {
        fetch('/alignment_data.json')  // JSON 파일의 경로 설정
            .then(response => response.json())
            .then(jsonData => {
                const { alignment_index } = jsonData;

                const data = Object.entries(alignment_index).map(([label, [start, end]]) => {
                    const value = end - start;  // 서열 길이 계산
                    return {
                        label,
                        value,
                        color: generatePastelColor(),
                        start,
                        end,
                    };
                });

                setChartData(data);
            })
            .catch(error => console.error('Error fetching sequence data:', error));
    }, []);

    return (
        <div>
            <div className="stacked-bar-chart">
                <StackedBar data={chartData} onBarClick={setSelectedRegion} />
            </div>
            <ProteinSeq selectedRegion={selectedRegion} setSelectedRegion={setSelectedRegion} />
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
      <div className='stackeb-bar-container'>
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
        </div></div>
    );
};


const ProteinSeq = ({ selectedRegion, setSelectedRegion }) => {  // setSelectedRegion 추가
  const [sequences, setSequences] = useState([]);
  const [selectedSequence, setSelectedSequence] = useState(null);
  const [alignmentIndex, setAlignmentIndex] = useState({});
  const [isModalOpen, setModalOpen] = useState(false);

  useEffect(() => {
    fetch('/alignment_data.json')
      .then(response => response.json())
      .then(jsonData => {
        const alignedSequences = jsonData.aligned_sequences;
        const sequenceData = Object.entries(alignedSequences).map(([label, sequence]) => ({
          label,
          sequence
        }));
        setSequences(sequenceData);
        setAlignmentIndex(jsonData.alignment_index);
      })
      .catch(error => console.error('Error fetching sequence data:', error));
  }, []);

  if (!sequences.length || !alignmentIndex[selectedRegion]) {
    return <p>Loading sequences...</p>;
  }

  const referenceSequence = sequences[0].sequence;
  const regionIndices = alignmentIndex[selectedRegion];
  const regionSequence = referenceSequence.slice(regionIndices[0], regionIndices[1]);

  const handleSequenceClick = (sequence) => {
    setSelectedSequence(sequence);
    setModalOpen(true);
  };

  const handleIndicesChange = (startIndex, endIndex) => {
    console.log('Start Index:', startIndex, 'End Index:', endIndex);
  };

  return (
    <div className="protein-sequence-container">
      <div className="region-selector">
        <label htmlFor="region-select">Select Protein Region: </label>
        <select id="region-select" value={selectedRegion} onChange={(e) => setSelectedRegion(e.target.value)}>
          {Object.keys(alignmentIndex).map(region => (
            <option key={region} value={region}>
              {region}
            </option>
          ))}
        </select>
      </div>
      <div>
        <SequenceDisplay 
          sequences={sequences} 
          referenceSequence={regionSequence}
          onSequenceClick={handleSequenceClick} 
          selectedSequence={selectedSequence}
          regionIndices={regionIndices}
        />
      </div>
      <Modal 
        isOpen={isModalOpen} 
        onClose={() => setModalOpen(false)} 
        selectedSequence={selectedSequence} 
        onIndicesChange={handleIndicesChange}
        sequences={sequences}  // sequences 전달
      />
    </div>
  );
};



const SequenceDisplay = ({ sequences, referenceSequence, onSequenceClick, selectedSequence, regionIndices }) => {
  useEffect(() => {
    const labels = document.querySelectorAll('.sequence-label');
    
    let maxWidth = 0;
    labels.forEach(label => {
      const labelWidth = label.offsetWidth;
      if (labelWidth > maxWidth) {
        maxWidth = labelWidth;
      }
    });
    maxWidth += 10;
    labels.forEach(label => {
      label.style.width = `${maxWidth}px`;
    });

    const indexesContainers = document.querySelectorAll('.sequence-indexes');
    indexesContainers.forEach((container) => {
      container.style.marginLeft = `${maxWidth}px`;
    });

  }, [sequences]);

  return (
    <div className="sequence-container">
      {splitSequences(sequences, referenceSequence, regionIndices).map((chunk, chunkIndex) => (
        <div key={chunkIndex} className="sequence-chunk">
          {chunk.map((seq, index) => (
            <div
              key={index}
              className={`sequence ${selectedSequence && seq.label === selectedSequence.label ? 'selected' : ''}`}
              onClick={() => onSequenceClick(seq)}
            >
              <div className="sequence-label">{seq.label}</div>
              <div className="sequence-boxes">
                {seq.lines.map((line, lineIndex) => (
                  <div key={lineIndex} className="sequence-line">
                    {line.map((charObj, charIndex) => (
                      <div 
                        key={charIndex} 
                        className={`sequence-box ${charObj.different ? 'different' : ''} ${charObj.char === '-' ? 'gap' : ''}`}
                      >
                        {charObj.char}
                      </div>
                    ))}
                  </div>
                ))}
              </div>
            </div>
          ))}
          <div className="sequence-indexes">
            {Array.from({ length: 5 }, (_, i) => {
              // Start indices from 1 relative to the start of the region
              const startPos = (chunkIndex * 50) + ((i + 1) * 10);
              return (
                <div key={i} className="sequence-index">
                  {startPos <= (regionIndices[1] - regionIndices[0]) ? startPos : '---'}
                </div>
              );
            })}
          </div>
        </div>
      ))}
    </div>
  );
};

// 서열을 나누는 함수
const splitSequences = (sequences, referenceSequence, regionIndices) => {
  const chunkSize = 50; // 각 청크의 크기 설정
  const result = [];
  const numChunks = Math.ceil((regionIndices[1] - regionIndices[0]) / chunkSize); // 청크 수 계산

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
              char: isInRange ? char : '', // 범위 밖일 경우 공백 출력
              different: isInRange && char !== referenceSequence[chunkIndex * chunkSize + index] // 참조 서열과 다른지 확인
            };
          })
      ]
    }));
    result.push(chunk);
  }

  return result;
};

