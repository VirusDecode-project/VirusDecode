import React from 'react';
import '../styles/Alignment.css';
import {Sequence, Chunk} from '../components/types';

interface SequenceDisplayProps {
  sequences: Sequence[];
  referenceSequence: string;
  onSequenceClick: (sequence: Sequence) => void;
  selectedSequence: Sequence | null;
  regionIndices: [number,number];
  selectedRegion: string;
}

const SequenceDisplay: React.FC<SequenceDisplayProps> = ({ sequences, referenceSequence, onSequenceClick, selectedSequence, regionIndices, selectedRegion }) => {
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

const splitSequences = (sequences: Sequence[], referenceSequence: string, regionIndices: [number, number]): Chunk[] => {
  const chunkSize = 50;
  const result: Chunk[] = [];
  const numChunks = Math.ceil((regionIndices[1] - regionIndices[0]) / chunkSize);

  for (let chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
    const chunk: Chunk = sequences.map((seq) => ({
      label: seq.label,
      lines: [seq.sequence
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

export default SequenceDisplay;