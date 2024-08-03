import React from 'react';
import { BrowserRouter as Router, Route, Routes, Link } from 'react-router-dom';
import './ProteinSeq.css';


//예시 단백질 시퀀스 데이터
const ProteinSeq = () => {
  const sequences = [
    {
      label: 'Reference',
      sequence: 'MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
    }, 
    {
      label: 'USA-WA1/2020',
      sequence: 'MESLHPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
    },
    {
      label: 'SARS-P2/2021',
      sequence: 'MESLHPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
    },
    {
      label: 'SH-P10-A-2-Shanghai/2020',
      sequence: 'MESLHPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
    }
  ];

  return (
    <Router>
      <div className="ProteinSeq">
        <SequenceDisplay sequences={sequences} />
        <Routes>
          {sequences.map((seq, index) => (
            <Route
              key={index}
              path={`/${seq.label.replace(/\s+/g, '-')}`}
              element={<SequenceDetails label={seq.label} sequence={seq.sequence} />}
            />
          ))}
        </Routes>
      </div>
    </Router>
  );
};

const SequenceDisplay = ({ sequences }) => {
  return (
    <div className="sequence-container">
      {splitSequences(sequences).map((chunk, chunkIndex) => (
        <div key={chunkIndex} className="sequence-chunk">
          {chunk.map((seq, index) => (
            <div key={index} className="sequence">
              <Link to={`/${seq.label.replace(/\s+/g, '-')}`}>
                <div className="sequence-label">{seq.label}</div>
              </Link>
              <div className="sequence-boxes">
                {seq.lines.map((line, lineIndex) => (
                  <div key={lineIndex} className="sequence-line">
                    {line.split('').map((char, charIndex) => (
                      <div key={charIndex} className={`sequence-box ${char === '-' ? 'gap' : ''}`}>
                        {char}
                      </div>
                    ))}
                  </div>
                ))}
              </div>
            </div>
          ))}
          <div className="sequence-indexes">
            {Array.from({ length: Math.ceil(chunk[0].lines[0].length / 10) }, (_, i) => (i + 1) * 10).map((num, i) => (
              <div key={i} className="sequence-index">
                {chunkIndex * 50 + num}
              </div>
            ))}
          </div>
        </div>
      ))}
    </div>
  );
};

const splitSequences = (sequences) => {
  const chunkSize = 50;
  const result = [];
  const numChunks = Math.ceil(sequences[0].sequence.length / chunkSize);

  for (let chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
    const chunk = sequences.map((seq) => ({
      label: seq.label,
      lines: [seq.sequence.slice(chunkIndex * chunkSize, (chunkIndex + 1) * chunkSize)]
    }));
    result.push(chunk);
  }

  return result;
};

const SequenceDetails = ({ label, sequence }) => {
  return (
    <div className="sequence-details">
      <h2>{label} 상세 정보 페이지 넘어가는 로직 연결하기 </h2>
    </div>
  );
};

export default ProteinSeq;
