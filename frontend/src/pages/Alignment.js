
import React, { useState, useRef, useEffect } from 'react';
import { Bar } from 'react-chartjs-2';
import 'chart.js/auto';
import { Link } from 'react-router-dom';
import './ProteinSeq.css';

function Alignment({ data }) {
    const chartRef = useRef(); // 차트에 대한 참조
    const [selectedData, setSelectedData] = useState(null);

    const chartData = {
        labels: ['Total'], // 하나의 레이블을 사용해 모든 데이터 표시
        datasets: data.map((item) => ({
            label: item.label,
            data: [item.value],
            backgroundColor: item.color,
            categoryPercentage: 0.1, // 카테고리 전체에서 바가 차지하는 비율 조정
        })),
    };

    const options = {
        indexAxis: 'y', // 차트를 가로 방향으로 설정
        layout: {
            padding: {
                top: 0, 
                right: 0,
                bottom: 0,
                left: 0,
            },
        },
        scales: {
            x: {
                stacked: true, // x축 데이터를 스택으로 쌓음
                display: false, // x축 지우기
            },
            y: {
                stacked: true, // y축 데이터를 스택으로 쌓음
                display: false, // y축 지우기
            },
        },
        plugins: {
            legend: {
                display: false,
            },
            tooltip: {
                enabled: false, // 툴팁 비활성화
            },
        },
        elements: {
            bar: {
                borderWidth: 0, // 바 테두리 지우기
            },
        },
        animation: {
            duration: 0, // 애니메이션 비활성화
        },
    };

    const handleClick = (event) => {
        const chart = chartRef.current; // 차트 참조를 가져옴
        if (!chart) return;

        const points = chart.getElementsAtEventForMode(event, 'nearest', { intersect: true }, false);
        if (points.length) {
            const firstPoint = points[0];
            const label = chart.data.datasets[firstPoint.datasetIndex].label;
            const value = chart.data.datasets[firstPoint.datasetIndex].data[firstPoint.index];
            setSelectedData({ label, value }); // 선택된 데이터를 상태에 저장
        }
    };

    const renderComponent = () => {
        if (!selectedData) return null; // 선택된 데이터가 없으면 아무것도 렌더링하지 않음

        switch (selectedData.label) {
            case 'Label 1':
                return <ComponentForLabel1 value={selectedData.value} />;
            case 'Label 2':
                return <ComponentForLabel2 value={selectedData.value} />;
            // 다른 레이블에 대한 케이스 추가하기
            default:
                return <ProteinSeq/>; // 디폴트는 ProteinSeq 컴포넌트 렌더링
        }
    };

    return (
      <div>
        <div className="chart-container">
            <Bar
                ref={chartRef}
                data={chartData}
                options={options}
                onClick={handleClick} // 클릭 이벤트 핸들러 추가
                plugins={[
                    {
                        id: 'custom-datalabels',
                        afterDatasetsDraw: (chart) => {
                            const ctx = chart.ctx; // 차트의 2D 캔버스 렌더링 컨텍스트 가져오기
                            const datasets = chart.data.datasets; // 차트의 데이터셋 가져오기
                            const meta = chart.getDatasetMeta(0); // 첫 번째 데이터셋의 메타데이터 가져오기
                            const bar = meta.data[0]; // 첫 번째 데이터셋의 바 메타데이터 가져오기

                            // 데이터셋의 총 너비 계산
                            let totalWidth = 0;
                            datasets.forEach((dataset, i) => {
                                const meta = chart.getDatasetMeta(i);
                                const bar = meta.data[0];
                                totalWidth += bar.width;
                            });

                            // 각 데이터셋의 중앙에 텍스트 배치
                            let startX = bar.x - totalWidth / 2; // 전체 바의 시작 X 좌표
                            datasets.forEach((dataset, i) => {
                                const value = dataset.data[0];
                                const meta = chart.getDatasetMeta(i);
                                const bar = meta.data[0];
                                const barWidth = bar.width; // 각 바의 너비
                                const centerX = startX + barWidth / 2; // 바의 중앙 X 좌표 계산
                                const centerY = bar.y; // 바의 중앙 Y 좌표 계산

                                // 텍스트 스타일 설정
                                ctx.fillStyle = 'black'; 
                                ctx.font = 'bold 12px Arial'; 
                                ctx.textAlign = 'center'; 
                                ctx.textBaseline = 'middle'; 

                                // 바의 중앙에 텍스트 그리기
                                ctx.fillText(dataset.label, centerX, centerY);

                                // 다음 바의 시작 X 좌표 업데이트
                                startX += barWidth;
                            });
                        },
                    },
                ]}
            />

                {selectedData && renderComponent(selectedData)}
            </div>

            <ProteinSeq />
        </div>
    );
}

const ComponentForLabel1 = ({ value }) => (
    <div>
        <h2>Label 1 Data</h2>
        <p>Value: {value}</p>
    </div>
);

const ComponentForLabel2 = ({ value }) => (
    <div>
        <h2>Label 2 Data</h2>
        <p>Value: {value}</p>
    </div>
);

export default Alignment;


const ProteinSeq = () => {
  const [sequences, setSequences] = useState([]);
  const [selectedSequence, setSelectedSequence] = useState(null);
  const [alignmentIndex, setAlignmentIndex] = useState({});
  const [selectedRegion, setSelectedRegion] = useState("ORF1ab");

  useEffect(() => {
    // JSON 데이터 가져오기
    fetch('../../public/alignment_data.json')  // JSON 파일의 경로 설정
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

  if (sequences.length === 0) return <p>Loading sequences...</p>; // 서열 데이터 로드 중 표시

  const referenceSequence = sequences[0].sequence;  // 첫 번째 서열을 참조 서열로 사용
  const regionIndices = alignmentIndex[selectedRegion]; // 선택된 지역의 시작 및 끝 인덱스
  const regionSequence = referenceSequence.slice(regionIndices[0], regionIndices[1]); // 참조 서열의 해당 영역

  const handleSequenceClick = (sequence) => {
    setSelectedSequence(sequence); // 선택된 서열 설정
  };

  const handleRegionChange = (event) => {
    setSelectedRegion(event.target.value); // 선택된 영역 변경
  };

  return (
    <div className="protein-sequence-container">
      <div className="region-selector">
        <label htmlFor="region-select">Select Protein Region: </label>
        <select id="region-select" value={selectedRegion} onChange={handleRegionChange}>
          {Object.keys(alignmentIndex).map(region => (
            <option key={region} value={region}>
              {region}
            </option>
          ))}
        </select>
      </div>
      <SequenceDisplay 
        sequences={sequences} 
        referenceSequence={regionSequence}
        onSequenceClick={handleSequenceClick} 
        selectedSequence={selectedSequence}
        regionIndices={regionIndices}
      />
      {selectedSequence && (
        <NextButton selectedSequence={selectedSequence} sequences={sequences} />
      )}
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
            {Array.from({ length: 5 }, (_, i) => (chunkIndex * 50) + ((i + 1) * 10)).map((num, i) => (
              <div key={i} className="sequence-index">
                {num <= referenceSequence.length ? num : '---'}
              </div>
            ))}
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
  const numChunks = Math.ceil(referenceSequence.length / chunkSize); // 청크 수 계산

  for (let chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
    const chunk = sequences.map((seq) => ({
      label: seq.label,
      lines: [
        seq.sequence
          .slice(regionIndices[0] + chunkIndex * chunkSize, regionIndices[0] + (chunkIndex + 1) * chunkSize)
          .split('')
          .map((char, index) => ({
            char,
            different: char !== referenceSequence[chunkIndex * chunkSize + index] // 참조 서열과 다른지 확인
          }))
      ]
    }));
    result.push(chunk);
  }

  return result;
};

const NextButton = ({ selectedSequence, sequences }) => {
  const currentIndex = sequences.findIndex(seq => seq.label === selectedSequence.label); // 현재 선택된 서열의 인덱스 찾기
  const nextIndex = (currentIndex + 1) % sequences.length; // 다음 서열 인덱스 계산
  const nextSequence = sequences[nextIndex];

  return (
    <div className="next-button">
      <button>
        <Link to={`/${nextSequence.label.replace(/\s+/g, '-')}`}>
          Convert {nextSequence.label}
        </Link>
      </button>
    </div>
  );
};
