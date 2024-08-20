import React, { useState, useRef, useEffect } from 'react';
import { Bar } from 'react-chartjs-2';
import 'chart.js/auto';
import './ProteinSeq.css';
import Modal from './Modal';
import { Link } from 'react-router-dom'; 

const generatePastelColor = () => {
  const hue = Math.floor(Math.random() * 360);
  const saturation = 70 + Math.floor(Math.random() * 10);
  const lightness = 85 + Math.floor(Math.random() * 10);
  return `hsl(${hue}, ${saturation}%, ${lightness}%)`;
};

function Alignment() {
    const chartRef = useRef(); // 차트에 대한 참조
    const [chartData, setChartData] = useState({ labels: [], datasets: [] });
    const [selectedRegion, setSelectedRegion] = useState("ORF1ab");

    useEffect(() => {
        fetch('/alignment_data.json')  // JSON 파일의 경로 설정
            .then(response => response.json())
            .then(jsonData => {
                const { alignment_index } = jsonData;

                // alignment_index 데이터를 사용하여 datasets 생성
                const datasets = Object.entries(alignment_index).map(([label, [start, end]]) => {
                    const value = end - start;  // 서열 길이 계산
                    let color = generatePastelColor();

                    return {
                        label,
                        data: [value],
                        backgroundColor: color,
                        categoryPercentage: 0.1, // 카테고리 전체에서 바가 차지하는 비율 조정
                        start,
                        end,
                    };
                });

                setChartData({
                    labels: ['Total'], // 하나의 레이블을 사용해 모든 데이터 표시
                    datasets: datasets
                });
            })
            .catch(error => console.error('Error fetching sequence data:', error));
    }, []);

    const options = {
        indexAxis: 'y', // 차트를 가로 방향으로 설정
        layout: {
            padding: {
                top: 0,
                right: 0,
                bottom: 30, // 스택바 아래에 텍스트를 표시할 공간 확보
                left: 0,
            },
        },
        scales: {
            x: {
                stacked: true, // x축 데이터를 스택으로 쌓음
                display: false, // x축 숨기기
            },
            y: {
                stacked: true, // y축 데이터를 스택으로 쌓음
                display: false, // y축 숨기기
            },
        },
        plugins: {
          legend: {
              display: false,
          },
          tooltip: {
              enabled: true, // 툴팁 활성화
              position: 'nearest', // 툴팁 위치 기본값 사용
              callbacks: {
                  title: () => '', // 타이틀을 빈 배열로 반환하여 제거
                  label: (tooltipItem) => {
                      const datasetLabel = tooltipItem.dataset.label || '';
                      const start_num = tooltipItem.dataset.start || 0;
                      const end_num = tooltipItem.dataset.end || 0;

                      return `${datasetLabel}: ${start_num}~${end_num}`;
                  }
              }
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
            setSelectedRegion(label); // 선택된 영역 상태 업데이트
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
                            const ctx = chart.ctx;
                            const datasets = chart.data.datasets;
              
                            if (!datasets || datasets.length === 0) return;

                            datasets.forEach((dataset, i) => {
                                const meta = chart.getDatasetMeta(i);
                                const bar = meta.data[0]; // 각 데이터셋의 바 메타데이터 가져오기
              
                                if (!bar) return;
                                const centerX = bar.x - bar.width / 2; // 바의 중앙 X 좌표 계산
                                const centerY = bar.y; // 바의 중앙 Y 좌표 계산
                                const barWidth = bar.width;

                                ctx.fillStyle = 'black'; 
                                ctx.font = 'bold 12px Arial'; 
                                ctx.textAlign = 'center'; 
                                ctx.textBaseline = 'middle'; 
                                if (barWidth > 30) {
                                    ctx.fillText(dataset.label, centerX, centerY);
                                }
                            });
                        },
                    },
                ]}
            />
        </div>

        <ProteinSeq selectedRegion={selectedRegion} setSelectedRegion={setSelectedRegion} />  {/* setSelectedRegion 전달 */}
      </div>
    );
}

export default Alignment;


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


const NextButton = ({ selectedSequence }) => {
  if (!selectedSequence) return null;

  return (
    <div className="next-button">
      <button>
        <Link to={`/${selectedSequence.label.replace(/\s+/g, '-')}`}>
          Convert {selectedSequence.label}
        </Link>
      </button>
    </div>
  );
};