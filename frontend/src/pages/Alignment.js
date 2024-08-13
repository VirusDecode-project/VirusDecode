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
            {renderComponent()} {/* 선택된 데이터에 따라 다른 컴포넌트 렌더링 */}
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


// 세부 단백질 시퀀스 데이터예시

const ProteinSeq = () => {
  const sequences = [
    {
      label: 'Reference',
      sequence: 'MESLNPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQMESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVATLEGIQ'
    },
    {
      label: 'USA-WA1/2020', 
      sequence: 'MESLHPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQMESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
    },
    {
      label: 'SARS-P2/2021',
      sequence: 'MESLHPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQMESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGoQ'
    },
    {
      label: 'SH-P10-A-2-Shanghai/2020', 
      sequence: 'MESLHPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQMESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
    }
  ];

  const [selectedSequence, setSelectedSequence] = useState(null); // 선택된 시퀀스를 상태로 관리

  const handleSequenceClick = (sequence) => {
    setSelectedSequence(sequence); // 선택된 시퀀스를 설정
  };

  const referenceSequence = sequences[0].sequence;  // 첫 번째 서열이 참조 서열

  return (
    <div className="ProteinSeq">
      <SequenceDisplay 
        sequences={sequences} 
        referenceSequence={referenceSequence}
        onSequenceClick={handleSequenceClick} 
        selectedSequence={selectedSequence} 
      />
      {/* NextButton이 selectedSequence가 설정될 때만 렌더링되도록 함 */}
      {selectedSequence && (
        <NextButton selectedSequence={selectedSequence} sequences={sequences} />
      )}
    </div>
  );
};

const SequenceDisplay = ({ sequences, referenceSequence, onSequenceClick, selectedSequence }) => {

  useEffect(() => {
    const labels = document.querySelectorAll('.sequence-label');
    
    let maxWidth = 0;
    labels.forEach(label => {
      const labelWidth = label.offsetWidth;
      if (labelWidth > maxWidth) {
        maxWidth = labelWidth + 10;
      }
    });

    labels.forEach(label => {
      label.style.width = `${maxWidth}px`;
    });

    const indexesContainers = document.querySelectorAll('.sequence-indexes');
    indexesContainers.forEach(container => {
      container.style.marginLeft = `${maxWidth}px`;
    });

  }, [sequences]);

  return (
    <div className="sequence-container">
      {splitSequences(sequences, referenceSequence).map((chunk, chunkIndex) => (
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
                {num <= referenceSequence.length ? num : ''}
              </div>
            ))}
          </div>
        </div>
      ))}
    </div>
  );
};

const splitSequences = (sequences, referenceSequence) => {
  const chunkSize = 50; // 서열을 50개씩 나눔
  const result = [];
  const numChunks = Math.ceil(referenceSequence.length / chunkSize); // 서열을 나눌 덩어리 수 계산

  for (let chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
    const chunk = sequences.map((seq) => ({
      label: seq.label,
      lines: [
        seq.sequence
          .slice(chunkIndex * chunkSize, (chunkIndex + 1) * chunkSize)
          .split('')
          .map((char, index) => ({
            char,
            different: char !== referenceSequence[chunkIndex * chunkSize + index] // 참조 서열과 다른지 여부
          }))
      ]
    }));
    result.push(chunk);
  }

  return result;
};


const NextButton = ({ selectedSequence, sequences }) => {
  const currentIndex = sequences.findIndex(seq => seq.label === selectedSequence.label); // 현재 선택된 시퀀스 인덱스 찾기
  const nextIndex = (currentIndex + 1) % sequences.length; // 다음 시퀀스 인덱스 계산
  const nextSequence = sequences[nextIndex]; // 다음 시퀀스 가져오기

  return (
    <div className="next-button">
      <button >
        <Link to={`/${nextSequence.label.replace(/\s+/g, '-')}`}>
          Convert {nextSequence.label} {/* 다음 시퀀스로 이동 */}
        </Link>
      </button>
    </div>
  );
};
