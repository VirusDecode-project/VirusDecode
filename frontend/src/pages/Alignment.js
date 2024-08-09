import React, { useState, useRef } from 'react';
import { Bar } from 'react-chartjs-2';
import 'chart.js/auto';
import { BrowserRouter as Router, Route, Routes, Link, useNavigate } from 'react-router-dom';
import './ProteinSeq.css';



function Alignment({ data }) {
    const chartRef = useRef();
    const [selectedData, setSelectedData] = useState(null);

    const chartData = {
        labels: ['Total'],
        datasets: data.map((item) => ({
            label: item.label,
            data: [item.value],
            backgroundColor: item.color,
            categoryPercentage: 0.1, // 카테고리 전체에서 바가 차지하는 비율 조정
        })),
    };

    const options = {
        indexAxis: 'y',
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
                stacked: true,
                display: false, // x축 지우기
            },
            y: {
                stacked: true,
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
            duration: 0, 
        },
    };

    const handleClick = (event) => {
        const chart = chartRef.current;
        if (!chart) return;

        const points = chart.getElementsAtEventForMode(event, 'nearest', { intersect: true }, false);
        if (points.length) {
            const firstPoint = points[0];
            const label = chart.data.datasets[firstPoint.datasetIndex].label;
            const value = chart.data.datasets[firstPoint.datasetIndex].data[firstPoint.index];
            setSelectedData({ label, value });
        }
    };

    const renderComponent = () => {
        if (!selectedData) return null;

        switch (selectedData.label) {
            case 'Label 1':
                return <ComponentForLabel1 value={selectedData.value} />;
            case 'Label 2':
                return <ComponentForLabel2 value={selectedData.value} />;
            // Add cases for other labels
            default:
                return <ProteinSeq/>;
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
                                ctx.fillStyle = 'black'; // 텍스트 색상 설정
                                ctx.font = 'bold 12px Arial'; // 텍스트 폰트 설정
                                ctx.textAlign = 'center'; // 텍스트 정렬 설정
                                ctx.textBaseline = 'middle'; // 텍스트 기준선 설정

                                // 바의 중앙에 텍스트 그리기
                                ctx.fillText(dataset.label, centerX, centerY);

                                // 다음 바의 시작 X 좌표 업데이트
                                startX += barWidth;
                            });
                        },
                    },
                ]}
            />
            {renderComponent()} {/* 선택된 데이터에 따라 다른 컴포넌트를 렌더링 */}
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





// 세부 단백질 시퀀스 표시


const ProteinSeq = () => {
    const sequences = [
        {
            label: 'Reference',
            sequence: 'MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQMESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
        },
        {
            label: 'USA-WA1/2020',
            sequence: 'MESLHPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQMESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
        },
        {
            label: 'SARS-P2/2021',
            sequence: 'MESLHPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQMESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
        },
        {
            label: 'SH-P10-A-2-Shanghai/2020',
            sequence: 'MESLHPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQMESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQ'
        }
    ];

    const [selectedSequence, setSelectedSequence] = useState(null);

    const handleSequenceClick = (sequence) => {
        setSelectedSequence(sequence);
    };

    return (
        <div className="ProteinSeq">
            <SequenceDisplay sequences={sequences} onSequenceClick={handleSequenceClick} />
            <Routes>
                {sequences.map((seq, index) => (
                    <Route
                        key={index}
                        path={`/${seq.label.replace(/\s+/g, '-')}`}
                        element={<SequenceDetails label={seq.label} sequence={seq.sequence} />}
                    />
                ))}
            </Routes>
            {selectedSequence && (
                <NextButton selectedSequence={selectedSequence} />
            )}
        </div>
    );
};

const SequenceDisplay = ({ sequences, onSequenceClick }) => {
    return (
        <div className="sequence-container">
            {splitSequences(sequences).map((chunk, chunkIndex) => (
                <div key={chunkIndex} className="sequence-chunk">
                    {chunk.map((seq, index) => (
                        <div key={index} className="sequence" onClick={() => onSequenceClick(seq)}>
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
            <h2>{label} 상세 정보 페이지</h2>
            <p>{sequence}</p> {/* 여기에서 시퀀스를 보여줌 */}
        </div>
    );
};

const NextButton = ({ selectedSequence }) => {
    const navigate = useNavigate();

    const handleNextClick = () => {
        navigate(`/${selectedSequence.label.replace(/\s+/g, '-')}`);
    };

    return (
        <div className="next-button">
            <button onClick={handleNextClick}>{selectedSequence.label} - Next</button>
        </div>
    );
};

