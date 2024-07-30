import React from 'react';
import { Bar } from 'react-chartjs-2';
import 'chart.js/auto';

function Alignment({ data }) {
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
                enabled: true,
                callbacks: {
                    label: function (context) {
                        const label = context.dataset.label || '';
                        const value = context.raw || 0;
                        return `${label}: ${value}`;
                    },
                },
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

return (
    <div className="chart-container">
    <Bar
        data={chartData}
        options={options}
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
                  const centerY = bar.y  ; // 바의 중앙 Y 좌표 계산

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
    </div>
);
};

export default Alignment;