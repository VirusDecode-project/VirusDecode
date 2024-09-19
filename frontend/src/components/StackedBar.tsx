import React, { useState, useEffect } from 'react';
import '../styles/Alignment.css';

interface DataItem {
  label: string;
  value: number;
  color: string;
  start: number; // Tooltip에 사용되는 start index
  end: number;   // Tooltip에 사용되는 end index
}

interface StackedBarProps {
  data: DataItem[];
  onBarClick: (label: string) => void;
}

const StackedBar: React.FC<StackedBarProps> = ({ data, onBarClick }) => {
  const [tooltip, setTooltip] = useState({ visible: false, text: '', x: 0, y: 0 });

  const handleMouseOver = (e: React.MouseEvent<HTMLDivElement>, item: DataItem) => {
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

export default StackedBar;



