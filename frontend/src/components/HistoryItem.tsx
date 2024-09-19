import React, {MouseEvent} from 'react';

interface HistoryItemProps{
  item: string;
  index: number;
  handleHistoryClick: (index: number) => void;
  handleEllipsisClick: (e: MouseEvent<HTMLButtonElement>, index: number) => void;
}

const HistoryItem: React.FC<HistoryItemProps> = ({ item, index, handleHistoryClick, handleEllipsisClick }) => {
  return (
    <div
      className="history-item"
      onClick={() => handleHistoryClick(index)}
      style={{ cursor: 'pointer' }}
    >
      <span>{item}</span>
      <div className="ellipsis-container">
        <button
          className="ellipsis-button"
          onClick={(e) => handleEllipsisClick(e, index)} // Pass event and index
        >
          ...
        </button>
      </div>
    </div>
  );
};

export default HistoryItem;
