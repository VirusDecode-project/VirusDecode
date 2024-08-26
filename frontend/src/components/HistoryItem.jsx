// HistoryItem.jsx
import React from 'react';

const HistoryItem = ({ item, index, handleHistoryClick, handleEllipsisClick }) => {
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
