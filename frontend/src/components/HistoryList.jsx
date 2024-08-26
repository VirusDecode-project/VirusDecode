// HistoryList.jsx
import React from 'react';
import HistoryItem from './HistoryItem';

const HistoryList = ({ history, handleHistoryClick, handleEllipsisClick }) => {
  return (
    <div className="history-list">
      {history.map((item, index) => (
        <HistoryItem
          key={index}
          item={item}
          index={index}
          handleHistoryClick={handleHistoryClick}
          handleEllipsisClick={handleEllipsisClick}
        />
      ))}
    </div>
  );
};

export default HistoryList;
