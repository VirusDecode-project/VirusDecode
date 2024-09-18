import React, {MouseEvent} from 'react';
import HistoryItem from './HistoryItem';

interface HistoryListProps {
  history: string[];
  handleHistoryClick: (index: number) => void;
  handleEllipsisClick: (e: MouseEvent<HTMLButtonElement>, index: number) => void;
}

const HistoryList: React.FC<HistoryListProps> = ({ history, handleHistoryClick, handleEllipsisClick }) => {
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
