import React, { useState } from 'react';
import './Tabs.css';
import Sequence from './tab_pages/Sequence';
import Gene from './tab_pages/Gene';
import Mutation from './tab_pages/Mutation';
import Pathogenecity from './tab_pages/Pathogenecity';
import Analyze from './tab_pages/Analyze';


const Tabs = () => {
  const [activeTab, setActiveTab] = useState('Sequence');

  const handleTabClick = (tab) => {
    setActiveTab(tab);
  };


  // 서버에서 annotation 결과 가져오기, 색상 랜덤 설정하기 
  const [data, setData] = useState([
    { label: 'ORF1ab', value: 50, color: 'rgba(144, 238, 144, 0.6)' },
    { label: 'S', value: 20, color: 'rgba(255, 99, 132, 0.6)' },
    { label: 'ORF3a', value: 10, color: 'rgba(255, 182, 193, 0.6)' },
    { label: 'E', value: 5, color: 'rgba(255, 255, 102, 0.6)' },
    { label: 'M', value: 10, color: 'rgba(135, 206, 235, 0.6)' },
    { label: 'ORF7b', value: 5, color: 'rgba(255, 160, 122, 0.6)' },
]);


  const renderContent = () => {
    switch (activeTab) {
      case 'Sequence':
        return <Sequence data={data}/>;
      case 'Gene':
        return <Gene/>;
      case 'Mutation':
        return <Mutation/>;
      case 'Pathogenecity':
        return <Pathogenecity/>;
      case 'Analyze':
        return <Analyze/>;
      default:
        return <div></div>;
    }
  };

  return (
    <div>
      <div className="tab-bar">
        {['Sequence', 'Gene', 'Mutation', 'Pathogenecity', 'Analyze'].map((tab) => (
          <div
            key={tab}
            className={`tab ${activeTab === tab ? 'active' : ''}`}
            onClick={() => handleTabClick(tab)}
          >
            {tab}
          </div>
        ))}
      </div>
      <div className="tab-content">
        {renderContent()}
      </div>
    </div>
  );
};

export default Tabs;
