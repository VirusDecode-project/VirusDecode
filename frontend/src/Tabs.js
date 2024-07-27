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

  const renderContent = () => {
    switch (activeTab) {
      case 'Sequence':
        return <Sequence/>;
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
