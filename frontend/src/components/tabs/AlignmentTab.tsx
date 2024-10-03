import React, { Dispatch, SetStateAction, useState, useEffect } from 'react';
import '../../styles/Alignment.css';
import Loading from '../Loading';
import StackedBar from '../StackedBar';
import ProteinSeq from '../ProteinSeq';
import {AlignmentData} from '../types';
import {MRNAData} from '../types';

interface ChartDataItem {
  label: string;
  value: number;
  color: string;
  start: number;
  end: number;
}

interface AlignmentProps {
  alignmentData: AlignmentData;
  setTab: Dispatch<SetStateAction<number>>;
  onRegionUpdate: (region: string) => void;
  setMRNAReceived: Dispatch<SetStateAction<boolean>>;
  setPDBReceived: Dispatch<SetStateAction<boolean>>;
  workingHistory: string;
  setLinearDesignData:Dispatch<SetStateAction<MRNAData | null>>;
  setPDBids: Dispatch<SetStateAction<string[]>>;
  setPDBInfo: Dispatch<SetStateAction<string[]>>;
  setSelectedPDBid: Dispatch<SetStateAction<string>>;
}

let lastHue = 0;

const generatePastelColor = () => {
  const minDifference = 60;
  let hue;
  do {
    hue = Math.floor(Math.random() * 360);
  } while (Math.abs(hue - lastHue) < minDifference);

  lastHue = hue;
  const saturation = 70 + Math.floor(Math.random() * 10);
  const lightness = 85 + Math.floor(Math.random() * 10);
  return `hsl(${hue}, ${saturation}%, ${lightness}%)`;
};

const Alignment: React.FC<AlignmentProps> = ({ alignmentData, setTab, onRegionUpdate, setMRNAReceived, setPDBReceived, workingHistory, setLinearDesignData, setPDBids, setPDBInfo, setSelectedPDBid }) => {
  const [chartData, setChartData] = useState<ChartDataItem[]>([]);
  const [selectedRegion, setSelectedRegion] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  
  useEffect(() => {
    if (alignmentData) {
      const firstRegion = Object.keys(alignmentData.alignment_index).length > 0 ? Object.keys(alignmentData.alignment_index)[0] : '';
      setSelectedRegion(firstRegion);
      try {
        const data = Object.entries(alignmentData.alignment_index).map(([label, [start, end]]) => {
          const value = end - start;
          return {
            label,
            value,
            color: generatePastelColor(),
            start,
            end,
          };
        });
        setChartData(data);
      } catch (error) {
        console.error('Error processing response data:', error);
      }
    }
  }, [alignmentData]);


  

  if (!alignmentData.alignment_index[selectedRegion]) {
    return (
      <div>
        <Loading text="Loading" />
      </div>
    );
  }
  

  return (
    <div>
      {isLoading ? (
        <Loading text="Converting" />
      ) : (
        <div>
          <div className="stacked-bar-chart">
            <StackedBar data={chartData} onBarClick={setSelectedRegion} />
          </div>
          {alignmentData && (
            <ProteinSeq
              selectedRegion={selectedRegion}
              setSelectedRegion={setSelectedRegion}
              alignmentData={alignmentData}
              setTab={setTab}
              setIsLoading={setIsLoading}
              onRegionUpdate={onRegionUpdate}
              setMRNAReceived={setMRNAReceived}
              setPDBReceived={setPDBReceived}
              workingHistory={workingHistory}
              setLinearDesignData={setLinearDesignData}
              setPDBids={setPDBids}
              setPDBInfo={setPDBInfo}
              setSelectedPDBid={setSelectedPDBid}
              isLoading={isLoading}
            />
          )}
        </div>
      )}
    </div>
  );
}

export default Alignment;