import React, { Dispatch, SetStateAction, useState, useEffect } from 'react';
import '../styles/Alignment.css';
import Loading from '../components/Loading';
import StackedBar from '../components/StackedBar';
import ProteinSeq from '../components/ProteinSeq';
import {ResponseData} from '../components/types';

interface ChartDataItem {
  label: string;
  value: number;
  color: string;
  start: number;
  end: number;
}

interface AlignmentProps {
  responseData: ResponseData;
  setTab: Dispatch<SetStateAction<number>>;
  onRegionUpdate: (region: string) => void;
  setMRNAReceived: Dispatch<SetStateAction<boolean>>;
  setPDBReceived: Dispatch<SetStateAction<boolean>>;
  workingHistory: string;
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

const Alignment: React.FC<AlignmentProps> = ({ responseData, setTab, onRegionUpdate, setMRNAReceived, setPDBReceived, workingHistory }) => {
  const [chartData, setChartData] = useState<ChartDataItem[]>([]);
  const [selectedRegion, setSelectedRegion] = useState('');
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    if (responseData) {
      const firstRegion = Object.keys(responseData.alignment_index)[0];
      setSelectedRegion(firstRegion);
      try {
        const data = Object.entries(responseData.alignment_index).map(([label, [start, end]]) => {
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
  }, [responseData]);

  if(chartData.length === 0) {
    return <div>Loading...</div>;
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
          {responseData && (
            <ProteinSeq
              selectedRegion={selectedRegion}
              setSelectedRegion={setSelectedRegion}
              responseData={responseData}
              setTab={setTab}
              setIsLoading={setIsLoading}
              onRegionUpdate={onRegionUpdate}
              setMRNAReceived={setMRNAReceived}
              setPDBReceived={setPDBReceived}
              workingHistory={workingHistory}
            />
          )}
        </div>
      )}
    </div>
  );
}

export default Alignment;