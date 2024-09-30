import React, { useEffect, Dispatch, SetStateAction, useState } from 'react';
import { Nav } from 'react-bootstrap';
import '../styles/Analysis.css';
import Alignment from '../components/tabs/AlignmentTab';
import MRNAdesign from '../components/tabs/MRNAdesignTab';
import Render3D from '../components/tabs/Render3DTab';
import { MRNAData, AlignmentData } from '../components/types';
import { useNavigate } from "react-router-dom";

interface AnalysisProps {
  tab: number;
  setTab: Dispatch<SetStateAction<number>>;
  mRNAReceived: boolean;
  setMRNAReceived: Dispatch<SetStateAction<boolean>>;
  PDBReceived: boolean;
  setPDBReceived: Dispatch<SetStateAction<boolean>>;
  workingHistory: string;
  setWorkingHistory: Dispatch<SetStateAction<string>>;
  linearDesignData: MRNAData | null;
  setLinearDesignData:Dispatch<SetStateAction<MRNAData | null>>;
  PDBids: string[];
  setPDBids: Dispatch<SetStateAction<string[]>>;
  PDBInfo: string[];
  setPDBInfo: Dispatch<SetStateAction<string[]>>;
  selectedPDBid: string;
  setSelectedPDBid: Dispatch<SetStateAction<string>>;
  alignmentData: AlignmentData;
  setHistory: Dispatch<SetStateAction<string[]>>;
  setAlignmentData: Dispatch<SetStateAction<AlignmentData>>;
}

const Analysis: React.FC<AnalysisProps> = ({ tab, setTab, mRNAReceived, setMRNAReceived, PDBReceived, setPDBReceived, workingHistory, setWorkingHistory, linearDesignData, setLinearDesignData, PDBids, setPDBids, PDBInfo, setPDBInfo, selectedPDBid, setSelectedPDBid, alignmentData, setHistory, setAlignmentData }) => {
  const [modalRegion, setModalRegion] = useState('');
  const handleModalRegion = (region: string) => {
    setModalRegion(region);
  };
  useEffect(() => {
    const fetchHistory = async () => {
      try {
        const serverResponse = await fetch(`${process.env.REACT_APP_BACKEND_URL}/history/list`, {
          method: 'GET',
          credentials: 'include',
        });
        if (!serverResponse.ok) {
          throw new Error("Failed to fetch history list");
        }
        const responseData = await serverResponse.json();
        setHistory(responseData);
      } catch (error) {
        console.error("Error fetching history:", error);
      }
    };

    fetchHistory();
  },[]);

  return (
    <div>
      <div className={"analysis-container"}>
        <>
          <Nav variant="tabs" defaultActiveKey="link0" className="justify-content-start">
            <Nav.Item>
              <Nav.Link eventKey="link0" active={tab === 0} onClick={() => setTab(0)}>Alignment</Nav.Link>
            </Nav.Item>
            <Nav.Item>
              <Nav.Link eventKey="link1" active={tab === 1} onClick={() => setTab(1)} disabled={!mRNAReceived}>mRNA design</Nav.Link>
            </Nav.Item>
            <Nav.Item>
              <Nav.Link eventKey="link2" active={tab === 2} onClick={() => setTab(2)} disabled={!PDBReceived}>3D viewer</Nav.Link>
            </Nav.Item>
          </Nav>
          <div>
            {tab === 0 && (
              <Alignment
                alignmentData={alignmentData}
                setMRNAReceived={setMRNAReceived}
                setPDBReceived={setPDBReceived}
                setTab={setTab}
                onRegionUpdate={handleModalRegion}
                workingHistory={workingHistory}
                setLinearDesignData={setLinearDesignData}
                setPDBids={setPDBids}
                setPDBInfo={setPDBInfo}
                setSelectedPDBid={setSelectedPDBid}
              />
            )}
            {tab === 1 && <MRNAdesign workingHistory={workingHistory} linearDesignData={linearDesignData}/>}
            {tab === 2 && <Render3D region={modalRegion} PDBids={PDBids} PDBInfo={PDBInfo} selectedPDBid={selectedPDBid} setSelectedPDBid={setSelectedPDBid} />}
          </div>
        </>
      </div>
    </div>
  );
}

export default Analysis;
