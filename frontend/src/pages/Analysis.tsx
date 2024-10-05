import React, { useEffect, Dispatch, SetStateAction, useState } from 'react';
import { Nav } from 'react-bootstrap';
import '../styles/Analysis.css';
import Alignment from '../components/tabs/AlignmentTab';
import MRNAdesign from '../components/tabs/MRNAdesignTab';
import Render3D from '../components/tabs/Render3DTab';
import { MRNAData, AlignmentData, PDBResponse } from '../components/types';
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
  setIsLoading: Dispatch<SetStateAction<boolean>>;
}

const Analysis: React.FC<AnalysisProps> = ({ tab, setTab, mRNAReceived, setMRNAReceived, PDBReceived, setPDBReceived, workingHistory, setWorkingHistory, linearDesignData, setLinearDesignData, PDBids, setPDBids, PDBInfo, setPDBInfo, selectedPDBid, setSelectedPDBid, alignmentData, setHistory, setAlignmentData, setIsLoading }) => {
  const [modalRegion, setModalRegion] = useState('');
  const handleModalRegion = (region: string) => {
    setModalRegion(region);
  };
  const checkPDBFileExists = async (url: string) => {
    try {
      const response = await fetch(url, { method: 'HEAD' });
      return response.ok;
    } catch (error) {
      console.error('Error checking PDB file existence:', error);
      return false;
    }
  };
  useEffect(() => {
    const fetchHistory = async () => {
      try {
        const historyListResponse = await fetch(`/api/history/list`, {
          method: 'GET',
          credentials: 'include',
        });
        if (!historyListResponse.ok) {
          throw new Error("Failed to fetch history list");
        }
        const responseHistoryListData = await historyListResponse.json();
        setHistory(responseHistoryListData);
        const requestData = { historyName: workingHistory };
        const serverResponse = await fetch(`/api/history/get`, {
          method: "POST",
          credentials: 'include',
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify(requestData),
        });
  
        if (!serverResponse.ok) {
          const errorMessage = await serverResponse.text();
          throw new Error(errorMessage);
        }
  
        const responseData = await serverResponse.json();
        console.log("History details fetched successfully: ", responseData);
  
        if (responseData.alignment) {
          const alignmentJson = JSON.parse(responseData.alignment);
          setAlignmentData(alignmentJson);
        } else {
          console.warn('Alignment data is null or undefined');
        }
  
        if (responseData.linearDesign) {
          const linearDesignJson = JSON.parse(responseData.linearDesign);
          setLinearDesignData(linearDesignJson);
          setMRNAReceived(true);
        } else {
          setMRNAReceived(false);
          setTab(0);
        }
  
        if (responseData.pdb) {
          const pdbJson: PDBResponse = JSON.parse(responseData.pdb);
          const keys = Object.keys(pdbJson);
          setPDBids(keys);
          const values = Object.values(pdbJson);
          setPDBInfo(values);
          if (keys.length > 0) {
            for (let i = 0; i < keys.length; i++) {
              const exist = await (checkPDBFileExists(`https://files.rcsb.org/download/${keys[i]}.pdb`))
              if (exist) {
                setSelectedPDBid(keys[i]);
                break;
              }
            }
          }
          setPDBReceived(true);
        } else {
          setPDBReceived(false);
        }


      } catch (error) {
        console.error("Error fetching history:", error);
      } finally{
        setIsLoading(false);
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
