import React, { MouseEvent, useEffect, useState, Dispatch, SetStateAction } from 'react';
import Viztein from 'viztein';
import "../../styles/Render3D.css";

interface Render3DProps {
  region: string;
  PDBids: string[];
  PDBInfo: string[];
  selectedPDBid: string;
  setSelectedPDBid: Dispatch<SetStateAction<string>>;
}



const Render3D: React.FC<Render3DProps> = ({ region, PDBids, PDBInfo, selectedPDBid, setSelectedPDBid }) => {
  const [representation, setRepresentation] = useState("default");
  const [error, setError] = useState("");
  const [tooltip, setTooltip] = useState({ visible: false, text: '', x: 0, y: 0 });
  const [viewportStyle, setViewportStyle] = useState({ width: '60vw', height: '70vh' });

  useEffect(() => {
    // 탭에 들어갈 때 스크롤 비활성화
    document.body.style.overflow = 'hidden';

    return () => {
      // 탭에서 나올 때 스크롤 다시 활성화
      document.body.style.overflow = 'auto';
    };
  }, []);
  // useEffect(() => {
  //   // 화면 크기 감지
  //   const updateViewportStyle = () => {
  //     if (window.innerWidth <= 480) {
  //       setViewportStyle({ width: '320px', height: '240px' }); // 모바일 크기
  //     } else {
  //       setViewportStyle({ width: '1080px', height: '720px' }); // 기본 크기
  //     }
  //   };

  //   // 처음 로드될 때와 창 크기 변경 시 업데이트
  //   updateViewportStyle();
  //   window.addEventListener('resize', updateViewportStyle);

  //   return () => {
  //     window.removeEventListener('resize', updateViewportStyle);
  //   };
  // }, []);

  const checkPDBFileExists = async (url: string) => {
    try {
      const response = await fetch(url, { method: 'HEAD' });
      return response.ok;
    } catch (error) {
      console.error('Error checking PDB file existence:', error);
      return false;
    }
  };

  if (PDBids.length === 0) {
    return <div>Loading...</div>;
  }
  // const viewportStyle = {
  //   width: '50vw',  // Viewport width의 100%로 설정
  //   height: '80vh',  // Viewport height의 70%로 설정
  // };
  const refData = {
    filename: `https://files.rcsb.org/download/${selectedPDBid}.pdb`,
    ...(representation !== "default" && {
      config: [{
        type: 'addRepresentation',
        input: representation
      }]
    })
  };


  const handlePDBSelection = async (id: string) => {
    const pdbUrl = `https://files.rcsb.org/download/${id}.pdb`;
    const exists = await checkPDBFileExists(pdbUrl);
    if (exists) {
      setSelectedPDBid(id);
      setError("");
    } else {
      setError(`PDB file ${id}.pdb does not exist.`);
    }
  };

  const handleMouseOver = (e: MouseEvent<HTMLDivElement>, info: string) => {
    setTooltip({
      visible: true,
      text: `${info}`,
      x: e.clientX,
      y: e.clientY,
    });
  };

  const handleMouseOut = () => {
    setTooltip({ visible: false, text: '', x: 0, y: 0 });
  };

  const vizteinKey = `${selectedPDBid}-${representation}`;

  return (
    <div className='reference3D'>
      <div className='left-column'>
          {/* <div className='renderTitle'>{region}</div> */}
        <div className='representation-container'>
          <select
            className="style"
            value={representation}
            onChange={(e) => setRepresentation(e.target.value)}
          >
            <option value="default"> Default </option>
            <option value="ball+stick"> Ball + Stick </option>
            <option value="cartoon"> Cartoon </option>
          </select>
          <Viztein key={vizteinKey} data={refData} viewportStyle={viewportStyle} />
          <p className='PDBname'>{selectedPDBid}.pdb</p>
        </div>
      </div>
      <div className='right-column'>
        <div className='list-header'>
          <span className='header-item'>PDB ID</span>
        </div>
        <div className='list-container'>
          {PDBids.map((id, index) => (
            <div
              key={id}
              className="list-row"
              onClick={() => handlePDBSelection(id)}
              onMouseOver={(e) => handleMouseOver(e, PDBInfo[index])}
              onMouseOut={handleMouseOut}
            >
              <div className="list-item">{id}</div>
            </div>
          ))}
        </div>
        {tooltip.visible && (
          <div
            className="custom-tooltip"
            style={{ top: tooltip.y + 10, left: tooltip.x + 10 }}
          >
            {tooltip.text}
          </div>
        )}
        {error && <div className="error-message">{error}</div>}


      </div>
    </div>
  );
}

export default Render3D;
