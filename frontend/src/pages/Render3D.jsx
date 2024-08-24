import React, { useEffect, useState } from 'react';
import Viztein from 'viztein';
import "./Render3D.css";
import helpIcon from '../image/help.png';
import helpImg1 from '../image/helpImg1.png'
import helpImg2 from '../image/helpImg2.png'
import helpImg3 from '../image/helpImg3.png'
import HelpModal from './HelpModal';

function Render3D() {
  //const [PDBids, setPDBids] = useState(['8VCI', '8UYS', '7O7Y', '7O7Z', '7O81', '7O80']); //test
  const [renderTitle, setRenderTitle] = useState("example");
  const [PDBids, setPDBids] = useState([]); //test2 (backend)
  /*
      GK
      1. selectedPDBid: PDBids 중 첫 번째 값으로 초기화 필요합니다.
      2. RCSB PDB API endpoint에서 받아온 json 값에 socre도 있던데, 함께 받아와서 아래에 표시해주면 좋을 것 같습니다.
      KY
      --> 완료
  */
  const [selectedPDBid, setSelectedPDBid] = useState("");
  const [PDBInfo, setPDBInfo] = useState([]);
  const [representation, setRepresentation] = useState("default");
  const [error, setError] = useState("");
  const [tooltip, setTooltip] = useState({ visible: false, text: '', x: 0, y: 0 });
  //------------------------------helpIcon-----------------------//
  const [isHelpModalOpen, setHelpModalOpen] = useState(false);
  //------------------------------helpIcon-----------------------//

  /*useEffect(() => {
    const fetchData = async () => {
      const url = 'url'; // backend url
      try {
        const response = await fetch(url);
        const data = await response.json();
        if (data.status === 'success') {
          setPDBids(data.PDBids); 
        }
      } catch (error) {
        console.error('Error fetching data:', error);
      }
    };
    fetchData();
  }, []);*/


  /* backend 추가 코드 시작 */
  useEffect(() => {
    // 서버로 GET 요청을 보내고, 응답을 받아와서 상태로 설정
    const fetchPDBids = async () => {
      try {
        const response = await fetch('http://localhost:8080/analysis/render3d'); // 서버의 엔드포인트로 요청
        if (!response.ok) {
          throw new Error(`Network response was not ok: ${response.statusText}`);
        }
        const PDBlist = await response.json(); // JSON 객체를 받아옴
        const keys = Object.keys(PDBlist);
        setPDBids(keys);
        // value 값만 리스트로 모아서 상태에 저장
        const values = Object.values(PDBlist);
        setPDBInfo(values); // value들만 저장
        // 리스트의 첫 번째 값으로 PDB id 초기화
        if (keys.length > 0) {
          setSelectedPDBid(keys[0]);
        }
        // 백엔드에서 가져온 값으로 title 지정
        // setRenderTitle(title);
      } catch (error) {
        console.error('Error fetching PDB IDs:', error);
      }
    };

    fetchPDBids(); // 데이터 가져오기
  }, []);
  /* backend 추가 코드 끝 */

  const viewportStyle = {
    width: '900px',
    height: '900px',
  };

  const refData = {
    filename: `https://files.rcsb.org/download/${selectedPDBid}.pdb`,
    ...(representation !== "default" && { // style 지정
      config: [{
        type: 'addRepresentation',
        input: representation
      }]
    })
  };

  const checkPDBFileExists = async (url) => {
    try {
      const response = await fetch(url, { method: 'HEAD' });
      return response.ok;
    } catch (error) {
      console.error('Error checking PDB file existence:', error);
      return false;
    }
  };

  const handlePDBSelection = async (id) => {
    const pdbUrl = `https://files.rcsb.org/download/${id}.pdb`;
    const exists = await checkPDBFileExists(pdbUrl);
    if (exists) {
      setSelectedPDBid(id);
      setError(""); // 파일이 존재하므로 에러 메시지를 지웁니다.
    } else {
      setError(`PDB file ${id}.pdb does not exist.`);
    }
  };

  const handleMouseOver = (e, info) => {
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

  const vizteinKey = `${selectedPDBid}-${representation}`; // PDBid 또는 style 변경 시 렌더링 업데이트에 필요한 key
  //------------------------------helpIcon-----------------------//
  const toggleHelpModal = () => {
    setHelpModalOpen(!isHelpModalOpen);
  };
  //------------------------------helpIcon-----------------------//
  return (
    <div className='reference3D'>
      <div className="help-icon-container">
        <img
          className="help-icon"
          src={helpIcon}
          onClick={toggleHelpModal}
          style={{ cursor: "pointer" }}
        />
      </div>
      <HelpModal
        isOpen={isHelpModalOpen}
        onClose={toggleHelpModal}
        content={
          <div>
            <p className='helpTitle'>mRNA Conversion Instructions</p>
            <p className='helpLevel'>1. Select CDS</p>
            <p className='helpContents'>Click on the CDS you wish to convert for mRNA design.</p>
            <img className='helpImageWide' src={helpImg1} />
            <p className='helpLevel'>2. Choose sublineage and input amino acid range</p>
            <p className='helpContents'>- Click on the appropriate sublineage from the initial input sequence.</p>
            <p className='helpContents'>- Enter the start and end positions for the amino acids.</p>
            <p className='helpWarnning'>*Note: Do not exceed 500 amino acids.</p>
            <img className='helpImage' src={helpImg2} />
            <p className='helpLevel'>3. Convert</p>
            <p className='helpContents'>Click the "Convert" button to proceed with the mRNA design.</p>
            <img className='helpImage' src={helpImg3} />
          </div>
        }
      />
      <div className='left-column'>
        <p className='renderTitle'>{renderTitle}</p>
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
        </div>
        <p className='PDBname'>{selectedPDBid}.pdb</p>
      </div>
      <div className='right-column'>
        <div className='list-header'>
          <span className='header-item'>PDB ID</span>
        </div>
        <ul>
          {PDBids && PDBids.map((id, index) => (
            <li
              key={id}
              onClick={() => handlePDBSelection(id)}
              onMouseOver={(e) => handleMouseOver(e, PDBInfo[index])} // 이곳에 PDB 정보가 들어감
              onMouseOut={handleMouseOut}
            >
              <span className='list-item'>{id}</span>
            </li>
          ))}
        </ul>
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