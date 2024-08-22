import Viztein from 'viztein';
import React, { useEffect, useState } from 'react';
import "./Render3D.css";

function Render3D() {
  //const [PDBids, setPDBids] = useState(['8VCI', '8UYS', '7O7Y', '7O7Z', '7O81', '7O80']); //test
  const [PDBids, setPDBids] = useState([]); //test2 (backend)
  const [PDBscores, setPDBscores] = useState([]);
  /*
      GK
      1. selectedPDBid: PDBids 중 첫 번째 값으로 초기화 필요합니다.
      2. RCSB PDB API endpoint에서 받아온 json 값에 socre도 있던데, 함께 받아와서 아래에 표시해주면 좋을 것 같습니다.
      KY
      --> 완료
  */
  const [selectedPDBid, setSelectedPDBid] = useState("");
  const [representation, setRepresentation] = useState("default");

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
        setPDBscores(values); // value들만 저장
        // 리스트의 첫 번째 값으로 PDB id 초기화
        if (keys.length > 0) {
          setSelectedPDBid(keys[0]);
        }
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
    filename: `https://files.rcsb.org/download/${selectedPDBid}.pdb`, // 에러처리 필요
    ...(representation !== "default" && { // style 지정
      config: [{
        type: 'addRepresentation',
        input: representation
      }]
    })
  };

  const vizteinKey = `${selectedPDBid}-${representation}`; // PDBid 또는 style 변경 시 렌더링 업데이트에 필요한 key

  return (
    <div className='reference3D'>
      <div className='left-column'>
        <h4>Reference</h4>
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
        <h5>{selectedPDBid}.pdb</h5>
      </div>
      <div className='right-column'>
        <div className='list-header'>
          <span className='header-item'>PDB ID</span>
          <span className='header-item'>Matching Score</span>
        </div>
        <ul>
          {PDBids && PDBids.map((id, index) => (
            <li key={id} onClick={() => setSelectedPDBid(id)}>
              <span className='list-item'>{id}</span>
              <span className='list-item'>{PDBscores[index]}</span>
            </li>
          ))}
        </ul>
      </div>
    </div>
  );
}
export default Render3D;
