import React, { MouseEvent, useEffect, useState } from 'react';
import Viztein from 'viztein';
import "../../styles/Render3D.css";

interface Render3DProps {
  region: string;
}

interface PDBResponse {
  [key: string]: string;
}

const Render3D: React.FC<Render3DProps> = ({ region }) => {
  const [PDBids, setPDBids] = useState<string[]>([]);
  const [selectedPDBid, setSelectedPDBid] = useState("");
  const [PDBInfo, setPDBInfo] = useState<string[]>([]);
  const [representation, setRepresentation] = useState("default");
  const [error, setError] = useState("");
  const [tooltip, setTooltip] = useState({ visible: false, text: '', x: 0, y: 0 });
  /*
  const [scrollIndex, setScrollIndex] = useState(0);
  const itemsToShow = 3;
  const listRef = useRef(null);
*/
  useEffect(() => {
    const fetchPDBids = async () => {
      try {
        const serverResponse = await fetch('http://localhost:8080/analysis/re-render3d');
        if (!serverResponse.ok) {
          const errorMessage = await serverResponse.text();
          throw new Error(errorMessage);
        }
        const responseData: PDBResponse = await serverResponse.json();

        const keys = Object.keys(responseData);
        setPDBids(keys);
        const values = Object.values(responseData);
        setPDBInfo(values);

        if (keys.length > 0) {
          setSelectedPDBid(keys[0]);
        }
      } catch (error) {
        if (error instanceof Error){
          console.error("An error occurred during the request: ", error.message);
        }
      }
    };
    fetchPDBids();
  }, []);

  if (PDBids.length === 0) {
    return <div>Loading...</div>;
  }
  const viewportStyle = {
    width: '900px',
    height: '900px',
  };

  const refData = {
    filename: `https://files.rcsb.org/download/${selectedPDBid}.pdb`,
    ...(representation !== "default" && {
      config: [{
        type: 'addRepresentation',
        input: representation
      }]
    })
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
  /*
    const handleArrowClick = (direction) => {
      if (direction === 'left' && scrollIndex > 0) {
        setScrollIndex(scrollIndex - 1);
      } else if (direction === 'right' && scrollIndex < PDBids.length - itemsToShow) {
        setScrollIndex(scrollIndex + 1);
      }
    };
  */
  const vizteinKey = `${selectedPDBid}-${representation}`;
  //const transformValue = -(scrollIndex * 150);

  return (
    <div className='reference3D'>
      <div className='left-column'>
        <p className='renderTitle'>{region}</p>
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
