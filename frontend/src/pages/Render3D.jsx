import React, { useEffect, useState, useRef } from 'react';
import Viztein from 'viztein';
import "./Render3D.css";

const Render3D = ({ region }) => {
  const [PDBids, setPDBids] = useState([]);
  const [selectedPDBid, setSelectedPDBid] = useState("");
  const [PDBInfo, setPDBInfo] = useState([]);
  const [representation, setRepresentation] = useState("default");
  const [error, setError] = useState("");
  const [tooltip, setTooltip] = useState({ visible: false, text: '', x: 0, y: 0 });
  const [scrollIndex, setScrollIndex] = useState(0);
  const itemsToShow = 3;
  const listRef = useRef(null);

  useEffect(() => {
    const fetchPDBids = async () => {
      try {
        const response = await fetch('http://localhost:8080/analysis/render3d');
        if (!response.ok) {
          throw new Error(`Network response was not ok: ${response.statusText}`);
        }
        const PDBlist = await response.json();
        const keys = Object.keys(PDBlist);
        setPDBids(keys);
        const values = Object.values(PDBlist);
        setPDBInfo(values);

        if (keys.length > 0) {
          setSelectedPDBid(keys[0]);
        }
      } catch (error) {
        console.error('Error fetching PDB IDs:', error);
      }
    };
    fetchPDBids();
  }, []);

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
      setError("");
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

  const handleArrowClick = (direction) => {
    if (direction === 'left' && scrollIndex > 0) {
      setScrollIndex(scrollIndex - 1);
    } else if (direction === 'right' && scrollIndex < PDBids.length - itemsToShow) {
      setScrollIndex(scrollIndex + 1);
    }
  };

  const vizteinKey = `${selectedPDBid}-${representation}`;
  const transformValue = -(scrollIndex * 150);

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
        <div className='arrow-container'>
          <div className='arrow arrow-left' onClick={() => handleArrowClick('left')}>◀</div>
          <div className='arrow arrow-right' onClick={() => handleArrowClick('right')}>▶</div>
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
