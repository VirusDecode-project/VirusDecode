import Viztein from 'viztein';
import "./Render3D.css";

function Render3D() {
  const viewportStyle = {
    width: '900px',
    height: '900px',
  };

  const pdbId = "7rn1";

  const refData = {
    filename: `https://files.rcsb.org/download/${pdbId}.pdb`,
    config: [{
      type: 'addRepresentation',
      input: 'ball+stick'
    }]
  };

  return (
    <div className='reference3D'>
      <h4>Reference</h4>
      <Viztein data={refData} viewportStyle={viewportStyle} />
    </div>
  );
}

export default Render3D;
