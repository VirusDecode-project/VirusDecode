import Viztein from 'viztein';
import "./Render3D.css";

function Render3D() {
  const viewportStyle = {
    width: '900px',
    height: '900px',
  };

  const pdbId = "5r7y"; // default, SARS-CoV-2

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
      <h5>{pdbId}.pdb</h5>
    </div>
  );
}

export default Render3D;
