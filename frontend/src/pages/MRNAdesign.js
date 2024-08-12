import React, { useEffect, useState } from "react";
import "./MRNAdesign.css";

function MRNAdesign() {
  const [data, setData] = useState(null);
  const [showFullSequence, setShowFullSequence] = useState(false);
  const [showFullStructure, setShowFullStructure] = useState(false);

  useEffect(() => {
    async function fetchData() {
      try {
        const response = await fetch("/genomic_data.json");
        const result = await response.json();
        setData(result);
      } catch (error) {
        console.error("Error fetching data:", error);
      }
    }

    fetchData();
  }, []);

  if (!data) {
    return <div>Loading...</div>;
  }
  const formatSequence = (sequence) => {
    const formatted = [];
    for (let i = 0; i < sequence.length; i += 50) {
      const line = sequence.slice(i, i + 50).match(/.{1,10}/g) || [];
      const lineNumber = i + 1;
      formatted.push(`${lineNumber} | ${line.join(" ")}`);
    }
    return formatted;
  };

  return (
    <div>
      <p className="detail">
        The detailed description of mRNA conversion can be found here.
      </p>

      <div className="mrna-container">
        <div className="mrna-column">
          <h2 className="mrna-title">mRNA Sequence</h2>
          <div className="mrna-sequence">
            {showFullSequence
              ? formatSequence(data.mRNA_sequence).map((seq, index) => (
                  <div key={index}>{seq}</div>
                ))
              : `${data.mRNA_sequence.substring(0, 50)} ...`}
            {!showFullSequence && (
              <span
                className="show-toggle"
                onClick={() => setShowFullSequence(true)}
              >
                show
              </span>
            )}
            {showFullSequence && (
              <span
                className="show-toggle"
                onClick={() => setShowFullSequence(false)}
              >
                {" "}
                ...hide
              </span>
            )}
          </div>

          <h2 className="mrna-title">mRNA Structure</h2>
          <p className="mrna-sequence">
            {showFullStructure
              ? data.mRNA_structure
              : `${data.mRNA_structure.substring(0, 50)} `}
            {!showFullStructure && (
              <span
                className="show-toggle"
                onClick={() => setShowFullStructure(true)}
              >
                show
              </span>
            )}
            {showFullStructure && (
              <span
                className="show-toggle"
                onClick={() => setShowFullStructure(false)}
              >
                {" "}
                hide
              </span>
            )}
          </p>

          <h3 className="mrna-subtitle">mRNA folding free energy</h3>
          <p className="mrna-value">{data.mRNA_folding_free_energy} kcal/mol</p>

          <h3 className="mrna-subtitle">mRNA CAI</h3>
          <p className="mrna-value">{data.mRNA_CAI}</p>

          <h3 className="mrna-subtitle">Molecular Weight</h3>
          <p className="mrna-value">{data.Molecular_Weight} Da</p>

          <h3 className="mrna-subtitle">Isoelectric Point(Pl)</h3>
          <p className="mrna-value">{data.Isoelectric_Point}</p>

          <h3 className="mrna-subtitle">Instability Index</h3>
          <p className="mrna-value">{data.Instability_Index}</p>

          <h3 className="mrna-subtitle">
            Secondary Structure Fraction (Helix, Turn, Sheet)
          </h3>
          <p className="mrna-value">
            ({data.Secondary_Structure_Fraction.Helix},{" "}
            {data.Secondary_Structure_Fraction.Turn},{" "}
            {data.Secondary_Structure_Fraction.Sheet})
          </p>

          <h3 className="mrna-subtitle">Gravy</h3>
          <p className="mrna-value">{data.Gravy}</p>

          <h3 className="mrna-subtitle">Aromaticity</h3>
          <p className="mrna-value">{data.Aromaticity} %</p>
        </div>

        <div className="mrna-column">
          <table className="amino-acid-table">
            <thead>
              <tr>
                <th> </th>
                <th>Amino Acid Count</th>
                <th>Amino Acid Percent (%)</th>
              </tr>
            </thead>
            <tbody>
              {data.Amino_Acid_Data.map((amino, index) => (
                <tr key={index}>
                  <td>{amino.Amino_Acid}</td>
                  <td>{amino.Count}</td>
                  <td>{amino.Percent}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}

export default MRNAdesign;
