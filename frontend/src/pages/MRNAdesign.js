import React, { useEffect, useState } from "react";
import "./MRNAdesign.css";
<<<<<<< HEAD
=======
import RNAVisualizer from './RNAVisualizer';
>>>>>>> fffe92e8325694fb1d1937c09790a2a3dc4c81ae

function MRNAdesign() {
  const [data, setData] = useState(null);
  const [showFullSequence, setShowFullSequence] = useState(false);
  const [showFullStructure, setShowFullStructure] = useState(false);
<<<<<<< HEAD
=======
  const zeroWidthSpace = "\u200B";
>>>>>>> fffe92e8325694fb1d1937c09790a2a3dc4c81ae

  useEffect(() => {
    async function fetchData() {
      try {
        const response = await fetch("/genomic_data.json"); //이부분이 백엔드 링크로 수정되어야함.
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
<<<<<<< HEAD
      formatted.push(`${lineNumber} | ${line.join(" ")}`);
=======
      // const zeroWidthSpace = "\u200B";
      if (lineNumber === 1) {
        formatted.push(
          `${zeroWidthSpace} ${zeroWidthSpace} ${lineNumber} | ${line.join(
            " "
          )}`
        );
      } else if (lineNumber === 51) {
        formatted.push(`${zeroWidthSpace} ${lineNumber} | ${line.join(" ")}`);
      } else {
        formatted.push(`${lineNumber} | ${line.join(" ")}`);
      }
    }
    return formatted;
  };

  const formatStructure = (structure) => {
    const formatted = [];
    for (let i = 0; i < structure.length; i += 60) {
      const line = structure.slice(i, i + 60).match(/.{1,10}/g) || [];
      formatted.push(`${line.join("")}`);
>>>>>>> fffe92e8325694fb1d1937c09790a2a3dc4c81ae
    }
    return formatted;
  };

  return (
    <div>
      <p className="detail">
        The detailed description of mRNA conversion can be found{" "}
<<<<<<< HEAD
        <a href="https://www.youtube.com/@%EC%98%A4%EC%B4%88%EB%AA%BD">here.</a>
      </p>

      <div className="mrna-container">
        <div className="mrna-column">
=======
        <a
          href="https://www.nature.com/articles/s41586-023-06127-z"
          target="_blank"
        >
          here.
        </a>
      </p>
      <div className="mrna-container">
        <div className="mrna-column">
          {/*visualization 추가*/}
          <h2 className="mrna-title">mRNA Visualization</h2>
          <RNAVisualizer sequence={data.mRNA_sequence} structure={data.mRNA_structure} />
          {/*visualization 끝*/}
>>>>>>> fffe92e8325694fb1d1937c09790a2a3dc4c81ae
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
<<<<<<< HEAD
=======
                
>>>>>>> fffe92e8325694fb1d1937c09790a2a3dc4c81ae
                show
              </span>
            )}
            {showFullSequence && (
              <span
                className="show-toggle"
                onClick={() => setShowFullSequence(false)}
              >
<<<<<<< HEAD
                {" "}
=======
>>>>>>> fffe92e8325694fb1d1937c09790a2a3dc4c81ae
                ...hide
              </span>
            )}
          </div>

          <h2 className="mrna-title">mRNA Structure</h2>
<<<<<<< HEAD
          <p className="mrna-sequence">
            {showFullStructure
              ? data.mRNA_structure
=======
          <p className="mrna-structure">
            {showFullStructure
              ? formatStructure(data.mRNA_structure).map((seq, index) => (
                  <div key={index}>{seq}</div>
                ))
>>>>>>> fffe92e8325694fb1d1937c09790a2a3dc4c81ae
              : `${data.mRNA_structure.substring(0, 50)} `}
            {!showFullStructure && (
              <span
                className="show-toggle"
                onClick={() => setShowFullStructure(true)}
              >
<<<<<<< HEAD
                show
=======
                {"   "}show
>>>>>>> fffe92e8325694fb1d1937c09790a2a3dc4c81ae
              </span>
            )}
            {showFullStructure && (
              <span
                className="show-toggle"
                onClick={() => setShowFullStructure(false)}
              >
<<<<<<< HEAD
                {" "}
                hide
=======
                ...hide
>>>>>>> fffe92e8325694fb1d1937c09790a2a3dc4c81ae
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
