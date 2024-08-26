import React, { useEffect, useState } from "react";
import "../styles/MRNAdesign.css";
import RNAVisualizer from '../components/MRNAVisualizer';

function MRNAdesign() {
  const [data, setData] = useState(null);
  const [showFullSequence, setShowFullSequence] = useState(false);
  const [showFullStructure, setShowFullStructure] = useState(false);
  const zeroWidthSpace = "\u200B";

  /* backend 수정 코드 시작 */
  useEffect(() => {
    const fetchJsonData = async () => {
      try {
        const serverResponse = await fetch('http://localhost:8080/analysis/re-mrnadesign');

        if (!serverResponse.ok) {
          const errorMessage = await serverResponse.text();
          throw new Error(errorMessage);
        }

        const responseData = await serverResponse.json();
        setData(responseData); // JSON 데이터를 상태로 설정
      } catch (error) {
        console.error("An error occurred during the request: ", error.message);
      }
    };
    fetchJsonData();
  }, []);


  if (!data) {
    return <div>Loading...</div>;
  }
  const formatSequence = (sequence) => {
    const formatted = [];
    for (let i = 0; i < sequence.length; i += 50) {
      const line = sequence.slice(i, i + 50).match(/.{1,10}/g) || [];
      const lineNumber = i + 1;
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
    }
    return formatted;
  };

  return (
    <div>
      <p className="detail">
        The detailed description of mRNA conversion can be found{" "}
        <a
          href="https://github.com/LinearDesignSoftware/LinearDesign"
          target="_blank"
          rel="noopener noreferrer"  // 'rel' 속성 추가
        >
          here.
        </a>
      </p>


      <div className="mrna-container">
        <div className="mrna-column">
          {/*visualization 추가*/}
          <h2 className="mrna-title">mRNA Parameters</h2>
          <h3 className="mrna-subtitle">mRNA Visualization</h3>
          <RNAVisualizer
            sequence={data.linearDesign.mRNA_sequence}
            structure={data.linearDesign.mRNA_structure}
          />
          {/*visualization 끝*/}

        </div>

        <div className="mrna-column">
          <h3 className="mrna-subtitle-margintop">mRNA Sequence</h3>
          <div className="mrna-sequence">
            {showFullSequence
              ? formatSequence(data.linearDesign.mRNA_sequence).map(
                (seq, index) => <div key={index}>{seq}</div>
              )
              : `${data.linearDesign.mRNA_sequence.substring(0, 50)} ...`}
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
                ...hide
              </span>
            )}
          </div>
          <h3 className="mrna-subtitle">mRNA Structure</h3>

          <p className="mrna-structure">
            {showFullStructure
              ? formatStructure(data.linearDesign.mRNA_structure).map(
                (seq, index) => <div key={index}>{seq}</div>
              )
              : `${data.linearDesign.mRNA_structure.substring(0, 50)} `}
            {!showFullStructure && (
              <span
                className="show-toggle"
                onClick={() => setShowFullStructure(true)}
              >
                {"   "}show
              </span>
            )}
            {showFullStructure && (
              <span
                className="show-toggle"
                onClick={() => setShowFullStructure(false)}
              >
                ...hide
              </span>
            )}
          </p>

          <h3 className="mrna-subtitle">mRNA folding free energy</h3>
          <p className="mrna-value">{data.linearDesign.free_energy} </p>

          <h3 className="mrna-subtitle">mRNA CAI</h3>
          <p className="mrna-value">{data.linearDesign.cai}</p>


        </div>
      </div>
      <div className="mrna-container">
        <div className="mrna-column">
          <h2 className="mrna-title">Protein Parameters</h2>

          <h3 className="mrna-subtitle">Molecular Weight</h3>
          <p className="mrna-value">{data.protParam.molecular_weight} Da</p>

          <h3 className="mrna-subtitle">Isoelectric Point(Pl)</h3>
          <p className="mrna-value">{data.protParam.isoelectric_point}</p>

          <h3 className="mrna-subtitle">Instability Index</h3>
          <p className="mrna-value">{data.protParam.instability_index}</p>

          <h3 className="mrna-subtitle">
            Secondary Structure Fraction (Helix, Turn, Sheet)
          </h3>
          <p className="mrna-value">
            ({data.protParam.secondary_structure_fraction[0]},{" "}
            {data.protParam.secondary_structure_fraction[1]},{" "}
            {data.protParam.secondary_structure_fraction[2]})
          </p>

          <h3 className="mrna-subtitle">Gravy</h3>
          <p className="mrna-value">{data.protParam.gravy}</p>

          <h3 className="mrna-subtitle">Aromaticity</h3>
          <p className="mrna-value">{data.protParam.aromaticity} %</p>
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
              {data && data.protParam.amino_acid_count ? (
                Object.keys(data.protParam.amino_acid_count).map(
                  (amino, index) => (
                    <tr key={index}>
                      <td>{amino}</td>
                      <td>{data.protParam.amino_acid_count[amino]}</td>
                      <td>
                        {(
                          data.protParam.amino_acid_percent[amino] * 100
                        ).toFixed(2)}
                      </td>
                    </tr>
                  )
                )
              ) : (
                <tr>
                  <td colSpan="2">Loading...</td>
                </tr>
              )}

            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}

export default MRNAdesign;
