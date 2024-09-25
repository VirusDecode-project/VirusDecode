import React, { useEffect, useState } from "react";
import "../../styles/MRNAdesign.css";
import RNAVisualizer from '../MRNAVisualizer';
import { MRNAData } from '../types';

interface MRNAdesignProps {
  workingHistory: string;
  linearDesignData: MRNAData | null;
}

const MRNAdesign: React.FC<MRNAdesignProps> = ({workingHistory, linearDesignData}) => {
  // const [linearDesignData, setLinearDesignData] = useState<Data | null>(null);
  const [showFullAminoAcidSequence, setShowFullAminoAcidSequence] = useState(false);
  const [showFullSequence, setShowFullSequence] = useState(false);
  const [showFullStructure, setShowFullStructure] = useState(false);
  const zeroWidthSpace = "\u200B";

  if (!linearDesignData) {
    return <div>Loading...</div>;
  }
  const formatSequence = (sequence: string) => {
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

  const formatStructure = (structure: string) => {
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
            sequence={linearDesignData.linearDesign.mRNA_sequence}
            structure={linearDesignData.linearDesign.mRNA_structure}
          />
          {/*visualization 끝*/}

        </div>

        <div className="mrna-column">
          <h3 className="mrna-subtitle-margintop">Amino Acid Sequence</h3>
          <div className="mrna-sequence">
            {showFullAminoAcidSequence
              ? formatSequence(linearDesignData.linearDesign.amino_acid_sequence).map(
                (seq, index) => <div key={index}>{seq}</div>
              )
              : `${linearDesignData.linearDesign.amino_acid_sequence.substring(0, 50)} ...`}
            {!showFullAminoAcidSequence && (
              <span
                className="show-toggle"
                onClick={() => setShowFullAminoAcidSequence(true)}
              >
                show
              </span>
            )}
            {showFullAminoAcidSequence && (
              <span
                className="show-toggle"
                onClick={() => setShowFullAminoAcidSequence(false)}
              >
                ...hide
              </span>
            )}
          </div>
          <h3 className="mrna-subtitle">mRNA Sequence</h3>

          <div className="mrna-sequence">
            {showFullSequence
              ? formatSequence(linearDesignData.linearDesign.mRNA_sequence).map(
                (seq, index) => <div key={index}>{seq}</div>
              )
              : `${linearDesignData.linearDesign.mRNA_sequence.substring(0, 50)} ...`}
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
              ? formatStructure(linearDesignData.linearDesign.mRNA_structure).map(
                (seq, index) => <div key={index}>{seq}</div>
              )
              : `${linearDesignData.linearDesign.mRNA_structure.substring(0, 50)} `}
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
          <p className="mrna-value">{linearDesignData.linearDesign.free_energy} </p>

          <h3 className="mrna-subtitle">mRNA CAI</h3>
          <p className="mrna-value">{linearDesignData.linearDesign.cai}</p>


        </div>
      </div>
      <div className="mrna-container">
        <div className="mrna-column">
          <h2 className="mrna-title">Protein Parameters</h2>

          <h3 className="mrna-subtitle">Molecular Weight</h3>
          <p className="mrna-value">{linearDesignData.protParam.molecular_weight} Da</p>

          <h3 className="mrna-subtitle">Isoelectric Point(Pl)</h3>
          <p className="mrna-value">{linearDesignData.protParam.isoelectric_point}</p>

          <h3 className="mrna-subtitle">Instability Index</h3>
          <p className="mrna-value">{linearDesignData.protParam.instability_index}</p>

          <h3 className="mrna-subtitle">
            Secondary Structure Fraction (Helix, Turn, Sheet)
          </h3>
          <p className="mrna-value">
            ({linearDesignData.protParam.secondary_structure_fraction[0]},{" "}
            {linearDesignData.protParam.secondary_structure_fraction[1]},{" "}
            {linearDesignData.protParam.secondary_structure_fraction[2]})
          </p>

          <h3 className="mrna-subtitle">Gravy</h3>
          <p className="mrna-value">{linearDesignData.protParam.gravy}</p>

          <h3 className="mrna-subtitle">Aromaticity</h3>
          <p className="mrna-value">{linearDesignData.protParam.aromaticity} %</p>
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
              {linearDesignData && linearDesignData.protParam.amino_acid_count ? (
                Object.keys(linearDesignData.protParam.amino_acid_count).map(
                  (amino, index) => (
                    <tr key={index}>
                      <td>{amino}</td>
                      <td>{linearDesignData.protParam.amino_acid_count[amino]}</td>
                      <td>
                        {(
                          linearDesignData.protParam.amino_acid_percent[amino] * 100
                        ).toFixed(2)}
                      </td>
                    </tr>
                  )
                )
              ) : (
                <tr>
                  <td colSpan={2}>Loading...</td>
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
