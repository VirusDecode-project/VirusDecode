export interface CharObj {
  char: string;
  different: boolean;
}

export interface Sequence {
  label: string;
  sequence: string;
  lines?: CharObj[][] | null;
}

export interface Chunk {
  label: string;
  lines: CharObj[][];
}

export interface AlignmentIndex {
  [key:string]: [number, number];
}
  
export interface ModalData {
  genome: string;
  protein: string;
}

export interface LinearDesign {
  mRNA_sequence: string;
  mRNA_structure: string;
  amino_acid_sequence: string;
  free_energy: number;
  cai: number;
}

export interface ProtParam {
  molecular_weight: number;
  isoelectric_point: number;
  instability_index: number;
  secondary_structure_fraction: [number, number, number];
  gravy: number;
  aromaticity: number;
  amino_acid_count: Record<string, number>;
  amino_acid_percent: Record<string, number>;
}

export interface MRNAData {
  linearDesign: LinearDesign;
  protParam: ProtParam;
}

export interface AlignmentData {
  aligned_sequences: Record<string, string>; // { label: sequence }
  alignment_index: Record<string, [number, number]>; // { region: [start, end] }
}

export interface PDBResponse {
  [key: string]: string;
}