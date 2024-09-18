export interface CharObj {
  char: string;
  different: boolean;
}

export interface Sequence {
  label: string;
  sequence: string;
  lines?: CharObj[][];
}

export interface ResponseData {
  aligned_sequences: Record<string, string>; // { label: sequence }
  alignment_index: Record<string, [number, number]>; // { region: [start, end] }
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
