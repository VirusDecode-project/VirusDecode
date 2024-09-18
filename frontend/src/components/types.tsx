export interface Sequence {
  label: string;
  sequence: string;
}

export interface ResponseData {
  aligned_sequences: Record<string, string>; // { label: sequence }
  alignment_index: Record<string, [number, number]>; // { region: [start, end] }
}

export interface Chunk {
  label: string;
  lines: {
    char: string;
    different: boolean;
  }[];
}

export interface AlignmentIndex {
  [key:string]: [number, number];
}
  
export interface ModalData {
  genome: string;
  protein: string;
}
