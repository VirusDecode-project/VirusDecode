declare module 'fornac' {
  export class FornaContainer {
      constructor(element: HTMLElement, options?: any);
      addRNA(structure: string, options?: any): void;
      clearNodes(): void;
  }
}
