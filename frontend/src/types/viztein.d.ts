declare module 'viztein' {
  import { ComponentType } from 'react';

  interface VizteinProps {
    data: {
      filename: string;
      config?: any[]; 
    };
    viewportStyle: React.CSSProperties;
  }

  const Viztein: ComponentType<VizteinProps>;

  export default Viztein;
}