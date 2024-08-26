import React, { useEffect } from 'react';
import './HelpModal.css';
import helpImg1 from '../image/helpImg1.png'
import helpImg2 from '../image/helpImg2.png'
import helpImg3 from '../image/helpImg3.png'

function HelpModal({ isOpen, onClose }) {
    useEffect(() => {
        const handleEsc = (event) => {
            if (event.key === 'Escape') {
                onClose();
            }
        };

        if (isOpen) {
            document.addEventListener('keydown', handleEsc);
        } else {
            document.removeEventListener('keydown', handleEsc);
        }

        return () => {
            document.removeEventListener('keydown', handleEsc);
        };
    }, [isOpen, onClose]);

    if (!isOpen) {
        return null;
    }

    return (
        <div className="help-modal">
            <button className="close-button" onClick={onClose}>×</button>
            {<div>
                <p className='helpTitle'>mRNA Conversion Instructions</p>
                <p className='helpLevel'>1. Select CDS</p>
                <p className='helpContents'>Click on the CDS you wish to convert for mRNA design.</p>
                <img className='helpImageWide' src={helpImg1} alt="Step 1: Select CDS"/>
                <p className='helpLevel'>2. Choose sublineage and input amino acid range</p>
                <p className='helpContents'>- Click on the appropriate sublineage from the initial input sequence.</p>
                <p className='helpContents'>- Enter the start and end positions for the amino acids.</p>
                <p className='helpWarnning'>*Note: Do not exceed 500 amino acids.</p>
                <img className='helpImage' src={helpImg2} alt="Step 2: Choose sublineage and input amino acid range" />
                <p className='helpLevel'>3. Convert</p>
                <p className='helpContents'>Click the "Convert" button to proceed with the mRNA design.</p>
                <img className='helpImage' src={helpImg3} alt="Step 3: Convert" />
            </div>
            }
        </div>
    );
}

export default HelpModal;
