import React, { useEffect, useRef } from 'react';
import '../styles/HelpModal.css';
import helpImg1 from '../assets/helpImg1.png';
import helpImg2 from '../assets/helpImg2.png';
import helpImg3 from '../assets/helpImg3.png';

function HelpModal({ isOpen, onClose }) {
    const modalRef = useRef(null); // 모달을 참조하기 위한 useRef

    useEffect(() => {
        const handleEsc = (event) => {
            if (event.key === 'Escape') {
                onClose();
            }
        };

        const handleClickOutside = (event) => {
            if (modalRef.current && !modalRef.current.contains(event.target)) {
                onClose();
            }
        };

        if (isOpen) {
            document.addEventListener('keydown', handleEsc);
            document.addEventListener('mousedown', handleClickOutside); // 모달 바깥 클릭 감지
        }

        // Clean up 이벤트 리스너
        return () => {
            document.removeEventListener('keydown', handleEsc);
            document.removeEventListener('mousedown', handleClickOutside);
        };
    }, [isOpen, onClose]);

    if (!isOpen) {
        return null;
    }

    return (
        <div className="help-modal-overlay">
            <div className="help-modal" ref={modalRef}>
                <button className="close-button" onClick={onClose}>×</button>
                <div>
                    <p className='helpTitle'>mRNA Conversion Instructions</p>
                    <p className='helpLevel'>1. Select CDS</p>
                    <p className='helpContents'>Click on the CDS you wish to convert for mRNA design.</p>
                    <img className='helpImageWide' src={helpImg1} alt="Step 1: Select CDS" />
                    <p className='helpLevel'>2. Choose sublineage and input amino acid range</p>
                    <p className='helpContents'>- Click on the appropriate sublineage from the initial input sequence.</p>
                    <p className='helpContents'>- Enter the start and end positions for the amino acids.</p>
                    <p className='helpWarnning'>*Note: Do not exceed 500 amino acids.</p>
                    <img className='helpImage' src={helpImg2} alt="Step 2: Choose sublineage and input amino acid range" />
                    <p className='helpLevel'>3. Convert</p>
                    <p className='helpContents'>Click the "Convert" button to proceed with the mRNA design.</p>
                    <img className='helpImage' src={helpImg3} alt="Step 3: Convert" />
                </div>
            </div>
        </div>
    );
}

export default HelpModal;
