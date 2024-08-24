import React, { useState, useEffect } from 'react';
import './HelpModal.css';

function HelpModal({ isOpen, onClose, content }) {
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
            {content}
        </div>
    );
}

export default HelpModal;
