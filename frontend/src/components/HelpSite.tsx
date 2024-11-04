import React, { useState, useEffect } from 'react';
import '../styles/HelpSite.css';
import guideImg1 from '../assets/guideImg1.png';
import guideImg2 from '../assets/guideImg2.png';
import guideImg3 from '../assets/guideImg3.png';
import guideImg4 from '../assets/guideImg4.png';

interface Page {
    id: number;
    image: string;
}

const pages: Page[] = [
    { id: 1, image: guideImg1 },
    { id: 2, image: guideImg2 },
    { id: 3, image: guideImg3 },
    { id: 4, image: guideImg4 },
];

interface HelpSiteProps {
    isOpen: boolean;
    onClose: () => void;
}

const HelpSite: React.FC<HelpSiteProps> = ({ isOpen, onClose }) => {
    const [currentPage, setCurrentPage] = useState<number>(0);

    // 도옴말 창 다시 열 떄 첫 페이지로 초기화
    useEffect(() => {
        if (isOpen) {
            setCurrentPage(0);
        }
    }, [isOpen]);

    // 이미지 미리 로드
    useEffect(() => {
        pages.forEach(page => {
            const img = new Image();
            img.src = page.image;
        });
    }, []);
    
    const goToNextPage = () => {
        setCurrentPage((prevPage) => (prevPage + 1) % pages.length);
    };

    const goToPreviousPage = () => {
        setCurrentPage((prevPage) => (prevPage - 1 + pages.length) % pages.length);
    };

    // 키보드 화살표 인식
    useEffect(() => {
        if (!isOpen) return; // 모달이 열려 있을 때만 실행

        const handleKeyDown = (event: KeyboardEvent) => {
            if (event.key === 'ArrowRight') {
                goToNextPage();
            } else if (event.key === 'ArrowLeft') {
                goToPreviousPage();
            }
        };

        window.addEventListener('keydown', handleKeyDown);

         // 컴포넌트가 언마운트되거나 모달이 닫힐 때 이벤트 리스너 제거
        return () => {
            window.removeEventListener('keydown', handleKeyDown);
        };
    }, [isOpen]);

    if (!isOpen) return null;

    return (
        <div className="modal-overlay" onClick={onClose}>
            <div className="modal-content" onClick={(e) => e.stopPropagation()}>
                <button className="close-button" onClick={onClose}>
                    &times;
                </button>
                <div className="page-content">
                    <img src={pages[currentPage].image} alt={`Page ${currentPage + 1}`} />
                </div>
                <div className="navigation-buttons">
                    <button onClick={goToPreviousPage}>&lt;</button>
                    <span>{currentPage + 1} / {pages.length}</span>
                    <button onClick={goToNextPage}>&gt;</button>
                </div>
            </div>
        </div>
    );
};

export default HelpSite;
