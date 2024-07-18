import { useState, useEffect } from 'react';
import Button from 'react-bootstrap/Button';
import Modal from 'react-bootstrap/Modal';
import './inputSeq.css';
import googleLoginBtn from './googleloginbtn.png';
import Offcanvas from 'react-bootstrap/Offcanvas';
import { useNavigate } from 'react-router-dom';
import GoogleLoginButton from '../GoogleLoginButton.js'; 

function InputSeq() {
    const [showModal, setShowModal] = useState(false);
    const [showOffcanvas, setShowOffcanvas] = useState(false);

    const navigate = useNavigate(); // useNavigate 초기화

    const handleCloseModal = () => setShowModal(false);
    const handleShowModal = () => setShowModal(true);

    const handleCloseOffcanvas = () => setShowOffcanvas(false);
    const handleShowOffcanvas = () => setShowOffcanvas(true);

    useEffect(() => {
        // 컴포넌트가 마운트될 때 모달을 표시합니다.
        setShowModal(true);
    }, []);

    return (
        <div className="next-page-container">
            <div>Reference Sequence ID<br/><input/><button>DONE</button></div>
            <h4>Variance Sequence<br/>추후 수정 예정</h4>

            <Modal show={showModal} onHide={handleCloseModal}>
                <Modal.Header closeButton>
                    <Modal.Title>Welcome to VirusDecode!</Modal.Title>
                </Modal.Header>
                <Modal.Body className="modal-body-centered">
                    Log in to get your<br />
                    virus analysis records.
                    <div className="google-login-button-container">
                        <GoogleLoginButton/>
                    </div>
                </Modal.Body>
                <Modal.Footer>
                    <Button variant="primary" onClick={handleCloseModal}>
                        Stay logged out
                    </Button>
                </Modal.Footer>
            </Modal>

            <Button variant="primary" onClick={handleShowOffcanvas}>
                사이드바
            </Button>

            <Offcanvas show={showOffcanvas} onHide={handleCloseOffcanvas}>
                <Offcanvas.Header closeButton>
                    <Offcanvas.Title>History</Offcanvas.Title>
                </Offcanvas.Header>
                <Offcanvas.Body>
                    <div>yesterday</div>
                    <div>Reference1</div>
                    <div>Reference2</div>
                    <br />
                    <div>previous 7 day</div>
                    <div>Reference1</div>
                    <div>Reference2</div>
                    <div>Reference3</div>
                </Offcanvas.Body>
            </Offcanvas>

            <button className="next-page" onClick={() => { navigate('/analysis') }}>NEXT</button>
        </div>
    );
}

export default InputSeq;
