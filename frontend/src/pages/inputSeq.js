import { useState, useEffect } from 'react';
import Button from 'react-bootstrap/Button';
import Modal from 'react-bootstrap/Modal';
import './inputSeq.css';
import googleLoginBtn from './googleloginbtn.png';
import Offcanvas from 'react-bootstrap/Offcanvas';

function InputSeq() {
    const [showModal, setShowModal] = useState(false);
    const [showOffcanvas, setShowOffcanvas] = useState(false);

    const handleCloseModal = () => setShowModal(false);
    const handleShowModal = () => setShowModal(true);

    const handleCloseOffcanvas = () => setShowOffcanvas(false);
    const handleShowOffcanvas = () => setShowOffcanvas(true);

    useEffect(() => {
        // 컴포넌트가 마운트될 때 모달을 표시합니다.
        setShowModal(true);
    }, []);

    return (
        <div>
            <div>시컨스를 입력하세요!</div>

            <Modal show={showModal} onHide={handleCloseModal}>
                <Modal.Header closeButton>
                    <Modal.Title>Welcome to VirusDecode!</Modal.Title>
                </Modal.Header>
                <Modal.Body className="modal-body-centered">
                    Log in to get your<br />
                    virus analysis records.
                    <div className="google-login-button-container">
                        <img
                            src={googleLoginBtn}
                            alt="Google Login Button"
                            className="google-login-button"
                            onClick={() => { /* 구글 로그인 함수 호출 */ }}
                        />
                    </div>
                </Modal.Body>
                <Modal.Footer>
                    <Button variant="primary" onClick={handleCloseModal}>
                        Stay logged out
                    </Button>
                </Modal.Footer>
            </Modal>

            <Button variant="primary" onClick={handleShowOffcanvas}>
                히스토리
            </Button>

            <Offcanvas show={showOffcanvas} onHide={handleCloseOffcanvas}>
                <Offcanvas.Header closeButton>
                    <Offcanvas.Title>Offcanvas</Offcanvas.Title>
                </Offcanvas.Header>
                <Offcanvas.Body>
                    Some text as placeholder. In real life you can have the elements you
                    have chosen. Like, text, images, lists, etc.
                </Offcanvas.Body>
            </Offcanvas>
        </div>
    );
}

export default InputSeq;
