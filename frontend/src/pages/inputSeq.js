import { useState, useEffect } from 'react';
import { Button, Modal, Form, Offcanvas } from 'react-bootstrap';
import { useNavigate } from 'react-router-dom';
import GoogleLoginButton from '../GoogleLoginButton.js'; // 경로 확인
import './inputSeq.css';
import historyIcon from './history.png';

function InputSeq() {
    const [showModal, setShowModal] = useState(false);
    const [showOffcanvas, setShowOffcanvas] = useState(false);

    const navigate = useNavigate();

    const handleCloseModal = () => setShowModal(false);
    const handleShowModal = () => setShowModal(true);

    const handleCloseOffcanvas = () => setShowOffcanvas(false);
    const handleShowOffcanvas = () => {
        setShowOffcanvas(true);
        document.body.style.overflow = 'hidden';  // 오버플로우 숨김
    };

    useEffect(() => {
        setShowModal(true);
        return () => {
            document.body.style.overflow = 'auto';  // 오버플로우 기본값으로 재설정
        };
    }, []);

    useEffect(() => {
        if (!showOffcanvas) {
            document.body.style.overflow = 'auto';  // 오버플로우 기본값으로 재설정
        }
    }, [showOffcanvas]);

    return (
        <div className={`next-page-container ${showOffcanvas ? 'shrink' : ''}`}>
            <div className="container mt-4">
                <Form>
                    <Form.Group controlId="referenceSequenceId">
                        <Form.Label>Reference Sequence ID</Form.Label>
                        <Form.Control type="text" placeholder="Enter sequence ID" />
                    </Form.Group>
                    <Button variant="primary" className="mt-3" >DONE</Button>

                    <h4 className="mt-5">Variant Sequence</h4>
                    <Form.Group controlId="formFile" className="mb-3">
                        <Form.Label>Upload File</Form.Label>
                        <div className="upload-box">
                            <Form.Control type="file" label="Drag your FASTA files here" custom />
                        </div>
                    </Form.Group>
                    <Form.Group controlId="pasteSequence">
                        <Form.Label>Paste Sequence</Form.Label>
                        <Form.Control as="textarea" rows={3} placeholder="sequence1" />
                    </Form.Group>
                    <Button variant="link" className="mt-3">+ Add Sequence</Button>
                </Form>
            </div>

            <Modal show={showModal} onHide={handleCloseModal}>
                <Modal.Header closeButton>
                    <Modal.Title>Welcome to VirusDecode!</Modal.Title>
                </Modal.Header>
                <Modal.Body className="modal-body-centered">
                    Log in to get your<br />
                    virus analysis records.
                    <div className="google-login-button-container">
                        <GoogleLoginButton />
                    </div>
                </Modal.Body>
                <Modal.Footer>
                    <Button variant="primary" onClick={handleCloseModal}>
                        Stay logged out
                    </Button>
                </Modal.Footer>
            </Modal>

            <Button variant="primary" onClick={handleShowOffcanvas} className="mt-4">
                사이드바
            </Button>

            <Offcanvas show={showOffcanvas} onHide={handleCloseOffcanvas} backdrop={false} style={{ width: '300px' }}>
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

            <h4 className="next-page" onClick={() => { navigate('/analysis') }}>{'Next ->'}</h4>
        </div>
    );
}

export default InputSeq;
