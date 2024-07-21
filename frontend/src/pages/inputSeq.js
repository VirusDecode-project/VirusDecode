import { useState, useEffect } from 'react';
import { Button, Modal, Form, Offcanvas, Row, Col } from 'react-bootstrap';
import { useNavigate } from 'react-router-dom';
import GoogleLoginButton from '../GoogleLoginButton.js'; // 경로 확인
import './inputSeq.css';
import historyIcon from './history.png';
import editIcon from './edit.png';

function InputSeq() {
    const [showModal, setShowModal] = useState(false);
    const [showOffcanvas, setShowOffcanvas] = useState(true);

    let navigate = useNavigate();

    const handleCloseModal = () => setShowModal(false);

    let handleCloseOffcanvas = () => setShowOffcanvas(false);
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
            <div className="header-bar">
                {!showOffcanvas && (
                    <>
                        <img onClick={handleShowOffcanvas} style={{ cursor: 'pointer' }} src={historyIcon} alt="History" className="history-icon" />
                        <img src={editIcon} alt="Edit" className="edit-icon" />
                    </>
                )}
                <span className='logo-text' onClick={() => { navigate('/') }} style={{ cursor: 'pointer' }}>VirusDecode</span>
            </div>


            <div className="container mt-4" style={{ marginLeft: '30px' }}>
                
                <Form>
                    
                    <h5 className="RS-id">Reference Sequence ID</h5>
                    <Row className="align-items-center">
                        <Col md={6}>
                            <Form.Group controlId="referenceSequenceId">
                                <Form.Control type="text" placeholder="Enter sequence ID" className="input-field" />
                            </Form.Group>
                        </Col>
                        <Col md={1} className="d-flex justify-content-end align-items-center">
                            <Button variant="primary" className="done-button">DONE</Button>
                        </Col>
                    </Row>

                </Form>

                <Form>

                    <div className="mb-5"></div>

                    <h5 className="RS-id">Variant Sequence</h5>


                    <Form.Group controlId="formFile" className="mb-3">
                        <Form.Label>Upload File</Form.Label>
                        <Col md={6}>
                            <div className="upload-box">
                                <Form.Control type="file" label="Drag your FASTA files here" custom />
                            </div>
                        </Col>
                    </Form.Group>


                    <Form.Group controlId="pasteSequence">
                        <Form.Label>Paste Sequence</Form.Label>
                        <Col md={9}>
                            <Form.Control as="textarea" rows={3} placeholder="sequence1" />
                        </Col>
                    </Form.Group>

                    <Row>
                        <Col md={9} className="d-flex justify-content-start">
                            <Button variant="link" className="mt-3">+ Add Sequence</Button>
                        </Col>
                    </Row>




                </Form>
            </div>

            <Modal show={showModal} onHide={handleCloseModal} centered>
                <Modal.Header className="modal-body-centered">
                    <Modal.Title>Welcome to VirusDecode!</Modal.Title>
                </Modal.Header>
                <Modal.Body className="modal-body-centered">
                    Log in to get your<br />
                    virus analysis records.
                    <div className="google-login-button-container">
                        <GoogleLoginButton />
                    </div>
                </Modal.Body>

                <Modal.Footer className="modal-body-centered">
                    <div>
                        New to VirusDecode?<br />
                        Sign up now!</div>
                    <div className="google-login-button-container">
                        <GoogleLoginButton />
                    </div>

                </Modal.Footer>
                <Modal.Footer>
                    <p className='logged-out' onClick={handleCloseModal}>
                        Stay logged out
                    </p>
                </Modal.Footer>
            </Modal>

            <Offcanvas className="custom-offcanvas" show={showOffcanvas} onHide={handleCloseOffcanvas} backdrop={false} style={{ width: '260px' }}>
                <Offcanvas.Header>
                    <Offcanvas.Title>
                        <img onClick={handleCloseOffcanvas} src={historyIcon} alt="History" className="history-icon" style={{ cursor: 'pointer' }} />
                    </Offcanvas.Title>
                    <img src={editIcon} alt="Edit" className="edit-icon" />
                </Offcanvas.Header>
                <Offcanvas.Body>
                    <div>Today</div>

                </Offcanvas.Body>
            </Offcanvas>

            <h4 className="next-page" onClick={() => { navigate('/analysis') }}>{'Next ->'}</h4>
        </div>
    );
}

export default InputSeq;
