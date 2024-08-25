import { useState } from "react";
import { Button, Form, Row, Col } from "react-bootstrap";
import { useNavigate } from "react-router-dom";
import "./inputSeq.css";
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import {
  faChevronDown,
  faChevronRight,
  faTrash,
} from "@fortawesome/free-solid-svg-icons";
import uploadIcon from "../image/upload_icon.png";
// GK - Loading 컴포넌트 추가
import Loading from '../components/Loading';

function InputSeq({ setTab }) {
  let navigate = useNavigate();

  /*-----------다솔님 코드 구현 함수, 변수---------------*/
  const [editingFileIndex, setEditingFileIndex] = useState(null);
  const [uploadedFiles, setUploadedFiles] = useState([]);
  const [sequences, setSequences] = useState([
    //GK name에서 공백 제거
    { id: 1, name: "Sequence1", value: "", visible: true },
  ]);
  const [editingId, setEditingId] = useState(null);
  const [nextId, setNextId] = useState(2);

  /*parkki */
  const [referenceSequenceId, setReferenceSequenceId] = useState("");
  const [responseMessage, setResponseMessage] = useState("");
  let [isLoading, setIsLoading] = useState(false);
  /*parkki */

  // next 비활성화
  const [responseReceived, setResponseReceived] = useState(false);

  /*parkki */
  const handleDoneSubmit = async (e) => {
    e.preventDefault(); // 폼의 기본 제출 동작 방지
    const requestData = { sequenceId: referenceSequenceId };

    try {
      const serverResponse = await fetch("http://localhost:8080/inputSeq/reference", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(requestData),
      });

      // 서버 응답이 올바르지 않으면 오류 메시지 생성
      if (!serverResponse.ok) {
        const errorMessage = await serverResponse.text();
        throw new Error(errorMessage);
      }

      // 서버에서 받은 응답을 텍스트로 처리한 후 JSON으로 파싱
      const responseData = await serverResponse.json();

      let formattedMessage = "";
      for (const [key, value] of Object.entries(responseData)) {
        formattedMessage += `<span class="key">${key}:</span> <span class="value">${value}</span><br />`;
      }

      // 서버 응답 메시지 설정
      setResponseMessage(formattedMessage);

      setResponseReceived(true);
    } catch (error) {
      console.error("An error occurred during the request: ", error.message);
      setResponseMessage("An error occurred during the request: " + error.message);
    }
  };

  // next button 클릭시 서버로 (Sequence ID, file, sequence)전송
  const handleFileUploadToServer = async (e) => {
    e.preventDefault(); // 폼의 기본 제출 동작 방지

    // sequences를 Map 형태로 변환
    const sequencesMap = new Map();

    sequences.forEach((seq) => {
      if (seq.value) {
        // seq.value가 존재하고 빈 문자열이 아닌 경우에만 저장
        sequencesMap.set(seq.name, seq.value);
      }
    });

    // 모든 파일을 문자열로 읽기
    const filesContent = await Promise.all(
      uploadedFiles.map(async (uploadedFile) => {
        const content = await uploadedFile.file.text(); // 파일을 텍스트로 읽기
        return { name: uploadedFile.name, content: content };
      })
    );

    // sequences와 files가 비어 있는지 확인
    const hasSequences = sequencesMap.size > 0;
    const hasFiles = filesContent.length > 0;

    // 빈 데이터에 대한 기본 처리
    const jsonData = {
      referenceSequenceId: referenceSequenceId || null,
      sequences: hasSequences ? Object.fromEntries(sequencesMap) : {}, // 비어 있을 경우 빈 객체
      files: hasFiles ? filesContent : [], // 비어 있을 경우 빈 배열
    };

    try {
      // GK - Loading 위치 이동
      setIsLoading(true);
      const serverResponse = await fetch("http://localhost:8080/inputSeq/alignment", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(jsonData),
      });

      if (!serverResponse.ok) {
        const errorMessage = await serverResponse.text();
        throw new Error(errorMessage);
      }

      const responseData = await serverResponse.json();
      setTab(0);
      setIsLoading(false);
      navigate("/analysis", { state: { responseData: responseData } });
    } catch (error) {
      console.error("An error occurred during the request: ", error.message);
    }
  };
  /* parkki */

  const handleFileUpload = (event) => {
    const files = Array.from(event.target.files);
    const newFiles = files.map((file) => ({ name: file.name, file }));
    setUploadedFiles([...uploadedFiles, ...newFiles]);
    setEditingFileIndex(null);
  };

  const handleFileNameChange = (index, name) => {
    const updatedFiles = [...uploadedFiles];
    updatedFiles[index] = { ...updatedFiles[index], name };
    setUploadedFiles(updatedFiles);
  };

  const deleteUploadedFile = (index) => {
    setUploadedFiles(uploadedFiles.filter((_, i) => i !== index));
  };
  const toggleVisibility = (id) => {
    setSequences(
      sequences.map((seq) =>
        seq.id === id ? { ...seq, visible: !seq.visible } : seq
      )
    );
  };

  const handleNameChange = (id, name) => {
    setSequences(
      sequences.map((seq) => (seq.id === id ? { ...seq, name } : seq))
    );
  };

  const deleteSequence = (id) => {
    setSequences(sequences.filter((seq) => seq.id !== id));
  };

  const handleSequenceChange = (id, value) => {
    setSequences(
      sequences.map((seq) => (seq.id === id ? { ...seq, value } : seq))
    );
  };
  const addSequence = (event) => {
    event.preventDefault();
    setSequences([
      ...sequences,
      //GK name에서 공백 제거
      { id: nextId, name: `Sequence${nextId}`, value: "", visible: true },
    ]);
    setNextId(nextId + 1);
  };

  return (
    <div>
      {isLoading ? (
        // GK - Loading 컴포넌트로 변경
        <Loading text="Analyzing" />
        // <div className="loading-container">
        //   <img src={loadingImage} alt="Loading..." className="loading-image" />
        //   <div className="loading-text">{loadingText}</div>
        // </div>
      ) : (
        <div>
          <div className="container mt-4" style={{ marginLeft: "75px" }}>
            <Form>
              <Row className="align-items-center">
                <Col md={6}>
                  <p className="RS-id">Reference Sequence ID</p>
                  <Form.Group controlId="referenceSequenceId">
                    {/*parkki */}
                    <Form.Control
                      type="text"
                      placeholder="Enter sequence ID"
                      className="input-field"
                      value={referenceSequenceId}
                      onChange={(e) => setReferenceSequenceId(e.target.value)}
                    />
                    {/*parkki */}
                  </Form.Group>
                </Col>
                <Col
                  md={1}
                  className="d-flex justify-content-end align-items-center"
                >
                  {/* <Button variant="primary" className="done-button">DONE</Button> */}
                  {/*parkki */}
                  <Button
                    type="submit"
                    variant="primary"
                    onClick={handleDoneSubmit}
                    className="done-button"
                  >
                    DONE
                  </Button>
                  {/*parkki */}
                </Col>
              </Row>
            </Form>
            {/*parkki */}
            {responseMessage && (
              <div className="response-message">
                <p className="metadata">Metadata</p>
                <div
                  className="metadata2"
                  dangerouslySetInnerHTML={{ __html: responseMessage }}
                />
              </div>
            )}
            {/*parkki */}

            <Form>
              <div className="mb-5"></div>

              <p className="RS-id">Variant Sequence</p>

              <Form.Group controlId="formFile" className="mb-3">
                <Form.Label>Upload File</Form.Label>
                {/* ------------다솔님 업로드 박스--------시작------- */}
                <Row className="align-items-center">
                  <Col md={6}>
                    <div className="upload-box">
                      <input
                        type="file"
                        className="file-input"
                        accept=".fasta"
                        multiple
                        onChange={handleFileUpload}
                      />
                      <div className="upload-text">
                        <img
                          src={uploadIcon}
                          alt="Upload Icon"
                          className="upload-icon"
                        />
                        <p>Drag your FASTA files here</p>
                      </div>
                    </div>
                    {uploadedFiles.map((uploadedFile, index) => (
                      <div key={index} className="uploaded-file">
                        {editingFileIndex === index ? (
                          <input
                            type="text"
                            value={uploadedFile.name}
                            onChange={(e) =>
                              handleFileNameChange(index, e.target.value)
                            }
                            onBlur={() => setEditingFileIndex(null)}
                            className="edit-file-name-input"
                            autoFocus
                          />
                        ) : (
                          <span onClick={() => setEditingFileIndex(index)}>
                            {uploadedFile.name}
                          </span>
                        )}
                        <FontAwesomeIcon
                          icon={faTrash}
                          className="delete-icon"
                          onClick={() => deleteUploadedFile(index)}
                        />
                      </div>
                    ))}
                  </Col>
                </Row>
              </Form.Group>
              {/* ------------다솔님 업로드 박스--------끝------- */}

              {/* -----------------다솔님 Paste Sequence ------------시작---- */}
              <Form.Group>
                <Form.Label>Paste Sequence</Form.Label>
                <Row>
                  <Col md={9} className="text-left">
                    {sequences.map((seq) => (
                      <div key={seq.id} className="form-group">
                        <div className="sequence-header d-flex align-items-center justify-content-start">
                          {" "}
                          {/* justify-content-start 클래스 추가 */}
                          <FontAwesomeIcon
                            icon={seq.visible ? faChevronDown : faChevronRight}
                            className="chevron-icon"
                            onClick={() => toggleVisibility(seq.id)}
                          />
                          {editingId === seq.id ? (
                            <input
                              type="text"
                              value={seq.name}
                              onChange={(e) =>
                                handleNameChange(seq.id, e.target.value)
                              }
                              onBlur={() => setEditingId(null)}
                              className="edit-name-input"
                              autoFocus
                            />
                          ) : (
                            <span onClick={() => setEditingId(seq.id)}>
                              {seq.name}
                            </span>
                          )}
                          <FontAwesomeIcon
                            icon={faTrash}
                            className="delete-icon"
                            onClick={() => deleteSequence(seq.id)}
                          />
                        </div>
                        {seq.visible && (
                          <textarea
                            placeholder="TAGCTAGCCGATCG....."
                            value={seq.value}
                            onChange={(e) =>
                              handleSequenceChange(seq.id, e.target.value)
                            }
                            className="w-100"
                          />
                        )}
                      </div>
                    ))}
                  </Col>
                </Row>
              </Form.Group>

              <button onClick={addSequence} className="add-sequence-button">
                + Add Sequence
              </button>
              {/* -----------------다솔님 Paste Sequence ------------끝---- */}
            </Form>
          </div>

          <Button
            className="next-button"
            onClick={handleFileUploadToServer}
            disabled={!responseReceived}
          >
            Next ➔
          </Button>

          {/* Backend 없이 실행할 때 바로위 코드 주석 처리, 해당 주석 제거 후 실행하시면 됩니당, 지우지 말아주세요
          <Button
            className="next-button"
            onClick={() => {
              navigate("/analysis");
            }}
          >
            Next ➔
          </Button> */}

        </div>
      )}
    </div>
  );
}

export default InputSeq;
