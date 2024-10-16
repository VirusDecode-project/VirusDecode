import React, { Dispatch, SetStateAction, MouseEvent, useState, useEffect } from "react";
import { Button, Form, Row, Col } from "react-bootstrap";
import { useNavigate } from "react-router-dom";
import "../styles/InputSeq.css";
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import {
  faChevronDown,
  faChevronRight,
  faTrash,
} from "@fortawesome/free-solid-svg-icons";
import uploadIcon from "../assets/upload_icon.png";
import Loading from '../components/Loading';
import { AlignmentData } from '../components/types';
import { setSelectionRange } from "@testing-library/user-event/dist/utils";

interface InputSeqProps {
  setTab: Dispatch<SetStateAction<number>>;
  setShow: Dispatch<SetStateAction<boolean>>;
  setWorkingHistory: Dispatch<SetStateAction<string>>;
  workingHistory: string;
  setMRNAReceived: Dispatch<SetStateAction<boolean>>;
  setPDBReceived: Dispatch<SetStateAction<boolean>>;
  setAlignmentData: Dispatch<SetStateAction<AlignmentData>>;
  setHistory: Dispatch<SetStateAction<string[]>>;
  setIsLoading: Dispatch<SetStateAction<boolean>>;
  isLoading: boolean;
  handleError: (message: string) => void;
}

interface UploadedFile {
  name: string;
  file: File;
}

const InputSeq: React.FC<InputSeqProps> = ({ setTab, setShow, setWorkingHistory, workingHistory, setMRNAReceived, setPDBReceived, setAlignmentData, setHistory, setIsLoading, isLoading, handleError }) => {
  let navigate = useNavigate();

  useEffect(() => {
    const fetchHistory = async () => {
      try {
        const serverResponse = await fetch(`/api/history/list`, {
          method: 'GET',
          credentials: 'include',
        });
        if (!serverResponse.ok) {
          throw new Error("Failed to fetch history list");
        }
        const responseData = await serverResponse.json();
        setHistory(responseData);
      } catch (error) {
        console.error("Error fetching history:", error);
      }
    };

    fetchHistory();
  }, []);


  const [editingFileIndex, setEditingFileIndex] = useState<number | null>(null);
  const [uploadedFiles, setUploadedFiles] = useState<UploadedFile[]>([]);
  const [sequences, setSequences] = useState([
    { id: 1, name: "Sequence1", value: "", visible: true },
  ]);
  const [editingId, setEditingId] = useState<number | null>(null);
  const [nextId, setNextId] = useState(2);
  const [referenceSequenceId, setReferenceSequenceId] = useState("");
  const [responseMessage, setResponseMessage] = useState("");
  const [metadata, setMetadata] = useState("");
  const [responseReceived, setResponseReceived] = useState(false);
  const [doneReceived, setDoneReceived] = useState(true);
  const [fileInputKey, setFileInputKey] = useState<number>(0); 

  const handleDoneSubmit = async (e: MouseEvent<HTMLButtonElement>) => {
    e.preventDefault(); // 폼의 기본 제출 동작 방지
    if (!referenceSequenceId.trim()) {  // 입력이 없거나 빈 문자열인 경우
      setResponseMessage("Please enter a valid sequence ID.");
      return; // 값이 없으므로 요청을 보내지 않음
    }
    const requestData = { sequenceId: referenceSequenceId };
    try {
      setDoneReceived(false);
      const serverResponse = await fetch(`/api/inputSeq/metadata`, {
        method: "POST",
        credentials: 'include',  // 세션 쿠키 포함
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
      setMetadata(responseData);
      let formattedMessage = "";
      for (const [key, value] of Object.entries(responseData)) {
        formattedMessage += `<span class="key">${key}:</span> <span class="value">${value}</span><br />`;
      }

      // 서버 응답 메시지 설정
      setResponseMessage(formattedMessage);
      setResponseReceived(true);
    } catch (error) {
      if (error instanceof Error) {
        console.error("An error occurred during the request: ", error.message);
        setResponseMessage(error.message);
      }
    } finally {
      setDoneReceived(true);
      setWorkingHistory(referenceSequenceId);
    }
  };

  // next button 클릭시 서버로 (Sequence ID, file, sequence)전송
  const handleFileUploadToServer = async (e: MouseEvent<HTMLButtonElement>) => {
    e.preventDefault(); // 폼의 기본 제출 동작 방지

    // sequences를 Map 형태로 변환
    const sequencesMap = new Map<string, string>();

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
      historyName: workingHistory,
      referenceId: JSON.parse(JSON.stringify(metadata))['Sequence ID'],
    };

    try {
      setIsLoading(true);
      const serverResponse = await fetch(`/api/inputSeq/alignment`, {
        method: "POST",
        credentials: 'include',  // 세션 쿠키 포함
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
      setAlignmentData(JSON.parse(responseData["alignment"]));
      const createdHistoryName = responseData["historyName"];
      setWorkingHistory(createdHistoryName);
      setTab(0);
      setMRNAReceived(false);
      setPDBReceived(false);
      navigate("/analysis");
    } catch (error) {
      if (error instanceof Error) {
        console.error("An error occurred during the request: ", error.message);
        // window.alert(error.message);
        handleError(error.message);
      }
    }
    finally {
      setIsLoading(false);
    }
  };


  // 파일 업로드 핸들러
  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const files = Array.from(event.target.files || []);
    const validExtensions = ['.fasta', '.fa'];

    const invalidFilesList = files.filter(file => {
      const fileExtension = file.name.slice(file.name.lastIndexOf('.'));
      return !validExtensions.includes(fileExtension);
    });

    const validFiles = files.filter(file => {
      const fileExtension = file.name.slice(file.name.lastIndexOf('.'));
      return validExtensions.includes(fileExtension);
    });

    if (invalidFilesList.length > 0) {
      handleError("FASTA 파일 형식만 지원됩니다.")
      setFileInputKey(prevKey => prevKey + 1);
    }

    const newFiles = validFiles.map((file) => ({ name: file.name, file }));
    setUploadedFiles([...uploadedFiles, ...newFiles]);

    setEditingFileIndex(null);
  };

  const handleFileNameChange = (index: number, name: string) => {
    const updatedFiles = [...uploadedFiles];
    updatedFiles[index] = { ...updatedFiles[index], name };
    setUploadedFiles(updatedFiles);
  };

  const deleteUploadedFile = (index: number) => {
    setUploadedFiles(uploadedFiles.filter((_, i) => i !== index));
  };
  const toggleVisibility = (id: number) => {
    setSequences(
      sequences.map((seq) =>
        seq.id === id ? { ...seq, visible: !seq.visible } : seq
      )
    );
  };

  const handleNameChange = (id: number, name: string) => {
    setSequences(
      sequences.map((seq) => (seq.id === id ? { ...seq, name } : seq))
    );
  };

  const deleteSequence = (id: number) => {
    setSequences(sequences.filter((seq) => seq.id !== id));
  };

  const handleSequenceChange = (id: number, value: string) => {
    setSequences(
      sequences.map((seq) => (seq.id === id ? { ...seq, value } : seq))
    );
  };
  const addSequence = (event: MouseEvent<HTMLButtonElement>) => {
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
        <Loading text="Converting" />
      ) : (

        <div>
          <div className="container mt-4 custom-container">
            <Form>
              <p className="RS-id">Reference Sequence ID</p>
              <div className="d-flex align-items-center justify-content-between">
                <Form.Group controlId="referenceSequenceId" className="flex-grow-1">
                  <Form.Control
                    type="text"
                    placeholder="Enter sequence ID"
                    value={referenceSequenceId}
                    onChange={(e) => setReferenceSequenceId(e.target.value)}
                  />
                </Form.Group>
                <Button
                  type="submit"
                  variant="primary"
                  onClick={handleDoneSubmit}
                  className="done-button"
                  disabled={!doneReceived}
                  style={{ marginLeft: '10px' }}
                >
                  DONE
                </Button>
              </div>
            </Form>




            {responseMessage && (
              <div className="response-message">
                <p className="metadata">Metadata</p>
                <div
                  className="metadata2"
                  dangerouslySetInnerHTML={{ __html: responseMessage }}
                />
              </div>
            )}

            <Form>
              <div className="mb-5"></div>

              <p className="RS-id">Variant Sequence</p>

              <Form.Group controlId="formFile" className="mb-3">
                <Form.Label>Upload File</Form.Label>
                  <Col md={10}>
                    <div className="upload-box">
                      <input
                        type="file"
                        className="file-input"
                        accept=".fasta"
                        multiple
                        key={fileInputKey}
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
              </Form.Group>

              <Form.Group className="mb-3">
                <Form.Label>Paste Sequence</Form.Label>
                  <Col md={10} className="text-left">
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
              </Form.Group>

              <button onClick={addSequence} className="add-sequence-button">
                + Add Sequence
              </button>
            </Form>
          </div>

          <Button
            className="next-button"
            onClick={handleFileUploadToServer}
            disabled={!responseReceived}
          >
            Next ➔
          </Button>

        </div>
      )}
    </div>
  );
}

export default InputSeq;
