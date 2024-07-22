import React, { useState } from 'react';
import './SeqInput.css';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faFile, faChevronDown, faChevronRight, faTrash } from '@fortawesome/free-solid-svg-icons';
import { useNavigate } from 'react-router-dom';

function SeqInput() {
  // 시퀀스를 저장하는 상태, 초기 상태로 하나의 시퀀스를 가지고 있음
  const [sequences, setSequences] = useState([{ id: 1, name: 'Sequence 1', value: '', visible: true }]);
  const [nextId, setNextId] = useState(2); // 다음 시퀀스 이름
  const [editingId, setEditingId] = useState(null); // 현재 편집 중인 시퀀스 이름
  const [uploadedFiles, setUploadedFiles] = useState([]); // 업로드된 파일들을 추적
  const [editingFileIndex, setEditingFileIndex] = useState(null); // 파일 이름을 편집 중인 인덱스
  const navigate = useNavigate();
  // 시퀀스 값이 변경될 때 호출되는 함수
  const handleSequenceChange = (id, value) => {
    setSequences(sequences.map(seq => seq.id === id ? { ...seq, value } : seq));
  };

  // 시퀀스 이름이 변경될 때 호출되는 함수
  const handleNameChange = (id, name) => {
    setSequences(sequences.map(seq => seq.id === id ? { ...seq, name } : seq));
  };

  // 파일 이름이 변경될 때 호출되는 함수
  const handleFileNameChange = (index, name) => {
    const updatedFiles = [...uploadedFiles];
    updatedFiles[index] = { ...updatedFiles[index], name };
    setUploadedFiles(updatedFiles);
  };

  // 새로운 시퀀스를 추가하는 함수
  const addSequence = () => {
    setSequences([...sequences, { id: nextId, name: `Sequence ${nextId}`, value: '', visible: true }]);
    setNextId(nextId + 1); // 다음 시퀀스 ID를 증가시킴
  };

  // 시퀀스 내용 토글하는 함수
  const toggleVisibility = (id) => {
    setSequences(sequences.map(seq => seq.id === id ? { ...seq, visible: !seq.visible } : seq));
  };

  // 시퀀스 삭제
  const deleteSequence = (id) => {
    setSequences(sequences.filter(seq => seq.id !== id));
  };

  // 파일 업로드
  const handleFileUpload = (event) => {
    const files = Array.from(event.target.files);
    const newFiles = files.map(file => ({ name: file.name, file }));
    setUploadedFiles([...uploadedFiles, ...newFiles]);
    setEditingFileIndex(null);
  };

  // 업로드된 파일 삭제
  const deleteUploadedFile = (index) => {
    setUploadedFiles(uploadedFiles.filter((_, i) => i !== index));
  };

    // Next 버튼 클릭 처리
    const handleNextClick = () => {
      // Next 버튼 클릭 후 로딩창을 띄우거나 분석 결과 페이지로 가는 과정 추가
      navigate('LoadingPage');
    };
    const test = () => {
      // 테스트용 함수
      navigate('Tabs');
    };

  return (
    <div className="container">
      <div className="form-group">
        <label>Reference Sequence ID</label>
        <input type="text" placeholder="Enter sequence ID" className="form-control" />
        <button className="done-button">DONE</button>
      </div>
      
      <div className="form-group">
        <label>Variant Sequence</label>
        <div className='underline'>Upload File</div>
        <div className="upload-box">
          <input type="file" className="file-input" accept=".fasta" multiple onChange={handleFileUpload} />
          <div className="upload-text"><FontAwesomeIcon className='file-icon' icon={faFile} /><p/>Drag your FASTA files here</div>
        </div>
        {uploadedFiles.map((uploadedFile, index) => (
          <div key={index} className="uploaded-file">
            {editingFileIndex === index ? (
              <input
                type="text"
                value={uploadedFile.name}
                onChange={(e) => handleFileNameChange(index, e.target.value)}
                onBlur={() => setEditingFileIndex(null)}
                className="edit-file-name-input"
                autoFocus
              />
            ) : (
              <span onClick={() => setEditingFileIndex(index)}>{uploadedFile.name}</span>
            )}
            <FontAwesomeIcon icon={faTrash} className="delete-icon" onClick={() => deleteUploadedFile(index)} />
          </div>
        ))}
      </div>
      <div className='underline'>Paste Sequence</div>
      
      {sequences.map(seq => (
        <div key={seq.id} className="form-group">
          <div className="sequence-header">
            <FontAwesomeIcon icon={seq.visible ? faChevronDown : faChevronRight} className="chevron-icon" onClick={() => toggleVisibility(seq.id)} />
            {editingId === seq.id ? (
              <input
                type="text"
                value={seq.name}
                onChange={(e) => handleNameChange(seq.id, e.target.value)}
                onBlur={() => setEditingId(null)}
                className="edit-name-input"
                autoFocus
              />
            ) : (
              <span onClick={() => setEditingId(seq.id)}>{seq.name}</span>
            )}
            <FontAwesomeIcon icon={faTrash} className="delete-icon" onClick={() => deleteSequence(seq.id)} />
          </div>
          {seq.visible && (
            <textarea
              placeholder="TAGCTAGCCGATCG....."
              value={seq.value}
              onChange={(e) => handleSequenceChange(seq.id, e.target.value)}
            />
          )}
        </div>
      ))}
      
      <button onClick={addSequence} className="add-sequence-button">+ Add Sequence</button>
      <button className="next-button" onClick={handleNextClick}>Next ➔</button>
      <button className="next-button" onClick={test}>시연용: 탭 페이지 버튼</button>
    </div>
  );
}

export default SeqInput;

{/*중복되는 파일이름 업로드 안됨, 같은 파일 드래그 업로드롸 클릭 업로드의 중복 업로드 문제,   */}