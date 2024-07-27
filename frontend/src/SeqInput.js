import React, { useState } from 'react';
import './SeqInput.css';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faFile, faChevronDown, faChevronRight, faTrash } from '@fortawesome/free-solid-svg-icons';
import { useNavigate } from 'react-router-dom';
import LoadingPage from './LoadingPage'; // Import the loading screen component

function SeqInput() {
  const [sequences, setSequences] = useState([{ id: 1, name: 'Sequence 1', value: '', visible: true }]);
  const [nextId, setNextId] = useState(2);
  const [editingId, setEditingId] = useState(null);
  const [uploadedFiles, setUploadedFiles] = useState([]);
  const [editingFileIndex, setEditingFileIndex] = useState(null);
  const [loading, setLoading] = useState(false); // Add loading state
  const [error, setError] = useState(null); // Add error state
  const navigate = useNavigate();

  const handleSequenceChange = (id, value) => {
    setSequences(sequences.map(seq => seq.id === id ? { ...seq, value } : seq));
  };

  const handleNameChange = (id, name) => {
    setSequences(sequences.map(seq => seq.id === id ? { ...seq, name } : seq));
  };

  const handleFileNameChange = (index, name) => {
    const updatedFiles = [...uploadedFiles];
    updatedFiles[index] = { ...updatedFiles[index], name };
    setUploadedFiles(updatedFiles);
  };

  const addSequence = () => {
    setSequences([...sequences, { id: nextId, name: `Sequence ${nextId}`, value: '', visible: true }]);
    setNextId(nextId + 1);
  };

  const toggleVisibility = (id) => {
    setSequences(sequences.map(seq => seq.id === id ? { ...seq, visible: !seq.visible } : seq));
  };

  const deleteSequence = (id) => {
    setSequences(sequences.filter(seq => seq.id !== id));
  };

  const handleFileUpload = (event) => {
    const files = Array.from(event.target.files);
    const newFiles = files.map(file => ({ name: file.name, file }));
    setUploadedFiles([...uploadedFiles, ...newFiles]);
    setEditingFileIndex(null);
  };

  const deleteUploadedFile = (index) => {
    setUploadedFiles(uploadedFiles.filter((_, i) => i !== index));
  };

  const handleNextClick = () => {
    setLoading(true); // Start loading
    setError(null); // Clear any previous errors

    // Simulate a delay for loading (e.g., replace with actual analysis process)
    setTimeout(() => {
      const isSuccess = Math.random() > 0.2; // Simulate success or error (80% success rate)
      setLoading(false); // Stop loading

      if (isSuccess) {
        navigate('Tabs'); // Navigate to tabs page on success
      } else {
        setError('An error occurred while analyzing the sequences. Please try again.'); // Set error message
      }
    }, 3000); // Adjust the delay as needed
  };

  const handleModalClose = () => {
    setError(null);
  };

  const test = () => {
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

      {loading && <LoadingPage />} {/* Show loading screen when loading */}

      {error && (
        <div className="modal-backdrop" onClick={handleModalClose}>
          <div className="modal-content" onClick={(e) => e.stopPropagation()}>
            <h2>Error</h2>
            <p>{error}</p>
            <button className="modal-button" onClick={handleModalClose}>Close</button>
          </div>
        </div>
      )}
    </div>
  );
}

export default SeqInput;
