import React, { useState } from "react";
import '../styles/Modal.css';

interface MessageModalProps {
  message: string;
  isOpen: boolean;  // 모달 열림/닫힘 상태
  onClose: () => void;  // 모달 닫기 함수
}

const MessageModal: React.FC<MessageModalProps> = ({ message, isOpen, onClose }) => {
  if (!isOpen) return null;  // 모달이 닫혀 있으면 아무것도 렌더링하지 않음

  return (
    <div className="message-modal-overlay">
      <div className="message-modal-content">
        <div></div>
        <p>{message}</p>
        <button onClick={onClose}>Close</button>
      </div>
    </div>
  );
};

export default MessageModal;
// function MessageModal({ message, isOpen, onClose }: { message: string; isOpen: boolean; onClose: () => void }) {
//   return (
//     <div className={isOpen ? "modal-overlay" : "hidden"}>
//       <div className="modal-content">
//         <p>{message}</p>
//         <button onClick={onClose}>닫기</button>
//       </div>
//     </div>
//   );
// }

// export default function Logout() {
//   const [isModalOpen, setIsModalOpen] = useState(false);
//   const [message, setMessage] = useState<string>("");

//   const handlePopupMessage = () => {
//     setIsModalOpen(true);
//     setMessage("아직은 안되지롱!");
//   };

//   const handleCloseModal = () => {
//     setIsModalOpen(false);
//   };

//   return (
//     <div>
//       <button onClick={handlePopupMessage}>로그아웃</button>
//       <PopupMessage message={message} isOpen={isModalOpen} onClose={handleCloseModal} />
//     </div>
//   );
// }
