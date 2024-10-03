import React, {  Dispatch, SetStateAction, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import '../styles/Home.css';
import LoginModal from '../components/LoginModal';
import { useRecoilValue } from "recoil";
import { authState } from "../state/authState";
interface HomeProps {
  setUserName:Dispatch<SetStateAction<string | null>>;
}

const Home: React.FC<HomeProps> = ({setUserName}) => {
  let navigate = useNavigate();
  const [isLoginModalOpen, setIsLoginModalOpen] = useState(false);
  const isLoggedIn = useRecoilValue(authState);

  const handleDecodeBtn = async(event: React.MouseEvent<HTMLButtonElement>) => {
    try {
      const nameResponse = await fetch(`${process.env.REACT_APP_BACKEND_URL}/api/auth/userinfo`, {
        method: "POST",
        credentials: 'include',
        headers: {
          "Content-Type": "application/json",
        },
      });
  
      if (nameResponse.ok) { // 로그인 상태라면 inputSeq로 넘어감
        navigate("/inputSeq");
      }else{ // 로그인 정보가 없으면 로그인 모달 오픈
        setIsLoginModalOpen(true);
      }
    } catch (error) {
      if(error instanceof Error){
        window.alert(error.message);
      }
    }
  }
  return (
    <div>
      <div className="main-bg"></div>
      <div className="text-box">
        <p>
          Decode the virus’s genetic code, analyze its mutations, and determine the vaccine sequence.
        </p>
      </div>
      <button className="decode-button" onClick={handleDecodeBtn}>
        Try Decoding
      </button>
          {isLoginModalOpen && (
            <LoginModal
              isOpen={isLoginModalOpen}
              onClose={() => {
                setIsLoginModalOpen(false);
              }
            }
            setUserName={setUserName}
            /> 
          )}
    </div>
  );
};

export default Home;
