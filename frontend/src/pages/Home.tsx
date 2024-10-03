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

  return (
    <div>
      <div className="main-bg"></div>
      <div className="text-box">
        <p>
          Decode the virus’s genetic code, analyze its mutations, and determine the vaccine sequence.
        </p>
      </div>
      <button className="decode-button" onClick={(e) => setIsLoginModalOpen(true)}>
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
