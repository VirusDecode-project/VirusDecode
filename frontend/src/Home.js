import React, { useState } from 'react';
import './Home.css';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faChevronRight } from '@fortawesome/free-solid-svg-icons';
import { useNavigate } from 'react-router-dom';
import { GoogleLogin } from '@react-oauth/google';
import { GoogleOAuthProvider } from '@react-oauth/google';

function Home() {
  const navigate = useNavigate();
  const clientId = '758755790796-5sf7i6gfss7m2tpvuju44tviakdghvtm.apps.googleusercontent.com';
  const [showLoginModal, setShowLoginModal] = useState(false);

  const handleLoginSuccess = (credentialResponse) => {
    console.log(credentialResponse);
    navigate('/app');
  };

  const handleLoginError = () => {
    console.log('Login Failed');
  };

  const handleTryButtonClick = () => {
    setShowLoginModal(true);
  };

  const handleModalClose = () => {
    setShowLoginModal(false);
  };

  return (
    <>
      <GoogleOAuthProvider clientId={clientId}>
        <div className="header">
          <img alt="VirusDecode Logo" src="/logo.png" className="logo" />
          <span className="title">VirusDecode</span>
        </div>
        <img alt="BackGroundImage" src="/mainImage.png" className="main-bg" />
        <div className="contents">
          <p className="text-box">
            Decode the virus’s genetic code,<br />
            analyze its mutations,<br />
            and determine the vaccine sequence.
          </p>
          <button className="image-button" onClick={handleTryButtonClick}>
            Try Decoding
            <FontAwesomeIcon icon={faChevronRight} className="icon" />
          </button>
        </div>

        {showLoginModal && (
          <div className="modal-backdrop" onClick={handleModalClose}>
            <div className="modal-content" onClick={(e) => e.stopPropagation()}>
              <h2>Welcome to VirusDecode</h2>
              <p>Log in to get your virus analysis records.</p>
              <div className="google-login-buttons">
                <GoogleLogin
                  text="sign in with"
                  onSuccess={handleLoginSuccess}
                  onError={handleLoginError}
                />
                <GoogleLogin
                  text="sign up with"
                  onSuccess={handleLoginSuccess}
                  onError={handleLoginError}
                  theme="filled_black"
                />
              </div>
              <button className="modal-link" onClick={() => navigate('/app')}>
                stay logged out
              </button>
            </div>
          </div>
        )}
      </GoogleOAuthProvider>
    </>
  );
}

export default Home;