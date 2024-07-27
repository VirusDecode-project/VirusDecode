import React from 'react';
import ReactDOM from 'react-dom/client';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import App from './App';
import Home from './Home';
import { GoogleOAuthProvider } from '@react-oauth/google';

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(
  <React.StrictMode>
    <Router>
    <GoogleOAuthProvider clientId="758755790796-5sf7i6gfss7m2tpvuju44tviakdghvtm.apps.googleusercontent.com">
      <Routes>
        <Route path="/" element={<Home />} />
        <Route path="/app/*" element={<App />} />
      </Routes>
      </GoogleOAuthProvider>
    </Router>
  </React.StrictMode>
);
