import React from 'react';
import ReactDOM from 'react-dom'; // version downgrade (18->16)
import './index.css';
import App from './App';
import { BrowserRouter as Router } from 'react-router-dom';


ReactDOM.render(     // version downgrade (18->16)
  <React.StrictMode>
    <Router>
      <App />
    </Router>
  </React.StrictMode>,

  document.getElementById('root'),
);