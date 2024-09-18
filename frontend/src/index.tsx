import React from 'react';
import ReactDOM from 'react-dom'; // version downgrade (18->16)
import './styles/index.css';
import App from './App';
import { BrowserRouter as Router } from 'react-router-dom';
import reportWebVitals from './reportWebVitals';


ReactDOM.render(     // version downgrade (18->16)
  <React.StrictMode>
    <Router>
      <App />
    </Router>
  </React.StrictMode>,

  document.getElementById('root'),
);

reportWebVitals();
