import React from 'react';

const LoginButtonComponent = () => {
    const handleLogin = () => {
        window.location.href = 'http://localhost:8080/oauth2/authorization/google';
    };

    return (
        <button className="image-button" onClick={handleLogin}>Login with Google</button>
    );
};

export default LoginButtonComponent;