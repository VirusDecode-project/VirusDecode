import { GoogleLogin } from "@react-oauth/google";
import { GoogleOAuthProvider } from "@react-oauth/google";
import { jwtDecode } from "jwt-decode";  // 수정된 import 구문

const GoogleLoginButton = ({setShowModal}) => {
    const clientId = '758755790796-5sf7i6gfss7m2tpvuju44tviakdghvtm.apps.googleusercontent.com';
    
    return (
        <>
            <GoogleOAuthProvider clientId={clientId}>
                <GoogleLogin
                    onSuccess={(credentialResponse) => {
                        console.log(credentialResponse);

                        // credential (JWT) 디코딩
                        const decoded = jwtDecode(credentialResponse.credential);
                        console.log("Decoded JWT:", decoded);

                        // 사용자 이름과 이메일 추출
                        const username = decoded.name;  // 사용자 이름
                        const email = decoded.email;    // 사용자 이메일
                        
                        console.log("User Name:", username);
                        console.log("User Email:", email);

                        setShowModal(false);  // 로그인 성공 후 모달 닫기
                    }}
                    onError={() => {
                        console.log('Login Failed');
                    }}
                />
            </GoogleOAuthProvider>
        </>
    );
};

export default GoogleLoginButton;
