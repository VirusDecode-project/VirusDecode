export const mockLogin = () => {
  cy.intercept('POST', '/api/auth/login', (req) => {
    const { loginId, password } = req.body;

    // 유효하지 않은 ID에 대한 응답
    if (loginId === 'InvalidId') {
      req.reply({
        statusCode: 400,
        body: '유효하지 않은 ID 입니다.',
      });
    } 
    // 잘못된 비밀번호에 대한 응답
    else if (loginId === 'testId' && password === 'wrongPw') {
      req.reply({
        statusCode: 401,
        body: '비밀번호가 틀렸습니다.',
      });
    } 
    // 성공적인 로그인 응답
    else if (loginId === 'testId' && password === 'testPw') {
      req.reply({
        statusCode: 200,
        body: { message: '로그인 성공', redirectUrl: '/inputSeq' },
      });
    } 
    // 그 외의 입력은 유효하지 않음
    else {
      req.reply({
        statusCode: 400,
        body: '유효하지 않은 ID 입니다.',
      });
    }
  }).as('loginRequest');
};

export const mockSignup = () => {
  cy.intercept('POST', '/api/auth/signup', (req) => {
    const { loginId } = req.body;

    // 중복된 ID에 대한 응답을 가정
    if (loginId === 'test_duplicated') {
      req.reply({
        statusCode: 400,
        body: '이미 존재하는 ID 입니다.',
      });
    } else {
      // 성공적인 응답 가정
      req.reply({
        statusCode: 200,
        body: { message: '회원가입이 완료되었습니다.' },
      });
    }
  }).as('signupRequest');  
}