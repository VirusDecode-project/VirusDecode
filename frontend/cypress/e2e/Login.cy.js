describe('Login Page Test', () => {

    // 각 테스트 실행 전마다 로그인 페이지로 이동
  beforeEach(() => {
    cy.visit('http://localhost:3000/login');

    cy.intercept('POST', '/api/auth/login', (req) => {
      const { loginId, password } = req.body;

      // 1. 유효하지 않은 ID에 대한 응답
      if (loginId === 'InvalidId') {
        req.reply({
          statusCode: 400,
          body: '유효하지 않은 ID 입니다.',
        });
      } 
      // 2. 잘못된 비밀번호에 대한 응답
      else if (loginId === 'testId' && password === 'wrongPw') {
        req.reply({
          statusCode: 401,
          body: '비밀번호가 틀렸습니다.',
        });
      } 
      // 3. 성공적인 로그인 응답
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
  });

  it('should successfully move to signup page', () => {
    // 1. 'Create a new account' 버튼을 누를 시
    cy.get('.gotoSignupBtn').click();

    // 회원가입 페이지로 이동
    cy.url().should('include', '/signup');
  });

  it('should successfully login with valid data', () => {
    // 2. 회원가입 완료된 ID와 Password 입력 시
    cy.get('input[name="loginId"]').type('testId');
    cy.get('input[name="password"]').type('testPw');
    cy.get('.loginBtn').click();

    cy.wait('@loginRequest').its('response.statusCode').should('eq', 200);
    // inputSeq 페이지로 이동
    cy.url().should('include', '/inputSeq');
  });

  it('should show error for invalid ID', () => {
    // 3. 가입하지 않은 ID를 입력 시
    cy.get('input[name="loginId"]').type('InvalidId');
    cy.get('input[name="password"]').type('testPw');
    cy.get('.loginBtn').click();

    cy.wait('@loginRequest').its('response.statusCode').should('eq', 400);
    // error alert
    cy.on('window:alert', (alertMessage) => {
      expect(alertMessage).to.equal('유효하지 않은 ID 입니다.'); 
    });
  });

  it('should show error for wrong ID', () => {
    // 4. 가입이 완료된 ID와 잘못된 Password를 입력 시
    cy.get('input[name="loginId"]').type('testId');
    cy.get('input[name="password"]').type('wrongPw');
    cy.get('.loginBtn').click();

    cy.wait('@loginRequest').its('response.statusCode').should('eq', 401);
    // error alert
    cy.on('window:alert', (alertMessage) => {
      expect(alertMessage).to.equal('비밀번호가 틀렸습니다.'); 
    });
  });

  it('should show error for empty form', () => {
    // 5.	아무것도 입력하지 않을 경우
    cy.get('.loginBtn').click();
    
    cy.wait('@loginRequest').its('response.statusCode').should('eq', 400);
    // error alert
    cy.on('window:alert', (alertMessage) => {
      expect(alertMessage).to.equal('유효하지 않은 ID 입니다.'); 
    });
  });
});
