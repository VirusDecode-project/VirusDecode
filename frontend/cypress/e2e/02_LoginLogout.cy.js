describe('1. 로그인', () => {
  beforeEach(() => {
    // 로그인 페이지로 이동
    cy.visit('http://localhost:3000/');
    cy.get('.decode-button').click();
    cy.get('.login-modal').should('be.visible');
    cy.get('.loginBtn_modal').click();
    cy.url().should('include', '/login');
    // mock
    cy.mockLogin();
  });
  
  it('1-1. 회원 아이디 검사', () => {
    cy.get('input[name="loginId"]').type('InvalidId');
    cy.get('input[name="password"]').type('testPw');
    cy.get('.loginBtn').click();

    cy.wait('@loginRequest').its('response.statusCode').should('eq', 400);
    // error alert
    cy.get('.message-modal-content') 
    .should('be.visible') 
    .and('contain', '유효하지 않은 ID 입니다.');
  });

  it('1-2. 회원 비밀번호 검사', () => {
    cy.get('input[name="loginId"]').type('testId');
    cy.get('input[name="password"]').type('wrongPw');
    cy.get('.loginBtn').click();

    cy.wait('@loginRequest').its('response.statusCode').should('eq', 401);
    // error alert
    cy.get('.message-modal-content') 
    .should('be.visible') 
    .and('contain', '비밀번호가 틀렸습니다.');
  });

  it('1-3. 로그인 성공', () => {
    cy.get('input[name="loginId"]').type('testId');
    cy.get('input[name="password"]').type('testPw');
    cy.get('.loginBtn').click();

    cy.wait('@loginRequest').its('response.statusCode').should('eq', 200);
    // inputSeq 페이지로 이동
    cy.url().should('include', '/inputSeq');
  });
  
  // it('should successfully move to signup page', () => {
  //   // #1 'Create a new account' 버튼을 누를 시
  //   cy.get('.gotoSignupBtn').click();

  //   // 회원가입 페이지로 이동
  //   cy.url().should('include', '/signup');
  // });

  // it('should show error for empty form', () => {
  //   // #5	아무것도 입력하지 않을 경우
  //   cy.get('.loginBtn').click();
    
  //   cy.wait('@loginRequest').its('response.statusCode').should('eq', 400);
  //   // error alert
  //   cy.get('.message-modal-content') 
  //   .should('be.visible') 
  //   .and('contain', '유효하지 않은 ID 입니다.');
  // });
});

describe('2. 로그아웃', () => {
  it('2-1. 로그아웃 및 메인페이지 이동', () => {
    cy.signupAndLoginIfDuplicate('testFName','testLName','testId','testPw','testPw');
    cy.intercept('POST', '/api/auth/userinfo').as('userinfoRequest');
    cy.get('.user-icon').click();
    cy.wait('@userinfoRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('.user-icon').click();
      cy.get('.userInfo-menu').invoke('show')
      .should('be.visible')
      cy.get('.logoutBtn').click();
      cy.get('.userInfo-menu').should('not.exist');
      cy.url().should('include', '/');
    });
  });
});