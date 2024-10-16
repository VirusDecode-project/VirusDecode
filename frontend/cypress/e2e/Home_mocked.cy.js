describe('Home Page Test', () => {
  let userInfoResponse;
  beforeEach(() => {
    cy.visit('http://localhost:3000/');

    cy.intercept('POST', '/api/auth/userinfo', (req) => {
      req.reply((res) => {
        res.send(userInfoResponse);
      });
    }).as('getUserInfo');

    cy.intercept('POST', '/api/auth/guest-login', {
      statusCode: 200,
      body: {
        success: true
      }
    }).as('guestLogin');
  });

  it('should successfully open and close LoginModal', () => {
    // #1 로그인 모달이 열리고 닫히는지 테스트
    userInfoResponse = { statusCode: 401 };
    cy.get('.decode-button').click();
    cy.get('.login-modal').should('be.visible');
    
    // close button
    cy.get('.close-button').click();
    cy.get('.login-modal').should('not.exist');

    cy.get('.decode-button').click();
    cy.get('.login-modal').should('be.visible');

    // ESC key
    cy.get('body').type('{esc}');
    cy.get('.login-modal').should('not.exist');

    cy.get('.decode-button').click();
    cy.get('.login-modal').should('be.visible');

    // 모달 바깥쪽 클릭
    cy.get('body').click(0, 0);
    cy.get('.login-modal').should('not.exist');
  });

  it('should successfully move to login page', () => {
    // #2 LoginModal에서 'Log in' 버튼을 누를 시
    userInfoResponse = { statusCode: 401 };
    cy.get('.decode-button').click();
    cy.get('.login-modal').should('be.visible');
    cy.get('.loginBtn_modal').click();

    // 로그인 페이지로 이동
    cy.url().should('include', '/login');
  });

  it('should successfully move to signup page', () => {
    // #3 LoginModal에서 'Sign up' 버튼을 누를 시
    userInfoResponse = { statusCode: 401 };
    cy.get('.decode-button').click();
    cy.get('.login-modal').should('be.visible');
    cy.get('.signupBtn_modal').click();
    
    // 회원가입 페이지로 이동
    cy.url().should('include', '/signup');
  });

  it('should successfully stay logged out and move to inputSeq', () => {
    // #4 LoginModal에서 'stay logged out' 버튼을 누를 시
    userInfoResponse = { statusCode: 401 };
    cy.get('.decode-button').click();
    cy.get('.login-modal').should('be.visible');

    cy.get('.stayLoggedOutBtn').click(); 
    cy.wait('@guestLogin').its('response.statusCode').should('eq', 200);
    
    // inputSeq 페이지로 이동
    cy.url().should('include', '/inputSeq');
  });

  it('should move to inputSeq without LoginModal if Logined', () => {
    // #5 홈 페이지에서 LoginModal없이 바로 inputSeq로 넘어가는지 테스트
    userInfoResponse = { statusCode: 200 };
    cy.get('.decode-button').click();
    cy.get('.login-modal').should('not.exist');
    cy.url().should('include', '/inputSeq');
  }); 
});