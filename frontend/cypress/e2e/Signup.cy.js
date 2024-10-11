describe('Signup Page Test', () => {
  
  // 각 테스트 실행 전마다
  beforeEach(() => {
    cy.visit('http://localhost:3000/signup');
    
    // mocking
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
  });

  it('should successfully move to login page', () => {
    // 1. 'Back to login' 버튼을 누를 시
    cy.get('.gotoLoginBtn').click();

    // 로그인 페이지로 이동
    cy.url().should('include', '/login');
  });

  it('should successfully signup with valid data', () => {
    // 2. 올바른 회원 정보 입력 시
    cy.get('input[name="firstName"]').type('testFirstName');
    cy.get('input[name="lastName"]').type('testLastName');
    cy.get('input[name="id"]').type('testId');
    cy.get('input[name="password"]').type('testPw');
    cy.get('input[name="cPassword"]').type('testPw');

    // 패스워드 일치 시 성공 아이콘 확인
    cy.get('.icon.success').should('be.visible');

    // 폼 제출 시 회원가입 성공 alert 후 로그인 페이지로 이동
    cy.get('.SignupBtn').click();
    cy.wait('@signupRequest');
    cy.on('window:alert', (alertMessage) => {
      expect(alertMessage).to.equal('회원가입이 완료되었습니다.'); 
    });
    cy.url().should('include', '/login');
  });

  it('should show error when passwords do not match', () => {
    // 3. 'Password'와 'Confirm Password'를 다르게 입력 시
    cy.get('input[name="firstName"]').type('testFirstName');
    cy.get('input[name="lastName"]').type('testLastName');
    cy.get('input[name="id"]').type('testId');
    cy.get('input[name="password"]').type('testPw');
    cy.get('input[name="cPassword"]').type('differentTestPw');

    // 패스워드 불일치 시 오류 아이콘 및 메시지 확인
    cy.get('.icon.error').should('be.visible');
    cy.get('.signupError').should('contain', 'Passwords do not match');

    // 폼 제출 시 오류 alert: !!!!기능 수정 필요!!!! - error-message가 이미 있으니 제출을 막기
    cy.get('.SignupBtn').click();
    cy.on('window:alert', (alertMessage) => {
      expect(alertMessage).to.equal('모든 필드를 올바르게 입력해 주세요.'); 
    });
  });

  it('should show error for duplicated ID', () => {
    // 4. 기존에 가입한 ID와 중복되는 ID 입력 시 
    cy.get('input[name="firstName"]').type('testFirstName');
    cy.get('input[name="lastName"]').type('testLastName');
    cy.get('input[name="id"]').type('test_duplicated');
    cy.get('input[name="password"]').type('testPw');
    cy.get('input[name="cPassword"]').type('testPw');

    // 폼 제출 시 오류 alert: !!!!기능 수정 필요!!!! - error-message: 이미 존재하는 ID 입니다.
    cy.get('.SignupBtn').click();
    cy.on('window:alert', (alertMessage) => {
      expect(alertMessage).to.equal('이미 존재하는 ID 입니다.'); 
    });
  });

  it('should show error when required fields are missing', () => {
    // 5. 입력하지 않고 버튼 클릭 시
    // !!!!기능 수정 필요!!!! - error-message: (빈 입력창 name)을 입력하세요.
    cy.get('.SignupBtn').click();
    cy.on('window:alert', (alertMessage) => {
      expect(alertMessage).to.equal('모든 필드를 올바르게 입력해 주세요.'); 
    });
  });
});
