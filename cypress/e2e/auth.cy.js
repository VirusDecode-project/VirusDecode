describe('VirusDecode Sign-Up E2E Tests', () => {

  let uniqueID; // 중복 아이디 테스트에서 사용할 변수를 정의

  beforeEach(() => {
    // 기본 URL로 애플리케이션에 접속
    cy.visit('https://virusdecode.com');
  });

  // 1. 회원가입 성공 시나리오
  it('should sign up successfully with valid information', () => {
    // 'Try Decoding' 버튼 클릭
    cy.contains('Try Decoding').click();

    // 'Sign up' 버튼 클릭
    cy.contains('Sign up').click();

    // 랜덤 유니크한 ID 생성 (타임스탬프 사용)
    uniqueID = `user_${Date.now()}`; // 변수에 저장

    // 회원가입 폼에 유효한 정보 입력
    cy.get('input[name="firstName"]').type('John');
    cy.get('input[name="lastName"]').type('Doe');
    cy.get('input[name="id"]').type(uniqueID);
    cy.get('input[name="password"]').type('TestPassword123');
    cy.get('input[name="cPassword"]').type('TestPassword123');

    // 브라우저 알림(alert) 메시지 확인
    cy.window().then((win) => {
      cy.stub(win, 'alert').as('alert');
    });

    // 'Signup' 버튼 클릭
    cy.contains('Signup').click();

    // 회원가입 성공 메시지 또는 페이지 전환 확인
    cy.get('@alert').should('have.been.calledWith', '회원가입이 완료되었습니다.'); // 회원가입 성공 메시지를 확인합니다.
  });

  // 2. 비밀번호 불일치 시나리오
  it('should show error for password mismatch', () => {
    cy.contains('Try Decoding').click();
    cy.contains('Sign up').click();

    // 회원가입 폼에 비밀번호 불일치 정보 입력
    cy.get('input[name="firstName"]').type('John');
    cy.get('input[name="lastName"]').type('Doe');
    cy.get('input[name="id"]').type('john_doe');
    cy.get('input[name="password"]').type('TestPassword123');
    cy.get('input[name="cPassword"]').type('DifferentPassword');



    // 브라우저 알림(alert) 메시지 확인
    cy.window().then((win) => {
      cy.stub(win, 'alert').as('alert');
    });

    // 비밀번호 불일치 오류 메시지 확인
    cy.get('.signupError.visible').should('contain', 'Passwords do not match.');

    // 'Signup' 버튼 클릭
    cy.contains('Signup').click();

    // 비밀번호 불일치 오류 메시지 확인
    cy.get('@alert').should('have.been.calledWith', '모든 필드를 올바르게 입력해 주세요.');
  });

  // 3. ID 중복 시나리오
  it('should show error for duplicate ID', () => {
    cy.contains('Try Decoding').click();
    cy.contains('Sign up').click();

    // 이미 사용 중인 ID 입력
    cy.get('input[name="firstName"]').type('Jane');
    cy.get('input[name="lastName"]').type('Doe');
    cy.get('input[name="id"]').type(uniqueID); // 이미 사용된 ID
    cy.get('input[name="password"]').type('AnotherPassword123');
    cy.get('input[name="cPassword"]').type('AnotherPassword123');

    // 브라우저 알림(alert) 메시지 확인
    cy.window().then((win) => {
      cy.stub(win, 'alert').as('alert');
    });

    // 'Signup' 버튼 클릭
    cy.contains('Signup').click();

    // 중복된 ID 오류 메시지 확인
    cy.get('@alert').should('have.been.calledWith', '이미 존재하는 ID 입니다.');
  });

  // 4. 필수 필드 누락 시나리오
  it('should show error when required fields are missing', () => {
    cy.contains('Try Decoding').click();
    cy.contains('Sign up').click();


    // 필드 중 일부만 입력 (예: ID와 비밀번호만 입력)
    cy.get('input[name="id"]').type('john_doe');
    cy.get('input[name="password"]').type('TestPassword123');

    // 브라우저 알림(alert) 메시지 확인
    cy.window().then((win) => {
      cy.stub(win, 'alert').as('alert');
    });

    // 'Signup' 버튼 클릭
    cy.contains('Signup').click();

    // 필수 필드 누락 alert 메시지 확인
    cy.get('@alert').should('have.been.calledWith', '모든 필드를 올바르게 입력해 주세요.');
  });

});

