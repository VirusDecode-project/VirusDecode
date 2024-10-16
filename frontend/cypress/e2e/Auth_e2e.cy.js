// 계정 생성 및 저장: 회원 가입 완료 #1
// 계정 생성 및 저장: 계정 데이터 저장 #1
// 로그인: 로그인 성공 #2

// 이 테스트 코드는 한 번 실행하면 데이터베이스에 계정 정보가 저장됩니다.
describe('Signup to Login e2e', () => {
  it('should successfully signup and login', () => {
    // 새 계정 생성
    cy.visit('http://localhost:3000/signup');
    cy.get('input[name="firstName"]').type('testFirstName');
    cy.get('input[name="lastName"]').type('testLastName');
    cy.get('input[name="loginId"]').type('testId_e2e');
    cy.get('input[name="password"]').type('testPw');
    cy.get('input[name="cPassword"]').type('testPw');

    // #1 폼 제출 시 회원가입 완료
    cy.get('.SignupBtn').click();
    cy.get('.message-modal-content')
    .should('be.visible')
    .and('contain', '회원가입이 완료되었습니다.');
    cy.get('.message-modal-content').contains('Close').click();

    // #2 로그인 페이지로 이동, 가입된 계정으로 로그인
    cy.url().should('include', '/login');
    cy.get('input[name="loginId"]').type('test_e2e');
    cy.get('input[name="password"]').type('testPw');
    cy.get('.loginBtn').click();

    // 로그인 성공 후 inputSeq 페이지로 이동
    cy.url().should('include', '/inputSeq');
    
    // 같은 정보로 회원가입 시도
    cy.visit('http://localhost:3000/signup');
    cy.get('input[name="firstName"]').type('testFirstName');
    cy.get('input[name="lastName"]').type('testLastName');
    cy.get('input[name="loginId"]').type('test_e2e');
    cy.get('input[name="password"]').type('testPw');
    cy.get('input[name="cPassword"]').type('testPw');

    // 아이디 중복 검사
    cy.get('.SignupBtn').click();
    cy.get('.message-modal-content') 
    .should('be.visible') 
    .and('contain', '이미 존재하는 ID 입니다.');
  });
});