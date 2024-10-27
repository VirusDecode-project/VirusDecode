describe('1. 로고 및 사이트명', () => {
  it('1-1. 메인페이지로 이동', () => {
    cy.signupAndLoginIfDuplicate('testFName','testLName','testId','testPw','testPw');
    cy.get('.logo-text').should('be.visible').click();
    cy.url().should('include', '/');
  });
});

describe('2. 회원 아이콘', () => {
  it('2-1. 회원 프로필 정보 확인', () => {
    cy.signupAndLoginIfDuplicate('testFName','testLName','testId','testPw','testPw');
    cy.intercept('POST', '/api/auth/userinfo').as('userinfoRequest');
    cy.get('.user-icon').click();
    cy.wait('@userinfoRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('.user-icon').click();
      cy.get('.userInfo-menu').invoke('show')
      .should('be.visible')
      .and('contain', 'testFName')
      .and('contain', 'testId');
      
      // 다시 클릭하여 닫힘
      cy.get('.user-icon').click();
      cy.get('.userInfo-menu').should('not.exist');
    });
  });
});

describe('3. 히스토리 아이콘', () => {
  it('3-1. 히스토리 여닫기', () => {
    cy.signupAndLoginIfDuplicate('testFName','testLName','testId','testPw','testPw');
    cy.get('.sidebar .history-icon').click();
    cy.get('.sidebar').should('not.have.class', 'show');
    cy.get('.header-bar .history-icon').click();
    cy.get('.sidebar').should('have.class', 'show');
  });
});

describe('4. 편집 아이콘', () => {
  it('4-1. 서열 입력 페이지로 이동', () => {
    cy.signupAndLoginIfDuplicate('testFName','testLName','testId','testPw','testPw');
    cy.visit('http://localhost:3000/analysis');
    // 사이드바가 닫혀 있을 때 헤더바의 edit-icon 클릭
    cy.get('.sidebar .history-icon').click();
    cy.get('.sidebar').should('not.have.class', 'show');

    // open history-modal
    cy.get('.header-bar .edit-icon').click();
    cy.get('.history-modal-content').should('be.visible');
    // cancel
    cy.get('.modal-close-button').click();
    cy.url().should('include', '/analysis');

    // open history-modal
    cy.get('.header-bar .edit-icon').click();
    cy.get('.history-modal-content').should('be.visible');
    // restart
    cy.get('.modal-next-button').click();
    cy.url().should('include', '/InputSeq');
  });
});