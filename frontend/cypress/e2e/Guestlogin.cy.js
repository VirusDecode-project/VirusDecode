describe('1. 비회원 이용', () => {
  const originHistoryName1 = 'SARS-CoV-2';
  const originHistoryName2 = 'RSV';

  it('1-1. 게스트 로그인 처리', () => {
    cy.visit('http://localhost:3000/');
    cy.get('.decode-button').click();
    cy.get('.stayLoggedOutBtn').click();
    
    // 회원 정보가 게스트인지 확인
    cy.intercept('POST', '/api/auth/userinfo').as('userinfoRequest');
    cy.get('.user-icon').click();
    cy.wait('@userinfoRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('.user-icon').click();
      cy.get('.userInfo-menu').invoke('show')
      .should('be.visible')
      .and('contain', 'Guest')
      
      cy.get('.user-icon').click();
      cy.get('.userInfo-menu').should('not.exist');
    });
  });

  it('1-2. 샘플 데이터 히스토리에 저장 ', () => {
    cy.guestlogin();
    cy.get('.history-item').first().should('be.visible').and('contain', originHistoryName1);
    cy.get('.history-item').eq(1).should('be.visible').and('contain', originHistoryName2);
  });
});