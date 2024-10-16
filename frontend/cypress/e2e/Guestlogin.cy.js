describe('1. 비회원 이용', () => {
  const originHistoryName1 = 'SARS-CoV-2';
  const originHistoryName2 = 'RSV';

  it('1-1. 게스트 로그인 처리', () => {
    cy.visit('http://localhost:3000/');
    cy.get('.decode-button').click();
    cy.get('.stayLoggedOutBtn').click();
  });

  it('1-2. 샘플 데이터 히스토리에 저장 ', () => {
    cy.guestlogin();
    cy.get('.history-item').first().should('be.visible').and('contain', originHistoryName1);
    cy.get('.history-item').eq(1).should('be.visible').and('contain', originHistoryName2);
  });
});