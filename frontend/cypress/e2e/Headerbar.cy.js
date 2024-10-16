// 로그아웃: 로그아웃 및 메인페이지 이동 #3
// 로고 및 사이트명: 메인페이지로 이동 #1
// 회원 아이콘: 회원 프로필 정보 확인 #2
// 히스토리 아이콘: 히스토리 여닫기 #4
// 편집 아이콘: 서열 입력 페이지로 이동 #5
describe('HeaderBar Test', () => {
  const guestlogin = () => {
    cy.visit('http://localhost:3000/');
    cy.get('.decode-button').click();
    cy.get('.stayLoggedOutBtn').click(); 
  }
  beforeEach(() => {
    guestlogin();
    cy.visit('http://localhost:3000/inputSeq');
    cy.intercept('POST', '/api/auth/userinfo').as('userinfoRequest');
  });

  it('navigates to home when logo is clicked', () => {
    // #1 home page로 이동
    cy.get('.logo-text').should('be.visible').click();
    cy.url().should('include', '/');
  });

  it('toggles user info menu when user icon is clicked (logged in)', () => {
    // #2 로그인 상태에서 userInfo 아이콘 클릭 시 메뉴 열림
    cy.get('.user-icon').click({ multiple: true });
    cy.wait('@userinfoRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('.userInfo-menu')
      .should('be.visible')
      .and('contain', 'Guest');
      
      // 다시 클릭하여 닫힘
      cy.get('.user-icon').click();
      cy.get('.userInfo-menu').should('not.exist');
    });
  });

  it('logs out and redirects to home after clicking logout button (logged in)', () => {
    // #3 로그아웃 버튼 클릭 시 로그아웃 처리 및 홈으로 이동
    cy.get('.user-icon').click({ multiple: true });
    cy.wait('@userinfoRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('.userInfo-menu').should('be.visible');
      cy.get('.logoutBtn').click();
      cy.get('.userInfo-menu').should('not.exist');
      cy.url().should('include', '/');
    });
  });

  it('open and close sidebar when history icon is clicked', () => {
    // #4 사이드바 여닫기
    cy.get('.sidebar .history-icon').click();
    cy.get('.sidebar').should('not.have.class', 'show');
    cy.get('.header-bar .history-icon').click();
    cy.get('.sidebar').should('have.class', 'show');
  });

  it('open history-modal when edit icon is clicked', () => {
    // #5 restart
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