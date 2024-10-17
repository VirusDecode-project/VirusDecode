describe('1. 히스토리 선택', () => {
  const originHistoryName1 = 'SARS-CoV-2';
  const originHistoryName2 = 'RSV';
  
  it('1-1. 저장된 히스토리 데이터 불러오기', () => {
    cy.guestlogin();
    cy.get('.history-item').first().should('be.visible').and('contain', originHistoryName1);
    cy.get('.history-item').eq(1).should('be.visible').and('contain', originHistoryName2);
    // 첫 번째 히스토리 불러오기
    cy.get('.history-item').first().click();
    cy.url().should('include', '/analysis');
    cy.get('.nav-tabs').contains('Alignment').should('not.have.class', 'disabled'); 
    cy.get('.nav-tabs').contains('mRNA').should('not.have.class', 'disabled'); 
    cy.get('.nav-tabs').contains('3D').should('not.have.class', 'disabled'); 
  });
});

describe('2. 히스토리 수정', () => {
  const originHistoryName1 = 'SARS-CoV-2';
  const originHistoryName2 = 'RSV';
  const openRenameModal = () => {
    cy.get('.history-item').first().should('be.visible').and('contain', originHistoryName1);
    cy.get('.ellipsis-button').first().invoke('show').should('be.visible').click();
    cy.get('.option-button').contains('Rename').click();
    cy.get('.history-modal-content').contains('Please enter a new name.').should('be.visible');
  }
  const openDeleteModal = () => {
    cy.get('.history-item').first().should('be.visible').and('contain', originHistoryName1);
    cy.get('.ellipsis-button').first().invoke('show').should('be.visible').click();
    cy.get('.option-button').contains('Delete').click();
    cy.get('.history-modal-content').contains('delete this history?').should('be.visible');
  }
  beforeEach(() => {
    cy.guestlogin();
  })
  it('2-1. 히스토리 이름 수정', () => {
    const newHistoryName = 'testRename';
    // cancel
    openRenameModal();
    cy.get('.modal-close-button').click();

    // empty input rename
    openRenameModal();
    cy.get('.modal-next-button').click();
    //error
    cy.get('.message-modal-content').should('be.visible').and('contain', 'Please enter correct name.');
    cy.get('.message-modal-content').contains('Close').click();

    // rename
    openRenameModal();
    cy.get('input[name="newName"]').type(newHistoryName);
    cy.get('.modal-next-button').click();
    // rename 확인
    cy.get('.history-item').first().should('contain', newHistoryName);
  });

  it('2-2. 히스토리 삭제', () => {
    // cancel
    openDeleteModal();
    cy.get('.modal-close-button').click();

    // delete
    openDeleteModal();
    cy.get('.modal-next-button').click();
    cy.get('.history-item').first().should('not.contain', originHistoryName1);
  });
});