// 비회원 이용: 게스트 로그인 처리 #1
// 비회원 이용: 샘플 데이터 불러오기 #2
// 히스토리 선택: 저장된 히스토리 데이터 불러오기 #3
// 히스토리 수정: 히스토리 이름 수정 #6
// 히스토리 수정: 히스토리 삭제 #7
describe('SideBar Test', () => {
  const referenceSeqId = 'NC_045512';
  const fileName1 = 'SARS_CoV_2/MT576556.1.spike.fasta';
  const originHistoryName1 = 'SARS-CoV-2';
  const originHistoryName2 = 'RSV';
  const guestlogin = () => {
    cy.visit('http://localhost:3000/');
    cy.get('.decode-button').click();
    cy.get('.stayLoggedOutBtn').click(); 
  }
  
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
    // #1 게스트 로그인
    guestlogin();
    cy.intercept('POST', '/api/inputSeq/metadata').as('metadataRequest');
    cy.intercept('POST', '/api/inputSeq/alignment').as('alignmentRequest');
  })

  it('get sample history list when guest login ', () => {
    // #2 게스트 로그인 시 샘플 히스토리 데이터 조회
    cy.get('.history-item').first().should('be.visible').and('contain', originHistoryName1);
    cy.get('.history-item').eq(1).should('be.visible').and('contain', originHistoryName2);
  });

  it('get data when history clicked', () => {
    // #3 히스토리 선택 시 데이터 불러오기
    cy.get('.history-item').first().should('be.visible').and('contain', originHistoryName1);
    cy.get('.history-item').eq(1).should('be.visible').and('contain', originHistoryName2);
    // 첫 번째 히스토리 불러오기
    cy.get('.history-item').first().click();
    cy.url().should('include', '/analysis');
    cy.get('.nav-tabs').contains('Alignment').should('not.have.class', 'disabled'); 
    cy.get('.nav-tabs').contains('mRNA').should('not.have.class', 'disabled'); 
    cy.get('.nav-tabs').contains('3D').should('not.have.class', 'disabled'); 
  });
  

  it('should rename a history item', () => {
    // #6 히스토리 이름 수정
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

  it('should delete a history item', () => {
    // #7 히스토리 삭제
    // cancel
    openDeleteModal();
    cy.get('.modal-close-button').click();

    // delete
    openDeleteModal();
    cy.get('.modal-next-button').click();
    cy.get('.history-item').first().should('not.contain', originHistoryName1);
  });
});