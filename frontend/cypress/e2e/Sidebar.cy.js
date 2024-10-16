// 비회원 이용: 게스트 로그인 처리 #1
// 비회원 이용: 샘플 데이터 불러오기 #2
// 히스토리 선택: 저장된 히스토리 데이터 불러오기 #3
// 히스토리 수정: 히스토리 이름 수정 #6
// 히스토리 수정: 히스토리 삭제 #7
// 히스토리 생성: 전송된 입력 데이터 저장 #4
// 히스토리 생성: 히스토리 이름 자동 생성 #4, #5
describe('SideBar Test', () => {
  const referenceSeqId = 'NC_045512';
  const fileName1 = 'MT576556.1.spike.fasta';
  const originHistoryName1 = 'SARS-CoV-2';
  const originHistoryName2 = 'RSV';
  const guestlogin = () => {
    cy.visit('http://localhost:3000/');
    cy.get('.decode-button').click();
    cy.get('.stayLoggedOutBtn').click(); 
  }
  const filesetup = () => {
    cy.get('input#referenceSequenceId') 
      .type(referenceSeqId)
      .should('have.value', referenceSeqId); 
    cy.get('button').contains('DONE').click(); 
    cy.wait('@metadataRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('input[type="file"]').attachFile([fileName1]);   
    cy.get('button.next-button').click();
    cy.url().should('include', '/analysis');
    });
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
  
  it('save new history', () => {
    // #4 사용자가 입력한 데이터를 새 히스토리에 저장
    filesetup();
    // 자동 생성 이름: 레퍼런스ID
    cy.get('.history-item').first().should('be.visible').and('contain', referenceSeqId);
  });

  it('save new history_duplicated referenceSeqId', () => {
    // #5 같은 레퍼런스ID 입력 후 분석 시작 시 히스토리 저장
    filesetup();
    cy.get('.sidebar .edit-icon').click();
    cy.get('.modal-next-button').click();
    filesetup();
    // 자동 생성 이름: 레퍼런스ID_1
    cy.get('.history-item').first().should('be.visible').and('contain', referenceSeqId + '_1');
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