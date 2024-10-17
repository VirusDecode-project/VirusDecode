describe('NCBI 레퍼런스 시퀀스 ID 입력 테스트', () => {

  beforeEach(() => {
    // 기본 URL로 애플리케이션에 접속
    cy.visit('http://localhost:3000');
    cy.contains('Try Decoding').click();
    cy.contains('stay logged out').click();
  });


  // 1. 시나리오 ID: TS_003 _1
  it('올바른 ID 입력 시: 시스템에 등록되고 "Sequence ID, Name, Description, Length”에 대한 정보가 나타남', () => {


    // NCBI 레퍼런스 시퀀스 ID 입력 필드에 NC_045512.2을 입력한다. ‘Done’ 버튼을 누른다.
    cy.get('input[id="referenceSequenceId"]').type('NC_045512.2');
    cy.get('button.done-button').click();


    // "Sequence ID", "Name", "Description", "Length"가 있는지 확인
    cy.contains('Sequence ID').should('be.visible');
    cy.contains('Name').should('be.visible');
    cy.contains('Description').should('be.visible');
    cy.contains('Length').should('be.visible');
  });

  // 2. 시나리오 ID: TS_003 _2
  it('잘못된 ID 입력 시: "NCBI에 요청한 nucleotide ID가 존재하지 않습니다."라는 메시지가 나타남.', () => {


    // 동일한 단계에서 잘못된 ID인 ATCG을 입력한다. ‘Done’ 버튼을 누른다.
    cy.get('input[id="referenceSequenceId"]').type('ATCG');
    cy.get('button.done-button').click();


    // 잘못된 ID 입력 시: "NCBI에 요청한 nucleotide ID가 존재하지 않습니다."라는 메시지가 나타남.
    cy.contains('NCBI에 요청한 nucleotide ID가 존재하지 않습니다.').should('be.visible');

  });

  // 3. 시나리오 ID: TS_003 _3
  it('10.	입력하지 않을 경우: “Please enter a valid sequence ID.”라는 메시지가 나타남.', () => {


    // 동일한 단계에서 잘못된 ID인 ATCG을 입력한다. ‘Done’ 버튼을 누른다.
    cy.get('button.done-button').click();


    // 잘못된 ID 입력 시: "NCBI에 요청한 nucleotide ID가 존재하지 않습니다."라는 메시지가 나타남.
    cy.contains('Please enter a valid sequence ID.').should('be.visible');

  });

});