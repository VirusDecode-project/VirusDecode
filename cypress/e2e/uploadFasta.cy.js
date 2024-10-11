describe('FASTA 파일 또는 시퀀스 업로드', () => {

  beforeEach(() => {
    // 기본 URL로 애플리케이션에 접속
    cy.visit('http://localhost:3000');
    cy.contains('Try Decoding').click();
    cy.contains('stay logged out').click();
  });

  // 시나리오 ID: TS_004_1.1. 한 개의 파일 업로드시
  it('사용자로부터 FASTA 파일 한개 업로드', () => {
    const filePath = 'fasta_example_3.fasta'; // fixtures 폴더에 있는 파일 경로

    // NCBI 레퍼런스 시퀀스 ID 입력 필드에 NC_045512.2을 입력한다. ‘Done’ 버튼을 누른다.
    cy.get('input[id="referenceSequenceId"]').type('NC_045512.2');
    cy.get('button.done-button').click();

    // "Sequence ID", "Name", "Description", "Length"가 있는지 확인
    cy.contains('Sequence ID').should('be.visible');
    cy.contains('Name').should('be.visible');
    cy.contains('Description').should('be.visible');
    cy.contains('Length').should('be.visible');

    // 파일을 input[type="file"] 요소에 업로드
    cy.get('input[type="file"]').attachFile(filePath);

    // 파일이 제대로 업로드 되었는지 확인 (예: 업로드 후 확인 메시지나 파일 이름이 표시되는지)
    cy.contains('fasta_example_3.fasta').should('be.visible');

    cy.get('button.next-button').click();

    // sequence-chunk 안의 sequence 클래스가 1 + 업로드한 FASTA 파일 개수(1개)인지 확인
    cy.get('.sequence-chunk').eq(0)  // 첫 번째 청크 선택
      .find('.sequence')  // 그 안에서 sequence 클래스 요소 찾기
      .should('have.length', 2);  // 예상 개수와 비교

  });

  // 시나리오 ID: TS_004_1.2.	두 개이상의 파일 업로드시
  it('사용자로부터 FASTA 파일 두개 이상 업로드', () => {
    const filePath1 = 'fasta_example_3.fasta'; // fixtures 폴더에 있는 파일 경로
    const filePath2 = 'fasta_example_2.fasta'; // fixtures 폴더에 있는 파일 경로

    // NCBI 레퍼런스 시퀀스 ID 입력 필드에 NC_045512.2을 입력한다. ‘Done’ 버튼을 누른다.
    cy.get('input[id="referenceSequenceId"]').type('NC_045512.2');
    cy.get('button.done-button').click();


    // 파일을 input[type="file"] 요소에 업로드
    cy.get('input[type="file"]').attachFile(filePath1);
    cy.get('input[type="file"]').attachFile(filePath2);

    // 파일이 제대로 업로드 되었는지 확인 (예: 업로드 후 확인 메시지나 파일 이름이 표시되는지)
    cy.contains('fasta_example_3.fasta').should('be.visible');
    cy.contains('fasta_example_2.fasta').should('be.visible');

    cy.get('button.next-button').click();

    // sequence-chunk 안의 sequence 클래스가 1 + 업로드한 FASTA 파일 개수(2개)인지 확인
    cy.get('.sequence-chunk').eq(0)  // 첫 번째 청크 선택
      .find('.sequence')  // 그 안에서 sequence 클래스 요소 찾기
      .should('have.length', 3);  // 예상 개수와 비교


  });

  // 시나리오 ID: TS_004_1.3.	한 개의 시퀀스 입력 시
  it('한 개의 시퀀스 입력', () => {

    // NCBI 레퍼런스 시퀀스 ID 입력 필드에 NC_045512.2을 입력한다. ‘Done’ 버튼을 누른다.
    cy.get('input[id="referenceSequenceId"]').type('NC_045512.2');
    cy.get('button.done-button').click();


    // 한 개의 시퀀스 입력
    cy.get('textarea[placeholder="TAGCTAGCCGATCG....."]')
      .type('ATGCGTCTGACGCGGACTTGCGCGAGTGGTGACGTGACGGTAGCCGGAGTGCGCTGCGGA');



    cy.get('button.next-button').click();

    // sequence-chunk 안의 sequence 클래스가 1 + 업로드한 FASTA 파일 개수(1개)인지 확인
    cy.get('.sequence-chunk').eq(0)  // 첫 번째 청크 선택
      .find('.sequence')  // 그 안에서 sequence 클래스 요소 찾기
      .should('have.length', 2);  // 예상 개수와 비교


  });

  // 시나리오 ID: TS_004_1.4.	두 개이상의 시퀀스 입력 시
  it('두 개이상의 시퀀스 입력', () => {

    // NCBI 레퍼런스 시퀀스 ID 입력 필드에 NC_045512.2을 입력한다. ‘Done’ 버튼을 누른다.
    cy.get('input[id="referenceSequenceId"]').type('NC_045512.2');
    cy.get('button.done-button').click();



    // 두 개의 시퀀스 입력
    cy.contains('div.sequence-header', 'Sequence1')  // 'Sequence1' 텍스트가 포함된 div를 찾음
      .parent()  // 부모 요소로 이동
      .find('textarea')  // 부모 요소 아래의 textarea를 찾음
      .type('ATGCGTCTGACGCGGACTTGCGCGAGTGGTGACGTGACGGTAGCCGGAGTGCGCTGCGGA');



    cy.get('button.add-sequence-button').click();

    cy.contains('div.sequence-header', 'Sequence2')  // 'Sequence2' 텍스트가 포함된 div를 찾음
      .parent()  // 부모 요소로 이동
      .find('textarea')  // 부모 요소 아래의 textarea를 찾음
      .type('ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT');




    cy.get('button.next-button').click();

    // sequence-chunk 안의 sequence 클래스가 1 + 업로드한 FASTA 파일 개수(2개)인지 확인
    cy.get('.sequence-chunk').eq(0)  // 첫 번째 청크 선택
      .find('.sequence')  // 그 안에서 sequence 클래스 요소 찾기
      .should('have.length', 3);  // 예상 개수와 비교

  });


  // 시나리오 ID: TS_004_2. 잘못된 형식의 파일을 업로드하거나 잘못된 시퀀스를 시퀀스 입력란에 입력하는 경우
  it('잘못된 형식의 파일을 업로드', () => {
    const filePath1 = 'sample_csv_file.csv'; // fixtures 폴더에 있는 파일 경로
    const filePath2 = 'sample_text_file.txt'; // fixtures 폴더에 있는 파일 경로

    // NCBI 레퍼런스 시퀀스 ID 입력 필드에 NC_045512.2을 입력한다. ‘Done’ 버튼을 누른다.
    cy.get('input[id="referenceSequenceId"]').type('NC_045512.2');
    cy.get('button.done-button').click();


    // 파일을 input[type="file"] 요소에 업로드
    cy.get('input[type="file"]').attachFile(filePath1);
    cy.get('input[type="file"]').attachFile(filePath2);

    // 현재 이 부분 에러 안나고 넘어갑니다....... 테스트 상으론 그래요. (alert)Failed to fetch 같기도 한데.. 일단은 넘어는 감.

    cy.get('button.next-button').click();

    // 넘어갔는지 확인하는 코드.
    cy.get('.sequence-chunk').eq(0)  // 첫 번째 청크 선택
      .find('.sequence')  // 그 안에서 sequence 클래스 요소 찾기
      .should('have.length', 1);  // 예상 개수와 비교


  });

  // 시나리오 ID: TS_004_3. 아무것도 입력하지 않는 경우
  it('아무것도 입력하지 않는 경우', () => {

    // NCBI 레퍼런스 시퀀스 ID 입력 필드에 NC_045512.2을 입력한다. ‘Done’ 버튼을 누른다.
    cy.get('input[id="referenceSequenceId"]').type('NC_045512.2');
    cy.get('button.done-button').click();


    cy.get('button.next-button').click();

    // sequence-chunk 안의 sequence 클래스가 1개인지 확인
    cy.get('.sequence-chunk').eq(0)  // 첫 번째 청크 선택
      .find('.sequence')  // 그 안에서 sequence 클래스 요소 찾기
      .should('have.length', 1);  // 예상 개수와 비교

  });


});
