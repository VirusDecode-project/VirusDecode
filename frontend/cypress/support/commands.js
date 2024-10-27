import { mockLogin, mockSignup } from "../mocks/mock";

Cypress.Commands.add('mockLogin', mockLogin);
Cypress.Commands.add('mockSignup', mockSignup);

Cypress.Commands.add('signupAndLoginIfDuplicate', (firstName, lastName, id, password, cPassword) => {
  cy.visit('http://localhost:3000/signup');
  cy.get('input[name="firstName"]').type(firstName);
  cy.get('input[name="lastName"]').type(lastName);
  cy.get('input[name="id"]').type(id);
  cy.get('input[name="password"]').type(password);
  cy.get('input[name="cPassword"]').type(cPassword);

  cy.intercept('POST', '/api/auth/signup').as('signupRequest');
  cy.get('.SignupBtn').click(); 
  cy.wait('@signupRequest').then((interception) => {
    const { response } = interception; 
    if (response.statusCode === 200) {
      cy.get('.message-modal-content')
        .should('be.visible')
        .and('contain', '회원가입이 완료되었습니다.');
      cy.get('.message-modal-content').contains('Close').click();
      cy.url().should('include', '/login');
      cy.login(id, password); 

    } else if (response.statusCode === 400 && response.body.includes("이미 존재하는 ID 입니다.")) {
      cy.get('.message-modal-content')
        .should('be.visible')
        .and('contain', '이미 존재하는 ID 입니다.');
      cy.get('.message-modal-content').contains('Close').click();
      cy.get('.gotoLoginBtn').click(); 
      cy.login(id, password); 

    } else {
      console.error(`Error response: ${JSON.stringify(response.body)}`);
      throw new Error(`Unexpected error: ${response.statusText}`);
    }
  });
});

Cypress.Commands.add('login', (loginId, password) => {
  cy.visit('http://localhost:3000/login');
  cy.get('input[name="loginId"]').type(loginId);
  cy.get('input[name="password"]').type(password);
  cy.get('.loginBtn').click();
});

Cypress.Commands.add('guestlogin', () => {
  cy.visit('http://localhost:3000/');
  cy.get('.decode-button').click();
  cy.get('.stayLoggedOutBtn').click(); 
});

Cypress.Commands.add('inputSeqSetup', () => {
  // NCBI로부터 ID 유효성 검사 및 메타데이터 가져오기
  cy.fixture("referenceId").then((referenceId) => {
    // 올바른 NCBI 레퍼런스 시퀀스 ID 입력 후 'Done' 버튼 클릭
    cy.get('input[id="referenceSequenceId"]').type(referenceId.SARS_CoV_2_ID);
    cy.get("button.done-button").click();
    cy.wait("@metadataRequest").then((interception) => {
      expect(interception.response.statusCode).to.eq(200);

      // "Sequence ID", "Name", "Description", "Length"가 있는지 확인
      cy.contains("Sequence ID").should("be.visible");
      cy.contains("Name").should("be.visible");
      cy.contains("Description").should("be.visible");
      cy.contains("Length").should("be.visible");
    });
    cy.contains("div.sequence-header", "Sequence1") // 'Sequence1' 텍스트가 포함된 div를 찾음
      .parent() // 부모 요소로 이동
      .find("textarea") // 부모 요소 아래의 textarea를 찾음
      .type(
        "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCC"
      );
    cy.get("button.next-button").click();
  });
});

Cypress.Commands.add('LinearDesignConvert', () => {
  cy.inputSeqSetup();
  let startValue = Math.floor(Math.random() * (1276 - 1 + 1)) + 1;
  let endValue = Math.min(
    startValue + Math.floor(Math.random() * 100),
    1276
  );
  cy.wait("@alignmentRequest").then((interception) => {
    expect(interception.response.statusCode).to.eq(200);
    cy.get(".sequence-boxes").eq(0).click();

    cy.contains(".modal-content label", "Select Coding Sequence:")
      .find("select")
      .select("S");

    cy.contains(".modal-content label", "Start Amino Acid Position:")
      .find('input[type="number"]')
      .type(startValue.toString());

    cy.contains(".modal-content label", "End Amino Acid Position:")
      .find('input[type="number"]')
      .type(endValue.toString());

    cy.get(".modal-next-button").click();
  });
});