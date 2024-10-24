describe('1. 아미노산 구간 설정 및 검증', () => {

    beforeEach(() => {
        cy.signupAndLoginIfDuplicate(
          "testFName",
          "testLName",
          "testId",
          "testPw",
          "testPw"
        );
        cy.intercept("POST", "/api/inputSeq/metadata").as("metadataRequest");
        cy.intercept("POST", "/api/inputSeq/alignment").as("alignmentRequest");

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

    it('1-1. 아미노산 구간 설정', () => {
        let startValue = Math.floor(Math.random() * (1273 - 1 + 1)) + 1;
        let endValue = Math.min(startValue + Math.floor(Math.random() * 500), 1273);

        cy.intercept('POST', '/api/analysis/linearDesign').as('linearDesignRequest');

        cy.get('.sequence-boxes').eq(0).click();

        cy.contains('.modal-content label', 'Select Coding Sequence:')
          .find('select')
          .select('S');

        cy.contains('.modal-content label', 'Start Amino Acid Position:')
          .find('input[type="number"]')
          .type(startValue.toString());

        cy.contains('.modal-content label', 'End Amino Acid Position:')
          .find('input[type="number"]')
          .type(endValue.toString());

        cy.get('.modal-next-button').click();

        cy.wait('@linearDesignRequest', { timeout: 90000 }).its('response.statusCode').should('eq', 200);
      });

    it('1-2. 범위 유효성 검증', () => {
        let startValue = Math.floor(Math.random() * (1273 - 1 + 1)) + 1274;
        let endValue = Math.floor(Math.random() * (1273 - 1 + 1)) + 1274;

        cy.get('.sequence-boxes').eq(0).click();

        cy.contains('.modal-content label', 'Select Coding Sequence:')
          .find('select')
          .select('S');

        cy.contains('.modal-content label', 'Start Amino Acid Position:')
          .find('input[type="number"]')
          .type(startValue.toString());

        cy.contains('.modal-content label', 'End Amino Acid Position:')
          .find('input[type="number"]')
          .type(endValue.toString());

        cy.get('.modal-next-button').click();

        cy.contains('index must be between').should('be.visible');
      });

});

describe("2. mRNA 서열 변환", () => {
  beforeEach(() => {
    // 기본 URL로 애플리케이션에 접속
    cy.signupAndLoginIfDuplicate(
      "testFName",
      "testLName",
      "testId",
      "testPw",
      "testPw"
    );
    cy.intercept("POST", "/api/inputSeq/metadata").as("metadataRequest");
    cy.intercept("POST", "/api/inputSeq/alignment").as("alignmentRequest");

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

      let startValue = Math.floor(Math.random() * (1273 - 1 + 1)) + 1;
      let endValue = Math.min(
        startValue + Math.floor(Math.random() * 100),
        1273
      );

      cy.intercept("POST", "/api/analysis/linearDesign").as(
        "linearDesignRequest"
      );

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

  it("2-1. Convert 버튼을 통해 데이터 전송 및 mRNA 디자인 탭으로 이동", () => {
    cy.wait(20000); 
    cy.get(".nav-tabs")
      .contains("mRNA design")
      .should("not.have.class", "disabled");
    cy.get(".nav-tabs").contains("mRNA design").click(); // mRNA 탭 클릭
    cy.get(".nav-tabs").contains("mRNA design").should("have.class", "active"); // mRNA 탭이 활성화되었는지 확인
    // mRNA 디자인 페이지의 고유 요소 확인 (예: mRNA Visualization)
    cy.contains("mRNA Visualization").should("be.visible"); // mRNA 시각화 요소가 표시되는지 확인
  });
});
