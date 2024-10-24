describe("1. mRNA 시각화 및 정보 제공", () => {
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

  it("1-1. mRNA 2D 구조 시각화", () => {
    cy.contains("mRNA Visualization").should("be.visible");
    cy.get("#plotting-area").should("be.visible");
  });

  it("1-2. RNA 및 단백질 파라미터 제공", () => {
    cy.contains("Protein Parameters").should("be.visible");
  });
});

describe("2. 히스토리 저장", () => {
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

  it("2-1. 생성된 mRNA 데이터 기록 저장", () => {
    cy.get(".history-list .history-item").first().click();

    cy.contains("mRNA Visualization").should("be.visible");

    cy.contains("Protein Parameters").should("be.visible");
  });
});
