describe("1. 아미노산 구간 설정 및 검증", () => {
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
    cy.intercept("POST", "/api/analysis/linearDesign").as("linearDesignRequest");
    cy.inputSeqSetup();
  });

  it("1-1. 아미노산 구간 설정", () => {
    let startValue = Math.floor(Math.random() * (1276 - 1 + 1)) + 1;
    let endValue = Math.min(startValue + Math.floor(Math.random() * 100), 1276);
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

    cy.wait('@linearDesignRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
    });
  });

  it("1-2. 범위 유효성 검증", () => {
    let startValue = Math.floor(Math.random() * (1276 - 1 + 1)) + 1277;
    let endValue = Math.floor(Math.random() * (1276 - 1 + 1)) + 1277;

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

    cy.contains("index must be between").should("be.visible");
  });
});

describe("2. mRNA 서열 변환", () => {
  it("2-1. Convert 버튼을 통해 데이터 전송 및 mRNA 디자인 탭으로 이동", () => {
    cy.signupAndLoginIfDuplicate(
      "testFName",
      "testLName",
      "testId",
      "testPw",
      "testPw"
    );
    cy.intercept("POST", "/api/inputSeq/metadata").as("metadataRequest");
    cy.intercept("POST", "/api/inputSeq/alignment").as("alignmentRequest");
    cy.intercept("POST", "/api/analysis/pdb").as("PDBrequest");
    cy.intercept("POST", "/api/analysis/linearDesign").as("linearDesignRequest");
    cy.LinearDesignConvert();
    cy.wait('@linearDesignRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get(".nav-tabs")
        .contains("mRNA design")
        .should("not.have.class", "disabled");
      
      cy.get(".nav-tabs").contains("mRNA design").click(); // mRNA 탭 클릭
      cy.get(".nav-tabs").contains("mRNA design").should("have.class", "active"); // mRNA 탭이 활성화되었는지 확인
      // mRNA 디자인 페이지의 고유 요소 확인 (예: mRNA Visualization)
      cy.contains("mRNA Visualization").should("be.visible"); // mRNA 시각화 요소가 표시되는지 확인
    });
  });
});
