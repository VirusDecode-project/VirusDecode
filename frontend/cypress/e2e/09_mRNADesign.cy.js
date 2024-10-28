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
    cy.intercept("POST", "/api/analysis/linearDesign").as("linearDesignRequest");
    cy.LinearDesignConvert();
  });

  it("1-1. mRNA 2D 구조 시각화", () => {
    cy.wait('@linearDesignRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.contains("mRNA Visualization").should("be.visible");
      cy.get("#plotting-area").should("be.visible");
    });
  });

  it("1-2. RNA 및 단백질 파라미터 제공", () => {
    cy.wait('@linearDesignRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.contains("Protein Parameters").should("be.visible");
    });
  });
});

describe("2. 히스토리 저장", () => {
  it("2-1. 생성된 mRNA 데이터 기록 저장", () => {
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
    cy.LinearDesignConvert();
    cy.wait('@linearDesignRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get(".history-list .history-item").first().click();
      cy.contains("mRNA Visualization").should("be.visible");
      cy.contains("Protein Parameters").should("be.visible");
    });
  });
});
