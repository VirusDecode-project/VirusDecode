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
      cy.get('.mrna-subtitle').eq(0).should('contain', 'mRNA Visualization');
      // 확대
      cy.get("#plotting-area")
        .trigger("wheel", { deltaY: -200 });
      // 축소
      cy.get("#plotting-area")
        .trigger("wheel", { deltaY: 200 });
    });
  });

  it("1-2. RNA 및 단백질 파라미터 제공", () => {
    cy.wait('@linearDesignRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('.mrna-title').eq(0).should('contain', 'mRNA Parameters');
      cy.get('.mrna-title').eq(1).should('contain', 'Protein Parameters');
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
    cy.intercept("POST", "/api/history/get").as("historyRequest");
    cy.LinearDesignConvert();
    cy.wait('@linearDesignRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      let displayedAminoAcidSeq;
      cy.get(".mrna-sequence").eq(0)
        .invoke("text")
        .then((text) => {
          displayedAminoAcidSeq = text.trim();
          cy.log("Displayed Amino: ", displayedAminoAcidSeq);
          cy.get(".history-list .history-item").eq(1).click();
          cy.wait("@historyRequest").then((interception) => {
            expect(interception.response.statusCode).to.eq(200);
            cy.get(".history-list .history-item").eq(0).click();
            cy.wait("@historyRequest").then((interception) => {
              expect(interception.response.statusCode).to.eq(200);
              cy.get(".mrna-sequence").eq(0).should('contain', displayedAminoAcidSeq);
            });
          });
        });
    });
  });
});
