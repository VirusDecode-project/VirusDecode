describe("1. PDB 관련 정보 제공", () => {
  it("1-1. RCSB로부터 PDB ID, 단백질 설명 가져오기", () => {
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
    cy.intercept("POST", "/api/analysis/pdb").as("PDBrequest");
    cy.LinearDesignConvert();
    cy.wait('@linearDesignRequest', {timeout: 20000}).then((interception) => {
      expect(interception.response.statusCode).to.eq(200);

      // PDB ID 입력 및 정보 요청

      // PDB 정보가 정상적으로 로드되는지 확인
      cy.wait("@PDBrequest", {timeout: 20000}).then((interception) => {
        expect(interception.response.statusCode).to.eq(200); // 요청 성공 여부 확인
        cy.get(".nav-tabs").contains("3D viewer").click();
        // PDB 관련 데이터 확인
        cy.contains("PDB ID").should("be.visible");
        // 첫 번째 list-row 클릭
        cy.get(".list-row").eq(1).click();

        // custom-tooltip이 나타나는지 확인
        cy.get(".custom-tooltip").should("be.visible");
      });
    });
  });
});

describe("2. 3D 이미지 생성 및 시각화", () => {
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
      cy.intercept("POST", "/api/analysis/pdb").as("PDBrequest");
      cy.LinearDesignConvert();
    });

    it("2-1. 선택한 PDB ID를 RCSB로부터 PDB 데이터 가져오기", () => {
      cy.wait('@linearDesignRequest', {timeout: 20000}).then((interception) => {
        expect(interception.response.statusCode).to.eq(200);

        // PDB ID 입력 및 정보 요청

        // PDB 정보가 정상적으로 로드되는지 확인
        cy.wait("@PDBrequest", {timeout: 20000}).then((interception) => {
          expect(interception.response.statusCode).to.eq(200); // 요청 성공 여부 확인
          cy.get(".nav-tabs").contains("3D viewer").click();
          // PDB 관련 데이터 확인
          cy.contains("PDB ID").should("be.visible");
          // 첫 번째 list-row 클릭
          cy.get(".list-row").eq(2).click();

          // custom-tooltip이 나타나는지 확인
          cy.get(".custom-tooltip").should("be.visible");
        });
      });
    });

    it("2-2. PDB 단백질 3D 구조 시각화", () => {
      cy.wait('@linearDesignRequest', {timeout: 20000}).then((interception) => {
        expect(interception.response.statusCode).to.eq(200);
      
        // PDB ID 입력 및 정보 요청

        // PDB 정보가 정상적으로 로드되는지 확인
        cy.wait("@PDBrequest", {timeout: 20000}).then((interception) => {
          expect(interception.response.statusCode).to.eq(200); // 요청 성공 여부 확인
          cy.get(".nav-tabs").contains("3D viewer").click();
          // PDB 관련 데이터 확인
          cy.contains("PDB ID").should("be.visible");
          // 첫 번째 list-row 클릭
          cy.get(".list-row").eq(3).click();

          cy.get('#viewport canvas').should('be.visible');
        });
      });
    });
  });

describe("3. 히스토리 저장", () => {
  it("3-1. 생성된 PDB 데이터 저장", () => {
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
    cy.wait('@linearDesignRequest', {timeout: 20000}).then((interception) => {
      expect(interception.response.statusCode).to.eq(200);

      // PDB 정보가 정상적으로 로드되는지 확인
      cy.wait("@PDBrequest", {timeout: 20000}).then((interception) => {
        expect(interception.response.statusCode).to.eq(200); // 요청 성공 여부 확인
        cy.get(".nav-tabs").contains("3D viewer").click();

        cy.get(".history-list .history-item").first().click();
        // PDB 관련 데이터 확인
        cy.contains("PDB ID").should("be.visible");
        // 첫 번째 list-row 클릭
        cy.get(".list-row").eq(2).click();

        // custom-tooltip이 나타나는지 확인
        cy.get(".custom-tooltip").should("be.visible");
      });
    });
  });
});
