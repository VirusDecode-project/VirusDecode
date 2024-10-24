describe("1. PDB 관련 정보 제공", () => {
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
    cy.intercept("POST", "/api/analysis/pdb").as("PDBrequest");

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

      let startValue = Math.floor(Math.random() * (1276 - 1 + 1)) + 1;
      let endValue = Math.min(
        startValue + Math.floor(Math.random() * 100),
        1276
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

  it("1-1. RCSB로부터 PDB ID, 단백질 설명 가져오기", () => {
    cy.wait(20000);
    cy.get(".nav-tabs").contains("3D viewer").click();

    // PDB ID 입력 및 정보 요청

    // PDB 정보가 정상적으로 로드되는지 확인
    cy.wait("@PDBrequest").then((interception) => {
      expect(interception.response.statusCode).to.eq(200); // 요청 성공 여부 확인

      // PDB 관련 데이터 확인
      cy.contains("PDB ID").should("be.visible");
      // 첫 번째 list-row 클릭
      cy.get(".list-row").eq(1).click();

      // custom-tooltip이 나타나는지 확인
      cy.get(".custom-tooltip").should("be.visible");
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
      cy.intercept("POST", "/api/analysis/pdb").as("PDBrequest");

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

        let startValue = Math.floor(Math.random() * (1276 - 1 + 1)) + 1;
        let endValue = Math.min(
          startValue + Math.floor(Math.random() * 100),
          1276
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

    it("2-1. 선택한 PDB ID를 RCSB로부터 PDB 데이터 가져오기", () => {
      cy.wait(20000);
      cy.get(".nav-tabs").contains("3D viewer").click();

      // PDB ID 입력 및 정보 요청

      // PDB 정보가 정상적으로 로드되는지 확인
      cy.wait("@PDBrequest").then((interception) => {
        expect(interception.response.statusCode).to.eq(200); // 요청 성공 여부 확인

        // PDB 관련 데이터 확인
        cy.contains("PDB ID").should("be.visible");
        // 첫 번째 list-row 클릭
        cy.get(".list-row").eq(2).click();

        // custom-tooltip이 나타나는지 확인
        cy.get(".custom-tooltip").should("be.visible");
      });
    });

    it("2-2. PDB 단백질 3D 구조 시각화", () => {
        cy.wait(20000);
        cy.get(".nav-tabs").contains("3D viewer").click();

        // PDB ID 입력 및 정보 요청

        // PDB 정보가 정상적으로 로드되는지 확인
        cy.wait("@PDBrequest").then((interception) => {
          expect(interception.response.statusCode).to.eq(200); // 요청 성공 여부 확인

          // PDB 관련 데이터 확인
          cy.contains("PDB ID").should("be.visible");
          // 첫 번째 list-row 클릭
          cy.get(".list-row").eq(3).click();

          cy.get('#viewport canvas').should('be.visible');
        });
      });

  });

describe("3. 히스토리 저장", () => {
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
    cy.intercept("POST", "/api/analysis/pdb").as("PDBrequest");

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

      let startValue = Math.floor(Math.random() * (1276 - 1 + 1)) + 1;
      let endValue = Math.min(
        startValue + Math.floor(Math.random() * 100),
        1276
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

  it("3-1. 생성된 PDB 데이터 저장", () => {
    cy.wait(20000);
    cy.get(".nav-tabs").contains("3D viewer").click();

    cy.get(".history-list .history-item").first().click();

    // PDB 정보가 정상적으로 로드되는지 확인
    cy.wait("@PDBrequest").then((interception) => {
      expect(interception.response.statusCode).to.eq(200); // 요청 성공 여부 확인

      // PDB 관련 데이터 확인
      cy.contains("PDB ID").should("be.visible");
      // 첫 번째 list-row 클릭
      cy.get(".list-row").eq(2).click();

      // custom-tooltip이 나타나는지 확인
      cy.get(".custom-tooltip").should("be.visible");
    });
  });
});
