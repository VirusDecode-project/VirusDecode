// 서열 분석 페이지 및 서열 정렬 탭
describe("1. 서열 정렬(alignment)", () => {
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

  // 1-1. 서열 정렬 알고리즘실행하여 정렬된 데이터 시각화
  it("1-1. 서열 정렬 알고리즘실행하여 정렬된 데이터 시각화", () => {
    cy.wait("@alignmentRequest", { timeout: 20000 }).then((interception) => {
      expect(interception.response.statusCode).to.eq(200);

      cy.get(".sequence-boxes").find(".sequence-box").should("exist");
    });
  });
});

describe("2. 불일치 구간(돌연변이) 색상 표시", () => {
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

  // 2-1. 불일치 아미노산 탐지 및 색상 변화 적용
  it("2-1. 불일치 아미노산 탐지 및 색상 변화 적용", () => {
    cy.get(".sequence-boxes").find(".sequence-box.gap.different").should("exist");
  });
});

describe("3. 유전체 바", () => {
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

  it("3-1. 유전체를 선택하여 해당하는 서열 표시", () => {
    // 변경 전 기본 옵션 값 확인
    cy.get('#region-select').should('have.value', 'ORF1ab');

    // 유전체바에서 새로운 옵션을 선택 (예: 'S')
    cy.get('.stacked-bar-label').contains('S').click();

    // 선택된 값이 'S'로 변경되었는지 확인
    cy.get('#region-select').should('have.value', 'S');

    // 유전체바에서 새로운 옵션을 선택 (예: 'ORF1ab')
    cy.get('.stacked-bar-label').contains('ORF1ab').click();
    cy.get('#region-select').should('have.value', 'ORF1ab');
    });
  });

describe("4. 도움말 아이콘", () => {
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

  // Help 버튼 클릭 시 모달 창이 뜨는지 확인
  it("4-1. 사용 설명서 안내", () => {
    // 도움말 버튼 클릭
    cy.get(".help-icon").click();

    // 모달 창이 열리는지 확인
    cy.get(".help-modal").should("be.visible");

    // 모달 창 제목이 'mRNA Conversion Instructions'인지 확인
    cy.get(".help-modal")
      .contains("mRNA Conversion Instructions")
      .should("be.visible");
  });
});

describe("5. mRNA 변환", () => {
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
  // Sequence Box 클릭 후 모달 창이 뜨는지 확인
  it("5-1. 서열을 클릭하여 서열 변환 창 실행", () => {
    // Sequence Box 갯수 가져와서 랜덤하게 하나 클릭
    cy.get('.sequence-boxes').its('length').then((len) => {
      const randomIndex = Math.floor(Math.random() * len);
      cy.get('.sequence-boxes').eq(randomIndex).click();
    });

    // 모달 창이 열리는지 확인
    cy.get(".modal-content").should("be.visible");
});
});


describe("6. 히스토리 저장", () => {
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
        cy.get('input[id="referenceSequenceId"]').type(referenceId.SARS_CoV_2_ID);
        cy.get("button.done-button").click();
        cy.wait("@metadataRequest").then((interception) => {
          expect(interception.response.statusCode).to.eq(200);
  
          // "Sequence ID", "Name", "Description", "Length" 확인
          cy.contains("Sequence ID").should("be.visible");
          cy.contains("Name").should("be.visible");
          cy.contains("Description").should("be.visible");
          cy.contains("Length").should("be.visible");
        });
        cy.contains("div.sequence-header", "Sequence1") // 'Sequence1' 텍스트가 포함된 div 찾음
          .parent() // 부모 요소로 이동
          .find("textarea") // 부모 요소 아래의 textarea 찾음
          .type(
            "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCC"
          );
        cy.get("button.next-button").click();
      });
    });
  
    // 히스토리 아이템 클릭 후 저장된 데이터 확인
    it("6-1. 생성된 정렬 데이터 확인", () => {
      // 히스토리 아이템 첫 번째 요소 클릭
      cy.get(".history-list .history-item").first().click();
  
      // 정렬된 데이터가 있는지 확인 (예: sequence-box 확인)
      cy.wait("@alignmentRequest", { timeout: 20000 }).then((interception) => {
        expect(interception.response.statusCode).to.eq(200);
  
        // 시퀀스 박스가 있는지 확인
        cy.get(".sequence-boxes").find(".sequence-box").should("exist");
      });
  
      // 불일치 구간이 있는지 확인
      cy.get(".sequence-boxes").find(".sequence-box.gap.different").should("exist");
    });
  });
  
