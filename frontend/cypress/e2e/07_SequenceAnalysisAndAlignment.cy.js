
// 서열 분석 페이지 및 서열 정렬 탭
describe("1. 서열 정렬(alignment)", () => {
  it("1-1. 서열 정렬 알고리즘 실행하여 정렬된 데이터 시각화", () => {
    cy.signupAndLoginIfDuplicate(
      "testFName",
      "testLName",
      "testId",
      "testPw",
      "testPw"
    );
    cy.intercept("POST", "/api/inputSeq/metadata").as("metadataRequest");
    cy.intercept("POST", "/api/inputSeq/alignment").as("alignmentRequest");
    cy.inputSeqSetupFile();

    cy.wait("@alignmentRequest").then((interception) => {
      expect(interception.response.statusCode).to.eq(200);

      cy.get('.stacked-bar-label', ).contains('S').click();
      cy.fixture('environment').then((environment) => {
        cy.get('.sequence-label').eq(0).should('contain', environment.SARS_CoV_2_ID);
        cy.get('.sequence-label').eq(1).should('contain', environment.VARIANT_1);
        cy.get('.sequence-label').eq(2).should('contain', environment.VARIANT_2);
      });

      let sequence1 = 'MFVFLVLLPLVSSQCVN';
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < sequence1.length; j++) {
          cy.get('.sequence-chunk').eq(0)
            .find('.sequence-boxes').eq(i)
            .find('.sequence-line')
            .find('.sequence-box').eq(j)
            .should('contain', sequence1[j]);
        }
      }
      let sequence2 = 'TQDLFLPFFSNVTWFH';
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < sequence2.length; j++) {
          cy.get('.sequence-chunk').eq(1)
            .find('.sequence-boxes').eq(i)
            .find('.sequence-line')
            .find('.sequence-box').eq(j)
            .should('contain', sequence2[j]);
        }
      }
    });
  });
});

describe("2. 불일치 구간(돌연변이) 색상 표시", () => {
  it("2-1. 불일치 아미노산 탐지 및 색상 변화 적용", () => {
    cy.signupAndLoginIfDuplicate(
      "testFName",
      "testLName",
      "testId",
      "testPw",
      "testPw"
    );
    cy.intercept("POST", "/api/inputSeq/metadata").as("metadataRequest");
    cy.intercept("POST", "/api/inputSeq/alignment").as("alignmentRequest");
    cy.inputSeqSetupFile();
    cy.wait("@alignmentRequest").then((interception) => {
      expect(interception.response.statusCode).to.eq(200);

      cy.get('.stacked-bar-label').contains('S').click();
      cy.fixture('environment').then((environment) => {
        cy.get('.sequence-label').eq(0).should('contain', environment.SARS_CoV_2_ID);
        cy.get('.sequence-label').eq(1).should('contain', environment.VARIANT_1);
        cy.get('.sequence-label').eq(2).should('contain', environment.VARIANT_2);
      });

      // 불일치 구간 1
      let sequence1 = ['LTT', 'FTN', 'LTT']
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < sequence1[i].length; j++) {
          cy.get('.sequence-chunk').eq(0)
            .find('.sequence-boxes').eq(i)
            .find('.sequence-line')
            .find('.sequence-box').eq(j + 17)
            .should('contain', sequence1[i][j]);
        }
      }

      // 불일치 구간 2
      let sequence2 = ['AIH', 'AIH', 'VI-']
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < sequence2[i].length; j++) {
          cy.get('.sequence-chunk').eq(1)
            .find('.sequence-boxes').eq(i)
            .find('.sequence-line')
            .find('.sequence-box').eq(j + 16)
            .should('contain', sequence2[i][j]);
        }
      }

      // 불일치 구간 3
      let sequence3 = ['---', '---', 'EPE']
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < sequence3[i].length; j++) {
          cy.get('.sequence-chunk').eq(4)
            .find('.sequence-boxes').eq(i)
            .find('.sequence-line')
            .find('.sequence-box').eq(j + 14)
            .should('contain', sequence3[i][j]);
        }
      }


    });
  });
});

describe("3. 유전체 바", () => {
  it("3-1. 유전체를 선택하여 해당하는 서열 표시", () => {
    cy.signupAndLoginIfDuplicate(
      "testFName",
      "testLName",
      "testId",
      "testPw",
      "testPw"
    );
    cy.intercept("POST", "/api/inputSeq/metadata").as("metadataRequest");
    cy.intercept("POST", "/api/inputSeq/alignment").as("alignmentRequest");
    cy.inputSeqSetup();
    cy.wait("@alignmentRequest").then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
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
});

describe("4. 도움말 아이콘", () => {
  // Help 버튼 클릭 시 모달 창이 뜨는지 확인
  it("4-1. 사용 설명서 안내", () => {
    cy.signupAndLoginIfDuplicate(
      "testFName",
      "testLName",
      "testId",
      "testPw",
      "testPw"
    );
    cy.intercept("POST", "/api/inputSeq/metadata").as("metadataRequest");
    cy.intercept("POST", "/api/inputSeq/alignment").as("alignmentRequest");
    cy.inputSeqSetup();
    cy.wait("@alignmentRequest").then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
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
});

describe("5. mRNA 변환", () => {
  // Sequence Box 클릭 후 모달 창이 뜨는지 확인
  it("5-1. 서열을 클릭하여 서열 변환 창 실행", () => {
    cy.signupAndLoginIfDuplicate(
      "testFName",
      "testLName",
      "testId",
      "testPw",
      "testPw"
    );
    cy.intercept("POST", "/api/inputSeq/metadata").as("metadataRequest");
    cy.intercept("POST", "/api/inputSeq/alignment").as("alignmentRequest");
    cy.inputSeqSetup();
    cy.wait("@alignmentRequest").then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      // Sequence Box 갯수 가져와서 랜덤하게 하나 클릭
      cy.get('.sequence-boxes').its('length').then((len) => {
        const randomIndex = Math.floor(Math.random() * len);
        cy.get('.sequence-boxes').eq(randomIndex).click();
        // 모달 창이 열리는지 확인
        cy.get(".modal-content").should("be.visible");
      });
    });
  });
});


describe("6. 히스토리 저장", () => {
  // 히스토리 아이템 클릭 후 저장된 데이터 확인
  it("6-1. 생성된 정렬 데이터 확인", () => {
    cy.signupAndLoginIfDuplicate(
      "testFName",
      "testLName",
      "testId",
      "testPw",
      "testPw"
    );
    cy.intercept("POST", "/api/inputSeq/metadata").as("metadataRequest");
    cy.intercept("POST", "/api/inputSeq/alignment").as("alignmentRequest");
    cy.inputSeqSetup();

    // 히스토리 아이템 첫 번째 요소 클릭
    cy.get(".history-list .history-item").first().click();

    // 정렬된 데이터가 있는지 확인 (예: sequence-box 확인)
    cy.wait("@alignmentRequest").then((interception) => {
      expect(interception.response.statusCode).to.eq(200);

      // 시퀀스 박스가 있는지 확인
      cy.get(".sequence-boxes").find(".sequence-box").should("exist");
      // 불일치 구간이 있는지 확인
      cy.get(".sequence-boxes").find(".sequence-box.gap.different").should("exist");
    });
  });
});

