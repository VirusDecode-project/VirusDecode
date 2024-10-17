// 서열 입력 페이지

describe('1. 레퍼런스 서열 입력 테스트', () => {
    beforeEach(() => {
        // 기본 URL로 애플리케이션에 접속
        cy.visit('http://localhost:3000');
        cy.contains('Try Decoding').click();
        cy.contains('stay logged out').click();
        cy.intercept('POST', '/api/inputSeq/metadata').as('metadataRequest');
        cy.intercept('POST', '/api/inputSeq/alignment').as('alignmentRequest');

    });

    // 1-1-1. Nucleotide ID 입력 실패(잘못된 ID 입력)
    it('잘못된 ID 입력 시: "NCBI에 요청한 nucleotide ID가 존재하지 않습니다."라는 메시지가 나타남.', () => {
        // 동일한 단계에서 잘못된 ID인 wrongReferenceId을 입력한다. ‘Done’ 버튼을 누른다.
        cy.get('input[id="referenceSequenceId"]').type('wrongReferenceId');
        cy.get('button.done-button').click();
        cy.wait('@metadataRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(500);

            // 잘못된 ID 입력 시: "NCBI에 요청한 nucleotide ID가 존재하지 않습니다."라는 메시지가 나타남.
            cy.contains('NCBI에 요청한 nucleotide ID가 존재하지 않습니다.').should('be.visible');
        })
    });

    // 1-1-2. Nucleotide ID 입력 실패(Null 입력)
    it('입력하지 않을 경우: “Please enter a valid sequence ID.”라는 메시지가 나타남.', () => {
        // 동일한 단계에서 입력하지 않고. ‘Done’ 버튼을 누른다.
        cy.get('button.done-button').click();

        // 잘못된 ID 입력 시: "NCBI에 요청한 nucleotide ID가 존재하지 않습니다."라는 메시지가 나타남.
        cy.contains('Please enter a valid sequence ID.').should('be.visible');

    });

    // 1-2. Nucleotide ID 입력 성공 및 NCBI로부터 ID 유효성 검사 및 메타데이터 가져오기
    it('올바른 ID 입력 시: 시스템에 등록되고 "Sequence ID, Name, Description, Length”에 대한 정보가 나타남', () => {
        // NCBI 레퍼런스 시퀀스 ID 입력 필드에 NC_045512.2을 입력한다. ‘Done’ 버튼을 누른다.
        cy.fixture('referenceId').then((referenceId) => {
            //   cy.get('input[id="referenceSequenceId"]').type('NC_045512.2');
            cy.get('input[id="referenceSequenceId"]').type(referenceId.SARS_CoV_2_ID);
            cy.get('button.done-button').click();
            cy.wait('@metadataRequest').then((interception) => {
                expect(interception.response.statusCode).to.eq(200);

                // "Sequence ID", "Name", "Description", "Length"가 있는지 확인
                cy.contains('Sequence ID').should('be.visible');
                cy.contains('Name').should('be.visible');
                cy.contains('Description').should('be.visible');
                cy.contains('Length').should('be.visible');
            });
        })
    });
});


describe('2. 분석 서열 입력', () => {
    beforeEach(() => {
        // 기본 URL로 애플리케이션에 접속
        cy.visit('http://localhost:3000');
        cy.contains('Try Decoding').click();
        cy.contains('stay logged out').click();
        cy.intercept('POST', '/api/inputSeq/metadata').as('metadataRequest');
        cy.intercept('POST', '/api/inputSeq/alignment').as('alignmentRequest');
    });

    it('Fasta Upload Test를 위한 reference id 입력', () => {
        cy.fixture('referenceId').then((referenceId) => {
            cy.get('input[id="referenceSequenceId"]').type(referenceId.SARS_CoV_2_ID);
            cy.get('button.done-button').click();
            cy.wait('@metadataRequest').then((interception) => {
                expect(interception.response.statusCode).to.eq(200);
                cy.contains('Sequence ID').should('be.visible');
                cy.contains('Name').should('be.visible');
                cy.contains('Description').should('be.visible');
                cy.contains('Length').should('be.visible');
            })
        })
    });

    it('잘못된 형식의 파일을 업로드', () => {
        const filePath = 'referenceId.json';
        cy.get('input[type="file"]').attachFile(filePath);

        cy.get('.message-modal-content')
            .should('be.visible')
            .and('contain', 'FASTA 파일 형식만 지원됩니다.');
        cy.get('.message-modal-content').contains('Close').click();
    });

    // 잘못된 염기서열 입력
    it('잘못된 형식의 시퀀스 입력', () => {
        cy.fixture('referenceId').then((referenceId) => {
            cy.get('input[id="referenceSequenceId"]').type(referenceId.SARS_CoV_2_ID);
            cy.get('button.done-button').click();
            cy.wait('@metadataRequest').then((interception) => {
                expect(interception.response.statusCode).to.eq(200);
                cy.contains('Sequence ID').should('be.visible');
                cy.contains('Name').should('be.visible');
                cy.contains('Description').should('be.visible');
                cy.contains('Length').should('be.visible');
            })
        })


        // 한 개의 시퀀스 입력
        cy.get('textarea[placeholder="TAGCTAGCCGATCG....."]')
            .type('안녕하세요');

        cy.get('button.next-button').click();
        cy.wait('@alignmentRequest', { timeout: 20000 }).then((interception) => {
            expect(interception.response.statusCode).to.eq(500);

            cy.get('.message-modal-content')
                .should('be.visible')
                .and('contain', '입력하신 서열 정보가 올바르지 않습니다. A, T, C, 그리고 G만 허용됩니다.');
            cy.get('.message-modal-content').contains('Close').click();

        })
    });

    // 시나리오 ID: TS_004_1.1. 한 개의 파일 업로드시
    it('사용자로부터 FASTA 파일 한개 업로드', () => {
        cy.fixture('referenceId').then((referenceId) => {

            cy.get('input[id="referenceSequenceId"]').type(referenceId.SARS_CoV_2_ID);
            cy.get('button.done-button').click();
            cy.wait('@metadataRequest').then((interception) => {
                expect(interception.response.statusCode).to.eq(200);
                // "Sequence ID", "Name", "Description", "Length"가 있는지 확인
                cy.contains('Sequence ID').should('be.visible');
                cy.contains('Name').should('be.visible');
                cy.contains('Description').should('be.visible');
                cy.contains('Length').should('be.visible');
            })


            const filePath = 'SARS_CoV_2/MT576556.1.spike.fasta'; // fixtures 폴더에 있는 파일 경로
            // 파일을 input[type="file"] 요소에 업로드
            cy.get('input[type="file"]').attachFile(filePath);

            // 파일이 제대로 업로드 되었는지 확인 (예: 업로드 후 확인 메시지나 파일 이름이 표시되는지)
            cy.contains('MT576556.1.spike.fasta').should('be.visible');

            cy.get('button.next-button').click();
            cy.wait('@alignmentRequest', { timeout: 20000 }).then((interception) => {
                expect(interception.response.statusCode).to.eq(200);

                // sequence-chunk 안의 sequence 클래스가 1 + 업로드한 FASTA 파일 개수(1개)인지 확인
                cy.get('.sequence-chunk').eq(0)  // 첫 번째 청크 선택
                    .find('.sequence')  // 그 안에서 sequence 클래스 요소 찾기
                    .should('have.length', 2);  // 예상 개수와 비교

            })
        })
    });


});