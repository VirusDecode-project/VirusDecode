// 서열 입력 페이지
describe('1. 레퍼런스 서열 입력 테스트', () => {
    let referenceSeqId;
    before(() => {
        cy.fixture('environment').then((environment) => {
            referenceSeqId = environment.SARS_CoV_2_ID;
        });
    });
    beforeEach(() => {
        cy.signupAndLoginIfDuplicate('testFName', 'testLName', 'testId', 'testPw', 'testPw');
        cy.intercept('POST', '/api/inputSeq/metadata').as('metadataRequest');
        cy.intercept('POST', '/api/inputSeq/alignment').as('alignmentRequest');
    });

    // 1-1-1. Nucleotide ID 입력 실패(잘못된 ID 입력)
    it('1-1-1. Nucleotide ID 입력 실패: 잘못된 ID 입력', () => {
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
    it('1-1-2. Nucleotide ID 입력 실패: Null 입력', () => {
        // 동일한 단계에서 입력하지 않고. ‘Done’ 버튼을 누른다.
        cy.get('button.done-button').click();

        // 잘못된 ID 입력 시: "NCBI에 요청한 nucleotide ID가 존재하지 않습니다."라는 메시지가 나타남.
        cy.contains('Please enter a valid sequence ID.').should('be.visible');

    });

    // 1-2. Nucleotide ID 입력 성공 및 NCBI로부터 ID 유효성 검사 및 메타데이터 가져오기
    it('1-2. NCBI로부터 ID 유효성 검사 및 메타데이터 가져오기', () => {

        //   cy.get('input[id="referenceSequenceId"]').type('NC_045512.2');
        cy.get('input[id="referenceSequenceId"]').type(referenceSeqId);
        cy.get('button.done-button').click();
        cy.wait('@metadataRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(200);

            // "Sequence ID", "Name", "Description", "Length"가 있는지 확인
            cy.contains('Sequence ID').should('be.visible');
            cy.contains('Name').should('be.visible');
            cy.contains('Description').should('be.visible');
            cy.contains('Length').should('be.visible');
        });
    });
});


describe('2. 분석 서열 입력', () => {
    let referenceSeqId;
    before(() => {
        cy.fixture('environment').then((environment) => {
            referenceSeqId = environment.SARS_CoV_2_ID;
        });
    });
    beforeEach(() => {
        cy.signupAndLoginIfDuplicate('testFName', 'testLName', 'testId', 'testPw', 'testPw');
        cy.intercept('POST', '/api/inputSeq/metadata').as('metadataRequest');
        cy.intercept('POST', '/api/inputSeq/alignment').as('alignmentRequest');
    });

    // 2-1-1 FASTA 파일 업로드 실패(잘못된 파일 업로드)
    it('2-1-1, 2-3 FASTA 파일 업로드 실패 및 입력 형식 검증: 잘못된 형식의 파일을 업로드', () => {
        const filePath = 'environment.json';
        cy.get('input[type="file"]').attachFile(filePath);

        cy.get('.message-modal-content')
            .should('be.visible')
            .and('contain', 'FASTA 파일 형식만 지원됩니다.');
        cy.get('.message-modal-content').contains('Close').click();
    });

    // 2-1-2 FASTA 파일 업로드 성공
    it('2-1-2 FASTA 파일 업로드 성공: 사용자로부터 하나의 FASTA 파일업로드', () => {
        cy.get('input[id="referenceSequenceId"]').type(referenceSeqId);
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
    });

    // 2-1-3 FASTA 파일 업로드 성공
    it('2-1-3 FASTA 파일 업로드 성공: 사용자로부터 여러 FASTA 파일업로드', () => {
        cy.get('input[id="referenceSequenceId"]').type(referenceSeqId);
        cy.get('button.done-button').click();
        cy.wait('@metadataRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(200);
            // "Sequence ID", "Name", "Description", "Length"가 있는지 확인
            cy.contains('Sequence ID').should('be.visible');
            cy.contains('Name').should('be.visible');
            cy.contains('Description').should('be.visible');
            cy.contains('Length').should('be.visible');
        })

        const filePath1 = 'SARS_CoV_2/MT576556.1.spike.fasta';
        const filePath2 = 'SARS_CoV_2/MW642250.1.spike.fasta';
        const filePath3 = 'SARS_CoV_2/OL672836.1.spike.fasta';
        const filePath4 = 'SARS_CoV_2/OM958567.1.spike.fasta';
        const filePath5 = 'SARS_CoV_2/OR240434.1.spike.fasta';
        const filePath6 = 'SARS_CoV_2/PP346415.1.spike.fasta';

        // 파일을 input[type="file"] 요소에 업로드
        cy.get('input[type="file"]').attachFile(filePath1);
        cy.get('input[type="file"]').attachFile(filePath2);
        cy.get('input[type="file"]').attachFile(filePath3);
        cy.get('input[type="file"]').attachFile(filePath4);
        cy.get('input[type="file"]').attachFile(filePath5);
        cy.get('input[type="file"]').attachFile(filePath6);

        // 파일이 제대로 업로드 되었는지 확인 (예: 업로드 후 확인 메시지나 파일 이름이 표시되는지)
        cy.contains('MT576556.1.spike.fasta').should('be.visible');
        cy.contains('MW642250.1.spike.fasta').should('be.visible');
        cy.contains('OL672836.1.spike.fasta').should('be.visible');
        cy.contains('OM958567.1.spike.fasta').should('be.visible');
        cy.contains('OR240434.1.spike.fasta').should('be.visible');
        cy.contains('PP346415.1.spike.fasta').should('be.visible');
    });

    // 2-2-1 염기서열(A,T,C,G) 입력 실패(잘못된 서열 입력)
    it('2-2-1, 2-3 염기서열(A,T,C,G) 입력 실패: 잘못된 형식의 시퀀스 입력', () => {
        cy.get('input[id="referenceSequenceId"]').type(referenceSeqId);
        cy.get('button.done-button').click();
        cy.wait('@metadataRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(200);
            cy.contains('Sequence ID').should('be.visible');
            cy.contains('Name').should('be.visible');
            cy.contains('Description').should('be.visible');
            cy.contains('Length').should('be.visible');
        })

        cy.get('.w-100')
            .type('wrongSequence');

        cy.get('button.next-button').click();
        cy.wait('@alignmentRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(500);

            cy.get('.message-modal-content')
                .should('be.visible')
                .and('contain', '입력하신 서열 정보가 올바르지 않습니다. A, T, C, 그리고 G만 허용됩니다.');
            cy.get('.message-modal-content').contains('Close').click();

        })
    });

    it('2-2-2 염기서열(A,T,C,G) 입력 성공', () => {
        cy.get('input[id="referenceSequenceId"]').type(referenceSeqId);
        cy.get('button.done-button').click();
        cy.wait('@metadataRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(200);
            // "Sequence ID", "Name", "Description", "Length"가 있는지 확인
            cy.contains('Sequence ID').should('be.visible');
            cy.contains('Name').should('be.visible');
            cy.contains('Description').should('be.visible');
            cy.contains('Length').should('be.visible');
        })

        cy.contains('div.sequence-header', 'Sequence1')  // 'Sequence1' 텍스트가 포함된 div를 찾음
            .parent()  // 부모 요소로 이동
            .find('textarea')  // 부모 요소 아래의 textarea를 찾음
            .type('ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCC');

        cy.get('button.add-sequence-button').click();

        cy.contains('div.sequence-header', 'Sequence2')  // 'Sequence2' 텍스트가 포함된 div를 찾음
            .parent()  // 부모 요소로 이동
            .find('textarea')  // 부모 요소 아래의 textarea를 찾음
            .type('ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCC');


        cy.get('button.next-button').click();
        cy.wait('@alignmentRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(200);

            // sequence-chunk 안의 sequence 클래스가 1 + 업로드한 FASTA 파일 개수(1개)인지 확인
            cy.get('.sequence-chunk').eq(0)  // 첫 번째 청크 선택
                .find('.sequence')  // 그 안에서 sequence 클래스 요소 찾기
                .should('have.length', 3);  // 예상 개수와 비교
        })
    });


});


describe('3. 유전체 분석 테스트', () => {
    let referenceSeqId;
    before(() => {
        cy.fixture('environment').then((environment) => {
            referenceSeqId = environment.SARS_CoV_2_ID;
        });
    });
    beforeEach(() => {
        cy.signupAndLoginIfDuplicate('testFName', 'testLName', 'testId', 'testPw', 'testPw');
        cy.intercept('POST', '/api/inputSeq/metadata').as('metadataRequest');
        cy.intercept('POST', '/api/inputSeq/alignment').as('alignmentRequest');
    });
    it('3-1 Next 버튼을 통해 데이터 전송 및 분석 페이지 이동', () => {
        cy.get('input[id="referenceSequenceId"]').type(referenceSeqId);
        cy.get('button.done-button').click();
        cy.wait('@metadataRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(200);
            // "Sequence ID", "Name", "Description", "Length"가 있는지 확인
            cy.contains('Sequence ID').should('be.visible');
            cy.contains('Name').should('be.visible');
            cy.contains('Description').should('be.visible');
            cy.contains('Length').should('be.visible');
        })

        const filePath1 = 'SARS_CoV_2/MT576556.1.spike.fasta';
        const filePath2 = 'SARS_CoV_2/MW642250.1.spike.fasta';

        // 파일을 input[type="file"] 요소에 업로드
        cy.get('input[type="file"]').attachFile(filePath1);
        cy.get('input[type="file"]').attachFile(filePath2);

        // 파일이 제대로 업로드 되었는지 확인 (예: 업로드 후 확인 메시지나 파일 이름이 표시되는지)
        cy.contains('MT576556.1.spike.fasta').should('be.visible');
        cy.contains('MW642250.1.spike.fasta').should('be.visible');

        // 파일이 제대로 업로드 되었는지 확인 (예: 업로드 후 확인 메시지나 파일 이름이 표시되는지)
        cy.contains('MT576556.1.spike.fasta').should('be.visible');

        // 두 개의 시퀀스 입력
        cy.contains('div.sequence-header', 'Sequence1')  // 'Sequence1' 텍스트가 포함된 div를 찾음
            .parent()  // 부모 요소로 이동
            .find('textarea')  // 부모 요소 아래의 textarea를 찾음
            .type('ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCC');

        cy.get('button.add-sequence-button').click();

        cy.contains('div.sequence-header', 'Sequence2')  // 'Sequence2' 텍스트가 포함된 div를 찾음
            .parent()  // 부모 요소로 이동
            .find('textarea')  // 부모 요소 아래의 textarea를 찾음
            .type('ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCC');


        cy.get('button.next-button').click();
        cy.wait('@alignmentRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(200);

            // sequence-chunk 안의 sequence 클래스가 1 + 업로드한 FASTA 파일 개수(1개)인지 확인
            cy.get('.sequence-chunk').eq(0)  // 첫 번째 청크 선택
                .find('.sequence')  // 그 안에서 sequence 클래스 요소 찾기
                .should('have.length', 5);  // 예상 개수와 비교
        })
    });
});



describe('4. 히스토리 생성 테스트', () => {
    let referenceSeqId;
    const fileName = 'SARS_CoV_2/MT576556.1.spike.fasta';
    const filesetup = () => {
        cy.get('input#referenceSequenceId')
            .type(referenceSeqId)
            .should('have.value', referenceSeqId);
        cy.get('button').contains('DONE').click();
        cy.wait('@metadataRequest').then((interception) => {
            expect(interception.response.statusCode).to.eq(200);
            cy.get('input[type="file"]').attachFile([fileName]);
            cy.get('button.next-button').click();
            cy.wait('@alignmentRequest').then((interception) => {
                expect(interception.response.statusCode).to.eq(200);
                cy.url().should('include', '/analysis');
            });
        });
    }
    before(() => {
        cy.fixture('environment').then((environment) => {
            referenceSeqId = environment.SARS_CoV_2_ID;
        });
    });
    beforeEach(() => {
        // 기본 URL로 애플리케이션에 접속
        cy.signupAndLoginIfDuplicate('testFName', 'testLName', 'testId', 'testPw', 'testPw');
        cy.intercept('POST', '/api/inputSeq/metadata').as('metadataRequest');
        cy.intercept('POST', '/api/inputSeq/alignment').as('alignmentRequest');
    });


    it('4-1 전송된 입력 데이터 저장', () => {
        // #4 사용자가 입력한 데이터를 새 히스토리에 저장
        filesetup();
        // 자동 생성 이름: 레퍼런스ID
        cy.get('.history-item').first().should('be.visible').and('contain', referenceSeqId);
    });

    it('4-2 히스토리 이름 자동 생성', () => {
        // #5 같은 레퍼런스ID 입력 후 분석 시작 시 히스토리 저장
        filesetup();
        cy.get('.sidebar .edit-icon').click();
        cy.get('.modal-next-button').click();
        filesetup();
        // 자동 생성 이름: 레퍼런스ID_1
        cy.get('.history-item').first().should('be.visible').and(($el) => {
            // 텍스트가 referenceSeqId로 시작하고 그 뒤에 _가 붙어있는지 확인
            const text = $el.text();
            expect(text).to.match(new RegExp(`^${referenceSeqId}_`));
        });

    });

});