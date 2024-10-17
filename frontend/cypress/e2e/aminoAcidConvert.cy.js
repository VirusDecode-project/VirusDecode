describe('아미노산 입력 및 변환', () => {

  beforeEach(() => {
    // 기본 URL로 애플리케이션에 접속
    cy.visit('http://localhost:3000');
    cy.contains('Try Decoding').click();
    cy.contains('stay logged out').click();
  });

  // 시나리오 ID: TS_005_1.	Select Coding sequence를 ‘S’로 선택한 다음, 
  //‘End Amino Acid Position’와 ‘Start Amino Acid Position’가 모두 1-1273 범위 이내이며, 
  //두 값의 차이가 <= 500인 범위를 입력 시: 성공적으로 변환됨.
  it('올바른 아미노산 시퀀스 범위를 입력하고 Convert버튼을 클릭', () => {

    // 1-1276 범위 내에서 랜덤한 시작 아미노산 위치 생성
    let startValue = Math.floor(Math.random() * (1276 - 1 + 1)) + 1;

    // 끝 아미노산 위치는 시작값 + 최대 500을 넘지 않도록 설정
    let endValue = Math.min(startValue + Math.floor(Math.random() * 500), 1276);

    // API 요청을 가로채고 별칭을 부여합니다.
    cy.intercept('POST', '/api/analysis/linearDesign').as('linearDesignRequest');

    cy.contains('SARS-CoV-2 유전체 분석').click();
    cy.get('.sequence-boxes').eq(0).click();

    cy.contains('.modal-content label', 'Select Coding Sequence:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('select')  // label 아래의 select 요소를 찾음
      .select('S');    // select에서 'S' 옵션 선택

    cy.contains('.modal-content label', 'Start Amino Acid Position:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('input[type="number"]')  // label 아래의 select 요소를 찾음
      .type(startValue.toString());  // 랜덤한 시작값 입력

    cy.contains('.modal-content label', 'End Amino Acid Position:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('input[type="number"]')  // label 아래의 select 요소를 찾음
      .type(endValue.toString());  // 랜덤한 끝값 입력

    cy.get('.modal-next-button')  // 클래스가 'modal-next-button'인 버튼을 선택
      .click();  // 클릭 이벤트 실행

    // 각 요청이 성공적으로 완료되었는지 상태 코드를 검증합니다.
    cy.wait('@linearDesignRequest').its('response.statusCode').should('eq', 200);

  });



  // 시나리오 ID: TS_005_2
  it('Select Coding Sequence를 S가 아닌 것을 선택한 경우', () => {
    cy.contains('SARS-CoV-2 유전체 분석').click();
    cy.get('.sequence').eq(1).click();

    cy.contains('.modal-content label', 'Select Coding Sequence:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('select')  // label 아래의 select 요소를 찾음
      .select('ORF1ab');    // select에서 S가 아닌 다른 옵션 선택

    cy.contains('.modal-content label', 'Start Amino Acid Position:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('input[type="number"]')  // label 아래의 select 요소를 찾음
      .type('1');  // 1

    cy.contains('.modal-content label', 'End Amino Acid Position:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('input[type="number"]')  // label 아래의 select 요소를 찾음
      .type('50');  // 50

    // 브라우저 알림(alert) 메시지 확인
    cy.window().then((win) => {
      cy.stub(win, 'alert').as('alert');
    });

    cy.get('.modal-next-button')  // 클래스가 'modal-next-button'인 버튼을 선택
      .click();  // 클릭 이벤트 실행


    cy.get('@alert').should('have.been.calledWith', '선택된 구간에 유효한 서열이 없습니다.');

  });






  // 시나리오 ID: TS_005_3
  it('500개 이상의 아미노산을 입력하는 경우', () => {

    // 1-1276 범위 내에서 랜덤한 시작 아미노산 위치 생성
    let startValue = Math.floor(Math.random() * (1276 - 1 + 1)) + 1;

    // 끝 아미노산 위치는 시작값 + 최소 501 이상 차이가 나도록 설정 (최대 1273까지)
    let endValue = Math.min(startValue + Math.floor(Math.random() * (1276 - startValue - 500)) + 501, 1276);

    // 해당 부분 시나리오와 다르게 무한 로딩 및 예기치못한 잉 메시지 나기도 함.
    cy.contains('SARS-CoV-2 유전체 분석').click();
    cy.get('.sequence').eq(1).click();

    cy.contains('.modal-content label', 'Select Coding Sequence:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('select')  // label 아래의 select 요소를 찾음
      .select('S');    // select에서 'S' 옵션 선택

    cy.contains('.modal-content label', 'Start Amino Acid Position:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('input[type="number"]')  // label 아래의 select 요소를 찾음
      .type(startValue.toString());  // 랜덤한 시작값 입력

    cy.contains('.modal-content label', 'End Amino Acid Position:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('input[type="number"]')  // label 아래의 select 요소를 찾음
      .type(endValue.toString());  // 랜덤한 끝값 입력

    // 브라우저 알림(alert) 메시지 확인
    cy.window().then((win) => {
      cy.stub(win, 'alert').as('alert');
    });

    cy.get('.modal-next-button')  // 클래스가 'modal-next-button'인 버튼을 선택
      .click();  // 클릭 이벤트 실행


    cy.get('@alert').should('have.been.calledWith', '필요한 파이썬 환경이 제대로 설치되지 않았습니다.');// 너무 오래 걸려서 확인 안댐...

  });




  // 시나리오 ID: TS_005_4
  it('변환 가능 범위에 벗어나는 amino acid position을 입력하였을 경우', () => {

    let startValue = Math.random() < 0.5
      ? Math.floor(Math.random() * (1276 - 1 + 1)) + 1 // 올바른 범위 내에서 생성
      : Math.floor(Math.random() * 100) + 1277; // 1276을 벗어나는 값 생성

    // 끝 아미노산 위치도 같은 방식으로 범위를 벗어나게 설정
    let endValue = Math.random() < 0.5
      ? Math.floor(Math.random() * (1276 - 1 + 1)) + 1 // 올바른 범위 내에서 생성
      : Math.floor(Math.random() * 100) + 1277; // 1276을 벗어나는 값 생성

    cy.contains('SARS-CoV-2 유전체 분석').click();
    cy.get('.sequence-boxes').eq(0).click();

    cy.contains('.modal-content label', 'Select Coding Sequence:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('select')  // label 아래의 select 요소를 찾음
      .select('S');    // select에서 'S' 옵션 선택

    cy.contains('.modal-content label', 'Start Amino Acid Position:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('input[type="number"]')  // label 아래의 select 요소를 찾음
      .type(startValue.toString());  // 랜덤한 시작값 입력

    cy.contains('.modal-content label', 'End Amino Acid Position:')  // 모달 안에서 'Select Coding Sequence:' 텍스트를 가진 label 찾기
      .find('input[type="number"]')  // label 아래의 select 요소를 찾음
      .type(endValue.toString());  // 랜덤한 끝값 입력

    cy.get('.modal-next-button')  // 클래스가 'modal-next-button'인 버튼을 선택
      .click();  // 클릭 이벤트 실행

    cy.contains('index must be between').should('be.visible');

  });

  // 시나리오 ID: TS_005_5
  it('아무것도 입력하지 않은 경우', () => {
    cy.contains('SARS-CoV-2 유전체 분석').click();
    cy.get('.sequence-boxes').eq(0).click();

    cy.get('.modal-next-button')  // 클래스가 'modal-next-button'인 버튼을 선택
      .click();  // 클릭 이벤트 실행

    cy.contains('Please enter valid indices.').should('be.visible');

  });



  // // 시나리오 ID: TS_005
  // it('올바른 아미노산 시퀀스 범위를 입력하고 Convert버튼을 클릭', () => {
  //   cy.contains('RSV 유전체 분석').click();

  // });


});