Cypress.Commands.add('signupAndLoginIfDuplicate', (firstName, lastName, id, password, cPassword) => {
  cy.visit('http://localhost:3000/signup');
  cy.get('input[name="firstName"]').type(firstName);
  cy.get('input[name="lastName"]').type(lastName);
  cy.get('input[name="id"]').type(id);
  cy.get('input[name="password"]').type(password);
  cy.get('input[name="cPassword"]').type(cPassword);

  cy.intercept('POST', '/api/auth/signup').as('signupRequest');
  cy.get('.SignupBtn').click(); 
  cy.wait('@signupRequest').then((interception) => {
    const { response } = interception; 
    if (response.statusCode === 200) {
      cy.get('.message-modal-content')
        .should('be.visible')
        .and('contain', '회원가입이 완료되었습니다.');
      cy.get('.message-modal-content').contains('Close').click();
      cy.url().should('include', '/login');
      cy.login(id, password); 

    } else if (response.statusCode === 400 && response.body.includes("이미 존재하는 ID 입니다.")) {
      cy.get('.message-modal-content')
        .should('be.visible')
        .and('contain', '이미 존재하는 ID 입니다.');
      cy.get('.message-modal-content').contains('Close').click();
      cy.get('.gotoLoginBtn').click(); 
      cy.login(id, password); 

    } else {
      console.error(`Error response: ${JSON.stringify(response.body)}`);
      throw new Error(`Unexpected error: ${response.statusText}`);
    }
  });
});

Cypress.Commands.add('login', (loginId, password) => {
  cy.visit('http://localhost:3000/login');
  cy.get('input[name="loginId"]').type(loginId);
  cy.get('input[name="password"]').type(password);
  cy.get('.loginBtn').click();
});

Cypress.Commands.add('guestlogin', () => {
  cy.visit('http://localhost:3000/');
  cy.get('.decode-button').click();
  cy.get('.stayLoggedOutBtn').click(); 
});
