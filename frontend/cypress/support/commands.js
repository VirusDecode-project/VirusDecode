Cypress.Commands.add('signup', (firstName, lastName, id, password, cPassword) => {
  cy.visit('http://localhost:3000/signup');
  cy.get('input[name="firstName"]').type(firstName);
  cy.get('input[name="lastName"]').type(lastName);
  cy.get('input[name="id"]').type(id);
  cy.get('input[name="password"]').type(password);
  cy.get('input[name="cPassword"]').type(cPassword);
  cy.get('.SignupBtn').click();
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

const checkDuplicateId = (firstName, lastName, loginId, password) => {
  return cy.request({
    method: 'POST',
    url: `/api/auth/signup`, 
    body: {
      firstName,
      lastName,
      loginId,
      password,
    },
    failOnStatusCode: false, 
  }).then((response) => {
    if (response.status === 200) {
      return false; 
    } else if (response.status === 400 && response.body.includes("이미 존재하는 ID 입니다.")) {
      return true; 
    } else {
      throw new Error(`Unexpected error: ${response.statusText}`);
    }
  });
};

Cypress.Commands.add('signupAndLoginIfDuplicate', (firstName, lastName, id, password, cPassword) => {
  cy.signup(firstName, lastName, id, password, cPassword); 
  checkDuplicateId(firstName, lastName, id, password).then((isDuplicate) => {
    if (isDuplicate) {
      cy.get('.message-modal-content').contains('Close').click();
      cy.get('.gotoLoginBtn').click();
      cy.login(id, password);
      cy.url().should('include', '/inputSeq');
    } else {
      cy.get('.message-modal-content')
      .should('be.visible')
      .and('contain', '회원가입이 완료되었습니다.');
      cy.get('.message-modal-content').contains('Close').click();
      cy.url().should('include', '/login');
      cy.login(id, password);
    }
  });
});