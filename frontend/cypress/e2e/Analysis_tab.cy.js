describe('Analysis component_tab', () => {
  const referenceSeqId = 'NC_045512';
  const fileName1 = 'MT576556.1.spike.fasta';
  const loginAndSetup = () => {
    cy.visit('http://localhost:3000/');
    cy.get('.decode-button').click();
    cy.get('.stayLoggedOutBtn').click(); 
    cy.get('input#referenceSequenceId') 
      .type(referenceSeqId)
      .should('have.value', referenceSeqId); 
    cy.get('button').contains('DONE').click(); 
    cy.wait('@metadataRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('input[type="file"]').attachFile([fileName1]);   
      cy.get('button.next-button').click();
      cy.wait('@alignmentRequest').then((interception) => {
        expect(interception.response.statusCode).to.eq(200);
        cy.url().should('include', '/analysis');
      });
    });
  }

  const convertmodal = () => {
    cy.get('.sequence-boxes').eq(0).click();
    cy.contains('.modal-content label', 'Select Coding Sequence:').find('select').select('S');    
    cy.contains('.modal-content label', 'Start Amino Acid Position:').find('input[type="number"]').type(1);
    cy.contains('.modal-content label', 'End Amino Acid Position:').find('input[type="number"]').type(10);  
    cy.get('.modal-next-button').click();
  }

  beforeEach(() => {
    cy.intercept('POST', '/api/inputSeq/metadata').as('metadataRequest');
    cy.intercept('POST', '/api/inputSeq/alignment').as('alignmentRequest');
    cy.intercept('POST', '/api/analysis/linearDesign').as('linearDesignRequest');
    cy.intercept('POST', '/api/analysis/pdb').as('PDBrequest');
    loginAndSetup();
  });

  it('enable tab step by step', () => {
    // 1. 분석 시작
    cy.get('.nav-tabs').contains('Alignment').should('have.class', 'active'); 
    cy.get('.nav-tabs').contains('mRNA').should('have.class','disabled'); 
    cy.get('.nav-tabs').contains('3D').should('have.class','disabled'); 

    // 2. mRNA Conversion
    convertmodal();
    cy.wait('@linearDesignRequest').then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('.nav-tabs').contains('mRNA').should('have.class', 'active'); 
      cy.get('.nav-tabs').contains('3D').should('have.class','disabled');
    });

    // 3. 3D Rendering
    cy.wait('@PDBrequest', { timeout: 20000 }).then((interception) => {
      expect(interception.response.statusCode).to.eq(200);
      cy.get('.nav-tabs').contains('mRNA').should('have.class', 'active'); 
      cy.get('.nav-tabs').contains('3D').should('not.have.class', 'disabled'); 
    });
    cy.get('.nav-tabs').contains('3D').click();
    cy.get('.nav-tabs').contains('3D').should('have.class', 'active'); 
  });
});