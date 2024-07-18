import { Nav } from 'react-bootstrap';

function Analysis() {
    return (
        <div>
            <Nav variant="tabs" defaultActiveKey="link0">
                <Nav.Item>
                    <Nav.Link eventKey="link0">Alignment</Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link eventKey="link1">Gene</Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link eventKey="link2">Mutation</Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link eventKey="link3">Pathogenecity</Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link eventKey="link4">Analyze</Nav.Link>
                </Nav.Item>
            </Nav>
        </div>

    );
}

export default Analysis;