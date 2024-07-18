import { Nav } from 'react-bootstrap';
import { useState } from 'react';

function Analysis() {
    let [tab, setTab] = useState(0)
    return (
        <div>
            <Nav variant="tabs" defaultActiveKey="link0 ">
                <Nav.Item>
                    <Nav.Link eventKey="link0" onClick={()=>{setTab(0)}}>Alignment</Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link eventKey="link1" onClick={()=>{setTab(1)}}>Gene</Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link eventKey="link2" onClick={()=>{setTab(2)}}>Mutation</Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link eventKey="link3" onClick={()=>{setTab(3)}}>Pathogenecity</Nav.Link>
                </Nav.Item>
                <Nav.Item>
                    <Nav.Link eventKey="link4" onClick={()=>{setTab(4)}}>Analyze</Nav.Link>
                </Nav.Item>
            </Nav>
            <Tab tab={tab}/>
        </div>

    );
}

function Tab(props) {
    if (props.tab == 0) {
        return <div>Alignment을 분석한 결과애 대한 내용입니다.<br/>(오픈예정)</div>
    } else if (props.tab == 1) {
        return <div>Gene에 관한 내용을 나타냅니다.<br/>(예정)</div>
    } else if (props.tab == 2) {
        return <div>Mutation에 관한 정보<br/>(베타서비스중)</div>
    } else if (props.tab == 3) {
        return <div>Pathogenecity를 분석한 내용입니다.<br/>(추후오픈예정)</div>
    } else if (props.tab == 4) {
        return <div>Analyze 결과 입니다.<br/>(유료서비스)</div>
    }
}

export default Analysis;