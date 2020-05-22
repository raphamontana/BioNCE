import { useEffect } from 'react';
import dynamic from 'next/dynamic';

const NonSSRHead = dynamic(
  () => import('next/head'),
  {
    ssr: false,
    loading: () => <p>Loading JSME...</p>
  }
);

// JSME Molecule Editor
function JSMEComponent(props) {

  useEffect(() => {
    function jsmeOnLoad(status) {
      const jsmeApplet = new JSApplet.JSME("jsme_container", "100%", "100%");
      jsmeApplet.setCallBack(
        "AfterStructureModified",
        () => props.callBack(jsmeApplet.smiles())
      );
    }
    window.jsmeOnLoad = jsmeOnLoad;
  });

  return(
    <>
      <NonSSRHead>
        <script type="text/javascript" src="https://peter-ertl.com/jsme/JSME_2017-02-26/jsme/jsme.nocache.js" />
      </NonSSRHead>
      <div id={ "jsme_container" } style={{height: '311px'}} />
    </>
  );
}

export default JSMEComponent;