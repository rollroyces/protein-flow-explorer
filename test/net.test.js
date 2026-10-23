const fs=require("fs"),vm=require("vm");
const noop=()=>{};
const el=()=>({value:"",textContent:"",hidden:false,disabled:false,style:{},innerHTML:"",dataset:{},
 classList:{add:noop,remove:noop,toggle:noop,contains:()=>false},addEventListener:noop,
 appendChild:noop,querySelector:()=>el(),querySelectorAll:()=>[],setAttribute:noop,
 removeAttribute:noop,focus:noop,click:noop,remove:noop});
global.document={getElementById:()=>el(),querySelectorAll:()=>[],querySelector:()=>({dataset:{filter:"all"}}),
 addEventListener:noop,readyState:"loading",createElement:()=>el()};
global.window={location:{search:"",pathname:"/"},addEventListener:noop,history:{replaceState:noop},prompt:noop};
global.navigator={clipboard:null};global.performance={now:()=>0};global.requestAnimationFrame=noop;
const realTimeout=setTimeout, realClear=clearTimeout;
global.setTimeout=()=>0;global.clearTimeout=noop;  // app must not rely on timers for correctness
global.$3Dmol={createViewer:()=>({setStyle:noop,removeAllModels:noop,removeAllShapes:noop,
 removeAllLabels:noop,removeAllSurfaces:noop,zoomTo:noop,render:noop,setBackgroundColor:noop,
 resize:noop,spin:noop,querySelector:()=>el(),addStyle:noop,addModel:()=>({}),superpose:()=>true,
 setClickable:noop,pngURI:()=>"",addSphere:noop,addCylinder:noop,addLabel:noop,addSurface:noop,
 removeModel:noop,getIntersection:()=>null}),SurfaceType:{VDW:1}};

// Stub fetch with programmable behaviour BEFORE loading app.js
let calls=[];
let behaviour=(url)=>({ok:true,status:200,text:async()=>"ATOM",json:async()=>({})});
global.fetch=async(url)=>{calls.push(url);return behaviour(url);};

const ROOT=require("path").join(__dirname,"..");
vm.runInThisContext(fs.readFileSync(ROOT+"/data/proteins.js","utf8"),{filename:"proteins.js"});
let src=fs.readFileSync(ROOT+"/app.js","utf8");
src=src.replace(/\n\}\)\(\);\s*$/,
 "\n  global.__T__={net,netGeneration,newNetworkGeneration,isAbort};\n})();");
vm.runInThisContext(src,{filename:"app.js"});
const T=global.__T__;
// The app uses setTimeout for backoff; give it a working one now that app.js is loaded.
global.setTimeout=realTimeout; global.clearTimeout=realClear;

(async()=>{
  // 1. happy path: no retry on 2xx
  calls=[];
  const r1=await T.net("https://x/ok");
    console.log("1. happy path        :", r1.status, "calls="+calls.length, calls.length===1?"PASS":"FAIL");

  // 2. 4xx is an answer, not a retry
  calls=[];
  behaviour=()=>({ok:false,status:404});
  const r2=await T.net("https://x/404");
  console.log("2. 4xx no retry      :", r2.status, "calls="+calls.length, calls.length===1?"PASS":"FAIL");

  // 3. 5xx retries once then succeeds
  calls=[];
  let n=0;
  behaviour=()=>{n++;return n===1?{ok:false,status:503}:{ok:true,status:200};};
  const r3=await T.net("https://x/503");
  console.log("3. 5xx retried       :", r3.status, "calls="+calls.length, (calls.length===2&&r3.status===200)?"PASS":"FAIL");

  // 4. 5xx exhausts retries then throws
  calls=[];
  behaviour=()=>({ok:false,status:500});
  let threw=false;
  try{await T.net("https://x/500");}catch(e){threw=true;}
  console.log("4. 5xx gives up      :", "threw="+threw, "calls="+calls.length, (threw&&calls.length===2)?"PASS":"FAIL");

  // 5. superseded generation is rejected before any network call
  calls=[];
  const g=T.newNetworkGeneration();
  let stale=false;
  try{await T.net("https://x/stale",undefined,g-1);}catch(e){stale=T.isAbort(e);}
  console.log("5. stale rejected    :", "isAbort="+stale, "calls="+calls.length, (stale&&calls.length===0)?"PASS":"FAIL");

  // 6. current generation still works
  const g2=T.newNetworkGeneration();
  behaviour=()=>({ok:true,status:200});
  calls=[];
  const r6=await T.net("https://x/fresh",undefined,g2);
  console.log("6. current gen works :", r6.status, calls.length===1?"PASS":"FAIL");
})();
