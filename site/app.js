"use strict";

const VBE = 0.7, VT = 0.026, VCESAT = 0.2;
const TAU = Math.PI*2;

const PARAMS = [
  { id:"VCC",  name:"V<sub>CC</sub>", min:5,    max:18,   step:0.5, val:12,   unit:"V",  fmt:"V",  grp:"Besleme & Öngerilim" },
  { id:"R1",   name:"R₁",            min:10000,max:100000,step:1000,val:47000,unit:"Ω", fmt:"R",  grp:"Besleme & Öngerilim" },
  { id:"R2",   name:"R₂",            min:2000, max:47000,step:500, val:10000,unit:"Ω", fmt:"R",  grp:"Besleme & Öngerilim" },
  { id:"RC",   name:"R<sub>C</sub>", min:500,  max:10000,step:100, val:2200, unit:"Ω", fmt:"R",  grp:"Yük & Dirençler" },
  { id:"RE",   name:"R<sub>E</sub>", min:100,  max:5000, step:100, val:1000, unit:"Ω", fmt:"R",  grp:"Yük & Dirençler" },
  { id:"RL",   name:"R<sub>L</sub>", min:1000, max:100000,step:1000,val:10000,unit:"Ω", fmt:"R",  grp:"Yük & Dirençler" },
  { id:"Rs",   name:"R<sub>s</sub>", min:50,   max:2000, step:50,  val:600,  unit:"Ω", fmt:"R",  grp:"Yük & Dirençler" },
  { id:"CB",   name:"C<sub>B</sub>", min:0.1,  max:100,  step:0.1, val:10,   unit:"µF", fmt:"uF", grp:"Kondansatörler" },
  { id:"CC",   name:"C<sub>C</sub>", min:0.1,  max:100,  step:0.1, val:10,   unit:"µF", fmt:"uF", grp:"Kondansatörler" },
  { id:"CE",   name:"C<sub>E</sub>", min:1,    max:470,  step:1,   val:100,  unit:"µF", fmt:"uF", grp:"Kondansatörler" },
  { id:"beta", name:"β (h<sub>FE</sub>)", min:50, max:400, step:10, val:150, unit:"", fmt:"int",  grp:"Transistör" },
  { id:"vs",   name:"v<sub>s</sub> genlik", min:1, max:50, step:1, val:10, unit:"mV", fmt:"mV",   grp:"Giriş İşareti" },
  { id:"freq", name:"frekans", min:100, max:5000, step:100, val:1000, unit:"Hz", fmt:"Hz",       grp:"Giriş İşareti" },
];
const S = {};
PARAMS.forEach(p => S[p.id] = p.val);

/* ===== biçimlendirme ===== */
const par = (a,b) => (a*b)/(a+b);
function fR(x){
  if(x>=1000){ const k=x/1000; const s=Number.isInteger(k)?String(k):k.toFixed(2).replace(/\.?0+$/,""); return s+" kΩ"; }
  return Math.round(x)+" Ω";
}
function fV(x){ return x.toFixed(2)+" V"; }
function fI(x){ // amper cinsinden
  const a=Math.abs(x);
  if(a>=1e-3) return (x*1e3).toFixed(2)+" mA";
  if(a>=1e-6) return (x*1e6).toFixed(2)+" µA";
  return (x*1e9).toFixed(1)+" nA";
}
function fmtParam(p,v){
  if(p.fmt==="R") return fR(v);
  if(p.fmt==="V") return v.toFixed(1)+" V";
  if(p.fmt==="int") return String(Math.round(v));
  if(p.fmt==="mV") return v+" mV";
  if(p.fmt==="uF") return (Number.isInteger(v)?v:v.toFixed(1))+" µF";
  if(p.fmt==="Hz") return v>=1000?(v/1000)+" kHz":v+" Hz";
  return v;
}

/* ===== karmaşık AC yardımcı fonksiyonları ===== */
function complex(real, imag){
  return { re:real, im:imag||0 };
}
function cAdd(a,b){ return complex(a.re+b.re,a.im+b.im); }
function cSub(a,b){ return complex(a.re-b.re,a.im-b.im); }
function cMul(a,b){ return complex(a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re); }
function cDiv(a,b){
  const denom=b.re*b.re+b.im*b.im;
  return complex((a.re*b.re+a.im*b.im)/denom,(a.im*b.re-a.re*b.im)/denom);
}
function cInv(a){ return cDiv(complex(1),a); }
function cAbs(a){ return Math.hypot(a.re,a.im); }
function cScale(a,scale){ return complex(a.re*scale,a.im*scale); }
function cParallel(a,b){ return cInv(cAdd(cInv(a),cInv(b))); }
function capZ(capacitanceUf, freq){
  const capacitance=Math.max(capacitanceUf,0.001)*1e-6;
  return complex(0,-1/(TAU*Math.max(freq,1)*capacitance));
}
function signedGain(value){
  const mag=cAbs(value);
  return value.re<0 ? -mag : mag;
}
function solveComplexLinear(matrix, vector){
  const size=vector.length;
  const a=matrix.map(row=>row.map(cell=>complex(cell.re,cell.im)));
  const b=vector.map(cell=>complex(cell.re,cell.im));

  for(let col=0; col<size; col++){
    let pivot=col;
    for(let row=col+1; row<size; row++){
      if(cAbs(a[row][col])>cAbs(a[pivot][col])) pivot=row;
    }
    if(pivot!==col){
      [a[col],a[pivot]]=[a[pivot],a[col]];
      [b[col],b[pivot]]=[b[pivot],b[col]];
    }
    for(let row=col+1; row<size; row++){
      const factor=cDiv(a[row][col],a[col][col]);
      for(let k=col; k<size; k++) a[row][k]=cSub(a[row][k],cMul(factor,a[col][k]));
      b[row]=cSub(b[row],cMul(factor,b[col]));
    }
  }

  const x=Array(size).fill(null).map(()=>complex(0));
  for(let row=size-1; row>=0; row--){
    let sum=complex(0);
    for(let col=row+1; col<size; col++) sum=cAdd(sum,cMul(a[row][col],x[col]));
    x[row]=cDiv(cSub(b[row],sum),a[row][row]);
  }
  return x;
}

function solveSmallSignal(params){
  const rb=par(params.R1,params.R2);
  const zSrc=cAdd(complex(params.Rs),capZ(params.CB,params.freq));
  const zEmitter=cParallel(complex(params.RE),capZ(params.CE,params.freq));
  const zOutCap=capZ(params.CC,params.freq);
  const ySrc=cInv(zSrc), yBias=cInv(complex(rb)), yPi=cInv(complex(params.rpi));
  const yEmitter=cInv(zEmitter), yRc=cInv(complex(params.RC));
  const yOutCap=cInv(zOutCap), yLoad=cInv(complex(params.RL));
  const gm=complex(params.gm);
  const source=complex(params.vsAmp);

  const matrix=[
    [cAdd(cAdd(ySrc,yBias),yPi), cScale(yPi,-1), complex(0), complex(0)],
    [cScale(cAdd(yPi,gm),-1), cAdd(cAdd(yPi,yEmitter),gm), complex(0), complex(0)],
    [gm, cScale(gm,-1), cAdd(yRc,yOutCap), cScale(yOutCap,-1)],
    [complex(0), complex(0), cScale(yOutCap,-1), cAdd(yOutCap,yLoad)]
  ];
  const vector=[cMul(source,ySrc), complex(0), complex(0), complex(0)];
  const [vb,ve,vc,vo]=solveComplexLinear(matrix,vector);
  const inputCurrent=cMul(cSub(source,vb),ySrc);
  return { vb, ve, vc, vo, inputCurrent, rb };
}

function solve(){
  const {VCC,R1,R2,RC,RE,RL,Rs,beta,CB,CC,CE,freq} = S;
  const vsAmp = S.vs/1000; // V

  // DC · Gerilim bölücü Thevenin eşdeğeri
  const VTh = VCC*R2/(R1+R2);
  const RTh = par(R1,R2);

  let region, IB, IC, IE, VCEq, VE, VB, VCB, ICsat;
  ICsat = (VCC-VCESAT)/(RC+RE);

  if(VTh <= 0.5){
    region="KESİM"; IB=0; IC=0; IE=0; VCEq=VCC; VE=0; VB=VTh; VCB=VCC-VBE;
  } else {
    IB = (VTh-VBE)/(RTh+(beta+1)*RE);
    let icActive = beta*IB;
    if(icActive >= ICsat){
      region="DOYMA"; IC=ICsat; IE=ICsat; VCEq=VCESAT;
      VE=IE*RE; VB=VE+VBE; VCB=VCEq-VBE; IB=icActive/beta;
    } else {
      region="AKTİF"; IC=icActive; IE=(beta+1)*IB; VE=IE*RE; VB=VE+VBE;
      VCEq=VCC-IC*RC-IE*RE; VCB=VCEq-VBE;
    }
  }

  const active = region==="AKTİF";
  const ICQ = IC;

  // AC · hibrit-π modeli (ro → sonsuz)
  let gm=0, rpi=0, Rg=0, Kv=0, Kv0=0, Ro=RC, Ki=0, Kvs=0, viPeak=0, voPeak=0;
  let viPhasor=complex(0), voPhasor=complex(0);
  if(active && ICQ>0){
    gm = ICQ/VT;
    rpi = beta/gm;
    const acParams={VCC,R1,R2,RC,RE,RL,Rs,beta,CB,CC,CE,freq,gm,rpi,vsAmp};
    const ac=solveSmallSignal(acParams);
    const acOpen=solveSmallSignal({...acParams,RL:1e12});
    const source=complex(vsAmp);
    const inputCurrentMag=cAbs(ac.inputCurrent);
    const kvComplex=cDiv(ac.vo,ac.vb);
    const kv0Complex=cDiv(acOpen.vo,acOpen.vb);
    const sourceToOutput=cDiv(ac.vo,source);
    const loadCurrent=cDiv(ac.vo,complex(RL));
    const currentGain=inputCurrentMag>1e-15 ? cDiv(loadCurrent,ac.inputCurrent) : complex(0);
    Rg = inputCurrentMag>1e-15 ? cAbs(cDiv(ac.vb,ac.inputCurrent)) : 0;
    Kv = signedGain(kvComplex);
    Kv0 = signedGain(kv0Complex);
    Ro = RC;
    Ki = signedGain(currentGain);
    Kvs = signedGain(sourceToOutput);
    viPeak = cAbs(ac.vb);
    voPeak = cAbs(ac.vo);
    viPhasor = ac.vb;
    voPhasor = ac.vo;
  }

 // Kırpılma öncesi çıkış salınım payı
  const headroom = active ? Math.min(VCEq-VCESAT, ICQ*par(RC,RL)) : 0;

  return {VTh,RTh,IB,IC,IE,ICQ,VCEq,VE,VB,VCB,ICsat,region,active,
          gm,rpi,Rg,Kv,Kv0,Ro,Ki,Kvs,viPeak,voPeak,viPhasor,voPhasor,headroom,vsAmp,freq};
}

 /* ===== SVG şema çizim elemanları ===== */
const SB = "var(--sch-bg)";
function wire(x1,y1,x2,y2){ return `<path class="wire" d="M${x1} ${y1} L${x2} ${y2}"/>`; }
function poly(pts){ return `<path class="wire" d="M${pts}"/>`; }
function node(x,y){ return `<circle class="node" cx="${x}" cy="${y}" r="3"/>`; }
function gnd(x,y){ return `<g stroke="var(--accent)" stroke-width="2" fill="none">
  <path d="M${x} ${y-10} L${x} ${y}"/><path d="M${x-11} ${y} L${x+11} ${y}"/>
  <path d="M${x-7} ${y+4} L${x+7} ${y+4}"/><path d="M${x-3} ${y+8} L${x+3} ${y+8}"/></g>`; }
// Direnç = teli örten kutu (tel alttan ayrı olarak çizilir)
function resV(x,y1,y2,label,vid,side){
  const m=(y1+y2)/2;
  const left=side==="left";
  const textX=left?x-14:x+15;
  const anchor=left?' text-anchor="end"':"";
  return `${wire(x,y1,x,y2)}
    <rect x="${x-9}" y="${m-22}" width="18" height="44" rx="2" fill="${SB}" class="comp"/>
    <text class="lbl" x="${textX}" y="${m-4}"${anchor}>${label}</text>
    ${vid?`<text class="lbl-v" id="${vid}" x="${textX}" y="${m+10}"${anchor}></text>`:""}`;
}
function resH(x1,x2,y,label,vid){
  const m=(x1+x2)/2;
  return `${wire(x1,y,x2,y)}
    <rect x="${m-22}" y="${y-9}" width="44" height="18" rx="2" fill="${SB}" class="comp"/>
    <text class="lbl" x="${m}" y="${y-15}" text-anchor="middle">${label}</text>
    ${vid?`<text class="lbl-v" id="${vid}" x="${m}" y="${y+24}" text-anchor="middle"></text>`:""}`;
}
function capV(x,y1,y2,label,vid){
  const m=(y1+y2)/2;
  return `${wire(x,y1,x,m-4)}${wire(x,m+4,x,y2)}
    <path class="comp" d="M${x-11} ${m-4} L${x+11} ${m-4}"/>
    <path class="comp" d="M${x-11} ${m+4} L${x+11} ${m+4}"/>
    <text class="lbl" x="${x+15}" y="${m+3}">${label}</text>${vid?`<text class="lbl-v" id="${vid}" x="${x+15}" y="${m+17}"></text>`:""}`;
}
function capH(x1,x2,y,label){
  const m=(x1+x2)/2;
  return `${wire(x1,y,m-4,y)}${wire(m+4,y,x2,y)}
    <path class="comp" d="M${m-4} ${y-11} L${m-4} ${y+11}"/>
    <path class="comp" d="M${m+4} ${y-11} L${m+4} ${y+11}"/>
    <text class="lbl" x="${m}" y="${y-15}" text-anchor="middle">${label}</text>`;
}
function source(x,y,label,vid){
  return `<circle cx="${x}" cy="${y}" r="16" fill="${SB}" class="comp"/>
    <path class="comp" d="M${x-8} ${y} q4 -8 8 0 q4 8 8 0" fill="none"/>
    <text class="lbl" x="${x-22}" y="${y+4}" text-anchor="end">${label}</text>${vid?`<text class="lbl-v" id="${vid}" x="${x-22}" y="${y+18}" text-anchor="end"></text>`:""}`;
}
function transistor(x,y){ // NPN, baz solda
  return `<circle cx="${x}" cy="${y}" r="22" fill="none" class="comp" opacity="0.55"/>
    <path class="comp" d="M${x-9} ${y-13} L${x-9} ${y+13}" stroke-width="2.5"/>
    <path class="comp" d="M${x-9} ${y-6} L${x+7} ${y-17}"/>
    <path class="comp" d="M${x-9} ${y+6} L${x+7} ${y+17}"/>
    <path class="comp-fill" d="M${x+1} ${y+8} L${x+8} ${y+18} L${x-2} ${y+15} Z"/>`;
}
function arrow(x,y,dir,color){ // yön: 'd' aşağı, 'u' yukarı, 'r' sağ
  color=color||"var(--accent-2)";
  let p;
  if(dir==="d") p=`${x-4},${y-5} ${x+4},${y-5} ${x},${y+3}`;
  else if(dir==="u") p=`${x-4},${y+5} ${x+4},${y+5} ${x},${y-3}`;
  else p=`${x-5},${y-4} ${x-5},${y+4} ${x+3},${y}`;
  return `<polygon points="${p}" fill="${color}"/>`;
}
function svg(vb,inner){ return `<svg viewBox="0 0 ${vb}" preserveAspectRatio="xMidYMid meet">${inner}</svg>`; }

/* ===== Şema 1: Tam Ortak Emitör (CE) Yükselteci ===== */
function drawAmp(){
  let g="";
  // VCC besleme hattı
  g+=wire(70,28,360,28);
  g+=`<text class="rail-t" x="40" y="32">+Vcc</text>`;
  g+=`<text class="lbl-v" id="ampVcc" x="62" y="20"></text>`;
  // R1 gerilim bölücü
  g+=wire(160,28,160,150)+resV(160,40,116,"R₁","ampR1");
  // Baz düğümü
  g+=node(160,150);
  // R2
  g+=wire(160,150,160,252)+resV(160,168,232,"R₂","ampR2")+gnd(160,262);
  // transistor
  g+=transistor(212,150);
  g+=wire(160,150,203,150); // Baz bağlantısı
  g+=`<text class="lbl-v" x="154" y="127" text-anchor="end">v<tspan dy='3'>i</tspan></text>`;
  g+=poly(`219 133 L232 122`); // Kollektör bağlantısı
  g+=poly(`219 167 L232 178`); // Emitör bağlantısı
  // RC ve kollektör
  g+=wire(232,28,232,122)+resV(232,40,96,"R<tspan dy='3'>C</tspan>","ampRC");
  g+=node(232,108);
  // Çıkış kuplajı
  g+=wire(232,108,290,108)+capH(290,322,108,"C<tspan dy='3'>C</tspan>");
  g+=wire(322,108,360,108)+node(360,108);
  g+=`<text class="lbl-v" x="360" y="97" text-anchor="middle">v<tspan dy='3'>o</tspan></text>`;
  // RL
  g+=wire(360,108,360,226)+resV(360,150,214,"R<tspan dy='3'>L</tspan>","ampRL")+gnd(360,236);
  // emitter / RE
  g+=node(232,178);
  g+=wire(232,178,232,290)+resV(232,210,276,"R<tspan dy='3'>E</tspan>","ampRE")+gnd(232,300);
  // CE bypass
  g+=node(232,208)+wire(232,208,300,208)+capV(300,208,260,"C<tspan dy='3'>E</tspan>")+gnd(300,270);
  // giriş kaynak zinciri
  g+=source(40,150,"v<tspan dy='3'>s</tspan>");
  g+=wire(40,166,40,196)+gnd(40,196);
  g+=wire(56,150,74,150)+resH(74,110,150,"R<tspan dy='3'>s</tspan>","ampRs");
  g+=wire(110,150,118,150)+capH(118,150,150,"C<tspan dy='3'>B</tspan>")+wire(150,150,160,150);
  return svg("480 340",g);
}

/* ===== Şema 2: DC Eşdeğeri ===== */
function drawDC(){
  let g="";
  // VCC rayı
  g+=wire(70,28,400,28);
  g+=`<text class="rail-t" x="40" y="32">+Vcc</text>`;
  // R1 / R2 bölücü (sol)
  g+=wire(170,28,170,150)+resV(170,40,116,"R₁","dcR1");
  g+=node(170,150);
  g+=wire(170,150,170,252)+resV(170,168,232,"R₂","dcR2")+gnd(170,262);
  // transistör (genişçe sağda)
  g+=transistor(270,150);
  g+=wire(170,150,259,150); // baz iletkeni
  g+=arrow(212,150,"r")+`<text class="lbl-v" id="dcIB" x="212" y="142" text-anchor="middle"></text><text class="lbl-m" x="212" y="167" text-anchor="middle">I<tspan dy='3'>B</tspan></text>`;
  // kollektör / emiter iletkenleri RC sütununa (x=310)
  g+=poly(`277 133 L310 122`);
  g+=poly(`277 167 L310 178`);
  // RC (kollektör)
  g+=wire(310,28,310,122)+resV(310,40,96,"R<tspan dy='3'>C</tspan>","dcRC");
  g+=arrow(310,112,"u","var(--accent-2)")+`<text class="lbl-v" id="dcIC" x="320" y="116"></text><text class="lbl-m" x="320" y="104">I<tspan dy='3'>C</tspan></text>`;
  // VCE köşeli ayraç (en sağ, ayrı sütun)
  g+=`<path class="comp" d="M362 122 L370 122 M366 122 L366 178 M362 178 L370 178" opacity="0.7"/>`;
  g+=`<text class="lbl-v" id="dcVCE" x="374" y="146"></text><text class="lbl-m" x="374" y="159">V<tspan dy='3'>CE</tspan></text>`;
  g+=node(310,178);
  // RE (emiter)
  g+=wire(310,178,310,290)+resV(310,210,276,"R<tspan dy='3'>E</tspan>","dcRE")+gnd(310,300);
  g+=arrow(310,235,"d")+`<text class="lbl-v" id="dcIE" x="276" y="238" text-anchor="end"></text><text class="lbl-m" x="300" y="253" text-anchor="end">I<tspan dy='3'>E</tspan></text>`;
  return svg("480 340",g);
}

/* ===== Şema 3: AC Eşdeğeri (Hibrit-π Modeli) ===== */
function drawAC(){
  let g="";
  // giriş kaynağı vs + Rs
  g+=source(40,150,"v<tspan dy='3'>s</tspan>");
  g+=wire(40,166,40,200)+gnd(40,200);
  g+=wire(56,150,72,150)+resH(72,108,150,"R<tspan dy='3'>s</tspan>",null);
  g+=wire(108,150,150,150);
  // baz düğümü B
  g+=node(150,150);
  g+=`<text class="lbl-m" x="150" y="141" text-anchor="middle">B</text>`;
  // RB = R1||R2  (baz → toprak)
  g+=wire(150,150,150,260)+resV(150,180,236,"R<tspan dy='3'>B</tspan>","acRB","left")+gnd(150,260);
  // rπ  (baz → toprak)
  g+=wire(150,150,205,150)+node(205,150);
  g+=wire(205,150,205,260)+resV(205,180,236,"r<tspan dy='3'>π</tspan>","acRpi")+gnd(205,260);
  // vbe (rπ üzerindeki gerilim)
  g+=`<path class="wire" d="M236 176 L236 232" stroke-dasharray="3 3" opacity="0.65"/>`;
  g+=arrow(236,178,"u")+`<text class="lbl-v" x="244" y="186">v<tspan dy='3'>be</tspan></text>`;
  // giriş | çıkış ayıracı
  g+=`<path d="M262 44 L262 282" stroke="#b9b8ae" stroke-width="1.4" stroke-dasharray="5 4"/>`;
  // kollektör rayı C
  g+=wire(300,120,440,120)+node(300,120);
  g+=`<text class="lbl-m" x="300" y="111" text-anchor="middle">C</text>`;
  // bağımlı akım kaynağı gm·vbe (eşkenar dörtgen → toprak)
  g+=wire(300,120,300,150);
  g+=`<path class="comp" d="M300 150 l15 17 l-15 17 l-15 -17 z" fill="${SB}"/>`;
  g+=arrow(300,159,"u")+`<text class="lbl" x="318" y="171">g<tspan dy='3'>m</tspan>v<tspan dy='3'>be</tspan></text>`;
  g+=wire(300,184,300,232)+gnd(300,232);
  // ro, RC, RL  (raydan → toprak)
  g+=node(345,120)+wire(345,120,345,232)+resV(345,150,206,"r<tspan dy='3'>o</tspan>","acRo")+gnd(345,232);
  g+=node(392,120)+wire(392,120,392,232)+resV(392,150,206,"R<tspan dy='3'>C</tspan>",null)+gnd(392,232);
  g+=wire(440,120,440,232)+resV(440,150,206,"R<tspan dy='3'>L</tspan>",null)+gnd(440,232);
  // çıkış düğümü vo
  g+=node(440,120);
  g+=`<text class="lbl-v" x="440" y="111" text-anchor="middle">v<tspan dy='3'>o</tspan></text>`;
  return svg("480 340",g);
}

function buildControls(){
  const c=document.getElementById("controls");
  let html="", curGrp="";
  PARAMS.forEach(p=>{
    if(p.grp!==curGrp){ curGrp=p.grp; html+=`<div class="ctrl-sec">${curGrp}</div>`; }
    html+=`
    <div class="ctrl">
      <span class="name">${p.name}</span>
      <span class="val" id="val_${p.id}">${fmtParam(p,p.val)}</span>
      <input type="range" id="rng_${p.id}" min="${p.min}" max="${p.max}" step="${p.step}" value="${p.val}">
    </div>`;
  });
  c.innerHTML=html;
  PARAMS.forEach(p=>{
    document.getElementById("rng_"+p.id).addEventListener("input",e=>{
      S[p.id]=parseFloat(e.target.value);
      document.getElementById("val_"+p.id).textContent=fmtParam(p,S[p.id]);
      update();
    });
  });
}

/* ===== sonuç tabloları ===== */
function row(k,v,invalid){ return `<tr${invalid?' class="invalid"':''}><td class="k">${k}</td><td class="v">${v}</td></tr>`; }

function renderResults(r){
  const inv=!r.active;
  // DC tablosu
  const dc=document.getElementById("resDC");
  dc.innerHTML =
    `<tr><th>Büyüklük</th><th>Değer</th></tr>`+
    row("I<sub>B</sub>", fI(r.IB)) +
    row("I<sub>C</sub> (I<sub>CQ</sub>)", fI(r.IC)) +
    row("I<sub>E</sub>", fI(r.IE)) +
    row("V<sub>CE(Q)</sub>", fV(r.VCEq)) +
    row("V<sub>E</sub>", fV(r.VE)) +
    row("V<sub>CB</sub>", fV(r.VCB)) +
    row("V<sub>BE</sub>", fV(r.region==="KESİM"?0:VBE)) +
    row("I<sub>C(sat)</sub>", fI(r.ICsat)) +
    row("V<sub>Th</sub>", fV(r.VTh)) +
    row("R<sub>Th</sub>", fR(r.RTh));
  // AC tablosu
  const ac=document.getElementById("resAC");
  ac.innerHTML =
    `<tr><th>Büyüklük</th><th>Değer</th></tr>`+
    row("K<sub>v</sub> (gerilim)", inv?"·":r.Kv.toFixed(1), inv) +
    row("K<sub>v0</sub> (yüksüz)", inv?"·":r.Kv0.toFixed(1), inv) +
    row("K<sub>i</sub> (akım)", inv?"·":r.Ki.toFixed(1), inv) +
    row("K<sub>vg</sub> (kaynak→yük)", inv?"·":r.Kvs.toFixed(1), inv) +
    row("R<sub>g</sub> (giriş)", inv?"·":fR(r.Rg), inv) +
    row("R<sub>o</sub> (çıkış)", inv?"·":fR(r.Ro), inv) +
    row("g<sub>m</sub>", inv?"·":(r.gm*1000).toFixed(2)+" mS", inv) +
    row("r<sub>π</sub>", inv?"·":fR(r.rpi), inv) +
    row("v<sub>o</sub> tepe", inv?"·":(Math.abs(r.voPeak)*1000).toFixed(1)+" mV", inv) +
    row("v<sub>i</sub> tepe", inv?"·":(r.viPeak*1000).toFixed(1)+" mV", inv);
  // bölge bandı
  const b=document.getElementById("regionBanner");
  let cls,msg;
  if(r.region==="AKTİF"){cls=""; msg=`AKTİF BÖLGE - yükselteç çalışıyor. I<sub>C(sat)</sub>=${fI(r.ICsat)}. Tablo 2 geçerli.`;}
  else if(r.region==="DOYMA"){cls="sat"; msg=`DOYMA BÖLGESİ - V<sub>CE</sub>≤${VCESAT} V. AC kazançları geçersiz.`;}
  else {cls="cut"; msg="KESİM BÖLGESİ - transistör iletmiyor. AC kazançları geçersiz.";}
  b.className="hg-desc "+cls; b.innerHTML=msg;
}

/* ===== Şemadaki değer etiketlerini güncelleme ===== */
function setTxt(id,v){ const e=document.getElementById(id); if(e) e.textContent=v; }
function updateLabels(r){
  setTxt("ampVcc",S.VCC.toFixed(1)+"V");
  setTxt("ampR1",fR(S.R1)); setTxt("ampR2",fR(S.R2)); setTxt("ampRC",fR(S.RC));
  setTxt("ampRE",fR(S.RE)); setTxt("ampRL",fR(S.RL)); setTxt("ampRs",fR(S.Rs));
  setTxt("dcR1",fR(S.R1)); setTxt("dcR2",fR(S.R2)); setTxt("dcRC",fR(S.RC)); setTxt("dcRE",fR(S.RE));
  setTxt("dcIB",fI(r.IB)); setTxt("dcIC",fI(r.IC)); setTxt("dcIE",fI(r.IE)); setTxt("dcVCE",fV(r.VCEq));
  setTxt("acRB", r.active?fR(par(S.R1,S.R2)):"·");
  setTxt("acRpi", r.active?fR(r.rpi):"·");
  setTxt("acRo","∞");
}

/* ===== region pill ===== */
function renderPill(r){
  const p=document.getElementById("regionPill");
  const regionText=document.getElementById("regionText");
  if(regionText) regionText.textContent=r.region;
  const cls=r.region==="DOYMA"?" sat":r.region==="KESİM"?" cut":"";
  if(p) p.className="region-pill"+cls;
  // durum çubuğu bölge göstergesi
  const sb=document.querySelector(".statusbar");
  if(sb){
    sb.className="statusbar"+cls;
    const mark=r.region==="DOYMA"?"▲":r.region==="KESİM"?"✕":"●";
    document.getElementById("sbRegion").textContent=mark+" "+r.region;
  }
}

/* ===== Dalga şekli grafiği ===== */
function plot(r){
  const cv=document.getElementById("plot"), ctx=cv.getContext("2d");
  const dpr=window.devicePixelRatio||1;
  const cssW=cv.clientWidth||760, cssH=cv.clientHeight||210;
  cv.width=Math.round(cssW*dpr); cv.height=Math.round(cssH*dpr);
  ctx.setTransform(dpr,0,0,dpr,0,0);
  const W=cssW,H=cssH, padL=44,padR=14,padT=14,padB=26;
  const x0=padL,x1=W-padR,y0=padT,y1=H-padB, mid=(y0+y1)/2;
  ctx.clearRect(0,0,W,H);
  ctx.fillStyle="#060606"; ctx.fillRect(0,0,W,H);
  const warn=document.getElementById("clipWarn");
  const FM="'IBM Plex Mono'";

  if(!r.active){
    ctx.fillStyle="#9aa0a6"; ctx.font="13px "+FM; ctx.textAlign="center";
    ctx.fillText("Transistör aktif bölgede değil · dalga biçimi yok.", W/2, mid);
    warn.textContent=""; return;
  }

  // dalga genliği
  const viA=r.viPeak, voA=Math.abs(r.voPeak), head=r.headroom;
  const clipped = voA>head;
  const scaleMax = Math.max(voA, head, viA)*1.18;
  const vy = v => mid - (v/scaleMax)*(mid-y0);
  const timeWindow = 0.005;
  const omega = TAU*r.freq;

  // beyaz noktalı ızgara
  ctx.save();
  ctx.setLineDash([1,3]); ctx.lineWidth=1; ctx.strokeStyle="rgba(255,255,255,0.18)";
  for(let i=0;i<=8;i++){const x=x0+(x1-x0)*i/8;ctx.beginPath();ctx.moveTo(x,y0);ctx.lineTo(x,y1);ctx.stroke();}
  for(let i=0;i<=4;i++){const y=y0+(y1-y0)*i/4;ctx.beginPath();ctx.moveTo(x0,y);ctx.lineTo(x1,y);ctx.stroke();}
  ctx.restore();
  // çerçeve + sıfır ekseni
  ctx.strokeStyle="rgba(255,255,255,0.32)";ctx.lineWidth=1;ctx.strokeRect(x0,y0,x1-x0,y1-y0);
  ctx.strokeStyle="rgba(255,255,255,0.45)";ctx.beginPath();ctx.moveTo(x0,mid);ctx.lineTo(x1,mid);ctx.stroke();

  // ekran etiketleri
  ctx.fillStyle="#b9bfc6"; ctx.font="10px "+FM; ctx.textAlign="right";
  ctx.fillText("+"+scaleMax.toFixed(2)+"V",x0-6,y0+8);
  ctx.fillText("0",x0-6,mid+3);
  ctx.fillText("-"+scaleMax.toFixed(2)+"V",x0-6,y1);
  ctx.textAlign="left"; ctx.fillStyle="#e02424"; ctx.fillText("v_o", x0+5, y0+12);
  ctx.textAlign="right"; ctx.fillStyle="#2f9bef"; ctx.fillText("v_i", x1-5, y0+12);
  const freqLabel=r.freq>=1000?(r.freq/1000).toFixed(r.freq%1000?1:0)+" kHz":r.freq+" Hz";
  ctx.textAlign="left"; ctx.fillStyle="#9aa0a6"; ctx.fillText("0 ms",x0,H-7);
  ctx.textAlign="center"; ctx.fillText(freqLabel,(x0+x1)/2,H-7);
  ctx.textAlign="right"; ctx.fillText((timeWindow*1000).toFixed(1)+" ms",x1,H-7);

  const N=900;
  function sample(phasor,time){
    return phasor.re*Math.cos(omega*time) - phasor.im*Math.sin(omega*time);
  }
  function draw(phasor,color,clamp){
    ctx.strokeStyle=color; ctx.lineWidth=1.8; ctx.lineJoin="round";
    ctx.shadowColor=color; ctx.shadowBlur=6; ctx.beginPath();
    for(let i=0;i<=N;i++){
      const t=i/N;
      let v=sample(phasor,t*timeWindow);
      if(clamp!=null){ if(v>clamp)v=clamp; if(v<-clamp)v=-clamp; }
      const X=x0+(x1-x0)*t, Y=vy(v);
      i?ctx.lineTo(X,Y):ctx.moveTo(X,Y);
    }
    ctx.stroke(); ctx.shadowBlur=0;
  }
  draw(r.viPhasor,"#2f9bef",null);
  draw(r.voPhasor,"#e02424",clipped?head:null);

  warn.textContent = clipped
    ? `Çıkış kırpılıyor: |v_o| tepe ${(voA*1000).toFixed(0)} mV > salınım sınırı ${(head*1000).toFixed(0)} mV. Giriş genliğini düşür.`
    : "";
}

let lastResult=null;
function update(){
  const r=solve();
  lastResult=r;
  renderPill(r); renderResults(r); updateLabels(r); plot(r);
}
let _rT;
window.addEventListener("resize",()=>{ clearTimeout(_rT); _rT=setTimeout(()=>{ if(lastResult) plot(lastResult); },120); });

/* ===== Başlatma ===== */
document.getElementById("schAmp").innerHTML=drawAmp();
document.getElementById("schDC").innerHTML=drawDC();
document.getElementById("schAC").innerHTML=drawAC();
buildControls();
update();

/* ===== EEschema durum çubuğu: imleç koordinatı + birim ===== */
(function statusbar(){
  const A4_W_MM=297, A4_H_MM=210, MM_PER_MIL=0.0254;
  let unit="mm";             
  const anchor={x:0,y:0};
  const $=id=>document.getElementById(id);

  function conv(mm){ return unit==="mm" ? mm : mm/MM_PER_MIL; }
  function fmt(mm){ return conv(mm).toFixed(unit==="mm"?3:1); }

  function applyUnit(){
    $("sbUnit").textContent=unit;
    const unitsToggle=$("unitsToggle");
    if(unitsToggle) unitsToggle.textContent=unit;
    $("sbGrid").textContent = unit==="mm" ? "1.270 mm" : "50.0 mil";
  }
  const unitsToggle=$("unitsToggle");
  if(unitsToggle){
    unitsToggle.addEventListener("click",()=>{
      unit = unit==="mm" ? "mil" : "mm"; applyUnit();
    });
  }

  const sheet=document.getElementById("sheet");
  function pos(e){
    const r=sheet.getBoundingClientRect();
    const fx=Math.min(1,Math.max(0,(e.clientX-r.left)/r.width));
    const fy=Math.min(1,Math.max(0,(e.clientY-r.top)/r.height));
    return {x:fx*A4_W_MM, y:fy*A4_H_MM};
  }
  sheet.addEventListener("mousemove",e=>{
    const p=pos(e);
    const dx=p.x-anchor.x, dy=p.y-anchor.y;
    $("sbX").textContent=fmt(p.x);
    $("sbY").textContent=fmt(p.y);
    $("sbDX").textContent=fmt(dx);
    $("sbDY").textContent=fmt(dy);
    $("sbDist").textContent=fmt(Math.hypot(dx,dy));
  });
  sheet.addEventListener("mousedown",e=>{ Object.assign(anchor,pos(e)); });

  applyUnit();
})();

/* ===== parametre diyalog penceresi ===== */
(function paramDialog(){
  const dlg=document.getElementById("paramDialog");
  const open=document.getElementById("openParams");
  const closeBtns=[document.getElementById("dlgClose"),document.getElementById("dlgOk")];
  const minBtn=document.querySelector(".dlg-min");
  const reset=document.getElementById("dlgReset");
  const bar=document.getElementById("dlgTitlebar");

  const DEFAULTS={}; PARAMS.forEach(p=>DEFAULTS[p.id]=p.val);

  const show=()=>dlg.hidden=false;
  const hide=()=>dlg.hidden=true;
  open.addEventListener("click",show);
  closeBtns.forEach(b=>b.addEventListener("click",hide));
  if(minBtn) minBtn.addEventListener("click",hide);
  document.querySelectorAll(".sheet-canvas").forEach(c=>c.addEventListener("dblclick",show));

  reset.addEventListener("click",()=>{
    PARAMS.forEach(p=>{
      S[p.id]=DEFAULTS[p.id];
      const rng=document.getElementById("rng_"+p.id);
      const val=document.getElementById("val_"+p.id);
      if(rng) rng.value=DEFAULTS[p.id];
      if(val) val.textContent=fmtParam(p,DEFAULTS[p.id]);
    });
    update();
  });

  // başlık çubuğundan sürükle
  let drag=null;
  bar.addEventListener("mousedown",e=>{
    if(e.target.closest(".dlg-winbtns")) return;
    const r=dlg.getBoundingClientRect();
    drag={dx:e.clientX-r.left,dy:e.clientY-r.top};
    e.preventDefault();
  });
  window.addEventListener("mousemove",e=>{
    if(!drag) return;
    let x=e.clientX-drag.dx, y=e.clientY-drag.dy;
    x=Math.max(0,Math.min(window.innerWidth-dlg.offsetWidth,x));
    y=Math.max(0,Math.min(window.innerHeight-dlg.offsetHeight,y));
    dlg.style.left=x+"px"; dlg.style.top=y+"px";
  });
  window.addEventListener("mouseup",()=>{ drag=null; });
})();
