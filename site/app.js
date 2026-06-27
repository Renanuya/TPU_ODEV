"use strict";

/* ===== constants ===== */
const VBE = 0.7, VT = 0.026, VCESAT = 0.2;

/* ===== slider params ===== */
const PARAMS = [
  { id:"VCC",  name:"V<sub>CC</sub>", min:5,    max:18,   step:0.5, val:12,   unit:"V",  fmt:"V" },
  { id:"R1",   name:"R₁",            min:10000,max:100000,step:1000,val:47000,unit:"Ω", fmt:"R" },
  { id:"R2",   name:"R₂",            min:2000, max:47000,step:500, val:10000,unit:"Ω", fmt:"R" },
  { id:"RC",   name:"R<sub>C</sub>", min:500,  max:10000,step:100, val:2200, unit:"Ω", fmt:"R" },
  { id:"RE",   name:"R<sub>E</sub>", min:100,  max:5000, step:100, val:1000, unit:"Ω", fmt:"R" },
  { id:"RL",   name:"R<sub>L</sub>", min:1000, max:100000,step:1000,val:10000,unit:"Ω", fmt:"R" },
  { id:"Rs",   name:"R<sub>s</sub>", min:50,   max:2000, step:50,  val:600,  unit:"Ω", fmt:"R" },
  { id:"beta", name:"β (h<sub>FE</sub>)", min:50, max:400, step:10, val:150, unit:"", fmt:"int" },
  { id:"vs",   name:"v<sub>s</sub> genlik", min:1, max:50, step:1, val:10, unit:"mV", fmt:"mV" },
  { id:"freq", name:"frekans", min:100, max:5000, step:100, val:1000, unit:"Hz", fmt:"Hz" },
];
const S = {};
PARAMS.forEach(p => S[p.id] = p.val);

/* ===== formatting ===== */
const par = (a,b) => (a*b)/(a+b);
function fR(x){
  if(x>=1000){ const k=x/1000; const s=Number.isInteger(k)?String(k):k.toFixed(2).replace(/\.?0+$/,""); return s+" kΩ"; }
  return Math.round(x)+" Ω";
}
function fV(x){ return x.toFixed(2)+" V"; }
function fI(x){ // amps
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
  if(p.fmt==="Hz") return v>=1000?(v/1000)+" kHz":v+" Hz";
  return v;
}

/* ===== core physics ===== */
function solve(){
  const {VCC,R1,R2,RC,RE,RL,Rs,beta} = S;
  const vsAmp = S.vs/1000; // V

  // DC — voltage divider Thevenin
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
      VE=IE*RE; VB=VE+VBE; VCB=VCEq-VBE; IB=icActive/beta; // keep IB as base drive
    } else {
      region="AKTİF"; IC=icActive; IE=(beta+1)*IB; VE=IE*RE; VB=VE+VBE;
      VCEq=VCC-IC*RC-IE*RE; VCB=VCEq-VBE;
    }
  }

  const active = region==="AKTİF";
  const ICQ = IC;

  // AC — hybrid-pi (ro -> inf)
  let gm=0, rpi=0, Rg=0, Kv=0, Kv0=0, Ro=RC, Ki=0, Kvs=0, viPeak=0, voPeak=0;
  if(active && ICQ>0){
    gm = ICQ/VT;
    rpi = beta/gm;
    Rg = par(par(R1,R2),rpi);
    const RCL = par(RC,RL);
    Kv = -gm*RCL;
    Kv0 = -gm*RC;
    Ro = RC;
    Ki = Kv*Rg/RL;
    Kvs = Kv*Rg/(Rg+Rs);
    viPeak = vsAmp*Rg/(Rg+Rs);
    voPeak = Kv*viPeak;
  }

  // output headroom for clipping
  const headroom = active ? Math.min(VCEq-VCESAT, ICQ*par(RC,RL)) : 0;

  return {VTh,RTh,IB,IC,IE,ICQ,VCEq,VE,VB,VCB,ICsat,region,active,
          gm,rpi,Rg,Kv,Kv0,Ro,Ki,Kvs,viPeak,voPeak,headroom,vsAmp};
}

/* ===== SVG schematic primitives ===== */
const SB = "var(--sch-bg)";
function wire(x1,y1,x2,y2){ return `<path class="wire" d="M${x1} ${y1} L${x2} ${y2}"/>`; }
function poly(pts){ return `<path class="wire" d="M${pts}"/>`; }
function node(x,y){ return `<circle class="node" cx="${x}" cy="${y}" r="3"/>`; }
function gnd(x,y){ return `<g stroke="var(--accent)" stroke-width="2" fill="none">
  <path d="M${x} ${y-10} L${x} ${y}"/><path d="M${x-11} ${y} L${x+11} ${y}"/>
  <path d="M${x-7} ${y+4} L${x+7} ${y+4}"/><path d="M${x-3} ${y+8} L${x+3} ${y+8}"/></g>`; }
// resistor = box masking the wire (wire drawn separately/under)
function resV(x,y1,y2,label,vid){
  const m=(y1+y2)/2;
  return `<rect x="${x-9}" y="${m-22}" width="18" height="44" rx="2" fill="${SB}" class="comp"/>
    <text class="lbl" x="${x+15}" y="${m-4}">${label}</text>
    <text class="lbl-v" id="${vid}" x="${x+15}" y="${m+10}"></text>`;
}
function resH(x1,x2,y,label,vid){
  const m=(x1+x2)/2;
  return `<rect x="${m-22}" y="${y-9}" width="44" height="18" rx="2" fill="${SB}" class="comp"/>
    <text class="lbl" x="${m}" y="${y-15}" text-anchor="middle">${label}</text>
    <text class="lbl-v" id="${vid}" x="${m}" y="${y+24}" text-anchor="middle"></text>`;
}
function capV(x,y1,y2,label,vid){ // wire passes; gap with plates
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
function source(x,y,label,vid){ // AC source circle
  return `<circle cx="${x}" cy="${y}" r="16" fill="${SB}" class="comp"/>
    <path class="comp" d="M${x-8} ${y} q4 -8 8 0 q4 8 8 0" fill="none"/>
    <text class="lbl" x="${x-22}" y="${y+4}" text-anchor="end">${label}</text>${vid?`<text class="lbl-v" id="${vid}" x="${x-22}" y="${y+18}" text-anchor="end"></text>`:""}`;
}
function transistor(x,y){ // NPN, base on left
  return `<circle cx="${x}" cy="${y}" r="22" fill="none" class="comp" opacity="0.55"/>
    <path class="comp" d="M${x-9} ${y-13} L${x-9} ${y+13}" stroke-width="2.5"/>
    <path class="comp" d="M${x-9} ${y-6} L${x+7} ${y-17}"/>
    <path class="comp" d="M${x-9} ${y+6} L${x+7} ${y+17}"/>
    <path class="comp-fill" d="M${x+1} ${y+8} L${x+8} ${y+18} L${x-2} ${y+15} Z"/>`;
}
function arrow(x,y,dir,color){ // dir: 'd' down,'u' up,'r' right
  color=color||"var(--accent-2)";
  let p;
  if(dir==="d") p=`${x-4},${y-5} ${x+4},${y-5} ${x},${y+3}`;
  else if(dir==="u") p=`${x-4},${y+5} ${x+4},${y+5} ${x},${y-3}`;
  else p=`${x-5},${y-4} ${x-5},${y+4} ${x+3},${y}`;
  return `<polygon points="${p}" fill="${color}"/>`;
}
function svg(vb,inner){ return `<svg viewBox="0 0 ${vb}" preserveAspectRatio="xMidYMid meet">${inner}</svg>`; }

/* ===== Schematic 1: full CE amplifier ===== */
function drawAmp(){
  let g="";
  // VCC rail
  g+=wire(70,28,360,28);
  g+=`<text class="rail-t" x="40" y="32">+Vcc</text>`;
  g+=`<text class="lbl-v" id="ampVcc" x="62" y="20"></text>`;
  // R1 divider
  g+=wire(160,28,160,150)+resV(160,40,116,"R₁","ampR1");
  // base node
  g+=node(160,150);
  // R2
  g+=wire(160,150,160,250)+resV(160,168,232,"R₂","ampR2")+gnd(160,262);
  // transistor
  g+=transistor(212,150);
  g+=wire(160,150,203,150); // base lead
  g+=poly(`212 144 L228 132 L232 122`); // collector lead
  g+=poly(`212 156 L228 168 L232 178`); // emitter lead
  // RC + collector
  g+=wire(232,28,232,122)+resV(232,40,96,"R<tspan dy='3'>C</tspan>","ampRC");
  g+=node(232,108);
  // output coupling
  g+=wire(232,108,290,108)+capH(290,322,108,"C<tspan dy='3'>C</tspan>");
  g+=wire(322,108,360,108)+node(360,108);
  g+=`<text class="lbl-v" x="366" y="104">v<tspan dy='3'>o</tspan></text>`;
  // RL
  g+=wire(360,108,360,150)+resV(360,150,214,"R<tspan dy='3'>L</tspan>","ampRL")+gnd(360,236);
  // emitter / RE
  g+=node(232,178);
  g+=wire(232,178,232,208)+resV(232,210,276,"R<tspan dy='3'>E</tspan>","ampRE")+gnd(232,300);
  // CE bypass
  g+=wire(232,208,300,208)+capV(300,208,260,"C<tspan dy='3'>E</tspan>")+gnd(300,272);
  // input source chain
  g+=source(40,150,"v<tspan dy='3'>s</tspan>");
  g+=wire(40,166,40,196)+gnd(40,196);
  g+=wire(56,150,74,150)+resH(74,110,150,"R<tspan dy='3'>s</tspan>","ampRs");
  g+=wire(110,150,118,150)+capH(118,150,150,"C<tspan dy='3'>B</tspan>");
  return svg("420 320",g);
}

/* ===== Schematic 2: DC equivalent ===== */
function drawDC(){
  let g="";
  g+=wire(70,28,300,28);
  g+=`<text class="rail-t" x="40" y="32">+Vcc</text>`;
  g+=wire(140,28,140,150)+resV(140,40,116,"R₁","dcR1");
  g+=node(140,150);
  g+=wire(140,150,140,250)+resV(140,168,232,"R₂","dcR2")+gnd(140,262);
  g+=transistor(192,150);
  g+=wire(140,150,183,150);
  g+=arrow(165,150,"r")+`<text class="lbl-v" id="dcIB" x="158" y="142" text-anchor="middle"></text><text class="lbl-m" x="150" y="166">I<tspan dy='3'>B</tspan></text>`;
  g+=poly(`192 144 L208 132 L212 122`);
  g+=poly(`192 156 L208 168 L212 178`);
  g+=wire(212,28,212,122)+resV(212,40,96,"R<tspan dy='3'>C</tspan>","dcRC");
  g+=arrow(212,112,"u","var(--accent-2)")+`<text class="lbl-v" id="dcIC" x="220" y="116"></text><text class="lbl-m" x="220" y="104">I<tspan dy='3'>C</tspan></text>`;
  // VCE bracket
  g+=`<path class="comp" d="M250 122 L258 122 M254 122 L254 178 M250 178 L258 178" opacity="0.7"/>`;
  g+=`<text class="lbl-v" id="dcVCE" x="262" y="146"></text><text class="lbl-m" x="262" y="158">V<tspan dy='3'>CE</tspan></text>`;
  g+=node(212,178);
  g+=wire(212,178,212,208)+resV(212,210,276,"R<tspan dy='3'>E</tspan>","dcRE")+gnd(212,300);
  g+=arrow(212,235,"d")+`<text class="lbl-v" id="dcIE" x="178" y="238" text-anchor="end"></text><text class="lbl-m" x="200" y="252" text-anchor="end">I<tspan dy='3'>E</tspan></text>`;
  return svg("360 320",g);
}

/* ===== Schematic 3: AC equivalent (hybrid-pi) ===== */
function drawAC(){
  let g="";
  // input
  g+=source(35,150,"v<tspan dy='3'>s</tspan>");
  g+=wire(35,166,35,196)+gnd(35,196);
  g+=wire(51,150,66,150)+resH(66,102,150,"R<tspan dy='3'>s</tspan>",null);
  g+=wire(102,150,140,150);
  // base node
  g+=node(140,150);
  g+=`<text class="lbl-m" x="140" y="142" text-anchor="middle">B</text>`;
  // RB = R1||R2
  g+=wire(140,150,140,176)+resV(140,178,234,"R<tspan dy='3'>B</tspan>","acRB")+gnd(140,266);
  // rpi
  g+=wire(140,150,195,150)+node(195,150);
  g+=wire(195,150,195,176)+resV(195,178,234,"r<tspan dy='3'>π</tspan>","acRpi")+gnd(195,266);
  // vbe arrow
  g+=`<path class="wire" d="M232 168 L232 232" stroke-dasharray="3 3" opacity="0.6"/>`;
  g+=arrow(232,170,"u")+`<text class="lbl-v" x="238" y="204">v<tspan dy='3'>be</tspan></text>`;
  // dashed divider
  g+=`<path d="M260 30 L260 270" stroke="var(--line)" stroke-width="1.5" stroke-dasharray="4 4"/>`;
  // collector rail
  g+=wire(300,120,425,120);
  g+=node(425,120);
  g+=`<text class="lbl-m" x="305" y="112">C</text>`;
  // dependent source (diamond) gm·vbe
  g+=wire(300,120,300,150);
  g+=`<path class="comp" d="M300 150 l15 16 l-15 16 l-15 -16 z" fill="${SB}"/>`;
  g+=arrow(300,158,"u")+`<text class="lbl" x="282" y="172" text-anchor="end">g<tspan dy='3'>m</tspan>v<tspan dy='3'>be</tspan></text>`;
  g+=wire(300,182,300,210)+gnd(300,210);
  // ro
  g+=wire(345,120,345,142)+resV(345,144,200,"r<tspan dy='3'>o</tspan>","acRo")+gnd(345,232);
  // RC
  g+=wire(385,120,385,142)+resV(385,144,200,"R<tspan dy='3'>C</tspan>",null)+gnd(385,232);
  // RL + vo
  g+=wire(425,120,425,142)+resV(425,144,200,"R<tspan dy='3'>L</tspan>",null)+gnd(425,232);
  g+=`<text class="lbl-v" x="425" y="112" text-anchor="middle">v<tspan dy='3'>o</tspan></text>`;
  return svg("470 290",g);
}

/* ===== build controls ===== */
function buildControls(){
  const c=document.getElementById("controls");
  c.innerHTML = PARAMS.map(p=>`
    <div class="ctrl">
      <span class="name">${p.name}</span>
      <span class="val" id="val_${p.id}">${fmtParam(p,p.val)}</span>
      <input type="range" id="rng_${p.id}" min="${p.min}" max="${p.max}" step="${p.step}" value="${p.val}">
    </div>`).join("");
  PARAMS.forEach(p=>{
    document.getElementById("rng_"+p.id).addEventListener("input",e=>{
      S[p.id]=parseFloat(e.target.value);
      document.getElementById("val_"+p.id).textContent=fmtParam(p,S[p.id]);
      update();
    });
  });
}

/* ===== metric card helper ===== */
function metric(k,v,invalid){ return `<div class="metric${invalid?" invalid":""}"><span class="k">${k}</span><span class="v">${v}</span></div>`; }

/* ===== render metrics ===== */
function renderDC(r){
  const el=document.getElementById("dcMetrics");
  el.innerHTML =
    metric("I<sub>B</sub>", fI(r.IB)) +
    metric("I<sub>C</sub> (I<sub>CQ</sub>)", fI(r.IC)) +
    metric("I<sub>E</sub>", fI(r.IE)) +
    metric("V<sub>CE</sub>", fV(r.VCEq)) +
    metric("V<sub>CB</sub>", fV(r.VCB)) +
    metric("V<sub>BE</sub>", fV(r.region==="KESİM"?0:VBE));
  const b=document.getElementById("regionBanner");
  let cls,msg;
  if(r.region==="AKTİF"){cls="active"; msg=`● AKTİF BÖLGE — yükselteç olarak çalışıyor. I<sub>C(sat)</sub>=${fI(r.ICsat)}`;}
  else if(r.region==="DOYMA"){cls="sat"; msg=`▲ DOYMA BÖLGESİ — V<sub>CE</sub>≤${VCESAT} V. AC kazançları geçersiz.`;}
  else {cls="cut"; msg="✕ KESİM BÖLGESİ — transistör iletmiyor. AC kazançları geçersiz.";}
  b.className="region-banner "+cls; b.innerHTML=msg;
}
function renderAC(r){
  const el=document.getElementById("acMetrics");
  const inv=!r.active;
  const g=x=> inv?"—":x;
  el.innerHTML =
    metric("K<sub>v</sub> (gerilim)", inv?"—":r.Kv.toFixed(1), inv) +
    metric("K<sub>v0</sub> (yüksüz)", inv?"—":r.Kv0.toFixed(1), inv) +
    metric("K<sub>i</sub> (akım)", inv?"—":r.Ki.toFixed(1), inv) +
    metric("K<sub>vg</sub> (kaynak→yük)", inv?"—":r.Kvs.toFixed(1), inv) +
    metric("R<sub>g</sub> (giriş dir.)", inv?"—":fR(r.Rg), inv) +
    metric("R<sub>o</sub> (çıkış dir.)", inv?"—":fR(r.Ro), inv) +
    metric("g<sub>m</sub>", inv?"—":(r.gm*1000).toFixed(2)+" mS", inv) +
    metric("r<sub>π</sub>", inv?"—":fR(r.rpi), inv) +
    metric("v<sub>o</sub> tepe", inv?"—":(Math.abs(r.voPeak)*1000).toFixed(1)+" mV", inv) +
    metric("v<sub>i</sub> tepe", inv?"—":(r.viPeak*1000).toFixed(1)+" mV", inv);
}

/* ===== update schematic value labels ===== */
function setTxt(id,v){ const e=document.getElementById(id); if(e) e.textContent=v; }
function updateLabels(r){
  setTxt("ampVcc",S.VCC.toFixed(1)+"V");
  setTxt("ampR1",fR(S.R1)); setTxt("ampR2",fR(S.R2)); setTxt("ampRC",fR(S.RC));
  setTxt("ampRE",fR(S.RE)); setTxt("ampRL",fR(S.RL)); setTxt("ampRs",fR(S.Rs));
  setTxt("dcR1",fR(S.R1)); setTxt("dcR2",fR(S.R2)); setTxt("dcRC",fR(S.RC)); setTxt("dcRE",fR(S.RE));
  setTxt("dcIB",fI(r.IB)); setTxt("dcIC",fI(r.IC)); setTxt("dcIE",fI(r.IE)); setTxt("dcVCE",fV(r.VCEq));
  setTxt("acRB", r.active?fR(par(S.R1,S.R2)):"—");
  setTxt("acRpi", r.active?fR(r.rpi):"—");
  setTxt("acRo","∞");
}

/* ===== region pill ===== */
function renderPill(r){
  const p=document.getElementById("regionPill");
  document.getElementById("regionText").textContent=r.region;
  p.className="region-pill"+(r.region==="DOYMA"?" sat":r.region==="KESİM"?" cut":"");
}

/* ===== waveform plot ===== */
function plot(r){
  const cv=document.getElementById("plot"), ctx=cv.getContext("2d");
  const W=cv.width,H=cv.height, padL=44,padR=14,padT=14,padB=26;
  const x0=padL,x1=W-padR,y0=padT,y1=H-padB, mid=(y0+y1)/2;
  ctx.clearRect(0,0,W,H);
  // bg
  ctx.fillStyle="#fbfaf5"; ctx.fillRect(0,0,W,H);
  const warn=document.getElementById("clipWarn");

  if(!r.active){
    ctx.fillStyle="#8a877d"; ctx.font="13px 'JetBrains Mono'"; ctx.textAlign="center";
    ctx.fillText("Transistör aktif bölgede değil — dalga biçimi yok.", W/2, mid);
    warn.textContent=""; return;
  }

  // amplitude scaling
  const viA=r.viPeak, voA=Math.abs(r.voPeak);
  const head=r.headroom;
  const clipped = voA>head;
  const scaleMax = Math.max(voA, head, viA)*1.18;
  const vy = v => mid - (v/scaleMax)*(mid-y0);

  // grid
  ctx.strokeStyle="#e3ded2"; ctx.lineWidth=1;
  for(let i=0;i<=8;i++){const x=x0+(x1-x0)*i/8;ctx.beginPath();ctx.moveTo(x,y0);ctx.lineTo(x,y1);ctx.stroke();}
  for(let i=0;i<=4;i++){const y=y0+(y1-y0)*i/4;ctx.beginPath();ctx.moveTo(x0,y);ctx.lineTo(x1,y);ctx.stroke();}
  // zero axis
  ctx.strokeStyle="#c3bdb0"; ctx.beginPath();ctx.moveTo(x0,mid);ctx.lineTo(x1,mid);ctx.stroke();
  // y labels
  ctx.fillStyle="#8a877d"; ctx.font="10px 'JetBrains Mono'"; ctx.textAlign="right";
  ctx.fillText("+"+scaleMax.toFixed(2)+"V",x0-6,y0+8);
  ctx.fillText("0",x0-6,mid+3);
  ctx.fillText("-"+scaleMax.toFixed(2)+"V",x0-6,y1);
  ctx.textAlign="center"; ctx.fillText("t →",(x0+x1)/2,H-7);

  const N=600, periods=2;
  // vo (red) — clamp to headroom if clipped
  function draw(amp,phase,color,clamp){
    ctx.strokeStyle=color; ctx.lineWidth=2; ctx.beginPath();
    for(let i=0;i<=N;i++){
      const t=i/N;
      let v=amp*Math.sin(2*Math.PI*periods*t+phase);
      if(clamp!=null){ if(v>clamp)v=clamp; if(v<-clamp)v=-clamp; }
      const X=x0+(x1-x0)*t, Y=vy(v);
      i?ctx.lineTo(X,Y):ctx.moveTo(X,Y);
    }
    ctx.stroke();
  }
  // vi (blue), vo inverted (phase pi), clamp vo to headroom
  draw(viA,0,"#1f5fd6",null);
  draw(voA,Math.PI,"#c0392b",clipped?head:null);

  warn.textContent = clipped
    ? `⚠ Çıkış kırpılıyor: |v_o| tepe ${(voA*1000).toFixed(0)} mV > salınım sınırı ${(head*1000).toFixed(0)} mV. Giriş genliğini düşür.`
    : "";
}

/* ===== master update ===== */
function update(){
  const r=solve();
  renderPill(r); renderDC(r); renderAC(r); updateLabels(r); plot(r);
}

/* ===== init ===== */
document.getElementById("schAmp").innerHTML=drawAmp();
document.getElementById("schDC").innerHTML=drawDC();
document.getElementById("schAC").innerHTML=drawAC();
buildControls();
update();
