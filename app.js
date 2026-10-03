(()=>{
  const $ = id => document.getElementById(id);
  let comparison = null;
  let playTimer = null;
  let engineReady = false;
  let workerEnv = null;
  let requestCounter = 0;
  const pending = new Map();
  const worker = new Worker('worker.js');

  function setStatus(message, isError=false) {
    $('status').textContent = message;
    $('status').classList.toggle('error', isError);
  }

  worker.onmessage = event => {
    const msg = event.data || {};
    if (msg.type === 'status') {
      if (!engineReady) setStatus(msg.message || 'Loading…');
      return;
    }
    if (msg.type === 'ready') {
      engineReady = true;
      workerEnv = msg.environment || {};
      $('runButton').disabled = false;
      const bits = [
        `Pyodide ${workerEnv.pyodide || ''}`,
        `Python ${workerEnv.python || ''}`,
        `NumPy ${workerEnv.numpy || ''}`,
        `SciPy ${workerEnv.scipy || ''}`
      ].filter(Boolean);
      setStatus(`Ready — browser engine loaded (${bits.join(' · ')}).`);
      return;
    }
    if (msg.type === 'fatal') {
      engineReady = false;
      $('runButton').disabled = true;
      setStatus(`Browser Python engine failed to load: ${msg.error || 'unknown error'}`, true);
      return;
    }
    if (msg.type === 'result' || msg.type === 'error') {
      const p = pending.get(msg.id);
      if (!p) return;
      pending.delete(msg.id);
      if (msg.type === 'result') p.resolve(msg.result);
      else p.reject(new Error(msg.error || 'Simulation failed'));
    }
  };
  worker.onerror = event => {
    engineReady = false;
    $('runButton').disabled = true;
    setStatus(`Worker error: ${event.message || 'unknown error'}`, true);
  };

  const cfg = window.SIMULATOR_CONFIG || {};
  const links = [
    ['scientificCodeLink', cfg.scientificCodeUrl],
    ['articleLink', cfg.articleUrl],
    ['sourceRepositoryLink', cfg.sourceRepositoryUrl]
  ];
  const shown = [];
  for (const [id,url] of links) {
    if (url) {
      const el = $(id);
      el.href = url;
      el.classList.remove('hidden');
      shown.push(id);
    }
  }
  if (shown.length) {
    $('resourceLinks').classList.remove('hidden');
    if (shown.length > 1) $('linkSep1').classList.remove('hidden');
    if (shown.length > 2) $('linkSep2').classList.remove('hidden');
  }

  function compareInBrowser(payload) {
    return new Promise((resolve,reject)=>{
      if (!engineReady) {
        reject(new Error('The browser Python engine is still loading.'));
        return;
      }
      const id = ++requestCounter;
      pending.set(id,{resolve,reject});
      worker.postMessage({id,action:'compare',payload});
    });
  }

  const modelInfo = {
    asexual_deterministic: {
      name: 'Asexual — deterministic', c: 0.2,
      text: `<strong>Asexual deterministic.</strong>
        <div class="equation">
          ∂<sub>t</sub>n(t,x,y) − (σ²/2)∂<sub>xx</sub>n(t,x,y) − (μ²/2)∂<sub>yy</sub>n(t,x,y)
          = [r<sub>max</sub> − (1/(2V<sub>s</sub>))(y − b(x−ct))² − (1/k)∫<sub>ℝ</sub>n(t,x,y′)dy′] n(t,x,y).
        </div>`
    },
    asexual_stochastic: {
      name: 'Asexual — stochastic', c: 0.2,
      text: `<strong>Asexual stochastic.</strong> An individual at position <span class="math">x∈ℝ</span> and phenotype <span class="math">y∈ℝ</span> undergoes:
        <ul>
          <li>death at rate <span class="math">(y−b(x−ct))²/(2V<sub>s</sub>) + I(t,x)/(kK)</span>;</li>
          <li>reproduction at rate <span class="math">r<sub>max</sub></span>;</li>
          <li>phenotypic change according to Brownian motion with diffusion parameter <span class="math">μ²</span>;</li>
          <li>spatial dispersal according to Brownian motion with diffusion parameter <span class="math">σ²</span>.</li>
        </ul>`
    },
    sexual: {
      name: 'Sexual', c: 0.2,
      text: `<strong>Sexual model.</strong>
        <div class="equation">
          ∂<sub>t</sub>n(t,x,y) − (σ²/2)Δ<sub>x</sub>n(t,x,y)<br>
          = [−(η−r<sub>max</sub>) − (1/(2V<sub>s</sub>))(y−b(x−ct))² − (1/k)∫<sub>ℝ</sub>n(t,x,w)dw] n(t,x,y)<br>
          &nbsp;&nbsp;+ η ∬<sub>ℝ²</sub> Γ<sub>V<sub>LE</sub>/2</sub>(y−(y<sub>*</sub>+y′<sub>*</sub>)/2)
          · n(t,x,y<sub>*</sub>)n(t,x,y′<sub>*</sub>) / ∫<sub>ℝ</sub>n(t,x,w)dw · dy<sub>*</sub>dy′<sub>*</sub>.<br>
          Γ<sub>V<sub>LE</sub>/2</sub>(y) = [1/√(πV<sub>LE</sub>)] exp(−y²/V<sub>LE</sub>).
        </div>`
    }
  };

  function updateModel(side, resetC=false) {
    const model = $(`model${side}`).value;
    const info = modelInfo[model];
    $(`description${side}`).innerHTML = info.text;
    if (resetC) $(`c${side}`).value = info.c;

    const specific = $(`specific${side}`);
    const specificName = $(`specificName${side}`);
    if (model === 'sexual') {
      specificName.textContent = 'VLE';
      specific.value = '0.3';
      specific.disabled = false;
      specific.min = '0.001'; specific.max = '5'; specific.step = '0.001';
    } else {
      specificName.textContent = 'μ²';
      specific.disabled = false;
      specific.min = '0.001'; specific.max = '10'; specific.step = '0.01';
      if (resetC) specific.value = '0.64';
    }

    const isStochastic = model === 'asexual_stochastic';
    $(`kLabel${side}`).classList.toggle('hidden', !isStochastic);
    $(`seedLabel${side}`).classList.toggle('hidden', !isStochastic);
  }

  $('modelA').addEventListener('change', () => updateModel('A', true));
  $('modelB').addEventListener('change', () => updateModel('B', true));
  updateModel('A');
  updateModel('B');

  function stopPlaying() {
    if (playTimer !== null) clearInterval(playTimer);
    playTimer = null;
    $('playButton').textContent = '▶ Play';
  }

  function fmt(value, digits=4) {
    if (value === null || value === undefined || !Number.isFinite(Number(value))) return '—';
    return Number(value).toPrecision(digits);
  }

  function gridLabel(p) {
    if (p.nx) return `${p.nx} × ${p.ny}`;
    if (p.m) return `${p.m} × ${p.n}`;
    return '—';
  }

  function renderMetrics(side, data) {
    const items = [
      ['c', fmt(data.c)],
      ['Right-front speed', fmt(data.speed_right)],
    ];
    if (data.speed_left !== null && data.speed_left !== undefined) items.push(['Left-front speed', fmt(data.speed_left)]);
    const p = data.parameters || {};
    if (data.model === 'sexual') {
      items.push(['Vs', fmt(p.Vs)]);
      items.push(['b', fmt(p.b)]);
      items.push(['σ', fmt(p.sigma)]);
      items.push(['VLE', fmt(p.VLE)]);
    } else {
      items.push(['Vs', fmt(p.Vs)]);
      items.push(['b', fmt(p.b ?? p.B)]);
      items.push(['σ', fmt(p.sigmax)]);
      items.push(['μ²', fmt(p.sigmay)]);
    }
    items.push(['rmax', fmt(p.rmax)]);
    items.push(['k', fmt(p.k)]);
    if (data.model === 'asexual_stochastic') {
      items.push(['K', fmt(p.K)]);
      items.push(['seed', fmt(p.seed)]);
    }
    items.push(['x₀', fmt(p.initial_x_center)]);
    items.push(['Varₓ(0)', fmt(p.initial_x_variance)]);
    items.push(['Grid', gridLabel(p)]);
    $(`metrics${side}`).innerHTML = items.map(([label,value]) =>
      `<div class="metric"><span>${label}</span><strong>${value}</strong></div>`
    ).join('');
    $(`title${side}`).textContent = `Simulation ${side} — ${modelInfo[data.model].name}`;
    $(`subtitle${side}`).textContent = `c = ${fmt(data.c)}`;
  }

  async function runBoth() {
    stopPlaying();
    $('runButton').disabled = true;
    $('results').classList.add('hidden');
    setStatus('Running both simulations in your browser…');

    const payload = {
      model_a: $('modelA').value,
      c_a: Number($('cA').value),
      model_b: $('modelB').value,
      c_b: Number($('cB').value),
      tfinal: Number($('tfinal').value),
      preset: $('preset').value,
      snapshot_count: Number($('snapshots').value),
      K_a: Number($('KValueA').value),
      K_b: Number($('KValueB').value),
      seed_a: Number($('seedValueA').value),
      seed_b: Number($('seedValueB').value),
      Vs: Number($('VsValue').value),
      rmax: Number($('rmaxValue').value),
      competition_k: Number($('competitionKValue').value),
      b_a: Number($('bA').value),
      b_b: Number($('bB').value),
      sigmax: Number($('sigmaxValue').value),
      sigmay_a: $('modelA').value === 'sexual' ? 0.64 : Number($('specificA').value),
      sigmay_b: $('modelB').value === 'sexual' ? 0.64 : Number($('specificB').value),
      VLE_a: $('modelA').value === 'sexual' ? Number($('specificA').value) : 0.3,
      VLE_b: $('modelB').value === 'sexual' ? Number($('specificB').value) : 0.3
    };

    try {
      const data = await compareInBrowser(payload);
      comparison = data;
      renderMetrics('A', data.simulation_a);
      renderMetrics('B', data.simulation_b);
      $('results').classList.remove('hidden');
      const count = Math.min(data.simulation_a.snapshots.length, data.simulation_b.snapshots.length);
      $('timeSlider').max = Math.max(0, count - 1);
      $('timeSlider').value = 0;
      drawPair(0);
      setStatus(`Done. ${count} synchronized snapshot pairs available. Everything was computed locally in this browser tab.`);
    } catch (err) {
      setStatus(err.message || String(err), true);
    } finally {
      if (engineReady) $('runButton').disabled = false;
    }
  }

  $('runButton').addEventListener('click', runBoth);
  $('timeSlider').addEventListener('input', () => { stopPlaying(); drawPair(Number($('timeSlider').value)); });
  $('playButton').addEventListener('click', () => {
    if (!comparison) return;
    if (playTimer !== null) { stopPlaying(); return; }
    $('playButton').textContent = '❚❚ Pause';
    playTimer = setInterval(() => {
      let i = Number($('timeSlider').value) + 1;
      const max = Number($('timeSlider').max);
      if (i > max) i = 0;
      $('timeSlider').value = i;
      drawPair(i);
    }, 600);
  });

  function canvasSetup(id) {
    const canvas = $(id), ctx = canvas.getContext('2d');
    ctx.clearRect(0,0,canvas.width,canvas.height);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0,0,canvas.width,canvas.height);
    return [canvas,ctx];
  }

  function palette(u) {
    u = Math.max(0, Math.min(1, u));
    const r = Math.round(255 * Math.min(1, 1.65*u));
    const g = Math.round(255 * Math.max(0, Math.min(1, 1.7*u - 0.22)));
    const b = Math.round(255 * Math.max(0, 1.15 - 1.4*u));
    return [r,g,b];
  }

  function tickText(value) {
    const v = Number(value), a = Math.abs(v);
    if (a >= 1000 || (a > 0 && a < 0.01)) return v.toExponential(1);
    if (a >= 100) return v.toFixed(0);
    if (a >= 10) return v.toFixed(1).replace(/\.0$/,'');
    return v.toFixed(2).replace(/0+$/,'').replace(/\.$/,'');
  }

  function sharedDomain() {
    const d = comparison?.plot_domain || {};
    return {
      xmin: Number.isFinite(Number(d.x_min)) ? Number(d.x_min) : 0,
      xmax: Number.isFinite(Number(d.x_max)) ? Number(d.x_max) : 1,
      zmin: Number.isFinite(Number(d.trait_min)) ? Number(d.trait_min) : 0,
      zmax: Number.isFinite(Number(d.trait_max)) ? Number(d.trait_max) : 1,
    };
  }

  function drawRangeMarkers(ctx, snap, X, top, bottom, xmin, xmax) {
    ctx.save();
    ctx.strokeStyle='#6b7280';
    ctx.lineWidth=1.6;
    ctx.setLineDash([6,4]);
    for (const value of [snap.range_left, snap.range_right]) {
      if(value===null || value===undefined || !Number.isFinite(Number(value))) continue;
      const xv=Number(value);
      if(xv<xmin || xv>xmax) continue;
      const px=X(xv);
      ctx.beginPath();
      ctx.moveTo(px,top);
      ctx.lineTo(px,bottom);
      ctx.stroke();
    }
    ctx.restore();
  }

  function drawHeatmap(side, snap) {
    const [canvas,ctx] = canvasSetup(`heatmap${side}`);
    const W=canvas.width,H=canvas.height,M={left:72,right:18,top:16,bottom:58};
    const pw=W-M.left-M.right, ph=H-M.top-M.bottom;
    const pop=snap.population, nx=pop.length, nz=nx?pop[0].length:0;
    if (!nx || !nz) return;
    const {xmin,xmax,zmin,zmax} = sharedDomain();
    const X=x=>M.left+(x-xmin)/Math.max(xmax-xmin,1e-12)*pw;
    const Y=z=>M.top+(zmax-z)/Math.max(zmax-zmin,1e-12)*ph;

    let pmax=0;
    for(let i=0;i<nx;i++) for(let j=0;j<nz;j++) pmax=Math.max(pmax,Math.max(0,Number(pop[i][j])));
    pmax=Math.max(pmax,1e-30);
    const image=ctx.createImageData(nx,nz), contrast=30, den=Math.log1p(contrast);
    for(let i=0;i<nx;i++) for(let j=0;j<nz;j++) {
      const rel=Math.max(0,Number(pop[i][j]))/pmax;
      const u=rel<1e-3?0:Math.log1p(contrast*rel)/den;
      const [r,g,b]=palette(u), yy=nz-1-j, q=4*(yy*nx+i);
      image.data[q]=r; image.data[q+1]=g; image.data[q+2]=b; image.data[q+3]=255;
    }
    const off=document.createElement('canvas'); off.width=nx; off.height=nz;
    off.getContext('2d').putImageData(image,0,0);

    ctx.save(); ctx.beginPath(); ctx.rect(M.left,M.top,pw,ph); ctx.clip();
    ctx.imageSmoothingEnabled=true; ctx.imageSmoothingQuality='high';
    const x0=Number(snap.x[0]), x1=Number(snap.x.at(-1));
    const z0=Number(snap.trait[0]), z1=Number(snap.trait.at(-1));
    ctx.drawImage(off, X(x0), Y(z1), X(x1)-X(x0), Y(z0)-Y(z1));

    const density=snap.density || [];
    let densityMax=0; for(const v of density) densityMax=Math.max(densityMax,Math.max(0,Number(v)));
    const cutoff=densityMax*1e-3;
    function curve(values,dash,stroke,useDensityMask=false) {
      ctx.save(); ctx.strokeStyle=stroke; ctx.lineWidth=2; ctx.setLineDash(dash);
      ctx.beginPath(); let pen=false;
      for(let i=0;i<snap.x.length;i++) {
        if(useDensityMask && Number(density[i]||0)<=cutoff){pen=false;continue;}
        const xv=Number(snap.x[i]), zv=Number(values[i]);
        if(!Number.isFinite(xv)||!Number.isFinite(zv)){pen=false;continue;}
        const px=X(xv), py=Y(zv);
        if(!pen){ctx.moveTo(px,py);pen=true;} else ctx.lineTo(px,py);
      }
      ctx.stroke(); ctx.restore();
    }
    curve(snap.mean_trait_heatmap || snap.mean_trait,[],'#111827',true);
    curve(snap.optimum_heatmap || snap.optimum,[7,5],'#ffffff',false);
    drawRangeMarkers(ctx,snap,X,M.top,M.top+ph,xmin,xmax);
    ctx.restore();

    ctx.strokeStyle='#4b5563'; ctx.lineWidth=1; ctx.strokeRect(M.left,M.top,pw,ph);
    ctx.fillStyle='#374151'; ctx.font='11px system-ui';
    for(let q=0;q<5;q++) {
      const f=q/4, xv=xmin+f*(xmax-xmin), px=X(xv), zv=zmin+f*(zmax-zmin), py=Y(zv);
      ctx.beginPath(); ctx.moveTo(px,M.top+ph); ctx.lineTo(px,M.top+ph+5); ctx.stroke();
      ctx.textAlign='center'; ctx.textBaseline='top'; ctx.fillText(tickText(xv),px,M.top+ph+8);
      ctx.beginPath(); ctx.moveTo(M.left-5,py); ctx.lineTo(M.left,py); ctx.stroke();
      ctx.textAlign='right'; ctx.textBaseline='middle'; ctx.fillText(tickText(zv),M.left-8,py);
    }
    ctx.font='12px system-ui'; ctx.textAlign='center'; ctx.textBaseline='alphabetic';
    ctx.fillText('space x',M.left+pw/2,H-8);
    ctx.save(); ctx.translate(16,M.top+ph/2); ctx.rotate(-Math.PI/2); ctx.fillText('trait / phenotype',0,0); ctx.restore();
  }

  function drawAxesAndCurve(canvasId, snap, values, yLabel, {mask=false, nonnegative=false, overlayValues=null, overlayMask=false}={}) {
    const [canvas,ctx] = canvasSetup(canvasId);
    const W=canvas.width,H=canvas.height,M={left:70,right:16,top:14,bottom:50};
    const pw=W-M.left-M.right, ph=H-M.top-M.bottom;
    const {xmin,xmax}=sharedDomain();
    const x=snap.x.map(Number), y=values.map(Number), dens=(snap.density||[]).map(Number);
    const overlay=overlayValues ? overlayValues.map(Number) : null;
    let dmax=0; for(const v of dens) dmax=Math.max(dmax,Math.max(0,v));
    const cutoff=dmax*1e-3;
    const valid=[];
    for(let i=0;i<y.length;i++) {
      if(Number.isFinite(x[i]) && Number.isFinite(y[i]) && !(mask && Number(dens[i]||0)<=cutoff)) valid.push(y[i]);
      if(overlay && Number.isFinite(x[i]) && Number.isFinite(overlay[i]) && !(overlayMask && Number(dens[i]||0)<=cutoff)) valid.push(overlay[i]);
    }
    if(!valid.length) valid.push(0,1);
    let ymin=nonnegative?0:Math.min(...valid), ymax=Math.max(...valid);
    if(nonnegative) ymax=Math.max(ymax,1e-8);
    if(Math.abs(ymax-ymin)<1e-10) {
      const pad=Math.max(Math.abs(ymax)*0.1,1);
      ymin=nonnegative?0:ymin-pad; ymax+=pad;
    } else {
      const pad=.07*(ymax-ymin);
      if(!nonnegative) ymin-=pad;
      ymax+=pad;
    }
    const X=v=>M.left+(Number(v)-xmin)/Math.max(xmax-xmin,1e-12)*pw;
    const Y=v=>M.top+ph-(Number(v)-ymin)/Math.max(ymax-ymin,1e-12)*ph;

    function drawCurve(curveValues, stroke, dash, useMask) {
      ctx.save();
      ctx.strokeStyle=stroke; ctx.lineWidth=2.2; ctx.setLineDash(dash);
      ctx.beginPath(); let pen=false;
      for(let i=0;i<curveValues.length;i++) {
        if(x[i]<xmin||x[i]>xmax || !Number.isFinite(curveValues[i]) || (useMask && Number(dens[i]||0)<=cutoff)) {pen=false;continue;}
        const px=X(x[i]),py=Y(curveValues[i]);
        if(!pen){ctx.moveTo(px,py);pen=true;} else ctx.lineTo(px,py);
      }
      ctx.stroke(); ctx.restore();
    }

    ctx.save(); ctx.beginPath(); ctx.rect(M.left,M.top,pw,ph); ctx.clip();
    drawCurve(y,'#315fbd',[],mask);
    if(overlay) drawCurve(overlay,'#c2410c',[8,5],overlayMask);
    drawRangeMarkers(ctx,snap,X,M.top,M.top+ph,xmin,xmax);
    ctx.restore();

    ctx.strokeStyle='#4b5563'; ctx.lineWidth=1;
    ctx.beginPath(); ctx.moveTo(M.left,M.top+ph); ctx.lineTo(M.left+pw,M.top+ph); ctx.moveTo(M.left,M.top); ctx.lineTo(M.left,M.top+ph); ctx.stroke();
    ctx.fillStyle='#374151'; ctx.font='11px system-ui';
    for(let q=0;q<5;q++) {
      const f=q/4, xv=xmin+f*(xmax-xmin), px=X(xv), yv=ymin+f*(ymax-ymin), py=Y(yv);
      ctx.beginPath();ctx.moveTo(px,M.top+ph);ctx.lineTo(px,M.top+ph+5);ctx.stroke();
      ctx.textAlign='center';ctx.textBaseline='top';ctx.fillText(tickText(xv),px,M.top+ph+8);
      ctx.beginPath();ctx.moveTo(M.left-5,py);ctx.lineTo(M.left,py);ctx.stroke();
      ctx.textAlign='right';ctx.textBaseline='middle';ctx.fillText(tickText(yv),M.left-8,py);
    }
    ctx.font='12px system-ui';ctx.textAlign='center';ctx.textBaseline='alphabetic';ctx.fillText('space x',M.left+pw/2,H-7);
    ctx.save();ctx.translate(16,M.top+ph/2);ctx.rotate(-Math.PI/2);ctx.fillText(yLabel,0,0);ctx.restore();
  }

  function drawProfile(side,snap) { drawAxesAndCurve(`profile${side}`,snap,snap.density,'population density',{nonnegative:true}); }
  function drawMeanTrait(side,snap) {
    drawAxesAndCurve(`meanTrait${side}`,snap,snap.mean_trait,'phenotypic trait',{mask:true,overlayValues:snap.optimum,overlayMask:true});
  }
  function drawVariance(side,snap) { drawAxesAndCurve(`variance${side}`,snap,snap.trait_variance,'variance V(x)',{mask:true,nonnegative:true}); }

  function drawPair(index) {
    if(!comparison) return;
    const a=comparison.simulation_a.snapshots[index], b=comparison.simulation_b.snapshots[index];
    if(!a||!b) return;
    $('timeA').textContent=`A: t = ${Number(a.time).toFixed(2)}`;
    $('timeB').textContent=`B: t = ${Number(b.time).toFixed(2)}`;
    drawHeatmap('A',a); drawProfile('A',a); drawMeanTrait('A',a); drawVariance('A',a);
    drawHeatmap('B',b); drawProfile('B',b); drawMeanTrait('B',b); drawVariance('B',b);
  }
})();
