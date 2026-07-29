import { project, divergence, toCentre } from './divb-project.mjs';

function rms(a){ let s=0; for (const v of a) s+=v*v; return Math.sqrt(s/a.length); }
function amax(a){ let m=0; for (const v of a) m=Math.max(m,Math.abs(v)); return m; }

const nx=152, ny=96, dx=1/ny, dy=1/ny;   // the engine's cells are square; dx=dy=1/h
const n=nx*ny;

// ---- test 1: a field that IS the curl of a known corner potential -------------
{
  const A=new Float64Array(n);
  for(let j=0;j<ny;j++) for(let i=0;i<nx;i++){
    let s=0;
    for(const [kx,ky,ph,amp] of [[1,2,0.3,1.0],[3,-1,1.1,0.6],[-2,4,2.2,0.35],[5,3,0.7,0.15]])
      s+=amp*Math.sin(2*Math.PI*(kx*i/nx+ky*j/ny)+ph);
    A[i+j*nx]=s;
  }
  const Aat=(i,j)=>A[((i%nx)+nx)%nx+(((j%ny)+ny)%ny)*nx];
  const bxf=new Float64Array(n), byf=new Float64Array(n);
  const B0=[0.31,0.63];
  for(let j=0;j<ny;j++) for(let i=0;i<nx;i++){
    bxf[i+j*nx]=(Aat(i,j+1)-Aat(i,j))/dy + B0[0];
    byf[i+j*nx]=-(Aat(i+1,j)-Aat(i,j))/dx + B0[1];
  }
  const c=toCentre(bxf,byf,nx,ny);
  const r=project(c.bx,c.by,nx,ny,dx,dy);
  const back=toCentre(r.bxf,r.byf,nx,ny);
  const ex=new Float64Array(n), ey=new Float64Array(n);
  for(let t=0;t<n;t++){ ex[t]=back.bx[t]-c.bx[t]; ey[t]=back.by[t]-c.by[t]; }
  const d=divergence(r.bxf,r.byf,nx,ny,dx,dy);
  const scale=Math.max(rms(c.bx),rms(c.by));
  console.log('T1 exactly-representable field');
  console.log('   round-trip cell-centred rel error  rms %s  max %s',
    (rms(ex)/scale).toExponential(2),(amax(ex)/scale).toExponential(2));
  console.log('   recovered mean field  %s %s   (true %s %s)',
    r.mean[0].toFixed(6), r.mean[1].toFixed(6), B0[0], B0[1]);
  console.log('   div(B_face) rms %s  max %s   [normalised by |B|/dx = %s]',
    (rms(d)).toExponential(2),(amax(d)).toExponential(2),(scale/dx).toExponential(2));
  console.log('   div normalised: rms %s  max %s',
    (rms(d)*dx/scale).toExponential(2),(amax(d)*dx/scale).toExponential(2));
  console.log('   T1 %s', (rms(ex)/scale<1e-10 && amax(d)*dx/scale<1e-12) ? 'PASS':'FAIL');
}

// ---- test 2: a field with real divergence in it (what Dedner actually hands us)
{
  let seed=12345; const rnd=()=>(seed=(seed*1103515245+12345)&0x7fffffff)/0x7fffffff-0.5;
  const bx=new Float64Array(n), by=new Float64Array(n);
  // smooth + a divergent component + grid-scale noise, mean field added
  for(let j=0;j<ny;j++) for(let i=0;i<nx;i++){
    const t=i+j*nx;
    bx[t]=0.30+0.5*Math.sin(2*Math.PI*(2*i/nx+1*j/ny))+0.10*Math.cos(2*Math.PI*(3*i/nx))+0.04*rnd();
    by[t]=0.63+0.5*Math.cos(2*Math.PI*(1*i/nx-2*j/ny))+0.10*Math.sin(2*Math.PI*(3*j/ny))+0.04*rnd();
  }
  const scale=Math.max(rms(bx),rms(by));
  // divergence of the naive face interpolation, for comparison
  const bxf0=new Float64Array(n), byf0=new Float64Array(n);
  const at=(a,i,j)=>a[((i%nx)+nx)%nx+(((j%ny)+ny)%ny)*nx];
  for(let j=0;j<ny;j++) for(let i=0;i<nx;i++){
    bxf0[i+j*nx]=0.5*(at(bx,i-1,j)+at(bx,i,j));
    byf0[i+j*nx]=0.5*(at(by,i,j-1)+at(by,i,j));
  }
  const d0=divergence(bxf0,byf0,nx,ny,dx,dy);
  const t0=performance.now();
  const r=project(bx,by,nx,ny,dx,dy);
  const ms=performance.now()-t0;
  const d=divergence(r.bxf,r.byf,nx,ny,dx,dy);
  const back=toCentre(r.bxf,r.byf,nx,ny);
  const ex=new Float64Array(n);
  for(let t=0;t<n;t++) ex[t]=Math.hypot(back.bx[t]-bx[t],back.by[t]-by[t]);
  console.log('\nT2 field with genuine divergence');
  console.log('   naive face interp:  div rms %s (normalised %s)',
    rms(d0).toExponential(2),(rms(d0)*dx/scale).toExponential(2));
  console.log('   after projection:   div rms %s (normalised %s)  max %s',
    rms(d).toExponential(2),(rms(d)*dx/scale).toExponential(2),(amax(d)*dx/scale).toExponential(2));
  console.log('   fidelity: |B_back - B_in| rms %s of |B|  (max %s)',
    (rms(ex)/scale).toExponential(2),(amax(ex)/scale).toExponential(2));
  console.log('   solve time %s ms at %dx%d', ms.toFixed(1), nx, ny);
  console.log('   T2 %s', amax(d)*dx/scale < 1e-12 ? 'PASS (divergence at round-off)':'FAIL');
}
