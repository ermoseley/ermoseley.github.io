const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const ok = await window.boot('mhd', window.__patch || []);
if (!ok) return { error: 'no field' };
await sleep(10000);
const out = { ms: {} };
for (const sv of [3, 2, 1, 0]) {
  window.F.set('solver', sv);
  await sleep(900);
  const a = []; for (let i = 0; i < 3; i++) a.push(Math.round(window.F.__bench(30)*100)/100);
  out.ms[['LLF','HLL','HLLC','HLLD'][sv]] = a;
}
window.F.set('solver', 3);
out.divb = window.F.readout().divbRms;
return out;
