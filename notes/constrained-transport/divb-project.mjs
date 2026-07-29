/* Cell-centred B  ->  divergence-free face-centred B, exactly, on a periodic box.
 *
 * In 2-D any in-plane divergence-free field is the curl of a scalar potential. Put A
 * at cell corners and define
 *
 *   Bx_f(i,j) [x-face i-1/2] = ( A(i,j+1) - A(i,j) ) / dy
 *   By_f(i,j) [y-face j-1/2] = -( A(i+1,j) - A(i,j) ) / dx
 *
 * and the discrete divergence telescopes the same four corner values with opposite
 * signs, so it is zero to round-off whatever A is. Convergence therefore buys
 * fidelity, never correctness -- which is what makes this safe to do mid-run.
 *
 * The cell-centred field the rest of the engine sees is the average of a cell's own
 * two faces, so the composite map A -> B_cell is linear and diagonal in Fourier
 * space. Solve it there in the least-squares sense and there is no iteration at all.
 *
 *   Bx_c(k) = A(k) * (ey - 1)(1 + ex) / (2 dy)   =  A(k) Cx
 *   By_c(k) = -A(k) * (ex - 1)(1 + ey) / (2 dx)  =  A(k) Cy
 *   A(k)    = ( conj(Cx) Bx_c(k) + conj(Cy) By_c(k) ) / (|Cx|^2 + |Cy|^2)
 *
 * k = 0 is the uniform field: not representable by a periodic A, and carried
 * separately as a mean, which is divergence-free on its own. The single
 * (Nyquist, Nyquist) checkerboard mode is also in the null space -- both Cx and Cy
 * vanish -- and is genuinely not expressible as the face average of any
 * divergence-free field, so it is dropped rather than fudged.
 */

// ---- a naive separable DFT with cached twiddles. Small grids, one-time cost. ----

const tw = new Map();
function twiddle(n) {
  let t = tw.get(n);
  if (t) return t;
  const c = new Float64Array(n * n), s = new Float64Array(n * n);
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      const a = (-2 * Math.PI * k * j) / n;
      c[k * n + j] = Math.cos(a);
      s[k * n + j] = Math.sin(a);
    }
  }
  t = { c, s };
  tw.set(n, t);
  return t;
}

// in-place-ish 1-D DFT of one strided complex line
function dft1(re, im, off, stride, n, inverse) {
  const { c, s } = twiddle(n);
  const or_ = new Float64Array(n), oi = new Float64Array(n);
  const sg = inverse ? -1 : 1;
  for (let k = 0; k < n; k++) {
    let ar = 0, ai = 0;
    const base = k * n;
    for (let j = 0; j < n; j++) {
      const xr = re[off + j * stride], xi = im[off + j * stride];
      const cr = c[base + j], ci = sg * s[base + j];
      ar += xr * cr - xi * ci;
      ai += xr * ci + xi * cr;
    }
    or_[k] = ar; oi[k] = ai;
  }
  const k = inverse ? 1 / n : 1;
  for (let j = 0; j < n; j++) {
    re[off + j * stride] = or_[j] * k;
    im[off + j * stride] = oi[j] * k;
  }
}

function dft2(re, im, nx, ny, inverse) {
  for (let j = 0; j < ny; j++) dft1(re, im, j * nx, 1, nx, inverse);
  for (let i = 0; i < nx; i++) dft1(re, im, i, nx, ny, inverse);
}

/* bxc, byc: Float64Array(nx*ny), cell-centred, row-major, index i + j*nx.
 * Returns { bxf, byf } face-centred: bxf(i,j) on the x-face at i-1/2,
 * byf(i,j) on the y-face at j-1/2. dx, dy in the same units the caller uses. */
export function project(bxc, byc, nx, ny, dx, dy) {
  const n = nx * ny;

  // the uniform part rides along untouched -- it is exactly divergence-free
  let mx = 0, my = 0;
  for (let t = 0; t < n; t++) { mx += bxc[t]; my += byc[t]; }
  mx /= n; my /= n;

  const xr = new Float64Array(n), xi = new Float64Array(n);
  const yr = new Float64Array(n), yi = new Float64Array(n);
  for (let t = 0; t < n; t++) { xr[t] = bxc[t] - mx; yr[t] = byc[t] - my; }

  dft2(xr, xi, nx, ny, false);
  dft2(yr, yi, nx, ny, false);

  const ar = new Float64Array(n), ai = new Float64Array(n);
  for (let jy = 0; jy < ny; jy++) {
    const py = (2 * Math.PI * jy) / ny;
    const eyr = Math.cos(py), eyi = Math.sin(py);
    for (let ix = 0; ix < nx; ix++) {
      const px = (2 * Math.PI * ix) / nx;
      const exr = Math.cos(px), exi = Math.sin(px);
      // Cx = (ey - 1)(1 + ex) / (2 dy)
      const a1r = eyr - 1, a1i = eyi;
      const b1r = 1 + exr, b1i = exi;
      const Cxr = (a1r * b1r - a1i * b1i) / (2 * dy);
      const Cxi = (a1r * b1i + a1i * b1r) / (2 * dy);
      // Cy = -(ex - 1)(1 + ey) / (2 dx)
      const a2r = exr - 1, a2i = exi;
      const b2r = 1 + eyr, b2i = eyi;
      const Cyr = -(a2r * b2r - a2i * b2i) / (2 * dx);
      const Cyi = -(a2r * b2i + a2i * b2r) / (2 * dx);

      const den = Cxr * Cxr + Cxi * Cxi + Cyr * Cyr + Cyi * Cyi;
      const t = ix + jy * nx;
      if (den < 1e-18) { ar[t] = 0; ai[t] = 0; continue; }   // k=0 and the checkerboard
      // conj(C) . B
      const nr = Cxr * xr[t] + Cxi * xi[t] + Cyr * yr[t] + Cyi * yi[t];
      const ni = Cxr * xi[t] - Cxi * xr[t] + Cyr * yi[t] - Cyi * yr[t];
      ar[t] = nr / den; ai[t] = ni / den;
    }
  }

  dft2(ar, ai, nx, ny, true);

  // exact discrete curl of the corner potential, plus the mean back
  const bxf = new Float64Array(n), byf = new Float64Array(n);
  const A = (i, j) => ar[((i % nx) + nx) % nx + (((j % ny) + ny) % ny) * nx];
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      bxf[i + j * nx] = (A(i, j + 1) - A(i, j)) / dy + mx;
      byf[i + j * nx] = -(A(i + 1, j) - A(i, j)) / dx + my;
    }
  }
  return { bxf, byf, a: ar, mean: [mx, my] };
}

/* Divergence of a face field, per cell, and the cell-centred reduction. */
export function divergence(bxf, byf, nx, ny, dx, dy) {
  const d = new Float64Array(nx * ny);
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const ip = (i + 1) % nx, jp = (j + 1) % ny;
      d[i + j * nx] = (bxf[ip + j * nx] - bxf[i + j * nx]) / dx
                    + (byf[i + jp * nx] - byf[i + j * nx]) / dy;
    }
  }
  return d;
}

export function toCentre(bxf, byf, nx, ny) {
  const bx = new Float64Array(nx * ny), by = new Float64Array(nx * ny);
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const ip = (i + 1) % nx, jp = (j + 1) % ny;
      bx[i + j * nx] = 0.5 * (bxf[i + j * nx] + bxf[ip + j * nx]);
      by[i + j * nx] = 0.5 * (byf[i + j * nx] + byf[i + jp * nx]);
    }
  }
  return { bx, by };
}
