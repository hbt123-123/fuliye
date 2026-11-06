// 数学核心自测：FFT + 级数重构 RMS
const M = 2048, KMAX = 1024, TWO_PI = Math.PI * 2;

function fft(re, im) {
  var n = re.length, i, j, bit, len, k;
  for (i = 1, j = 0; i < n; i++) {
    bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j |= bit;
    if (i < j) {
      var tr = re[i]; re[i] = re[j]; re[j] = tr;
      var ti = im[i]; im[i] = im[j]; im[j] = ti;
    }
  }
  for (len = 2; len <= n; len <<= 1) {
    var ang = -TWO_PI / len, wr = Math.cos(ang), wi = Math.sin(ang), half = len >> 1;
    for (i = 0; i < n; i += len) {
      var cwr = 1, cwi = 0;
      for (k = 0; k < half; k++) {
        var ur = re[i + k], ui = im[i + k];
        var xr = re[i + k + half], xi = im[i + k + half];
        var vr = xr * cwr - xi * cwi, vi = xr * cwi + xi * cwr;
        re[i + k] = ur + vr; im[i + k] = ui + vi;
        re[i + k + half] = ur - vr; im[i + k + half] = ui - vi;
        var nwr = cwr * wr - cwi * wi;
        cwi = cwr * wi + cwi * wr;
        cwr = nwr;
      }
    }
  }
}

// 测试1：FFT 提取已知正弦系数
{
  const re = new Float64Array(M), im = new Float64Array(M);
  for (let i = 0; i < M; i++) re[i] = Math.sin(TWO_PI * 3 * i / M) * 2;
  fft(re, im);
  const a3 = 2 * re[3] / M, b3 = -2 * im[3] / M;
  console.log('sin(6πt)·2  ->  a3 =', a3.toFixed(4), ' b3 =', b3.toFixed(4), '(期望 0, 2)');
  console.log('频谱泄漏 k=2:', (2 * re[2] / M).toFixed(6), ' k=4:', (2 * re[4] / M).toFixed(6));
}

// 测试2：方波两级数重构的 RMS 收敛（应为 1/k 收敛）
{
  const sq = new Float64Array(M * 2);
  for (let i = 0; i < M; i++) {
    const th = i / M * TWO_PI;
    sq[i * 2] = Math.cos(th) > 0 ? 1 : -1;
    sq[i * 2 + 1] = Math.sin(th) > 0 ? 0.5 : -0.5;
  }
  let mx = 0, my = 0;
  for (let i = 0; i < M; i++) { mx += sq[i * 2]; my += sq[i * 2 + 1]; }
  mx /= M; my /= M;
  let mr = 0;
  for (let i = 0; i < M; i++) {
    sq[i * 2] -= mx; sq[i * 2 + 1] -= my;
    mr = Math.max(mr, Math.hypot(sq[i * 2], sq[i * 2 + 1]));
  }
  for (let i = 0; i < M * 2; i++) sq[i] /= mr;

  const xr = new Float64Array(M), xi = new Float64Array(M),
        yr = new Float64Array(M), yi = new Float64Array(M);
  for (let i = 0; i < M; i++) { xr[i] = sq[i * 2]; yr[i] = sq[i * 2 + 1]; }
  fft(xr, xi); fft(yr, yi);
  const A = new Float64Array(KMAX + 1), B = new Float64Array(KMAX + 1),
        D = new Float64Array(KMAX + 1), E = new Float64Array(KMAX + 1);
  const a0 = xr[0] / M, b0 = yr[0] / M;
  for (let k = 1; k <= KMAX; k++) {
    A[k] = 2 * xr[k] / M; B[k] = -2 * xi[k] / M;
    D[k] = 2 * yr[k] / M; E[k] = -2 * yi[k] / M;
  }
  A[KMAX] = xr[KMAX] / M; D[KMAX] = yr[KMAX] / M;

  function evalXY(t, K) {
    let x = a0, y = b0;
    for (let k = 1; k <= K; k++) {
      const th = TWO_PI * k * t, c = Math.cos(th), s = Math.sin(th);
      x += A[k] * c + B[k] * s;
      y += D[k] * c + E[k] * s;
    }
    return [x, y];
  }
  for (const K of [1, 4, 64, 1024]) {
    let num = 0, den = 0;
    for (let i = 0; i < M; i += 3) {
      const [x, y] = evalXY(i / M, K);
      const ex = x - sq[i * 2], ey = y - sq[i * 2 + 1];
      num += ex * ex + ey * ey;
      den += sq[i * 2] ** 2 + sq[i * 2 + 1] ** 2;
    }
    console.log('方波 K=' + String(K).padEnd(4), 'RMS =',
      (Math.sqrt(num / (M / 3)) / Math.sqrt(den / (M / 3)) * 100).toFixed(4) + '%');
  }
}
