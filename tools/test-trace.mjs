// 图片轮廓提取管线自测：合成二值图 → 最大连通域 → Moore 边界追踪 → DP 简化 → 归一化
// 与 index.html 中的实现保持一致（复制自同版本代码）

function largestComponent(bin, w, h) {
  var label = new Int32Array(w * h), stack = new Int32Array(w * h);
  var best = 0, bestId = 0, id = 0, i;
  for (i = 0; i < w * h; i++) {
    if (bin[i] && !label[i]) {
      id++; var sp = 0, cnt = 0;
      stack[sp++] = i; label[i] = id;
      while (sp) {
        var p = stack[--sp]; cnt++;
        var px = p % w;
        if (px > 0 && bin[p - 1] && !label[p - 1]) { label[p - 1] = id; stack[sp++] = p - 1; }
        if (px < w - 1 && bin[p + 1] && !label[p + 1]) { label[p + 1] = id; stack[sp++] = p + 1; }
        if (p >= w && bin[p - w] && !label[p - w]) { label[p - w] = id; stack[sp++] = p - w; }
        if (p < w * (h - 1) && bin[p + w] && !label[p + w]) { label[p + w] = id; stack[sp++] = p + w; }
      }
      if (cnt > best) { best = cnt; bestId = id; }
    }
  }
  if (!bestId) return null;
  var out = new Uint8Array(w * h);
  for (i = 0; i < w * h; i++) out[i] = label[i] === bestId ? 1 : 0;
  return out;
}

function traceBoundary(mask, w, h) {
  var start = -1, x, y;
  for (y = 0; y < h && start < 0; y++)
    for (x = 0; x < w; x++) { if (mask[y * w + x]) { start = y * w + x; break; } }
  if (start < 0) return null;
  var dirs = [[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1]];
  function fg(x2, y2) { return x2 >= 0 && y2 >= 0 && x2 < w && y2 < h && mask[y2 * w + x2] === 1; }
  var sx = start % w, sy = (start / w) | 0, cx = sx, cy = sy;
  var bx = cx - 1, by = cy;
  var pts = [], guard = 0;
  do {
    pts.push([cx, cy]);
    var d0 = 0, ddx = bx - cx, ddy = by - cy;
    for (var d = 0; d < 8; d++) { if (dirs[d][0] === ddx && dirs[d][1] === ddy) { d0 = d; break; } }
    var found = false;
    for (var i = 1; i <= 8; i++) {
      var dd = (d0 + i) % 8;
      var nx = cx + dirs[dd][0], ny = cy + dirs[dd][1];
      if (fg(nx, ny)) { bx = cx; by = cy; cx = nx; cy = ny; found = true; break; }
    }
    if (!found) break;
  } while (!(cx === sx && cy === sy) && ++guard < 40000);
  return pts.length > 8 ? pts : null;
}

function simplifyRing(ring, cap) {
  var n = ring.length;
  if (n <= cap) return ring;
  var minx = 1e9, miny = 1e9, maxx = -1e9, maxy = -1e9, i;
  for (i = 0; i < n; i++) {
    if (ring[i][0] < minx) minx = ring[i][0]; if (ring[i][0] > maxx) maxx = ring[i][0];
    if (ring[i][1] < miny) miny = ring[i][1]; if (ring[i][1] > maxy) maxy = ring[i][1];
  }
  var size = Math.max(maxx - minx, maxy - miny) || 1;
  var pts = ring.map(p => [(p[0] - minx) / size, (p[1] - miny) / size]);
  var iA = 0, iB = (n / 2) | 0, best = -1;
  for (i = 0; i < n; i += Math.max(1, n >> 8)) {
    for (var j = i + 1; j < n; j += Math.max(1, n >> 8)) {
      var d = (pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2;
      if (d > best) { best = d; iA = i; iB = j; }
    }
  }
  if (iA > iB) { var tmp = iA; iA = iB; iB = tmp; }
  function dp(i0, i1, keep, eps) {
    if (i1 <= i0 + 1) return;
    var ax = pts[i0][0], ay = pts[i0][1], bx2 = pts[i1][0], by2 = pts[i1][1];
    var dx = bx2 - ax, dy = by2 - ay, len = Math.hypot(dx, dy) || 1e-12;
    var maxD = -1, maxI = -1;
    for (var k2 = i0 + 1; k2 < i1; k2++) {
      var dist = Math.abs((pts[k2][0] - ax) * dy - (pts[k2][1] - ay) * dx) / len;
      if (dist > maxD) { maxD = dist; maxI = k2; }
    }
    if (maxD > eps) {
      keep.add(maxI);
      dp(i0, maxI, keep, eps);
      dp(maxI, i1, keep, eps);
    }
  }
  var eps = 0.002, out = pts;
  for (;;) {
    var keep = new Set([0, n - 1, iA, iB]);
    dp(iA, iB, keep, eps);
    dp(iB, n - 1, keep, eps);
    dp(0, iA, keep, eps);
    var idx = Array.from(keep).sort((a, b) => a - b);
    out = idx.map(i2 => pts[i2]);
    if (out.length <= cap || eps > 0.08) break;
    eps *= 1.4;
  }
  if (out.length > cap) {
    var step2 = out.length / cap, arr = [];
    for (i = 0; i < cap; i++) arr.push(out[(i * step2) | 0]);
    out = arr;
  }
  return out;
}

/* ---------- 合成场景 1：大圆 + 小方块噪点（应只取大圆） ---------- */
{
  const w = 200, h = 200, bin = new Uint8Array(w * h);
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    const d = Math.hypot(x - 100, y - 100);
    if (d < 70) bin[y * w + x] = 1;
  }
  for (let y = 10; y < 20; y++) for (let x = 10; x < 20; x++) bin[y * w + x] = 1;  // 噪点块
  const comp = largestComponent(bin, w, h);
  let cnt = 0; for (let i = 0; i < w * h; i++) cnt += comp[i];
  const ring = traceBoundary(comp, w, h);
  const simp = simplifyRing(ring, 520);
  // 归一化后测到理想圆的相对误差
  let cx = 0, cy = 0;
  for (const p of simp) { cx += p[0]; cy += p[1]; }
  cx /= simp.length; cy /= simp.length;
  let mr = 0;
  for (const p of simp) mr = Math.max(mr, Math.hypot(p[0] - cx, p[1] - cy));
  let err = 0;
  for (const p of simp) err += Math.abs(Math.hypot(p[0] - cx, p[1] - cy) / mr - 70 / 70.0);
  console.log('[圆形+噪点] 前景', cnt, '(期望 ~15386)  边界点', ring.length, '→ 简化', simp.length,
    '  归一化半径误差', (err / simp.length * 100).toFixed(2) + '%');
}

/* ---------- 合成场景 2：双连通域取大者 ---------- */
{
  const w = 120, h = 80, bin = new Uint8Array(w * h);
  for (let y = 10; y < 40; y++) for (let x = 10; x < 40; x++) bin[y * w + x] = 1;    // 30×30
  for (let y = 30; y < 70; y++) for (let x = 60; x < 110; x++) bin[y * w + x] = 1;  // 50×40 更大
  const comp = largestComponent(bin, w, h);
  let c1 = 0, c2 = 0;
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    if (comp[y * w + x]) (x < 50 ? c1++ : c2++);
  }
  const ring = traceBoundary(comp, w, h);
  console.log('[双域取大] 左域像素', c1, '(期望 0)  右域像素', c2, '(期望 2000)  边界闭合点数', ring.length);
}

/* ---------- 合成场景 3：环形（只取外边界） ---------- */
{
  const w = 160, h = 160, bin = new Uint8Array(w * h);
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    const d = Math.hypot(x - 80, y - 80);
    if (d < 60 && d > 30) bin[y * w + x] = 1;
  }
  const ring = traceBoundary(largestComponent(bin, w, h), w, h);
  const simp = simplifyRing(ring, 400);
  let cx = 0, cy = 0;
  for (const p of simp) { cx += p[0]; cy += p[1]; }
  cx /= simp.length; cy /= simp.length;
  let rMax = 0, rMin = 1e9;
  for (const p of simp) {
    const r = Math.hypot(p[0] - cx, p[1] - cy);
    if (r > rMax) rMax = r; if (r < rMin) rMin = r;
  }
  console.log('[环形] 外边界点', ring.length, '→ 简化', simp.length,
    '  半径比 rMin/rMax', (rMin / rMax).toFixed(3), '(接近 1 → 只含外圈 ~60，未混入 30 的内圈)');
}
