// 内置轮廓数据构建脚本
// 用法: node tools/build-shapes.mjs
// 读取 .ref/geo/*.json（阿里云 DataV GeoAtlas 边界），
// 取最大外环 → Douglas-Peucker 简化 → 归一化(外接圆半径=1, y 向下) → 输出 shapes.js
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const SRC = join(root, '.ref', 'geo');
const OUT = join(root, 'shapes.js');

// [key, 名称, admcode, 保留点数上限]
const SHAPES = [
  ['beijing', '北京', '110000', 360],
  ['tianjin', '天津', '120000', 360],
  ['hebei', '河北', '130000', 440],
  ['shanxi1', '山西', '140000', 440],
  ['neimenggu', '内蒙古', '150000', 480],
  ['liaoning', '辽宁', '210000', 440],
  ['jilin', '吉林', '220000', 440],
  ['heilongjiang', '黑龙江', '230000', 480],
  ['shanghai', '上海', '310000', 360],
  ['jiangsu', '江苏', '320000', 440],
  ['zhejiang', '浙江', '330000', 440],
  ['anhui', '安徽', '340000', 440],
  ['fujian', '福建', '350000', 440],
  ['jiangxi', '江西', '360000', 440],
  ['shandong', '山东', '370000', 440],
  ['henan', '河南', '410000', 440],
  ['hubei', '湖北', '420000', 440],
  ['hunan', '湖南', '430000', 440],
  ['guangdong', '广东', '440000', 440],
  ['guangxi', '广西', '450000', 440],
  ['hainan', '海南', '460000', 400],
  ['sichuan', '四川', '510000', 480],
  ['guizhou', '贵州', '520000', 440],
  ['yunnan', '云南', '530000', 480],
  ['xizang', '西藏', '540000', 440],
  ['shaanxi', '陕西', '610000', 440],
  ['gansu', '甘肃', '620000', 560],
  ['qinghai', '青海', '630000', 440],
  ['ningxia', '宁夏', '640000', 440],
  ['xinjiang', '新疆', '650000', 520],
  ['taiwan', '台湾', '710000', 440],
  ['xianggang', '香港', '810000', 400],
  ['aomen', '澳门', '820000', 400],
  ['china', '中国', '100000', 720],
];

// 取最大外环
function largestRing(geom) {
  let best = null, bestLen = 0;
  const polys = geom.type === 'MultiPolygon' ? geom.coordinates
    : geom.type === 'Polygon' ? [geom.coordinates] : [];
  for (const poly of polys) {
    const ring = poly[0];
    if (ring && ring.length > bestLen) { best = ring; bestLen = ring.length; }
  }
  return best;
}

function perim(pts, closed = true) {
  let L = 0;
  for (let i = 1; i < pts.length; i++) L += Math.hypot(pts[i][0] - pts[i - 1][0], pts[i][1] - pts[i - 1][1]);
  if (closed && pts.length > 1) {
    const a = pts[pts.length - 1], b = pts[0];
    L += Math.hypot(b[0] - a[0], b[1] - a[1]);
  }
  return L;
}

// Douglas-Peucker（闭合环：以最长直径两端为分段锚点，避免切割点随意）
function dpBetween(pts, i0, i1, keep, eps) {
  if (i1 <= i0 + 1) return;
  const ax = pts[i0][0], ay = pts[i0][1], bx = pts[i1][0], by = pts[i1][1];
  const dx = bx - ax, dy = by - ay;
  const len = Math.hypot(dx, dy) || 1e-12;
  let maxD = -1, maxI = -1;
  for (let i = i0 + 1; i < i1; i++) {
    const d = Math.abs((pts[i][0] - ax) * dy - (pts[i][1] - ay) * dx) / len;
    if (d > maxD) { maxD = d; maxI = i; }
  }
  if (maxD > eps) {
    keep.add(maxI);
    dpBetween(pts, i0, maxI, keep, eps);
    dpBetween(pts, maxI, i1, keep, eps);
  }
}

function simplify(ring, cap) {
  // 归一化到 0..1 便于选择 epsilon
  let minx = 1e9, miny = 1e9, maxx = -1e9, maxy = -1e9;
  for (const p of ring) {
    if (p[0] < minx) minx = p[0]; if (p[0] > maxx) maxx = p[0];
    if (p[1] < miny) miny = p[1]; if (p[1] > maxy) maxy = p[1];
  }
  const size = Math.max(maxx - minx, maxy - miny) || 1;
  const pts = ring.map(p => [(p[0] - minx) / size, (p[1] - miny) / size]);

  // 闭合环上找相距最远的两点作锚
  let n = pts.length, iA = 0, iB = Math.floor(n / 2), best = -1;
  for (let i = 0; i < n; i += Math.max(1, n >> 8)) {
    for (let j = i + 1; j < n; j += Math.max(1, n >> 8)) {
      const d = (pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2;
      if (d > best) { best = d; iA = i; iB = j; }
    }
  }
  if (iA > iB) [iA, iB] = [iB, iA];

  let eps = 0.0008, out = pts;
  for (;;) {
    const keep = new Set([iA, iB, 0, n - 1]);
    dpBetween(pts, iA, iB, keep, eps);
    dpBetween(pts, iB, n - 1, keep, eps);
    dpBetween(pts, 0, iA, keep, eps);
    // 环是闭合的，最后一段 n-1 -> 0 已含在 0..iA / iB..n-1 之外，需补上 n-1→0 段
    dpBetween(pts, n - 1, n - 1, keep, eps); // no-op 占位
    const idx = [...keep].sort((a, b) => a - b);
    out = idx.map(i => pts[i]);
    if (out.length <= cap || eps > 0.08) break;
    eps *= 1.35;
  }
  // 若仍超上限则等距抽稀
  if (out.length > cap) {
    const step = out.length / cap, arr = [];
    for (let i = 0; i < cap; i++) arr.push(out[Math.floor(i * step)]);
    out = arr;
  }
  return out;
}

function normalize(pts) {
  // 中心=外接框中心，缩放=最大半径 1，y 翻转（屏幕坐标向下）
  let minx = 1e9, miny = 1e9, maxx = -1e9, maxy = -1e9;
  for (const p of pts) {
    if (p[0] < minx) minx = p[0]; if (p[0] > maxx) maxx = p[0];
    if (p[1] < miny) miny = p[1]; if (p[1] > maxy) maxy = p[1];
  }
  const cx = (minx + maxx) / 2, cy = (miny + maxy) / 2;
  let maxR = 0;
  for (const p of pts) maxR = Math.max(maxR, Math.hypot(p[0] - cx, p[1] - cy));
  const s = maxR || 1;
  return pts.map(p => [
    +((p[0] - cx) / s).toFixed(4),
    +(-(p[1] - cy) / s).toFixed(4),
  ]);
}

const result = {};
for (const [key, name, code, cap] of SHAPES) {
  const gj = JSON.parse(readFileSync(join(SRC, `${code}.json`), 'utf8'));
  const feat = Array.isArray(gj.features) ? gj.features[0] : gj;
  const ring = largestRing(feat.geometry);
  if (!ring) { console.error(`${name}: 未找到轮廓`); continue; }
  const simp = simplify(ring, cap);
  const norm = normalize(simp);
  result[key] = { name, pts: norm.flat() };
  console.log(`${name}: ${ring.length} -> ${simp.length} 点, 周长比例 ${perim(simp).toFixed(3)}`);
}

const js = '// 由 tools/build-shapes.mjs 生成：内置轮廓（归一化坐标, y 向下）\n'
  + '// 数据来源：阿里云 DataV GeoAtlas\n'
  + 'window.FOURIER_SHAPES = ' + JSON.stringify(result) + ';\n';
mkdirSync(dirname(OUT), { recursive: true });
writeFileSync(OUT, js);
console.log('已写出', OUT, (js.length / 1024).toFixed(1) + 'KB');
