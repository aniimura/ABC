#!/usr/bin/env node
/**
 * 持ち場に書いた「見当」が当たったかを**数える**(backlog M45 / メタ第 14 回)。
 *
 * ★なぜ要るか(M45 の負の結果)
 *   メタ第 13 回は「持ち場に未確認の断りを書くと誤りが覆る」という本体の見立てを
 *   検証しようとして、**原理的に測れない**という結論に達した:
 *     - 断りを書いた持ち場 21 件のうち、当否が判定できた 14 件・覆った 8 件(57%)。
 *     - ★しかし**断りを書かなかった**持ち場でも覆っている(#16 の `pdfPage := 8`)。
 *     - ★★**覆らなかった誤りは記録に現れない**ので、**分母が作れない**。
 *
 * ★★この道具が変えるのはただ 1 点: **分母を作る。**
 *   ⇒ **見当は「自信がないもの」だけでなく、断定して渡すものも含めて全部 `GUESS:` に書く。**
 *     `確度=高` の行が無ければ「断定した見当が何件あって、そのうち何件外れたか」が
 *     永遠に分からない。★**自信のあるものを書かないと、この道具は M45 と同じ壁に当たる。**
 *
 * ★書式(`ResearchPaper/autonomy-policy.md` §4.7 が正典)
 *
 *     GUESS[<持ち場>]: <見当を 1 行で> | 確度=<高|中|低> | 検算=<型|実機|索引|原文|なし>
 *     VERDICT[<持ち場>]: <当たり|外れ|半分|未確認> — <一言>
 *
 *   - `<持ち場>` は本体が既に使っている呼び名(`Y19e` / `#16` / `第1075波-3` など)。
 *     ★**GUESS と VERDICT で同じ字面**にすること。ここだけが対応の鍵である。
 *   - `GUESS` は**配るとき**に、`VERDICT` は**報告が返ったとき**に書く。
 *   - 1 つの持ち場に見当が複数あるなら `Y19e-a` `Y19e-b` のように枝番を付ける。
 *
 * ★これは**ゲートではない**。`check.mjs` に入れてはならない ——
 *   落とすようにすると「見当を書かない」ことで逃げられ、正直さを罰する
 *   (G7 の基準・G1 の `Found/` 非対称と同じ理由)。
 *
 * 使い方:
 *   node tools/unverified.mjs                # 集計
 *   node tools/unverified.mjs --open         # まだ VERDICT が無いものだけ
 *   node tools/unverified.mjs --json
 *   node tools/unverified.mjs --selftest     # 器具の較正(書式**と統計**を変えたら必ず走らせる)
 *   node tools/unverified.mjs --file <path>  # 別の記録を読む
 */

import { readFileSync, existsSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const args = process.argv.slice(2);
const flag = (f) => args.includes(f);
const opt = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };

/** 既定で読む記録。★増やすときはここだけを触る。 */
const SOURCES = [
  join(ROOT, 'ResearchPaper', 'decisions-pending.md'),
  join(ROOT, 'ResearchPaper', 'meta-backlog.md'),
];

const GUESS_RE = /^\s*(?:[★☆*\-\s]*)GUESS\[([^\]]+)\]\s*[:：]\s*(.*)$/;
const VERDICT_RE = /^\s*(?:[★☆*\-\s]*)VERDICT\[([^\]]+)\]\s*[:：]\s*(.*)$/;
const CONF_RE = /確度\s*=\s*([高中低])/;
const HOW_RE = /検算\s*=\s*(型|実機|索引|原文|なし)/;
const OUT_RE = /^\s*(当たり|外れ|半分|未確認)/;
const CONFS = ['高', '中', '低', '(無記入)'];

/** 記録 1 本を読んで GUESS / VERDICT を拾う。 */
function parse(path) {
  const guesses = new Map();   // 持ち場 → { conf, how, text, file, line }
  const verdicts = new Map();  // 持ち場 → { outcome, text, file, line }
  if (!existsSync(path)) return { guesses, verdicts, missing: true };
  const lines = readFileSync(path, 'utf8').split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    let m = GUESS_RE.exec(lines[i]);
    // ★★書式の見本(`GUESS[<持ち場>]`)を数えない(メタ第 19 回 §7 の指摘)。
    //   見出しの件数と表の合計が 1 だけ合わなかった原因。
    //   雛形は `確度=<高|中|低>` で列挙値に合わず表からは落ちていた。
    if (m && /[<>]/.test(m[1])) continue;
    if (m) {
      guesses.set(m[1].trim(), {
        conf: CONF_RE.exec(m[2])?.[1] ?? '(無記入)',
        how: HOW_RE.exec(m[2])?.[1] ?? '(無記入)',
        text: m[2].split('|')[0].trim(),
        file: path, line: i + 1,
      });
      continue;
    }
    m = VERDICT_RE.exec(lines[i]);
    if (m && /[<>]/.test(m[1])) continue;
    if (m) {
      verdicts.set(m[1].trim(), {
        outcome: OUT_RE.exec(m[2])?.[1] ?? '(読めない)',
        text: m[2].replace(OUT_RE, '').replace(/^\s*[—-]\s*/, '').trim(),
        file: path, line: i + 1,
      });
    }
  }
  return { guesses, verdicts, missing: false };
}

function collect(paths) {
  const guesses = new Map(); const verdicts = new Map(); const missing = [];
  for (const p of paths) {
    const r = parse(p);
    if (r.missing) { missing.push(p); continue; }
    for (const [k, v] of r.guesses) guesses.set(k, v);
    for (const [k, v] of r.verdicts) verdicts.set(k, v);
  }
  return { guesses, verdicts, missing };
}

/* ★★★統計 —— 「向きの断定」をやめるために足した部分(メタ第 22 回 / M90)
 *
 * ★なぜ要ったか(★この道具自身の罪)
 *   この道具は 33 件の標本に対して
 *     「断定(確度=高)の覆った率 33% 対 断り付きの 46% —— 断り付きの方がよく覆っている
 *       (断りが実装者の検算を誘発している可能性)」
 *   と**向きを断定して**印字していた。★同じ行の下で「向きを断定するな」と書きながら、である。
 *   ★実際の Fisher 両側は p = 0.6982(高 対 中低)/ p = 1.0000(高 対 低)。★単調ですらない。
 *   ★★台帳 M88 は「データを見てから選んだ 2 つの兆しが標本 +8 件で両方とも平均へ戻った」
 *   ことを実演している(look-elsewhere)。★道具が同じ罠を踏んでいては見張りにならない。
 *
 * ★方針(ここが変えた点)
 *   1. 比較の族(FAMILY)を**データを見る前にコードへ固定する**。増やせば m が増える。
 *   2. Fisher の正確検定(両側・確率法)を BigInt で**厳密に**計算する(丸め誤差なし)。
 *   3. Bonferroni(α = 0.05 / m)で **「言える / 言えない」だけ**を出す。★向きは書かない。
 *   4. 「あと何件で言えるか」も**多重比較込み**で出す(M78 の「あと 1 件」の自己訂正の再発防止)。
 *   5. ★交絡(確度 × 検算)を必ず出す。★片方だけの表を主役にしない。
 *   6. ★束(持ち場)ごとのばらつきを出す。33 件は 8 束から来ており独立ではない。
 */

/** 表示幅(全角と ★ を 2 と数える)。★padEnd は符号点数なので日本語の表が崩れる。 */
const WIDE = /[\u2605\u2606\u2190-\u21FF\u2500-\u257F\u3000-\u303F\u3040-\u30FF\u3400-\u4DBF\u4E00-\u9FFF\uF900-\uFAFF\uFF01-\uFF60\uFFE0-\uFFE6]/;
const dw = (t) => [...String(t)].reduce((n, ch) => n + (WIDE.test(ch) ? 2 : 1), 0);
const padE = (t, n) => String(t) + ' '.repeat(Math.max(0, n - dw(t)));
const padS = (t, n) => ' '.repeat(Math.max(0, n - dw(t))) + String(t);

/** ★厳密な二項係数(BigInt)。r*(n-i)/(i+1) は各段で必ず割り切れる。 */
function binom(n, k) {
  if (k < 0 || k > n) return 0n;
  const K = BigInt(Math.min(k, n - k)); const N = BigInt(n);
  let r = 1n;
  for (let i = 0n; i < K; i++) r = (r * (N - i)) / (i + 1n);
  return r;
}
/** BigInt の比を double に落とす(桁あふれを避けるため先に 10^15 倍する)。 */
function bigRatio(a, b) { const S = 10n ** 15n; return Number((a * S) / b) / 1e15; }

/**
 * Fisher の正確検定(両側・確率法)。表は [[a,b],[c,d]]。
 * ★確率の比較は BigInt の分子どうしで行う(共通分母 C(n,k))ので**厳密**。
 *   ★浮動小数の「p <= p_obs」判定でよくある取りこぼしが原理的に起きない。
 */
function fisher2x2(a, b, c, d) {
  const r1 = a + b, r2 = c + d, k = a + c, n = a + b + c + d;
  if (!r1 || !r2 || !k || k === n) return 1;
  const num = (x) => binom(r1, x) * binom(r2, k - x);
  const obs = num(a);
  let tot = 0n, hit = 0n;
  for (let x = Math.max(0, k - r2); x <= Math.min(k, r1); x++) {
    const p = num(x); tot += p; if (p <= obs) hit += p;
  }
  return bigRatio(hit, tot);
}

/** ★対数ガンマ(Lanczos)。★外挿(needMore)で N が数千になると BigInt では遅すぎるため。 */
const LGC = [676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
  12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
function lgamma(z) {
  if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
  z -= 1; let x = 0.99999999999980993;
  for (let i = 0; i < 8; i++) x += LGC[i] / (z + i + 1);
  const t = z + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}
const lfact = (n) => lgamma(n + 1);
/** Fisher 両側(対数版)。★厳密版と食い違わないことを selftest で確かめてある。 */
function fisherLog(a, b, c, d) {
  const r1 = a + b, r2 = c + d, k = a + c, n = a + b + c + d;
  if (!r1 || !r2 || !k || k === n) return 1;
  const base = lfact(r1) + lfact(r2) + lfact(k) + lfact(n - k) - lfact(n);
  const lp = (x) => base - lfact(x) - lfact(r1 - x) - lfact(k - x) - lfact(r2 - k + x);
  const obs = lp(a); let sum = 0;
  for (let x = Math.max(0, k - r2); x <= Math.min(k, r1); x++) {
    const v = lp(x);
    if (v <= obs + 1e-9) sum += Math.exp(v);
  }
  return Math.min(1, sum);
}
/** ★表が小さいうちは厳密(BigInt)、大きくなったら対数版。 */
const fisherAny = (a, b, c, d) => (a + b + c + d <= 300 ? fisher2x2(a, b, c, d) : fisherLog(a, b, c, d));

/**
 * ★★比較の族 —— **データを見る前に**ここに固定する。
 *   ★項目を足したら m が増え、全部の判定が厳しくなる。「効いていそうな欄」を
 *   見てから足すのは M88 が実演した罠そのものなので、足すときは台帳に理由を書くこと。
 */
const FAMILY = [
  { id: 'F1', label: '確度: 高 対 中・低',
    A: (r) => r.conf === '高', B: (r) => r.conf === '中' || r.conf === '低' },
  { id: 'F2', label: '確度: 高 対 低(単調性)',
    A: (r) => r.conf === '高', B: (r) => r.conf === '低' },
  { id: 'F3', label: '検算: なし 対 あり',
    A: (r) => r.how === 'なし', B: (r) => r.how !== 'なし' && r.how !== '(無記入)' },
  { id: 'F4', label: '検算: なし 対 あり(★確度=中 に限る。交絡を外した層)',
    A: (r) => r.how === 'なし' && r.conf === '中',
    B: (r) => r.how !== 'なし' && r.how !== '(無記入)' && r.conf === '中' },
];
const ALPHA = 0.05;

/** 今の内訳のまま t 倍に増えたとき、Bonferroni 補正後に α を割るか。★届かないなら null。 */
function needMore(t, m, cap = 400) {
  const [[a, b], [c, d]] = t;
  const n = a + b + c + d;
  if (!n) return null;
  const rA = (a + b) ? a / (a + b) : 0, rB = (c + d) ? c / (c + d) : 0;
  if (rA === rB) return { impossible: true };
  for (let mult = 1; mult <= cap; mult++) {
    if (fisherAny(a * mult, b * mult, c * mult, d * mult) * m < ALPHA) {
      return { mult, total: n * mult, more: n * (mult - 1) };
    }
  }
  return null;
}

/** 持ち場の鍵から束(枝番を落としたもの)を作る。`B5-a` → `B5` / `第1075波-3` → `第1075波`。 */
function bundleOf(key) { return key.replace(/-[^-]*$/, ''); }

/** 束内相関(ANOVA 推定)と設計効果。★33 件が独立でないことを数字にする。 */
function clusterStats(judged, over) {
  const by = new Map();
  for (const r of judged) { const k = bundleOf(r.k); by.set(k, [...(by.get(k) ?? []), over(r) ? 1 : 0]); }
  const cl = [...by].map(([k, v]) => ({ k, n: v.length, o: v.reduce((s, x) => s + x, 0) }));
  const N = judged.length, K = cl.length;
  if (K < 2 || N === K) return { cl, N, K, icc: null };
  const pbar = cl.reduce((s, c) => s + c.o, 0) / N;
  const n0 = (N - cl.reduce((s, c) => s + c.n * c.n, 0) / N) / (K - 1);
  const MSB = cl.reduce((s, c) => s + c.n * Math.pow(c.o / c.n - pbar, 2), 0) / (K - 1);
  const MSW = cl.reduce((s, c) => s + c.n * (c.o / c.n) * (1 - c.o / c.n), 0) / (N - K);
  const icc = MSW === 0 && MSB === 0 ? 0 : (MSB - MSW) / (MSB + (n0 - 1) * MSW);
  const nbar = N / K;
  const deff = Math.max(1, 1 + (nbar - 1) * Math.max(0, icc));
  return { cl, N, K, pbar, icc, deff, neff: N / deff };
}

function report({ guesses, verdicts }) {
  const rows = [...guesses].map(([k, g]) => ({ k, ...g, v: verdicts.get(k) ?? null }));
  const orphan = [...verdicts].filter(([k]) => !guesses.has(k));
  const judged = rows.filter((r) => r.v && r.v.outcome !== '未確認' && r.v.outcome !== '(読めない)');
  const over = (r) => r.v.outcome === '外れ' || r.v.outcome === '半分';

  const nOver = judged.filter((r) => r.v.outcome === '外れ' || r.v.outcome === '半分').length;
  console.log(`# 見当と当否 —— GUESS ${rows.length} 件 / VERDICT ${verdicts.size} 件` +
    (judged.length ? ` / 判定済 ${judged.length} 件 / 覆った ${nOver} 件 = ${(nOver * 100 / judged.length).toFixed(0)}%` : ''));
  if (!rows.length) {
    console.log('');
    console.log('★GUESS が 1 件も無い。書式は autonomy-policy.md §4.7。');
    console.log('★★M45 の壁(分母が作れない)は、**断定して渡す見当も `確度=高` で書く**ことでしか越えられない。');
    return 0;
  }
  console.log('');
  console.log('## 内訳 1: 確度 × 当否(★これは結論ではない。下の「言えるか」を見ること)');
  console.log('');
  console.log('  確度      判定済   覆った   当たり   覆った率');
  for (const c of CONFS) {
    const g = judged.filter((r) => r.conf === c);
    if (!g.length) continue;
    const o = g.filter(over).length;
    console.log(`  ${c.padEnd(8)}  ${String(g.length).padStart(6)}  ${String(o).padStart(6)}  ` +
      `${String(g.length - o).padStart(6)}  ${(o * 100 / g.length).toFixed(0).padStart(7)}%`);
  }
  // ---- 内訳 2: 検算 × 当否 ------------------------------------------------
  console.log('');
  console.log('## 内訳 2: 検算の指定ごと(★これも結論ではない)');
  const HOWS = ['なし', '型', '実機', '索引', '原文', '(無記入)'];
  const hows = HOWS.filter((h) => judged.some((r) => r.how === h));
  for (const h of hows) {
    const g = judged.filter((r) => r.how === h);
    console.log(`  ${h.padEnd(8)} ${String(g.length).padStart(3)} 件 / 覆った ${g.filter(over).length}`);
  }

  // ---- 交絡 ---------------------------------------------------------------
  console.log('');
  console.log('## ★★交絡 —— 確度 と 検算 は独立でない(★片方だけの表を主役にしない)');
  console.log('');
  console.log(`  ${padE('確度＼検算', 12)}${hows.map((h) => padS(h, 8)).join('')}${padS('計', 9)}   (件数/覆り)`);
  for (const c of CONFS) {
    const g = judged.filter((r) => r.conf === c);
    if (!g.length) continue;
    const cells = hows.map((h) => {
      const q = g.filter((r) => r.how === h);
      return padS(q.length ? `${q.length}/${q.filter(over).length}` : '·', 8);
    });
    console.log(`  ${padE(c, 12)}${cells.join('')}${padS(`${g.length}/${g.filter(over).length}`, 9)}`);
  }
  {
    const hiNone = judged.filter((r) => r.conf === '高' && r.how === 'なし').length;
    const hiSome = judged.filter((r) => r.conf === '高' && r.how !== 'なし').length;
    const loNone = judged.filter((r) => r.conf !== '高' && r.how === 'なし').length;
    const loSome = judged.filter((r) => r.conf !== '高' && r.how !== 'なし').length;
    const pc = fisher2x2(hiNone, hiSome, loNone, loSome);
    console.log('');
    console.log(`  ★2×2 に畳むと 確度=高 は 検算なし ${hiNone} / あり ${hiSome}、` +
      `それ以外は なし ${loNone} / あり ${loSome} —— Fisher 両側 p = ${pc.toFixed(4)}`);
    console.log('  ★★これは**当否についての検定ではない**(族の外。2 つの欄が独立かを見ているだけ)。');
    console.log('  ⇒ p が小さいほど 2 欄は絡んでおり、「検算が効くか」を見ているつもりで');
    console.log('    「確度」を見ている恐れが強い。★上の 2 つの内訳表だけで結論を出さないこと。');
  }

  // ---- 言えるか / 言えないか ----------------------------------------------
  console.log('');
  console.log('## ★★★言えるか / 言えないか(★向きは断定しない)');
  console.log('');
  const tests = [];
  for (const f of FAMILY) {
    const A = judged.filter(f.A), B = judged.filter(f.B);
    const t = [[A.filter(over).length, A.length - A.filter(over).length],
               [B.filter(over).length, B.length - B.filter(over).length]];
    tests.push({ f, A, B, t, ok: A.length > 0 && B.length > 0 });
  }
  const run = tests.filter((x) => x.ok);
  const m = run.length || 1;
  console.log(`  ★比較の族は**データを見る前に**コードへ固定してある(FAMILY)。` +
    `いま検定できるのは ${m} 本 ⇒ Bonferroni の閾は α = 0.05/${m} = ${(ALPHA / m).toFixed(4)}`);
  console.log('');
  const LW = Math.max(28, ...tests.map((x) => dw(x.f.label)));
  console.log(`  ${padE('比較', LW)}  ${padE('件数(覆り)', 18)}${padS('Fisher両側', 12)}${padS('×m', 9)}  判定`);
  let said = 0;
  for (const x of tests) {
    if (!x.ok) {
      console.log(`  ${padE(x.f.label, LW)}  ${padE(`${x.A.length} 対 ${x.B.length}`, 18)}${padS('—', 12)}${padS('—', 9)}  ★件数不足(検定していない)`);
      continue;
    }
    const p = fisher2x2(x.t[0][0], x.t[0][1], x.t[1][0], x.t[1][1]);
    x.p = p;
    const adj = Math.min(1, p * m);
    const verdict = adj < ALPHA ? '★言える(補正後)' : '言えない';
    if (adj < ALPHA) said++;
    const counts = `${x.t[0][0]}/${x.A.length} 対 ${x.t[1][0]}/${x.B.length}`;
    console.log(`  ${padE(x.f.label, LW)}  ${padE(counts, 18)}${padS(p.toFixed(4), 12)}${padS(adj.toFixed(4), 9)}  ${verdict}`);
  }
  console.log('');
  if (!said) {
    console.log('  ⇒ ★★**どれも閾を割らない。今の件数では、どの欄にも識別力があるとは言えない。**');
    console.log('     ★「言えない」は「差が無い」ではない。★そして**向きも書かない**' +
      '(率の大小を読み上げた瞬間に M88 の罠に落ちる)。');
  } else {
    console.log(`  ⇒ ${said} 本が閾を割った。★それでも「向き」は次の標本で確かめること` +
      '(M88: 兆し 2 つが標本 +8 件で両方消えた)。');
  }

  // ---- あと何件で言えるか(多重比較込み) ----------------------------------
  console.log('');
  console.log('## あと何件で言えるか(★多重比較を織り込んだ数)');
  const LW2 = Math.max(28, ...run.map((x) => dw(x.f.label)));
  for (const x of run) {
    const n = x.A.length + x.B.length;
    const need = needMore(x.t, m);
    if (need && need.impossible) {
      console.log(`  ${padE(x.f.label, LW2)}  ★覆った率が両群で同じ ⇒ この比のままなら何件集めても届かない`);
    } else if (!need) {
      console.log(`  ${padE(x.f.label, LW2)}  ★400 倍(N>${n * 400})でも届かない`);
    } else {
      if (need.mult === 1) {
        console.log(`  ${padE(x.f.label, LW2)}  ★もう届いている(あと 0 件)`);
        continue;
      }
      console.log(`  ${padE(x.f.label, LW2)}  今 N=${String(n).padStart(3)} ⇒ ` +
        `この比のまま ★あと ${need.more} 件(合計 ${need.total}、${need.mult} 倍)`);
    }
  }
  console.log('  ★これは「今見えている差が本物なら」の話であって、差がある証拠ではない。');
  console.log('  ★★M78 は「あと 1 件で決まる」と書き、Bonferroni で「あと 3 件」に自己訂正し、');
  console.log('     ★その 3 件目が当たりだったので兆しごと消えた(M88)。★同じ書き方をしないこと。');

  // ---- 束(独立でない) ----------------------------------------------------
  {
    const cs = clusterStats(judged, over);
    console.log('');
    console.log(`## ★独立でない —— 判定済 ${cs.N} 件は **${cs.K} 束**(持ち場)から来ている`);
    console.log('  ' + cs.cl.sort((a, b) => b.o / b.n - a.o / a.n)
      .map((c) => `${c.k} ${c.o}/${c.n}`).join(' ・ '));
    if (cs.icc !== null) {
      console.log(`  束内相関 ICC = ${cs.icc.toFixed(3)} / 設計効果 DEFF = ${cs.deff.toFixed(2)} ` +
        `⇒ ★有効件数は ${cs.N} ではなく **約 ${cs.neff.toFixed(0)} 件**`);
      console.log('  ⇒ 上の p は「1 件ずつ独立」として計算している。★束ごと覆る癖があるぶん');
      console.log('    p は**小さめに出る**(甘い)。★「言えない」はより強く言えるが、');
      console.log('    「言える」が出たときは束を疑うこと。');
    }
  }

  const open = rows.filter((r) => !r.v);
  console.log('');
  console.log(`## VERDICT 待ち ${open.length} 件`);
  for (const r of open) console.log(`  ${r.k.padEnd(14)} 確度=${r.conf} 検算=${r.how}  ${r.text.slice(0, 60)}`);
  if (orphan.length) {
    console.log('');
    console.log(`## ★対応する GUESS が無い VERDICT ${orphan.length} 件(鍵の字面が食い違っている)`);
    for (const [k, v] of orphan) console.log(`  ${k.padEnd(14)} ${v.outcome}`);
  }
  return rows.length;
}

/** 器具の較正。★書式(正規表現)を変えたら必ずここを走らせること。 */
function selftest() {
  const fx = [
    'GUESS[Y19e]: N(-α) は最小多項式の定数項に落ちる | 確度=高 | 検算=型',
    '★GUESS[#16-a]: pdfPage は 8 | 確度=中 | 検算=原文',
    '- GUESS[Z1]: 何か | 確度=低 | 検算=実機',
    'GUESS[Z2]: 判定待ちのもの | 確度=高 | 検算=索引',
    'VERDICT[Y19e]: 当たり — 型が合った',
    '★VERDICT[#16-a]: 外れ — 実装者が 7 を採った',
    'VERDICT[Z1]: 半分 — 片方だけ',
    'VERDICT[ZZ]: 当たり — 対応する GUESS が無い',
    'この行は GUESS[ を含むが書式ではない',
  ].join('\n');
  const path = join(ROOT, '.cache', '_unverified-selftest.md');
  // ★一時ファイルを作らずに済むよう、parse を文字列で叩けるようにしてある。
  const lines = fx.split('\n');
  const guesses = new Map(); const verdicts = new Map();
  for (let i = 0; i < lines.length; i++) {
    let m = GUESS_RE.exec(lines[i]);
    if (m) { guesses.set(m[1].trim(), { conf: CONF_RE.exec(m[2])?.[1] ?? '(無記入)', how: HOW_RE.exec(m[2])?.[1] ?? '(無記入)', text: m[2].split('|')[0].trim(), file: path, line: i + 1 }); continue; }
    m = VERDICT_RE.exec(lines[i]);
    if (m) verdicts.set(m[1].trim(), { outcome: OUT_RE.exec(m[2])?.[1] ?? '(読めない)', text: '', file: path, line: i + 1 });
  }
  const checks = [
    ['GUESS を 4 件拾う', guesses.size === 4],
    ['VERDICT を 4 件拾う', verdicts.size === 4],
    ['★ を前置しても拾う', guesses.has('#16-a') && verdicts.has('#16-a')],
    ['`- ` を前置しても拾う', guesses.has('Z1')],
    ['確度を読む', guesses.get('Y19e').conf === '高' && guesses.get('Z1').conf === '低'],
    ['検算を読む', guesses.get('Y19e').how === '型' && guesses.get('#16-a').how === '原文'],
    ['当否を読む', verdicts.get('#16-a').outcome === '外れ' && verdicts.get('Z1').outcome === '半分'],
    ['書式でない行を拾わない', !guesses.has('を含むが書式ではない')],
    ['VERDICT 待ちを 1 件と数える', [...guesses.keys()].filter((k) => !verdicts.has(k)).length === 1],
    ['孤児 VERDICT を 1 件と数える', [...verdicts.keys()].filter((k) => !guesses.has(k)).length === 1],
  ];
  // ---- ★★★足した器具(統計)の較正 ---------------------------------------
  // ★既知値は scipy 1.17.1 の fisher_exact と突き合わせてある(メタ第 22 回)。
  const near = (x, y, e = 5e-7) => Math.abs(x - y) < e;
  const F = fisher2x2;
  checks.push(['Fisher 厳密: 紅茶 [[3,1],[1,3]] = 0.485714', near(F(3, 1, 1, 3), 0.4857142857)]);
  checks.push(['Fisher 厳密: [[3,9],[11,3]] = 0.016234', near(F(3, 9, 11, 3), 0.0162336788)]);
  checks.push(['Fisher 厳密: [[3,6],[11,13]] = 0.698198', near(F(3, 6, 11, 13), 0.6981979978)]);
  checks.push(['Fisher: 転置しても同じ', near(F(3, 6, 11, 13), F(3, 11, 6, 13))]);
  checks.push(['Fisher: 群を入れ替えても同じ', near(F(3, 6, 11, 13), F(11, 13, 3, 6))]);
  checks.push(['Fisher: 率が同じなら 1', near(F(2, 2, 3, 3), 1)]);
  {
    // ★対数版(外挿用)が厳密版と食い違わないこと。★ここが崩れると「あと何件」が嘘になる。
    let worst = 0, seed = 12345;
    const rnd = (n) => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed % n; };
    for (let i = 0; i < 300; i++) {
      const a = rnd(12), b = rnd(12), c = rnd(12), d = rnd(12);
      if (a + b === 0 || c + d === 0 || a + c === 0 || b + d === 0) continue;
      worst = Math.max(worst, Math.abs(F(a, b, c, d) - fisherLog(a, b, c, d)));
    }
    checks.push([`対数版が厳密版と一致(300 表、最大差 ${worst.toExponential(1)})`, worst < 1e-9]);
  }
  checks.push(['Bonferroni: p=0.02 は m=4 で言えない', !(0.02 * 4 < ALPHA)]);
  checks.push(['Bonferroni: p=0.02 は m=1 なら言える', 0.02 * 1 < ALPHA]);
  checks.push(['needMore: 率が同じなら「届かない」', needMore([[2, 2], [3, 3]], 4).impossible === true]);
  checks.push(['needMore: もう有意なら 1 倍', needMore([[20, 0], [0, 20]], 4).mult === 1]);
  checks.push(['needMore: 弱い差は倍数が要る', needMore([[3, 6], [11, 13]], 4).mult === 16]);
  checks.push(['束の鍵: B5-a → B5', bundleOf('B5-a') === 'B5']);
  checks.push(['束の鍵: 第1075波-3 → 第1075波', bundleOf('第1075波-3') === '第1075波']);
  {
    const mk = (n, o) => Array.from({ length: n }, (_, i) => ({ k: `X${i}-a`, v: { outcome: i < o ? '外れ' : '当たり' } }));
    const ov = (r) => r.v.outcome === '外れ' || r.v.outcome === '半分';
    // 束ごとに完全に揃っている(4 件束が 0/4 か 4/4)⇒ ICC は高い
    const tight = [];
    for (let b = 0; b < 6; b++) for (let i = 0; i < 4; i++) {
      tight.push({ k: `B${b}-${i}`, v: { outcome: b % 2 ? '外れ' : '当たり' } });
    }
    const cs = clusterStats(tight, ov);
    checks.push([`束: 6 束 24 件を数える`, cs.K === 6 && cs.N === 24]);
    checks.push([`束: 束ごと揃っていれば ICC≈1(${cs.icc === null ? "算出できず" : cs.icc.toFixed(2)})`, cs.icc !== null && cs.icc > 0.9]);
    checks.push(['束: 有効件数が件数より小さい', cs.neff != null && cs.neff < cs.N]);
    void mk;
  }
  {
    // ★★★回帰: **向きを断定する文を二度と出さない**(メタ第 22 回の主題)。
    //   ★確度=高 が全部当たり・低が全部外れという「向きが露骨な」標本を与えても、
    //   ★出力に「どちらがよく覆るか」を書いてはならない。
    const guesses = new Map(); const verdicts = new Map();
    for (let i = 0; i < 20; i++) {
      guesses.set(`H${i}-a`, { conf: '高', how: '型', text: '', file: '', line: 0 });
      verdicts.set(`H${i}-a`, { outcome: '当たり', text: '', file: '', line: 0 });
      guesses.set(`L${i}-a`, { conf: '低', how: 'なし', text: '', file: '', line: 0 });
      verdicts.set(`L${i}-a`, { outcome: '外れ', text: '', file: '', line: 0 });
    }
    const buf = []; const real = console.log;
    console.log = (...a) => buf.push(a.join(' '));
    try { report({ guesses, verdicts }); } finally { console.log = real; }
    const out = buf.join('\n');
    const banned = ['よく覆っている', '断定の方が', '断り付きの方が', '誘発している可能性'];
    const hit = banned.filter((w) => out.includes(w));
    checks.push([`★向きを断定する文を出さない${hit.length ? '(出た: ' + hit.join('/') + ')' : ''}`, hit.length === 0]);
    checks.push(['★露骨な差なら「言える」と出る', out.includes('言える(補正後)')]);
    checks.push(['★交絡の表を必ず出す', out.includes('交絡')]);
    checks.push(['★束(独立でない)を必ず出す', out.includes('独立でない')]);
  }

  {
    // ★★回帰: **多重比較の補正が効いていること**。
    //   ★確度=高 4 件が全部覆り・低 4 件が全部当たり ⇒ Fisher 両側 p = 0.0286。
    //   ★補正なし(α=0.05)なら「言える」。★補正あり(m=2 ⇒ 0.025)なら「言えない」。
    //   ⇒ この 1 件が、`* m` を落としたときに鳴る。
    const guesses = new Map(); const verdicts = new Map();
    for (let i = 0; i < 4; i++) {
      guesses.set(`H${i}-a`, { conf: '高', how: 'なし', text: '', file: '', line: 0 });
      verdicts.set(`H${i}-a`, { outcome: '外れ', text: '', file: '', line: 0 });
      guesses.set(`L${i}-a`, { conf: '低', how: 'なし', text: '', file: '', line: 0 });
      verdicts.set(`L${i}-a`, { outcome: '当たり', text: '', file: '', line: 0 });
    }
    const buf = []; const real = console.log;
    console.log = (...a) => buf.push(a.join(' '));
    try { report({ guesses, verdicts }); } finally { console.log = real; }
    const out = buf.join('\n');
    checks.push(['★補正: p=0.0286 は m=2 で「言えない」', !out.includes('言える(補正後)') && out.includes('言えない')]);
    checks.push(['★検定できない比較は m に数えない(m=2)', out.includes('0.05/2')]);
    checks.push(['★件数不足の比較はそう書く', out.includes('件数不足')]);
  }

  let ok = 0;
  for (const [label, pass] of checks) { console.log(`  ${pass ? 'ok ' : 'NG '} ${label}`); if (pass) ok++; }
  console.log(`\n  selftest: ${ok}/${checks.length} PASS`);
  return ok === checks.length;
}

if (flag('--selftest')) process.exit(selftest() ? 0 : 1);

const paths = opt('--file') ? [opt('--file')] : SOURCES;
const data = collect(paths);
for (const m of data.missing) console.log(`(記録が無い: ${m})`);
if (flag('--json')) {
  console.log(JSON.stringify({
    guesses: [...data.guesses].map(([k, v]) => ({ key: k, ...v })),
    verdicts: [...data.verdicts].map(([k, v]) => ({ key: k, ...v })),
  }, null, 1));
} else if (flag('--open')) {
  for (const [k, g] of data.guesses) {
    if (!data.verdicts.has(k)) console.log(`${k}\t確度=${g.conf}\t検算=${g.how}\t${g.text}`);
  }
} else {
  report(data);
}
