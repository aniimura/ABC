#!/usr/bin/env node
// agent-timing.mjs —— M6(セッションの実挙動が観測できない)の観測点。
//
// ★発見(メタ第 23 回): agent の完了通知 `<usage><subagent_tokens>/<tool_uses>/<duration_ms></usage>` は
//   会話ログ `~/.claude/projects/<slug>/<session>.jsonl` に **そのまま残っている**。
//   さらに `<session>/subagents/agent-*.jsonl` に **子 agent の全 tool 呼び出しと brief** が残っている。
//   ⇒ 本体が通知から手で書き写す必要は無い。この道具が読み出す。
//
// 使い方:
//   node tools/agent-timing.mjs --list                  … 通知を新しい順に出す
//   node tools/agent-timing.mjs --stats                 … 1 ノードあたりの時間の分布
//   node tools/agent-timing.mjs --explain               … 何が時間を説明するか(族は固定)
//   node tools/agent-timing.mjs --estimate              … brief の見積 行数 と 実測 行数
//   node tools/agent-timing.mjs --cost                  … 自己申告 COST[…] を実測に突き合わせる
//   node tools/agent-timing.mjs --selftest              … 自己検査
//   共通: --type lean-prover  --day 2026-09-07  --name <正規表現>  --limit N  --json
//         --ms 919271,1086819,…  … duration_ms を並べて過去の表をそのまま再現する
//
// ★統計の作法(M90 に倣う):
//   - 比較の族 FAMILY を **データを見る前に**コードへ固定する。
//   - 判定は「言える / 言えない」だけ。★率や相関の**向きを断定しない**。
//   - 多重比較(Holm)で補正する。
//   - 交絡を必ず出す。標本が独立でないとき(束)は ICC / DEFF を出す。

// ★★★自己申告の書式(メタ第 24 回。M96 が「今の記録では測れない」と書いた穴を塞ぐ)
//
//     COST[<持ち場>]: <安|並|高>  — <一言>
//     COST[<鍵>]: <安|並|高> | 持ち場=<agent に配った呼び名>  — <一言>     ★鍵が違うときの橋
//
//   - `<持ち場>` は **agent に配ったときの呼び名**(`Y19e` `B2+B3` `Λ12` など)。
//     ★`GUESS/VERDICT`(autonomy-policy §4.7)と**同じ字面**にすること。ここだけが対応の鍵。
//   - ★★**枝番を付けない**(`LT-a` ではなく `LT`)。費用は**持ち場 1 つに 1 行**である
//     (見当は 1 持ち場に何本あってもよいが、かかった時間は 1 つしかない)。
//   - ★★**実測は「鍵が持ち場名と一致すること」に依存する。**★第 24 回の実測:
//     いま `decisions-pending.md` にある鍵 13 個のうち、**agent の呼び名に当たるのは 2 個だけ**
//     (`段1` `C49` `L9` `TL` `L12` `HS` `AE` `GC` `LT` `IR` `SR` は当たらない。略号だから)。
//     ⇒ 鍵を揃えられないときは `| 持ち場=局所Tate双対性の道` と**明示の橋**を書くこと。
//     ⇒ ★当たらなかった申告は `--cost` が**名指しで**出す(黙って落とさない)。
//   - 意味は「**見積に対して**思ったより安かった / 見積どおり / 思ったより高かった」。
//     ★`安` は「行数が少なかった」ではない(それは `--estimate` が既に測っている)。
//   - **報告が返ったときに 1 行**書く(`VERDICT` と同じ間合い)。書かなくてもゲートは落ちない。
//   - 書き場所は `ResearchPaper/decisions-pending.md`(`GUESS/VERDICT` と同じ本)。
//
// ★なぜ語の検索ではだめか(M96 の実測): 「安い」系の語を含む節は 695 節中 21 節あるのに、
//   測れた 21 件の持ち場に**紐づくのは 1 件だけ**だった。★節は持ち場に対応していない。
//   ⇒ **鍵付きの 1 行**でなければ分母が立たない。
//
// ★なぜ `unverified.mjs` ではなく**ここ**に置いたか(3 つ)
//   1. ★M96 が測れなかったのは「申告と**実測**の一致」である。実測(所要時間 / 行数 / 見積)を
//      持っているのは**この道具だけ**。unverified 側に置くと「申告が N 件ある」までしか言えず、
//      穴は塞がらない。
//   2. ★unverified.mjs は **Bonferroni の族を固定してある**(M90)。同じ本の中に別種の欄を
//      足すと m の意味が濁り、**既に出ている判定が動く**。★族を守るために分ける。
//   3. ★unverified.mjs は import すると即座に集計が走る作り(main の番人が無い)。
//      読み手として使えないので、写すか作り替えるかになる。★どちらも今回の目的より高い。
//   ⇒ ★**書式は GUESS/VERDICT と揃え、置き場所も同じ本にし、数える道具だけ分ける。**

import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';

// ════════════════════════════════════════════════════════════════════
// 0. 比較の族 —— ★データを見る前に固定する(本体の brief が挙げた候補そのもの)
// ════════════════════════════════════════════════════════════════════
export const FAMILY = [
  // ★★族の正典は `lines`(いまの wc -l)。`linesWritten`(当時の Write 本文)は**族に入れない**。
  //   ── メタ第 25 回が測って決めた。理由は ρ の大小ではない(★下記のとおり ρ では決められない)。
  //   (1) ★ρ では決められない: n=191 で lines ρ=0.260(n=131) / linesWritten ρ=0.263(n=81)。
  //       ★メタ第 24 回の「0.758 対 0.511」(n=50)の差は **n を増やしたら消えた**。
  //       ⇒ ρ を根拠に選ぶのは「結果を見てから選ぶ」ことでもあり、そもそも差が無い。
  //   (2) ★どちらが何件間違うかで決めた。両方取れる 81 件のうち値が違うのは 53 件。
  //       ・そのうち **40 件は その agent 自身が最後の Write の後に Edit した**ぶん。
  //         ⇒ ここでは `linesWritten` のほうが**間違い**(自分の仕事を数え落とす)。
  //       ・自分は Write 以後さわっていないのにずれるのは 13 件。うち 9 件は **2 行以下**の差。
  //         ⇒ 「別の agent に後から書き換えられた」で実質的に効くのは **4 件 / 81 = 5%**。
  //       ⇒ 誤差 5% の `lines` を採り、誤差 49% の `linesWritten` は採らない。
  //   (3) ★取りこぼしも lines が有利: `linesWritten` は Write を使わなかった agent で欠測し、
  //       n が 131 → 81 に落ちる。★さらに族を差し替えると Holm の階段が早く止まり、
  //       **cores と checkFails の判定が「言える」→「言えない」に落ちる**(副作用が大きい)。
  //   (4) ★本当の正解はどちらでもない:「完了時点の行数」= 最後の Write 本文 + その後の自分の
  //       Edit を当てた行数。★これは識別もでき取りこぼしも無い。★未実装(meta-backlog M111)。
  //       ★実装したら**測ってから**差し替えること。今の判定を根拠に先に差し替えないこと。
  { key: 'lines',      label: '行数(wc -l)' },
  { key: 'toolUses',   label: 'tool_uses' },
  { key: 'tokens',     label: 'subagent_tokens' },
  { key: 'cores',      label: '抽象核の本数' },
  { key: 'leanChecks', label: 'lean_check の回数' },
  { key: 'checkFails', label: 'lean_check が失敗した回数(≒ 一発で通らなかった)' },
  // ★★事前登録(メタ第 23 回 M96 が「族の外」と自己申告し、★次の波で確かめると書いた欄)。
  //   ★第 24 回が**データを見る前に**ここへ移した。★向きは印字しない。
  { key: 'estMid',     label: '見積中点(brief の「見積 N–M 行」)' },
];

/** ★欄ごとの最小件数。これ未満なら検定しない(★データを見る前に決める)。 */
export const MIN_N = 8;

/**
 * ★★事前登録(メタ第 25 回)—— 自己申告 COST を FAMILY に足してよいか。
 * ★n=5 の時点では「言えない」ですらなく **判定を出さない**(--cost が既にそう振る舞う)。
 * ★足す条件を**データを見る前に**ここへ固定する:
 *   (a) 突き合わせが MIN_N(8)件以上、かつ
 *   (b) ★水準が偏っていないこと —— 各水準に 3 件以上。
 * ★(b) が要る理由: 2026-09-07 時点の 5 件は 安:4 / 並:1 / 高:0 で、
 *   件数だけ 8 に達しても「全部 安」なら順位相関は定義できず、**分散ゼロの欄**になる。
 * ★足すときは m が 7 → 8 になり、既存の判定が Holm で罰される。★それを承知で 1 度だけ足すこと。
 */
export function costFamilyReady(levels) {
  const ls = Array.isArray(levels) ? levels : [];
  if (ls.length < MIN_N) return false;
  return COST_LEVELS.every(l => ls.filter(x => x === l).length >= 3);
}

/**
 * 族の 1 欄について「検定に使える行」を返す。
 * ★★欄ごとの complete-case(listwise ではない)。理由: `estMid` は brief に見積を
 *   書かなかった持ち場で欠測する。listwise だと **1 件の欠測で欄ごと落ちる**ので、
 *   事前登録した欄が「測れなかった」に化けてしまう。
 * ★欠測が欄によって違うので、**n を欄ごとに印字する**(そうしないと比較が読めない)。
 */
export function usableRows(rows, key) {
  return rows.filter((r) => Number.isFinite(r[key]));
}

// ════════════════════════════════════════════════════════════════════
// 1. 純関数 —— 通知の解析
// ════════════════════════════════════════════════════════════════════

/** task-notification の本文 1 件を解析する。usage が無ければ null。 */
export function parseNotification(content) {
  if (typeof content !== 'string') return null;
  const g = (re) => re.exec(content)?.[1];
  const d = g(/<duration_ms>(\d+)<\/duration_ms>/);
  if (d === undefined) return null;
  const summary = g(/<summary>([\s\S]*?)<\/summary>/) ?? '';
  return {
    taskId: g(/<task-id>([^<]*)<\/task-id>/) ?? '',
    toolUseId: g(/<tool-use-id>([^<]*)<\/tool-use-id>/) ?? '',
    status: g(/<status>([^<]*)<\/status>/) ?? '',
    summary,
    // `Agent "…" finished` の中身。二重引用符は content 側で素のまま入っている。
    name: /^Agent "([\s\S]*)" finished\s*$/.exec(summary.trim())?.[1] ?? summary.trim(),
    durationMs: Number(d),
    toolUses: Number(g(/<tool_uses>(\d+)<\/tool_uses>/) ?? NaN),
    tokens: Number(g(/<subagent_tokens>(\d+)<\/subagent_tokens>/) ?? NaN),
  };
}

/** brief(子 agent の最初の user message)から見積 行数 を読む。 */
export function parseEstimate(brief) {
  if (typeof brief !== 'string') return null;
  // 「見積 **400–700 行**」「見積 120–250 行（Y4）」「見積 **約 300 行**」等
  const m = /見積[^0-9]{0,12}(\d{2,4})\s*[-–—~〜]\s*(\d{2,4})\s*行/.exec(brief);
  if (m) return { lo: +m[1], hi: +m[2] };
  const s = /見積[^0-9]{0,12}(\d{2,4})\s*行/.exec(brief);
  if (s) return { lo: +s[1], hi: +s[1] };
  return null;
}

/** .lean の本文から「抽象核」の節を数える(節見出し `/-! ## … 抽象核 …`)。 */
export function countAbstractCores(src) {
  if (typeof src !== 'string') return 0;
  let n = 0;
  for (const line of src.split(/\r?\n/)) {
    if (/^\s*\/-!\s*#+\s.*抽象核/.test(line)) n++;
  }
  return n;
}

/** `COST[<持ち場>]: 安` を 1 行から読む。書式の**見本**(鍵に <> を含む)は数えない。 */
export const COST_RE = /^\s*(?:[★☆*\-\s]*)COST\[([^\]]+)\]\s*[:：]\s*(.*)$/;
export const COST_LEVELS = ['安', '並', '高'];
export function parseCostLine(line) {
  const m = COST_RE.exec(String(line ?? ''));
  if (!m) return null;
  const key = m[1].trim();
  if (/[<>]/.test(key)) return null;              // ★書式の見本を数えない(M71/unverified と同じ穴)
  const lv = /^\s*(安|並|高)(?![^\s—\-|]) ?/.exec(m[2]) ?? /^\s*(安|並|高)/.exec(m[2]);
  if (!lv) return null;
  const rest = m[2].slice(lv[0].length);
  // ★明示の橋: `| 持ち場=<agent の呼び名>`(GUESS の `確度=` `検算=` と同じ流儀)
  const slot = /持ち場\s*=\s*([^|—]+)/.exec(rest)?.[1]?.trim() || null;
  return {
    key,
    level: lv[1],
    slot,                                  // ★null なら key で突き合わせる
    note: rest.replace(/\|?\s*持ち場\s*=\s*[^|—]+/, '').replace(/^\s*[—\-|]\s*/, '').trim(),
  };
}

/** 記録(複数)から COST を集める。★同じ鍵は後の記述で上書きする。 */
export function readCosts(paths) {
  const out = new Map(); const missing = [];
  for (const p of paths) {
    if (!fs.existsSync(p)) { missing.push(p); continue; }
    const lines = fs.readFileSync(p, 'utf8').split(/\r?\n/);
    for (let i = 0; i < lines.length; i++) {
      const c = parseCostLine(lines[i]);
      if (c) out.set(c.key, { ...c, file: p, line: i + 1 });
    }
  }
  return { costs: out, missing };
}

/**
 * 申告の鍵と agent の持ち場名を突き合わせる。
 * ★前方一致だけだと `Y19` が `Y19e …` を拾ってしまうので、**鍵の直後は英数字でない**ことを要る。
 */
export function matchesKey(description, key) {
  const d = String(description ?? '').trim();
  const k = String(key ?? '').trim();
  if (!k) return false;
  if (d === k) return true;
  if (!d.startsWith(k)) return false;
  return !/[0-9A-Za-z]/.test(d.slice(k.length, k.length + 1));
}

// ════════════════════════════════════════════════════════════════════
// 2. 純関数 —— 統計
// ════════════════════════════════════════════════════════════════════

/** type-7 分位点(numpy / R の既定)。xs は昇順でなくてよい。 */
export function quantile(xs, q) {
  const a = [...xs].sort((x, y) => x - y);
  if (a.length === 0) return NaN;
  if (a.length === 1) return a[0];
  const h = (a.length - 1) * q;
  const lo = Math.floor(h), hi = Math.ceil(h);
  return a[lo] + (h - lo) * (a[hi] - a[lo]);
}

/** 平均順位(同順位は平均)。 */
export function ranks(xs) {
  const idx = xs.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
  const r = new Array(xs.length);
  let i = 0;
  while (i < idx.length) {
    let j = i;
    while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
    const avg = (i + j) / 2 + 1;
    for (let k = i; k <= j; k++) r[idx[k][1]] = avg;
    i = j + 1;
  }
  return r;
}

export function pearson(x, y) {
  const n = x.length;
  if (n < 2) return NaN;
  let mx = 0, my = 0;
  for (let i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
  mx /= n; my /= n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    const a = x[i] - mx, b = y[i] - my;
    sxy += a * b; sxx += a * a; syy += b * b;
  }
  if (sxx === 0 || syy === 0) return NaN;
  return sxy / Math.sqrt(sxx * syy);
}

export function spearman(x, y) { return pearson(ranks(x), ranks(y)); }

/** 決定的な RNG(mulberry32)。seed を固定するので再現する。 */
export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * Spearman の両側 p を並べ替え検定で出す(|rho| 基準)。
 * ★決定的(seed 固定)。★p は (超えた数 + 1)/(反復 + 1) で 0 にならない。
 */
export function permP(x, y, iters = 20000, seed = 20260907) {
  const rx = ranks(x), ry = ranks(y);
  const r0 = Math.abs(pearson(rx, ry));
  if (!Number.isFinite(r0)) return NaN;
  const rnd = mulberry32(seed);
  const arr = ry.slice();
  let ge = 0;
  for (let it = 0; it < iters; it++) {
    for (let i = arr.length - 1; i > 0; i--) {
      const j = Math.floor(rnd() * (i + 1));
      const t = arr[i]; arr[i] = arr[j]; arr[j] = t;
    }
    if (Math.abs(pearson(rx, arr)) >= r0 - 1e-12) ge++;
  }
  return (ge + 1) / (iters + 1);
}

/** Holm-Bonferroni の調整後 p(単調化して 1 で頭打ち)。 */
export function holm(ps) {
  const m = ps.length;
  const ord = ps.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0]);
  const adj = new Array(m);
  let run = 0;
  for (let k = 0; k < m; k++) {
    const v = (m - k) * ord[k][0];
    run = Math.max(run, v);
    adj[ord[k][1]] = Math.min(1, run);
  }
  return adj;
}

/** 偏 Spearman(順位に対する 1 変量の偏相関)。 */
export function partial(rxy, rxz, ryz) {
  const d = Math.sqrt((1 - rxz * rxz) * (1 - ryz * ryz));
  if (d === 0) return NaN;
  return (rxy - rxz * ryz) / d;
}

/** 一元配置の ICC(1) と設計効果。groups: 値の配列の配列。 */
export function icc1(groups) {
  const gs = groups.filter(g => g.length > 0);
  const k = gs.length;
  const N = gs.reduce((s, g) => s + g.length, 0);
  if (k < 2 || N <= k) return { icc: NaN, deff: NaN, nEff: N, k, N };
  const grand = gs.flat().reduce((s, v) => s + v, 0) / N;
  let msb = 0, msw = 0;
  for (const g of gs) {
    const m = g.reduce((s, v) => s + v, 0) / g.length;
    msb += g.length * (m - grand) ** 2;
    for (const v of g) msw += (v - m) ** 2;
  }
  msb /= (k - 1);
  msw /= (N - k);
  const n0 = (N - gs.reduce((s, g) => s + g.length ** 2, 0) / N) / (k - 1);
  const icc = (msb - msw) / (msb + (n0 - 1) * msw);
  const mbar = N / k;
  const deff = 1 + (mbar - 1) * Math.max(0, icc);
  return { icc, deff, nEff: N / deff, k, N, n0 };
}

// ════════════════════════════════════════════════════════════════════
// 3. 読み出し
// ════════════════════════════════════════════════════════════════════

export function projectRoots(slugPrefix = 'D--Math-ABC3') {
  const base = path.join(os.homedir(), '.claude', 'projects');
  if (!fs.existsSync(base)) return [];
  return fs.readdirSync(base)
    .filter(d => d === slugPrefix || d.startsWith(slugPrefix + '-'))
    .map(d => path.join(base, d));
}

/** 通知 + 子 agent の meta を突き合わせて 1 件 = 1 agent にする。 */
export function collect(opts = {}) {
  const roots = opts.roots ?? projectRoots();
  const byToolUse = new Map();   // toolUseId -> record
  const metaByToolUse = new Map();
  const subagentFile = new Map();

  for (const root of roots) {
    if (!fs.existsSync(root)) continue;
    // (a) 子 agent の meta / 本文
    for (const sess of fs.readdirSync(root)) {
      const dir = path.join(root, sess, 'subagents');
      if (!fs.existsSync(dir)) continue;
      for (const f of fs.readdirSync(dir)) {
        if (!f.endsWith('.meta.json')) continue;
        let m; try { m = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8')); } catch { continue; }
        if (!m.toolUseId) continue;
        metaByToolUse.set(m.toolUseId, m);
        const jf = path.join(dir, f.replace(/\.meta\.json$/, '.jsonl'));
        if (fs.existsSync(jf)) subagentFile.set(m.toolUseId, jf);
      }
    }
    // (b) 親セッションの通知
    for (const f of fs.readdirSync(root)) {
      if (!f.endsWith('.jsonl')) continue;
      const text = fs.readFileSync(path.join(root, f), 'utf8');
      for (const line of text.split('\n')) {
        if (!line.includes('duration_ms')) continue;
        let o; try { o = JSON.parse(line); } catch { continue; }
        const content = typeof o.content === 'string' ? o.content : null;
        if (!content) continue;
        const p = parseNotification(content);
        if (!p) continue;
        const prev = byToolUse.get(p.toolUseId);
        // 同じ agent が 2 回以上通知しうる(resume)。★最後の(= 最大の)duration を採る。
        if (!prev || p.durationMs > prev.durationMs) {
          byToolUse.set(p.toolUseId, { ...p, ts: o.timestamp ?? null });
        }
      }
    }
  }

  const out = [];
  for (const [tu, rec] of byToolUse) {
    const m = metaByToolUse.get(tu) ?? {};
    out.push({
      ...rec,
      agentType: m.agentType ?? '?',
      description: m.description ?? rec.name,
      spawnDepth: m.spawnDepth ?? null,
      subagentFile: subagentFile.get(tu) ?? null,
      day: (rec.ts ?? '').slice(0, 10),
    });
  }
  out.sort((a, b) => (a.ts ?? '') < (b.ts ?? '') ? -1 : 1);
  return out;
}

/** 子 agent の本文から共変量を読む(重いので必要なときだけ)。 */
export function covariates(rec, repoRoot) {
  const zero = { leanChecks: 0, checkFails: 0, writes: 0, edits: 0, bashes: 0,
                 briefChars: 0, est: null, estMid: NaN, files: [], lines: NaN, cores: NaN, file: null,
                 linesWritten: NaN, coresWritten: NaN };  // ★linesWritten は「最後の Write の本文」であって完了時の行数ではない(その後の Edit を含まない)
  if (!rec.subagentFile || !fs.existsSync(rec.subagentFile)) return zero;
  const text = fs.readFileSync(rec.subagentFile, 'utf8');
  const lines = text.split('\n');
  const c = { ...zero, files: [] };
  const fileHits = new Map();
  const written = new Map();
  let pendingCheck = new Set();
  for (const line of lines) {
    if (!line) continue;
    let o; try { o = JSON.parse(line); } catch { continue; }
    const msg = o.message;
    if (!msg) continue;
    const content = Array.isArray(msg.content) ? msg.content : [];
    for (const b of content) {
      if (b.type === 'tool_use') {
        const n = b.name || '';
        if (n === 'Write') c.writes++;
        else if (n === 'Edit' || n === 'MultiEdit') c.edits++;
        else if (n === 'Bash') c.bashes++;
        if (/lean_check$/.test(n)) { c.leanChecks++; pendingCheck.add(b.id); }
        const fp = b.input && (b.input.file_path || b.input.filePath);
        if (typeof fp === 'string' && /\.lean$/i.test(fp) && /Found[\\/]/.test(fp)) {
          fileHits.set(fp, (fileHits.get(fp) || 0) + 1);
          // ★完了時点の行数。★ファイルは後から別の agent に書き換えられるので、
          //   「いま wc -l した行数」ではなく **その agent が書いた最後の本文** を採る。
          if (n === 'Write' && typeof b.input.content === 'string') {
            written.set(fp, b.input.content);
          }
        }
      } else if (b.type === 'tool_result') {
        if (pendingCheck.has(b.tool_use_id)) {
          pendingCheck.delete(b.tool_use_id);
          const s = typeof b.content === 'string' ? b.content : JSON.stringify(b.content ?? '');
          if (/error|✗|failed/i.test(s)) c.checkFails++;
        }
      }
    }
    if (c.briefChars === 0 && msg.role === 'user' && typeof msg.content === 'string') {
      c.briefChars = msg.content.length;
      c.est = parseEstimate(msg.content);
      c.estMid = c.est ? (c.est.lo + c.est.hi) / 2 : NaN;
    }
  }
  c.files = [...fileHits.entries()].sort((a, b) => b[1] - a[1]).map(e => e[0]);
  if (c.files.length) {
    c.file = c.files[0];
    const w = written.get(c.file);
    if (typeof w === 'string') {
      c.linesWritten = w.split('\n').length - (w.endsWith('\n') ? 1 : 0);
      c.coresWritten = countAbstractCores(w);
    }
    const abs = path.isAbsolute(c.file) ? c.file : path.join(repoRoot, c.file);
    const cand = [abs, path.join(repoRoot, c.file.replace(/^.*?lean[\\/]/, 'lean/'))];
    for (const p2 of cand) {
      if (fs.existsSync(p2)) {
        const src = fs.readFileSync(p2, 'utf8');
        c.lines = src.split('\n').length - (src.endsWith('\n') ? 1 : 0);
        c.cores = countAbstractCores(src);
        break;
      }
    }
  }
  return c;
}

// ════════════════════════════════════════════════════════════════════
// 4. 印字
// ════════════════════════════════════════════════════════════════════
const min1 = (ms) => (ms / 60000).toFixed(1);
const pad = (s, n) => String(s).padStart(n);
const padr = (s, n) => { s = String(s); return s + ' '.repeat(Math.max(0, n - width(s))); };
function width(s) { let w = 0; for (const ch of String(s)) w += /[　-鿿＀-￯]/.test(ch) ? 2 : 1; return w; }

function describe(vals) {
  return {
    n: vals.length,
    min: Math.min(...vals), q1: quantile(vals, 0.25), med: quantile(vals, 0.5),
    q3: quantile(vals, 0.75), max: Math.max(...vals),
    mean: vals.reduce((s, v) => s + v, 0) / vals.length,
    sum: vals.reduce((s, v) => s + v, 0),
  };
}

function cmdStats(recs) {
  const d = describe(recs.map(r => r.durationMs));
  console.log(`## 1 ノードあたりの時間の分布 —— n = ${d.n}`);
  console.log('');
  console.log('  統計量        ミリ秒        分');
  const row = (k, v) => console.log(`  ${padr(k, 12)}${pad(Math.round(v), 10)}   ${pad(min1(v), 7)}`);
  row('最小', d.min); row('第1四分位', d.q1); row('中央値', d.med);
  row('第3四分位', d.q3); row('最大', d.max); row('平均', d.mean);
  console.log(`  ${padr('四分位範囲', 12)}${pad(Math.round(d.q3 - d.q1), 10)}   ${pad(min1(d.q3 - d.q1), 7)}`);
  console.log(`  ${padr('合計', 12)}${pad(Math.round(d.sum), 10)}   ${pad(min1(d.sum), 7)}`);
  console.log('');
  // 実時間の窓(重なりを含む) —— 直列/並行の差を出す
  const ts = recs.map(r => r.ts).filter(Boolean).sort();
  if (ts.length >= 2) {
    const span = new Date(ts[ts.length - 1]) - new Date(ts[0]);
    console.log(`  最初の完了 ${ts[0]} / 最後の完了 ${ts[ts.length - 1]}`);
    console.log(`  ★完了時刻の幅 ${min1(span)} 分 に対し agent 時間の合計は ${min1(d.sum)} 分`);
    console.log(`    ⇒ 重なり率 = 合計 / 幅 = ${(d.sum / span).toFixed(2)} 倍(1.0 なら直列)`);
    console.log(`    ⇒ ★「幅 / 件数」= ${min1(span / d.n)} 分/ノード(★これが 37.2 分と同じ定義の量)`);
  }
  return d;
}

function fmtP(p) { return Number.isFinite(p) ? p.toFixed(4) : '  ——  '; }

function cmdExplain(rows) {
  console.log('## ★何が時間を説明するか(★向きは断定しない)');
  console.log('');
  console.log('  ★比較の族は**データを見る前に**コードへ固定してある(FAMILY)。');
  const usable = [], dropped = [];
  for (const f of FAMILY) {
    const sub = usableRows(rows, f.key);
    if (sub.length >= MIN_N && new Set(sub.map(r => r[f.key])).size > 1) usable.push({ f, sub });
    else dropped.push({ f, sub });
  }
  console.log(`  検定できるのは ${usable.length} 本 / 族 ${FAMILY.length} 本 ⇒ Holm で補正する(m = ${usable.length})`);
  for (const { f, sub } of dropped) {
    const miss = rows.length - sub.length;
    const why = sub.length < MIN_N ? `件数不足 n=${sub.length} < ${MIN_N}` : '全件同値';
    console.log(`  ★測れなかった: ${f.label}(欠測 ${miss} / ${rows.length}、${why})`);
  }
  console.log('');
  const raw = usable.map(({ f, sub }) => {
    const yy = sub.map(r => r.durationMs);
    const xx = sub.map(r => r[f.key]);
    return { f, n: sub.length, rho: spearman(xx, yy), p: permP(xx, yy) };
  });
  const adj = holm(raw.map(r => r.p));
  console.log('  説明変数                                        n   Spearman ρ    並べ替え p   Holm 後   判定');
  raw.forEach((r, i) => {
    const verdict = adj[i] < 0.05 ? '★言える(補正後)' : '言えない';
    console.log(`  ${padr(r.f.label, 42)} ${pad(r.n, 4)} ${pad(r.rho.toFixed(3), 8)}   ${pad(fmtP(r.p), 10)}  ${pad(fmtP(adj[i]), 8)}   ${verdict}`);
  });
  console.log('');
  console.log('  ★n が欄ごとに違うのは **欄ごとの complete-case** だから(欠測の理由が欄で違う)。');
  console.log('');
  console.log('  ★「言えない」は「差が無い」ではない。★ρ の符号を根拠に向きを主張しないこと。');
  console.log('');
  // 交絡
  console.log('## ★交絡(説明変数どうしの Spearman ρ)');
  console.log('');
  const keys = usable.map(u => u.f.key);
  const pairRho = (a, b) => {
    const s2 = rows.filter(r => Number.isFinite(r[a]) && Number.isFinite(r[b]));
    return spearman(s2.map(r => r[a]), s2.map(r => r[b]));
  };
  console.log('  ' + padr('', 12) + keys.map(k => pad(k.slice(0, 10), 11)).join(''));
  for (const a of keys) {
    let line = '  ' + padr(a.slice(0, 12), 12);
    for (const b of keys) {
      const v = a === b ? 1 : pairRho(a, b);
      line += pad(Number.isFinite(v) ? v.toFixed(3) : '——', 11);
    }
    console.log(line);
  }
  console.log('');
  // 偏相関(族の外。断定しない)
  if (keys.includes('toolUses') && keys.includes('tokens')) {
    const s3 = rows.filter(r => Number.isFinite(r.toolUses) && Number.isFinite(r.tokens));
    const y = s3.map(r => r.durationMs);
    const t = s3.map(r => r.toolUses), k = s3.map(r => r.tokens);
    const rtk = spearman(t, k), ryt = spearman(y, t), ryk = spearman(y, k);
    console.log(`  ★tool_uses と tokens の ρ = ${rtk.toFixed(3)}`);
    console.log(`    duration~tool_uses を tokens で偏らせる  ρ = ${partial(ryt, ryk, rtk).toFixed(3)}`);
    console.log(`    duration~tokens を tool_uses で偏らせる  ρ = ${partial(ryk, ryt, rtk).toFixed(3)}`);
    console.log('    ★これは族の外(検定していない)。★交絡の大きさを見るためだけに出す。');
  }
  console.log('');
  return { raw, adj };
}

function bundleKey(desc) {
  // 「B2+B3 …」「段1c …」「Λ9 …」等、先頭の記号列を束の鍵にする(枝番は落とす)
  const head = String(desc).trim().split(/[\s—]/)[0].split('+')[0];
  return head.replace(/[0-9]+[a-z]*$/, '') || head;
}

function cmdIcc(rows) {
  const g = new Map();
  for (const r of rows) {
    const k = bundleKey(r.description);
    if (!g.has(k)) g.set(k, []);
    g.get(k).push(r.durationMs);
  }
  const groups = [...g.values()];
  const s = icc1(groups);
  console.log('## ★標本は独立か —— 束(持ち場の頭)ごとのばらつき');
  console.log('');
  console.log('  ' + [...g.entries()].sort((a, b) => b[1].length - a[1].length)
    .map(([k, v]) => `${k}:${v.length}`).join(' '));
  console.log(`  束 ${s.k} / 件数 ${s.N} / ICC = ${Number.isFinite(s.icc) ? s.icc.toFixed(3) : '——'}` +
              ` / DEFF = ${Number.isFinite(s.deff) ? s.deff.toFixed(2) : '——'}` +
              ` ⇒ 有効件数 ≒ ${Number.isFinite(s.nEff) ? s.nEff.toFixed(0) : '——'}`);
  console.log('  ★束が 1 件ずつなら ICC は意味を持たない(そのときは「粗い」と読むこと)。');
  console.log('');
  return s;
}

function cmdEstimate(rows) {
  const withEst = rows.filter(r => r.est && Number.isFinite(r.lines));
  console.log('## ★見積(brief の「見積 N–M 行」)と実測(wc -l)');
  console.log('');
  if (withEst.length === 0) { console.log('  ★見積を書いた brief が 0 件。測れない。'); return null; }
  console.log('  持ち場                                     見積      実測   実測/見積中点   範囲内');
  let inRange = 0;
  for (const r of withEst) {
    const mid = (r.est.lo + r.est.hi) / 2;
    const ok = r.lines >= r.est.lo && r.lines <= r.est.hi;
    if (ok) inRange++;
    console.log(`  ${padr(String(r.description).slice(0, 34), 36)} ${pad(r.est.lo + '-' + r.est.hi, 9)} ${pad(r.lines, 7)}   ${pad((r.lines / mid).toFixed(2), 10)}    ${ok ? '○' : '×'}`);
  }
  const ratio = withEst.map(r => r.lines / ((r.est.lo + r.est.hi) / 2));
  console.log('');
  console.log(`  n = ${withEst.length} / 範囲内 ${inRange} 件(${(100 * inRange / withEst.length).toFixed(0)}%)`);
  console.log(`  実測/見積中点: 中央値 ${quantile(ratio, 0.5).toFixed(2)} / 四分位 ${quantile(ratio, 0.25).toFixed(2)}–${quantile(ratio, 0.75).toFixed(2)} / 範囲 ${Math.min(...ratio).toFixed(2)}–${Math.max(...ratio).toFixed(2)}`);
  console.log('  ★これは「行数の見積」であって「かかった時間の見積」ではない。');
  return { n: withEst.length, inRange, ratio };
}

function cmdCost(rows, paths) {
  const { costs, missing } = readCosts(paths);
  console.log('## ★自己申告(COST[<持ち場>]: 安|並|高)と実測');
  console.log('');
  for (const m of missing) console.log(`  (記録が無い: ${m})`);
  console.log(`  読んだ記録: ${paths.filter(p => !missing.includes(p)).join(' / ') || '(無し)'}`);
  console.log(`  申告 ${costs.size} 件 / 実測できる持ち場 ${rows.length} 件`);
  if (costs.size === 0) {
    console.log('');
    console.log('  ★申告が 0 件。**分母が立たないので何も言わない。**');
    console.log('  ★書式:  COST[<持ち場>]: <安|並|高>  — <一言>   (ResearchPaper/decisions-pending.md)');
    console.log('  ★「安い/高い」の語を本文から拾う数え方は M96 で**測れない**と分かっている');
    console.log('    (695 節中 21 節が語を含むのに、持ち場に紐づいたのは 1 件だけだった)。');
    return { n: 0, joined: 0, costs };
  }
  const joined = [];
  const unjoined = [];
  for (const [k, c] of costs) {
    const hits = rows.filter(r => matchesKey(r.description, c.slot ?? k));
    if (hits.length === 0) { unjoined.push([k, c]); continue; }
    for (const r of hits) joined.push({ k, c, r });
  }
  console.log('');
  console.log('  持ち場                              申告   所要(分)   実測行   見積中点  実測/見積');
  for (const j of joined) {
    const mid = j.r.est ? (j.r.est.lo + j.r.est.hi) / 2 : NaN;
    console.log(`  ${padr(String(j.r.description).slice(0, 30), 32)} ${padr(j.c.level, 5)} ${pad(min1(j.r.durationMs), 8)} ${pad(Number.isFinite(j.r.lines) ? j.r.lines : '——', 8)} ${pad(Number.isFinite(mid) ? mid : '——', 10)} ${pad(Number.isFinite(mid) && Number.isFinite(j.r.lines) ? (j.r.lines / mid).toFixed(2) : '——', 9)}`);
  }
  console.log('');
  const tally = COST_LEVELS.map(l => `${l}:${joined.filter(j => j.c.level === l).length}`).join(' / ');
  console.log(`  突き合わせ ${joined.length} 件(${tally})/ 実測に当たらなかった申告 ${unjoined.length} 件`);
  for (const [k, c] of unjoined) console.log(`    ・${k}(${c.level}) —— ${c.file}:${c.line}`);
  if (unjoined.length) {
    console.log('    ⇒ ★鍵が agent の呼び名と違うだけなら `| 持ち場=<呼び名>` を 1 つ足せば繋がる。');
  }
  const usable = joined.filter(j => j.r.est && Number.isFinite(j.r.lines));
  console.log('');
  console.log(`  ★見積が読めて実測もある組 ${usable.length} 件。`);
  const levels = usable.map(j => j.c.level);
  const tally2 = COST_LEVELS.map(l => `${l}:${levels.filter(x => x === l).length}`).join(' / ');
  if (!costFamilyReady(levels)) {
    console.log(`  ★★族に入れる条件を満たさない(内訳 ${tally2})。**判定を出さない。★向きも書かない。**`);
    console.log(`    ⇒ 事前登録した条件は 2 つ: (a) ${MIN_N} 件以上 (b) ★各水準に 3 件以上。`);
    console.log('    ⇒ ★(b) が要るのは、件数だけ足りても「全部 安」なら分散ゼロで順位相関が定義できないため。');
    console.log('    ⇒ 両方満たしたときに **1 度だけ**「申告 × 実測/見積」を族に足すこと(m が 7 → 8 になる)。');
  } else {
    const lv = usable.map(j => COST_LEVELS.indexOf(j.c.level));
    const rt = usable.map(j => j.r.lines / ((j.r.est.lo + j.r.est.hi) / 2));
    const p = permP(lv, rt);
    console.log(`  申告(安=0/並=1/高=2)と 実測/見積 の Spearman ρ = ${spearman(lv, rt).toFixed(3)} / 並べ替え p = ${fmtP(p)}`);
    console.log('  ★これは**族の外**(この道具が初めて出す欄)。★次の波で FAMILY に入れて確かめること。');
  }
  return { n: costs.size, joined: joined.length, costs };
}

// ════════════════════════════════════════════════════════════════════
// 5. selftest
// ════════════════════════════════════════════════════════════════════
function selftest() {
  let ok = 0, ng = 0;
  const t = (name, cond) => { if (cond) { ok++; } else { ng++; console.log('  NG ' + name); } };
  const close = (a, b, e = 1e-9) => Math.abs(a - b) <= e;

  // --- parseNotification
  const N = '<task-id>abc</task-id>\n<tool-use-id>toolu_1</tool-use-id>\n<status>completed</status>\n' +
            '<summary>Agent "B1 テスト" finished</summary>\n<result>x</result>\n' +
            '<usage><subagent_tokens>123</subagent_tokens><tool_uses>7</tool_uses><duration_ms>4500</duration_ms></usage>';
  const p = parseNotification(N);
  t('parse: duration', p && p.durationMs === 4500);
  t('parse: tool_uses', p && p.toolUses === 7);
  t('parse: tokens', p && p.tokens === 123);
  t('parse: name から Agent "…" finished を剥がす', p && p.name === 'B1 テスト');
  t('parse: toolUseId', p && p.toolUseId === 'toolu_1');
  t('parse: usage が無ければ null', parseNotification('<summary>x</summary>') === null);
  t('parse: 文字列でなければ null', parseNotification(null) === null);

  // --- parseEstimate
  t('est: 見積 **400–700 行**', JSON.stringify(parseEstimate('見積 **400–700 行**。')) === '{"lo":400,"hi":700}');
  t('est: 見積 120-250 行', JSON.stringify(parseEstimate('見積 120-250 行（Y4）')) === '{"lo":120,"hi":250}');
  t('est: 見積 約 300 行', JSON.stringify(parseEstimate('見積 約 300 行')) === '{"lo":300,"hi":300}');
  t('est: 無ければ null', parseEstimate('行数は 400 行') === null);

  // --- countAbstractCores
  const src = ['/-! ## 0. 抽象核——一般の Galois 拡大', '/-- **抽象核 1**: … -/', '/-! ## 1. 具体層', '/-! ## 2. ★抽象核 B —— x'].join('\n');
  t('cores: 節見出しだけ数える(2)', countAbstractCores(src) === 2);
  t('cores: docstring は数えない', countAbstractCores('/-- 抽象核 1 -/') === 0);

  // --- quantile(numpy type 7 と一致)
  const q = [1, 2, 3, 4];
  t('quantile: median 2.5', close(quantile(q, 0.5), 2.5));
  t('quantile: q1 1.75', close(quantile(q, 0.25), 1.75));
  t('quantile: q3 3.25', close(quantile(q, 0.75), 3.25));
  t('quantile: 単一要素', quantile([5], 0.5) === 5);
  t('quantile: 空は NaN', Number.isNaN(quantile([], 0.5)));

  // --- ranks
  t('ranks: 同順位は平均', JSON.stringify(ranks([10, 20, 20, 30])) === '[1,2.5,2.5,4]');
  t('ranks: 逆順', JSON.stringify(ranks([3, 2, 1])) === '[3,2,1]');

  // --- pearson / spearman(既知値)
  t('pearson: 完全一致 1', close(pearson([1, 2, 3], [1, 2, 3]), 1));
  t('pearson: 完全逆 -1', close(pearson([1, 2, 3], [3, 2, 1]), -1));
  // scipy.stats.pearsonr([1,2,3,4,5],[2,1,4,3,5]) = 0.8 (scipy 1.17.1 で実測)
  t('pearson: scipy 既知値 0.8', close(pearson([1, 2, 3, 4, 5], [2, 1, 4, 3, 5]), 0.8, 1e-12));
  // scipy.stats.spearmanr([1,2,3,4,5],[5,6,7,8,7]) = 0.8207826816681233 (scipy 1.17.1 で実測)
  t('spearman: scipy 既知値 0.8207826816681233',
    close(spearman([1, 2, 3, 4, 5], [5, 6, 7, 8, 7]), 0.8207826816681233, 1e-12));
  // scipy.stats.spearmanr([1,2,3,4,5,6],[6,5,4,3,2,1]) = -1
  t('spearman: scipy 既知値 -1', close(spearman([1, 2, 3, 4, 5, 6], [6, 5, 4, 3, 2, 1]), -1, 1e-12));
  t('spearman: 定数は NaN', Number.isNaN(spearman([1, 1, 1], [1, 2, 3])));

  // --- mulberry32 の決定性
  const r1 = mulberry32(42), r2 = mulberry32(42);
  t('rng: seed が同じなら同じ列', r1() === r2() && r1() === r2());
  t('rng: seed が違えば違う列', mulberry32(1)() !== mulberry32(2)());

  // --- permP
  const pa = permP([1, 2, 3, 4, 5, 6], [6, 5, 4, 3, 2, 1], 5000, 7);
  t('permP: 完全逆順は最小 p(2/(n!)=1/360 付近)', pa > 0 && pa < 0.01);
  t('permP: 決定的(2 回同じ)',
    permP([1, 2, 3, 4, 5], [2, 1, 4, 3, 5], 3000, 11) === permP([1, 2, 3, 4, 5], [2, 1, 4, 3, 5], 3000, 11));
  t('permP: p は 0 にならない', permP([1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6, 7], 2000, 3) > 0);
  {
    // 無相関に近い列では p が大きいこと(p ≧ 0.2)
    const pb = permP([1, 2, 3, 4, 5, 6, 7, 8], [3, 1, 4, 8, 2, 7, 5, 6], 5000, 5);
    t('permP: ばらばらな列は p が大きい', pb > 0.2);
  }

  // --- holm
  const h = holm([0.01, 0.04, 0.03]);
  // 昇順 0.01(×3=0.03), 0.03(×2=0.06), 0.04(×1=0.04→単調化で 0.06)
  t('holm: 最小 p の調整', close(h[0], 0.03, 1e-12));
  t('holm: 単調化(0.04 が 0.06 に持ち上がる)', close(h[1], 0.06, 1e-12));
  t('holm: 中間', close(h[2], 0.06, 1e-12));
  t('holm: 1 で頭打ち', holm([0.9, 0.9, 0.9]).every(v => v === 1));
  t('holm: m=1 は素通し', close(holm([0.02])[0], 0.02, 1e-12));

  // --- partial
  t('partial: 交絡ゼロなら素通し', close(partial(0.5, 0, 0), 0.5, 1e-12));
  t('partial: 完全交絡で消える', close(partial(0.5, 0.5, 1.0), 0, 1e-9) || Number.isNaN(partial(0.5, 0.5, 1.0)));

  // --- icc1
  const same = icc1([[1, 1], [5, 5], [9, 9]]);   // 束内分散 0 ⇒ ICC = 1
  t('icc: 束内が同一なら ICC = 1', close(same.icc, 1, 1e-9));
  const none = icc1([[1, 9], [1, 9], [1, 9]]);   // 束間分散 0 ⇒ ICC ≦ 0
  t('icc: 束間に差が無ければ ICC ≦ 0', none.icc <= 1e-9);
  t('icc: DEFF は 1 以上', same.deff >= 1 && none.deff >= 1);
  t('icc: 束 1 つでは NaN', Number.isNaN(icc1([[1, 2, 3]]).icc));

  // --- bundleKey
  t('bundle: 枝番を落とす(段1c → 段)', bundleKey('段1c 形式群則の混合項と hρ') === '段');
  t('bundle: Λ9 → Λ', bundleKey('Λ9 円分子を副有限捩れから') === 'Λ');
  t('bundle: B1 → B', bundleKey('B1 アーベル閉包の定義') === 'B');
  t('bundle: 数字が無ければそのまま', bundleKey('hρ を外す帳簿 3 手') === 'hρ');
  t('bundle: B2+B3 は先頭側で束ねる', bundleKey('B2+B3 E_σ と K^ab の分解') === 'B');

  // --- ★COST(自己申告)の書式
  t('cost: 素の 1 行', JSON.stringify(parseCostLine('COST[Y19e]: 安')) === '{"key":"Y19e","level":"安","slot":null,"note":""}');
  t('cost: ★や - の飾りを許す', parseCostLine('- ★COST[B2+B3]: 高 — 抽象核が 2 本増えた')?.level === '高');
  t('cost: 一言を切り出す', parseCostLine('COST[Λ12]: 並 — 見積どおり')?.note === '見積どおり');
  t('cost: 全角コロン', parseCostLine('COST[Y9]：安')?.key === 'Y9');
  t('cost: 書式の見本は数えない(値も見本)', parseCostLine('COST[<持ち場>]: <安|並|高>') === null);
  // ★★値が正しくて**鍵だけが見本**の行。★これが無いと `[<>]` の番人を外しても鳴らない
  //   (第 24 回がわざと壊して気づいた。値の側の列挙で弾かれていて空虚だった)。
  t('cost: ★鍵だけが見本でも数えない', parseCostLine('COST[<持ち場>]: 安 — 見本') === null);
  t('cost: 値が列挙外なら読まない', parseCostLine('COST[X]: そこそこ') === null);
  t('cost: COST を含むだけの散文は読まない', parseCostLine('この行は COST[ を含むが書式ではない') === null);
  t('cost: 鍵の前後の空白を落とす', parseCostLine('COST[ Y19e ]: 高')?.key === 'Y19e');
  t('cost: 橋が無ければ slot は null', parseCostLine('COST[Y19e]: 安 — x')?.slot === null);
  t('cost: ★明示の橋を読む', parseCostLine('COST[LT]: 安 | 持ち場=局所Tate双対性の道')?.slot === '局所Tate双対性の道');
  t('cost: 橋を一言から取り除く',
    parseCostLine('COST[LT]: 安 | 持ち場=局所Tate双対性の道 — 早かった')?.note === '早かった');
  {
    const tmp3 = (process.env.TEMP || '.') + '/_agent-timing-cost-selftest3.md';
    fs.writeFileSync(tmp3, ['COST[LT]: 安 | 持ち場=局所Tate双対性の道'].join('\n'));
    const rows4 = [{ description: '局所Tate双対性の道', durationMs: 60000, lines: 100, est: { lo: 80, hi: 120 } }];
    const buf4 = []; const real4 = console.log;
    console.log = (...a) => buf4.push(a.join(' '));
    try { cmdCost(rows4, [tmp3]); } finally { console.log = real4; }
    t('cost: ★橋があれば鍵が違っても当たる', buf4.join('\n').includes('突き合わせ 1 件'));
    fs.unlinkSync(tmp3);
  }
  // ★突き合わせの境界(ここが緩いと Y19 が Y19e を食う)
  t('join: 完全一致', matchesKey('Y19e', 'Y19e'));
  t('join: 鍵 + 空白 + 説明', matchesKey('Y19e compat の心臓', 'Y19e'));
  t('join: ★Y19 は Y19e を拾わない', !matchesKey('Y19e compat の心臓', 'Y19'));
  t('join: Y19 は Y19 の持ち場を拾う', matchesKey('Y19 絶対Galois群の分岐', 'Y19'));
  t('join: 記号を含む鍵', matchesKey('B2+B3 E_σ と K^ab の分解', 'B2+B3'));
  t('join: 空の鍵は拾わない', !matchesKey('Y19e', ''));
  {
    const tmp = process.env.TEMP || '.';
    const p = tmp + '/_agent-timing-cost-selftest.md';
    fs.writeFileSync(p, ['# 見本', 'COST[<持ち場>]: <安|並|高>', 'COST[Y19e]: 安 — 早かった',
                         '★COST[Λ12]: 高', 'COST[Y19e]: 並 — 書き直した'].join('\n'));
    const { costs } = readCosts([p]);
    t('readCosts: 見本を除いて 2 鍵', costs.size === 2);
    t('readCosts: 同じ鍵は後の記述で上書き', costs.get('Y19e').level === '並');
    t('readCosts: 行番号を持つ', costs.get('Λ12').line === 4);
    const { missing } = readCosts([p + '.nope']);
    t('readCosts: 無い記録は missing に', missing.length === 1);
    fs.unlinkSync(p);
  }
  {
    // ★★空虚さの見張り: 申告が 0 件のとき **判定めいた文を出さない**
    const grab2 = (fn) => { const buf = []; const real = console.log;
      console.log = (...a) => buf.push(a.join(' ')); try { fn(); } finally { console.log = real; }
      return buf.join('\n'); };
    const out0 = grab2(() => cmdCost([{ description: 'Y19e x', durationMs: 1, lines: 1, est: null }], ['/nope/none.md']));
    t('cost: 申告 0 件なら「何も言わない」と書く', out0.includes('分母が立たない'));
    t('cost: 申告 0 件で ρ を出さない', !out0.includes('Spearman'));
    const tmp2 = (process.env.TEMP || '.') + '/_agent-timing-cost-selftest2.md';
    fs.writeFileSync(tmp2, ['COST[Y1]: 安', 'COST[Y2]: 高', 'COST[ZZ]: 並'].join('\n'));
    const rows0 = [{ description: 'Y1 なにか', durationMs: 60000, lines: 100, est: { lo: 80, hi: 120 } },
                   { description: 'Y2 べつ', durationMs: 60000, lines: 300, est: { lo: 80, hi: 120 } }];
    const out1 = grab2(() => cmdCost(rows0, [tmp2]));
    t('cost: 突き合わせ件数を出す', out1.includes('突き合わせ 2 件'));
    t('cost: 当たらなかった申告を名指しする', out1.includes('ZZ'));
    t('cost: ★件数不足なら判定を出さない', out1.includes('判定を出さない') && !out1.includes('Spearman'));
    // ★★水準が偏っていたら、件数が足りていても判定を出さない(事前登録した (b))
    {
      const tmp4 = (process.env.TEMP || '.') + '/_agent-timing-cost-selftest4.md';
      const keys = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'B1'];
      fs.writeFileSync(tmp4, keys.map((k, i) => `COST[${k}]: ${i < 7 ? '安' : '並'}`).join('\n'));
      const rows4 = keys.map(k => ({ description: k + ' x', durationMs: 60000, lines: 100, est: { lo: 80, hi: 120 } }));
      const out4 = grab2(() => cmdCost(rows4, [tmp4]));
      t('cost: ★8 件あっても水準が偏れば判定を出さない', out4.includes('突き合わせ 8 件')
        && out4.includes('判定を出さない') && !out4.includes('Spearman'));
    }
    fs.unlinkSync(tmp2);
  }

  // --- ★族(事前登録)—— 欄が黙って消えたら鳴る
  t('family: 事前登録した欄が族に居る(estMid)', FAMILY.some(f => f.key === 'estMid'));
  t('family: 族は 7 本', FAMILY.length === 7);
  t('family: 鍵が重複しない', new Set(FAMILY.map(f => f.key)).size === FAMILY.length);
  // ★★メタ第 25 回が決めた正典を固定する。黙って差し替わったら鳴る(理由は FAMILY の頭に書いた)。
  t('family: ★行数の正典は lines(linesWritten ではない)', FAMILY.some(f => f.key === 'lines'));
  t('family: ★linesWritten は族に入れない(m を増やさない)', !FAMILY.some(f => f.key === 'linesWritten'));
  // ★★COST(自己申告)は n が足りるまで族に入れない —— 事前登録(メタ第 25 回)。
  //    条件は 2 つとも満たすこと: (a) 突き合わせ 8 件以上 (b) 各水準に 3 件以上。
  t('family: ★COST はまだ族に入れない', !FAMILY.some(f => f.key === 'costLevel'));
  t('cost: ★族に入れる条件(8 件かつ各水準 3 件)', costFamilyReady([]) === false
     && costFamilyReady(['安', '安', '安', '並', '並', '並', '高', '高']) === false
     && costFamilyReady(['安', '安', '安', '並', '並', '並', '高', '高', '高']) === true);

  // --- ★欄ごとの complete-case
  {
    const rs = [{ a: 1, b: NaN }, { a: 2, b: 5 }, { a: 3, b: NaN }];
    t('usableRows: NaN の行を落とす', usableRows(rs, 'b').length === 1 && usableRows(rs, 'a').length === 3);
    t('MIN_N は 8', MIN_N === 8);
  }
  {
    // ★★空虚さの見張り: 欠測のある欄を足しても **他の欄の ρ が動かない** こと。
    //   ★ここが壊れる(listwise に戻す/欠測行を 0 で埋める)と、この 2 件が鳴る。
    const mk = (n) => Array.from({ length: n }, (_, i) => ({
      durationMs: (i + 1) * 1000, lines: (i + 1) * 10, toolUses: i + 1, tokens: (i + 1) * 7,
      cores: i % 3, leanChecks: i % 4, checkFails: i % 5,
      // ★半分だけ見積がある(欠測が欄で違う状況を作る)
      estMid: i % 2 === 0 ? 100 * (n - i) : NaN,
      description: 'X' + i,
    }));
    const grab = (fn) => { const buf = []; const real = console.log;
      console.log = (...a) => buf.push(a.join(' ')); try { fn(); } finally { console.log = real; }
      return buf.join('\n'); };
    const rows = mk(20);
    const out = grab(() => cmdExplain(rows));
    t('explain: 欠測のある欄も検定に入る(m = 7)', out.includes('m = 7'));
    t('explain: 欄ごとの n を印字する(10 と 20 が並ぶ)', /\s10\s/.test(out) && /\s20\s/.test(out));
    // 同じ盤面から estMid を丸ごと外したものと、lines の ρ を比べる
    const rows2 = mk(20).map(r => ({ ...r, estMid: NaN }));
    const rho1 = spearman(rows.map(r => r.lines), rows.map(r => r.durationMs));
    const rho2 = spearman(rows2.map(r => r.lines), rows2.map(r => r.durationMs));
    t('explain: 欠測欄を足しても他の欄の ρ は動かない', close(rho1, rho2, 1e-12));
    const out2 = grab(() => cmdExplain(rows2));
    t('explain: 全欠測の欄は「測れなかった」に落ちる', out2.includes('測れなかった') && out2.includes('m = 6'));
    // ★件数不足(MIN_N 未満)の欄はそう書く
    const rows3 = mk(20).map((r, i) => ({ ...r, estMid: i < 5 ? r.estMid : NaN }));
    const out3 = grab(() => cmdExplain(rows3));
    t('explain: MIN_N 未満は件数不足と書いて m から外す', out3.includes('件数不足') && out3.includes('m = 6'));
    t('explain: 向きを断定する文を出さない', !/(長いほど|多いほど|ほど遅い|ほど速い)/.test(out));
  }

  // --- describe
  const d = describe([1, 2, 3, 4]);
  t('describe: 中央値 2.5', close(d.med, 2.5));
  t('describe: 合計 10', d.sum === 10);

  // --- 5b. --diag(診断の族 × 費用。★メタ第 27 回)
  {
    // Fisher —— 既知の値で較正する
    t('diag: Fisher お茶の実験 [[3,1],[1,3]] = 0.4857', close(fisher(3, 1, 1, 3), 0.4857142857, 1e-6));
    t('diag: Fisher 完全分離 [[10,0],[0,10]]', close(fisher(10, 0, 0, 10), 2 / 184756, 1e-9));
    t('diag: Fisher 独立な表は 1', close(fisher(5, 5, 5, 5), 1, 1e-9));
    t('diag: Fisher 空の行は 1', fisher(0, 0, 3, 3) === 1);
    t('diag: Fisher は対称', close(fisher(7, 2, 3, 8), fisher(2, 7, 8, 3), 1e-12));
    // needBalancedN —— 割合が同じなら null、違えば α を割る n を返す
    t('diag: あと何件 割合が同じなら null', needBalancedN(5, 5, 5, 5, 0.005) === null);
    const nb = needBalancedN(8, 2, 2, 8, 0.005);
    t('diag: あと何件 は MIN_N 以上の整数', Number.isInteger(nb) && nb >= MIN_N);
    t('diag: あと何件 は本当に α を割る', (() => { const A = Math.round(0.8 * nb), C = Math.round(0.2 * nb); return fisher(A, nb - A, C, nb - C) <= 0.005; })());
    t('diag: あと何件 その 1 つ手前では割らない', (() => { const n = nb - 1; if (n < MIN_N) return true; const A = Math.round(0.8 * n), C = Math.round(0.2 * n); return fisher(A, n - A, C, n - C) > 0.005; })());
    // tallyDiagnostics —— doc / meta / 本体セッション(ag 空)は入れない
    const errs = [
      { ag: 'T1', src: 'tree', msg: 'failed to synthesize instance Foo' },
      { ag: 'T1', src: 'scratch', msg: 'unsolved goals x' },
      { ag: 'T1', src: 'doc', msg: 'failed to synthesize (文書を印字しただけ)' },
      { ag: 'T1', src: 'meta', msg: 'failed to synthesize (改善係の probe)' },
      { ag: '', src: 'tree', msg: 'failed to synthesize (本体セッション)' },
      { ag: 'T2', src: 'scratch', msg: 'Type mismatch here' },
    ];
    const tal = tallyDiagnostics(errs);
    t('diag: 本体セッション(ag 空)は数えない', !tal.has(''));
    t('diag: doc / meta は数えない', tal.get('T1').B.total === 2);
    // ★スコープが挙げても doc / meta は数えない(★門番が二重であることを確かめる)
    const talX = tallyDiagnostics(errs, DIAG_FAMILIES, [{ id: 'X', srcs: ['tree', 'scratch', 'doc', 'meta'] }]);
    t('diag: スコープが挙げても doc / meta は数えない', talX.get('T1').X.total === 2);
    t('diag: スコープ A は tree だけ', tal.get('T1').A.total === 1 && tal.get('T1').A.fam[0] === 1 && tal.get('T1').A.fam[1] === 0);
    t('diag: スコープ B は tree + scratch', tal.get('T1').B.fam[0] === 1 && tal.get('T1').B.fam[1] === 1);
    t('diag: 族は部分文字列で当てる', tal.get('T2').B.fam[4] === 1);
    // medianSplit —— 同値は low 側(★事前に決めた向き)
    const ms2 = medianSplit([1, 1, 1, 5]);
    t('diag: 中央値 2 値化 同値は low', ms2.isHigh(1) === false && ms2.isHigh(5) === true);
    // cmdDiag —— 件数不足なら検定しない / 向きを断定しない / m を印字する
    const gr = (fn) => { const buf = []; const real = console.log;
      console.log = (...a) => buf.push(a.join(' ')); try { fn(); } finally { console.log = real; }
      return buf.join('\n'); };
    const mkRec = (i, fam) => ({ toolUseId: 'A' + i, agentType: DIAG_TYPE, toolUses: 10, durationMs: (i + 1) * 60000, day: '2026-09-0' + (1 + (i % 3)), fam });
    const recsA = Array.from({ length: 24 }, (_, i) => mkRec(i, i % 2 === 0));
    const digA = { v: 3, errors: recsA.filter(r => r.fam).map(r => ({ ag: r.toolUseId, src: 'tree', msg: 'failed to synthesize instance' })) };
    const outD = gr(() => cmdDiag(recsA, digA));
    t('diag: m = 10 を印字する', outD.includes('m = 10'));
    t('diag: 曝露が MIN_N を満たせば検定する', /言える|言えない/.test(outD));
    t('diag: 曝露ゼロの族は件数不足と書く', outD.includes('件数不足'));
    t('diag: 向きを断定する文を出さない', !/(ほど遅い|ほど速い|のほうが高い|のほうが遅い|を食っている)/.test(outD));
    t('diag: 「族が時間を食う」は否定形でしか書かない', /が時間を食う」とは言えない/.test(outD));
    // ★★多重比較の補正が本当に効いているか —— 素の p は 0.05 を割るが、m = 10 を掛けると割らない盤面。
    //   ★ここが素通りすると「言えない」が「★言える」に化ける。
    const hi = new Set([...Array.from({ length: 9 }, (_, k) => 12 + k), 0, 1, 2]);   // 曝露 12 件(うち高 9)
    const recsC = Array.from({ length: 24 }, (_, i) => mkRec(i, hi.has(i)));
    const digC = { v: 3, errors: recsC.filter(r => r.fam).map(r => ({ ag: r.toolUseId, src: 'tree', msg: 'failed to synthesize instance' })) };
    const outC = gr(() => cmdDiag(recsC, digC));
    t('diag: 盤面が意図どおり(素の p ≒ 0.039)', /0\.039/.test(outC));
    // ★「★言えるのは…」という地の文にも同じ字面が入るので、★行末の判定欄だけを見る。
    t('diag: 補正すると言えない(補正を外すと言えるになる盤面)', !/★言える\s*$/m.test(outC));
    const digEmpty = { v: 3, errors: [] };
    const outE = gr(() => cmdDiag(recsA, digEmpty));
    t('diag: 診断ゼロでも落ちない(全部 件数不足)', (outE.match(/件数不足/g) || []).length === 10);
    t('diag: 突き合わせられない件数を印字する', outD.includes('本体セッション'));
  }

  console.log(`\nselftest: ${ok}/${ok + ng}`);
  return ng === 0;
}

// ════════════════════════════════════════════════════════════════════
// 5b. 診断の族 × 費用(`--diag`)—— ★★事前登録(メタ第 27 回)
// ════════════════════════════════════════════════════════════════════
//
// ★何を測るか / ★★何は測れないか(先に書く)
// ------------------------------------------------
// ★測れない: **「族 → 時間」の因果**。難しいノードは長くもあり診断も多い(交絡)。
//   ★測れるのはせいぜい「族の混合比が違う agent の間で、**1 往復あたりの時間**が違うか」まで。
//   ⇒ 出す判定は「言える / 言えない」だけ。★向き(どちらが高い)は断定しない。
// ★★さらに強い制約(メタ第 27 回の実測): `tree` の診断 4,379 件のうち **4,176 件(95%)は
//   本体セッションが出したもの**で、★**本体には duration_ms が無い**(agent ではないため)。
//   ⇒ 費用と突き合わせられるのは残り 203 件だけ。★これは道具の不備ではなく**観測点の構造**である。
//
// ★事前登録(★データを見る前にここへ焼く。後から欄を足さない)
// ------------------------------------------------
//   母集団 : `agentType === 'lean-prover'`(= 実装者)かつ `toolUses >= 1`。
//            ★役割で切る理由: math-planner / Explore は Lean を叩かないので
//              「診断ゼロ かつ 速い」が構造的に決まり、族の効果と役割が完全に交絡する。
//   結果 Y : `duration_ms / tool_uses`(1 往復あたりのミリ秒)。★中央値で 2 値化(high = Y > 中央値)。
//   曝露 X : その agent が族 f の診断を **1 件以上**出したか。
//   族     : DIAG_FAMILIES(5 つ。★M115 の tree 上位 5 族をそのまま固定)。
//   スコープ: A = `tree` のみ(主) / B = `tree` + `scratch`(副)。★`doc` と `meta` は常に外す。
//   検定   : Fisher 正確検定(両側)。★多重比較は **m = 5 族 × 2 スコープ = 10** の Bonferroni。
//   最小件数: 曝露あり群・曝露なし群ともに MIN_N(8)件以上。★未満なら**検定しない**(M112 と同じ作法)。
export const DIAG_FAMILIES = [
  { id: 'D1', key: 'failed to synthesize' },
  { id: 'D2', key: 'unsolved goals' },
  { id: 'D3', key: 'Did not find an occurrence of the pattern' },
  { id: 'D4', key: 'Application type mismatch' },
  { id: 'D5', key: 'Type mismatch' },
];
export const DIAG_SCOPES = [
  { id: 'A', label: 'tree のみ(木の本を検査して出た失敗)', srcs: ['tree'] },
  { id: 'B', label: 'tree + scratch(REPL の断片も入れる)', srcs: ['tree', 'scratch'] },
];
export const DIAG_TYPE = 'lean-prover';
export const DIAG_ALPHA = 0.05;
export const DIAG_M = DIAG_FAMILIES.length * DIAG_SCOPES.length;   // = 10

/** ★対数ガンマ(Lanczos)。Fisher を n が大きくても引けるようにするため。 */
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
/**
 * Fisher の正確検定(両側・確率法)。★`unverified.mjs` の `fisherLog` と同じ式。
 * ★あちらは厳密版(BigInt)と食い違わないことを selftest で確かめてある。ここでは
 *   ★既知の値(お茶の実験 [[3,1],[1,3]] = 0.4857 など)で較正する。
 */
export function fisher(a, b, c, d) {
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

/**
 * ★「あと何件あれば言えるか」—— 観測した**割合を固定**したまま、均衡した群サイズを増やして
 * Fisher p が α' を割る最小の 1 群あたり件数を返す。★多重比較込みの α' を渡すこと。
 * ★これは検出力計算ではない(割合が本当にそうなら、という条件付きの最小 n)。★見つからなければ null。
 */
export function needBalancedN(a, b, c, d, alpha, cap = 4000) {
  const n1 = a + b, n2 = c + d;
  if (!n1 || !n2) return null;
  const p1 = a / n1, p2 = c / n2;
  if (Math.abs(p1 - p2) < 1e-12) return null;      // 割合が同じなら何件増やしても言えない
  for (let n = MIN_N; n <= cap; n++) {
    const A = Math.round(p1 * n), C = Math.round(p2 * n);
    if (fisher(A, n - A, C, n - C) <= alpha) return n;
  }
  return null;
}

/** digest(idiom-recur v3)を読み、agent(toolUseId)ごとの族の件数にする。★純関数。 */
export function tallyDiagnostics(errors, families = DIAG_FAMILIES, scopes = DIAG_SCOPES) {
  const byAgent = new Map();
  for (const e of errors || []) {
    const ag = e.ag || '';
    if (!ag) continue;                              // ★本体セッション(agent ではない)は突き合わせられない
    const src = e.src || 'tree';
    if (src === 'doc' || src === 'meta') continue;  // ★Lean のエラーではない
    const msg = String(e.msg || '');
    let rec = byAgent.get(ag);
    if (!rec) { rec = {}; for (const s of scopes) rec[s.id] = { total: 0, fam: families.map(() => 0) }; byAgent.set(ag, rec); }
    for (const s of scopes) {
      if (!s.srcs.includes(src)) continue;
      rec[s.id].total++;
      families.forEach((f, i) => { if (msg.includes(f.key)) rec[s.id].fam[i]++; });
    }
  }
  return byAgent;
}

/** 中央値で 2 値化する。★high = 値 > 中央値(同値は low 側)。★事前に決めた向き。 */
export function medianSplit(vals) {
  const med = describe(vals).med;
  return { med, isHigh: (v) => v > med };
}

function cmdDiag(recs, digest) {
  const errs = digest?.errors || [];
  const byAgent = tallyDiagnostics(errs);
  const pop = recs.filter(r => r.agentType === DIAG_TYPE && r.toolUses >= 1 && r.durationMs > 0 && r.toolUseId);
  const rows = pop.map(r => ({
    id: r.toolUseId, day: r.day, ms: r.durationMs, tu: r.toolUses,
    y: r.durationMs / r.toolUses,
    diag: byAgent.get(r.toolUseId) || null,
  }));
  const { med, isHigh } = medianSplit(rows.map(r => r.y));
  const dY = rows.length ? describe(rows.map(r => r.y)) : null;

  console.log('== --diag: 診断の族 × 1 往復あたりの費用(★事前登録。向きは断定しない) ==\n');
  console.log(`  母集団      : agentType=${DIAG_TYPE} かつ tool_uses>=1 …… ${rows.length} 件`);
  console.log(`  結果 Y      : duration_ms / tool_uses。中央値 ${(med / 1000).toFixed(1)} 秒/往復(high = これより上)`);
  if (dY) console.log(`  ★Y の散らばり: 四分位 ${(dY.q1 / 1000).toFixed(1)} — ${(dY.q3 / 1000).toFixed(1)} 秒/往復`
    + `(最小 ${(dY.min / 1000).toFixed(1)} / 最大 ${(dY.max / 1000).toFixed(1)})。★下の族ごとの中央値はこの幅と比べて読むこと。`);
  console.log(`  多重比較    : Bonferroni m = ${DIAG_M}(族 ${DIAG_FAMILIES.length} × スコープ ${DIAG_SCOPES.length})  α' = ${(DIAG_ALPHA / DIAG_M).toFixed(4)}`);
  const attributable = errs.filter(e => e.ag && e.src === 'tree').length;
  const treeAll = errs.filter(e => e.src === 'tree').length;
  console.log(`  ★突き合わせられる tree 診断: ${attributable} / ${treeAll} 件`
    + `(残り ${treeAll - attributable} 件は**本体セッション**が出したもので duration_ms が無い)`);

  const P = (n, w) => String(n).padStart(w);
  for (const s of DIAG_SCOPES) {
    console.log(`\n  -- スコープ ${s.id}: ${s.label} --`);
    console.log('   族                              曝露あり  うち高  曝露なし  うち高  Fisher p   補正後   判定');
    for (let i = 0; i < DIAG_FAMILIES.length; i++) {
      const f = DIAG_FAMILIES[i];
      const ex = rows.filter(r => r.diag && r.diag[s.id].fam[i] > 0);
      const un = rows.filter(r => !(r.diag && r.diag[s.id].fam[i] > 0));
      const a = ex.filter(r => isHigh(r.y)).length, b = ex.length - a;
      const c = un.filter(r => isHigh(r.y)).length, d = un.length - c;
      let p = NaN, adj = NaN, verdict;
      if (ex.length < MIN_N || un.length < MIN_N) {
        verdict = `件数不足(要 ${MIN_N})`;
      } else {
        p = fisher(a, b, c, d); adj = Math.min(1, p * DIAG_M);
        verdict = adj <= DIAG_ALPHA ? '★言える' : '言えない';
      }
      console.log(`   ${f.key.slice(0, 30).padEnd(32)}${P(ex.length, 6)}${P(a, 8)}${P(un.length, 9)}${P(c, 8)}`
        + `${isFinite(p) ? P(p.toFixed(4), 10) : P('—', 10)}${isFinite(adj) ? P(adj.toFixed(3), 9) : P('—', 9)}  ${verdict}`);
      // ★記述(★判定ではない。向きの主張でもない)
      const mEx = ex.length ? describe(ex.map(r => r.y)).med / 1000 : NaN;
      const mUn = un.length ? describe(un.map(r => r.y)).med / 1000 : NaN;
      const need = (ex.length >= MIN_N && un.length >= MIN_N) ? needBalancedN(a, b, c, d, DIAG_ALPHA / DIAG_M) : null;
      // ★所要時間と往復数そのものも並べる(★本体の問いは「時間・往復数はどう違うか」だった)。
      //   ★これは記述であって検定ではない。★検定は Y(秒/往復)だけに事前登録してある。
      const mdEx = ex.length ? describe(ex.map(r => r.ms)).med / 60000 : NaN;
      const mdUn = un.length ? describe(un.map(r => r.ms)).med / 60000 : NaN;
      const tuEx = ex.length ? describe(ex.map(r => r.tu)).med : NaN;
      const tuUn = un.length ? describe(un.map(r => r.tu)).med : NaN;
      console.log(`     記述: 秒/往復 中央値 ${isFinite(mEx) ? mEx.toFixed(1) : '—'} 対 ${isFinite(mUn) ? mUn.toFixed(1) : '—'}`
        + `  ・所要 分 ${isFinite(mdEx) ? mdEx.toFixed(1) : '—'} 対 ${isFinite(mdUn) ? mdUn.toFixed(1) : '—'}`
        + `  ・往復 ${isFinite(tuEx) ? tuEx : '—'} 対 ${isFinite(tuUn) ? tuUn : '—'}`
        + `  ・交絡: 診断総数 ${ex.length ? describe(ex.map(r => r.diag[s.id].total)).med : '—'} 対 ${un.length ? describe(un.map(r => (r.diag ? r.diag[s.id].total : 0))).med : '—'}`
        + `  ・★あと何件(1 群): ${need ?? '—'}`);
    }
  }
  // ★標本は独立でない —— 日ごとの束で ICC を出す(M90 / M99 と同じ作法)
  const byDay = new Map();
  for (const r of rows) { const k = r.day || '?'; if (!byDay.has(k)) byDay.set(k, []); byDay.get(k).push(r.y); }
  const ic = icc1([...byDay.values()]);
  console.log(`\n  ★束(日)= ${ic.k} / N = ${ic.N} / ICC ${isFinite(ic.icc) ? ic.icc.toFixed(3) : '—'}`
    + ` / DEFF ${isFinite(ic.deff) ? ic.deff.toFixed(2) : '—'} / 実効 n ${isFinite(ic.nEff) ? ic.nEff.toFixed(1) : '—'}`);
  console.log('  ★★上の「あと何件」は**実効 n ではなく素の n** の話である(束を無視している。実際にはもっと要る)。');
  console.log('  ★★この表から「族 f が時間を食う」とは言えない。★言えるのは「族 f を出した agent と');
  console.log('    出さなかった agent とで、1 往復あたりの時間の高低の分かれ方が違う(違わない)」まで。');
}

// ════════════════════════════════════════════════════════════════════
// 6. CLI
// ════════════════════════════════════════════════════════════════════
function main() {
  const argv = process.argv.slice(2);
  const has = (f) => argv.includes(f);
  const val = (f, d) => { const i = argv.indexOf(f); return i >= 0 && argv[i + 1] ? argv[i + 1] : d; };

  if (has('--selftest')) { process.exit(selftest() ? 0 : 1); }

  const repoRoot = path.resolve(path.join(new URL('.', import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, '$1'), '..'));
  let recs = collect();
  const type = val('--type', null);
  const day = val('--day', null);
  const name = val('--name', null);
  if (type) recs = recs.filter(r => r.agentType === type);
  if (day) recs = recs.filter(r => r.day === day);
  if (name) recs = recs.filter(r => new RegExp(name).test(r.description));
  const ms = val('--ms', null);           // ★duration_ms のリストで拾う(過去の表を再現するため)
  if (ms) {
    const want = new Set(ms.split(/[,\s]+/).filter(Boolean).map(Number));
    recs = recs.filter(r => want.has(r.durationMs));
  }

  const wantCov = has('--stats') || has('--explain') || has('--estimate') || has('--cov') || has('--cost');
  if (wantCov) {
    for (const r of recs) Object.assign(r, covariates(r, repoRoot));
  }

  if (has('--json')) { console.log(JSON.stringify(recs, null, 1)); return; }

  const filt = `type=${type ?? '全部'} day=${day ?? '全部'}${name ? ' name=/' + name + '/' : ''}`;
  console.log(`== agent-timing —— 通知 ${recs.length} 件(${filt}) ==\n`);

  if (has('--stats')) { cmdStats(recs); console.log(''); cmdIcc(recs); return; }
  if (has('--explain')) { cmdExplain(recs); cmdIcc(recs); return; }
  if (has('--estimate')) { cmdEstimate(recs); return; }
  if (has('--diag')) {
    // ★digest は `idiom-recur.mjs --rescan` が作る。v3 でないと `ag`(toolUseId)が無い。
    const dp = val('--digest', path.join(repoRoot, '.cache', 'idiom-recur-digest.json'));
    if (!fs.existsSync(dp)) { console.error(`digest がない: ${dp}\n  まず \`node tools/idiom-recur.mjs --rescan\``); process.exit(2); }
    const dg = JSON.parse(fs.readFileSync(dp, 'utf8'));
    if ((dg.v || 1) < 3) { console.error(`digest が古い(v${dg.v || 1}、要 v3 以上)。\`node tools/idiom-recur.mjs --rescan\` で作り直すこと。`); process.exit(2); }
    cmdDiag(recs, dg);
    return;
  }
  if (has('--cost')) {
    // ★worktree の decisions-pending.md は本体より古い(M99)。★--costs で本体の版を指せる。
    const cp = val('--costs', path.join(repoRoot, 'ResearchPaper', 'decisions-pending.md'));
    cmdCost(recs, [cp]);
    return;
  }

  // 既定: 一覧
  const lim = Number(val('--limit', '30'));
  const show = lim > 0 ? recs.slice(-lim) : recs;
  console.log('  完了時刻(UTC)          分   tool  tokens  type            持ち場');
  for (const r of show) {
    console.log(`  ${padr((r.ts ?? '').replace('T', ' ').slice(0, 19), 20)} ${pad(min1(r.durationMs), 6)} ${pad(r.toolUses, 5)} ${pad(r.tokens, 8)}  ${padr(r.agentType.slice(0, 14), 15)} ${String(r.description).slice(0, 44)}`);
  }
  if (lim > 0 && recs.length > lim) console.log(`  … 残り ${recs.length - lim} 件(--limit 0 で全部)`);
}

if (import.meta.url.endsWith(process.argv[1].replace(/\\/g, '/')) || process.argv[1].endsWith('agent-timing.mjs')) {
  main();
}
