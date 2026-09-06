#!/usr/bin/env node
/**
 * 前線 —— 「次にどのノードをやると一番効くか」を依存グラフから計算する。
 *
 * ★この道具の位置づけ
 *   `tools/graph.mjs` は「グラフがどうなっているか」を答える。
 *   こちらは **「だから次に何をすべきか」** を答える。
 *   Orchestrator（本体セッション）が sub-agent に仕事を配る前に、必ずここを引く。
 *
 * ★何を測るか（`ResearchPaper/orchestration.md` §2）
 *   sorry を持つノード v について:
 *     downstream  … v を推移的に import しているノード数（v が解けると解放される数）
 *     dsItems     … その下流が持つ `.src` 項目の総数（＝原典の主張いくつ分か）
 *     blockers    … v が推移的に import している **他の** sorry ノード
 *     startable   … blockers が空。すなわち **今すぐ着手できる**
 *
 *   ★`startable` でない節点に人手（agent）を割いても、上流の sorry に当たって
 *     止まる。だから配る順序は「startable のうち downstream が大きいもの」から。
 *
 * ★同時に起動する agent は **5 個まで**（`ResearchPaper/orchestration.md` §0）。
 *   だから既定の出力も 5 件で切る。足りなければ**次の波**にする。
 *
 * ★★印（2026-09-06、メタ第 9 回。backlog M18）
 *   `startable` の判定は**変えない**。持ち場を配る前に見えるべき事実を 3 つ**足すだけ**である。
 *     配線?N … その `sorry` に `Found/` の既存宣言で解ける候補がある（`unwired.mjs` の
 *              「裸 かつ 一致率 80% 以上 かつ 情報量 15 以上」。2026-09-06 に 6/9 が当たりだった）
 *     空撃ち／消費者なし … 誰もその語を使っていない（`unwired.mjs --dead`）
 *     保留 Dn … `ResearchPaper/decisions-pending.md` に**保留として載っている**
 *              （`決定` / `採用` / `解決` だけの節は数えない）
 *   ★実測（2026-09-06）: 前線 19 ノードのうち **11 が保留**、**5 が配線候補**、
 *   **8 が空撃ち／消費者なし**。すなわち **agent を配る前に人の判断が要るものが 6 割**ある。
 *   印を出さないと、これは `frontier.mjs` の画面からは 1 つも見えなかった。
 *
 * 使い方:
 *   node tools/frontier.mjs                 # 着手可能なものを効果の大きい順に（上位 5 件）
 *   node tools/frontier.mjs --all           # 着手不可のものも出す（blockers 付き）
 *   node tools/frontier.mjs --owner pGC     # 所属で絞る
 *   node tools/frontier.mjs --limit 0       # 件数の上限を外す（俯瞰したいときだけ）
 *   node tools/frontier.mjs --json          # Orchestrator が食う形
 *   node tools/frontier.mjs --no-marks      # 印を付けない（1.1 秒。付けると 3.2 秒）
 */

import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));

const args = process.argv.slice(2);
const flag = (f) => args.includes(f);
const opt = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };

/** `graph.mjs --json` を唯一の真実として使う（ここで木を再走査しない）。 */
function loadGraph() {
  const out = execFileSync('node', [join(ROOT, 'tools', 'graph.mjs'), '--json'], {
    encoding: 'utf8', maxBuffer: 1 << 28,
  });
  return JSON.parse(out);
}

const { nodes } = loadGraph();
const byMod = new Map(nodes.map((n) => [n.mod, n]));

/* ★★`import` 依存と数学的依存は一致しない(2026-09-05)
 *
 *   実例: `Skeleton/PGC/Section1Cor13.lean` は 2026-09-05 まで
 *   **着手可能・下流 27・第 1 位**と出ていた。しかし中身の `inertia_recoverable` は
 *   Prop 1.2 に完全に帰着することが証明済みで(`Found/PGC/InertiaTransport.lean`)、
 *   その Prop 1.2 は `Skeleton/PGC/Section1.lean` の **sorry のまま**である。
 *   `Section1Cor13.lean` は Prop 1.2 を型として要らないので import していない。
 *   ★つまりスケジューラの第 1 位が、数学的には塞がっているノードだった。
 *
 *   `graph.mjs` が `.needs`(`.otherPaper` / `.derivation` / `.implicitStep` / `.folklore`)
 *   から拾った辺を `mathEdges` として渡してくる。ここではそれを import と同じ扱いにする。
 *   ★`--no-math` で切れる(2026-09-05 以前の挙動に戻す)。 */
const useMath = !flag('--no-math');
const mathUp = new Map();                    // mod → 数学的な上流 mod[]
const mathVia = new Map();                   // `${from}→${to}` → 理由
for (const n of nodes) {
  const es = useMath ? (n.mathEdges ?? []) : [];
  const ms = [...new Set(es.flatMap((e) => e.mods))].filter((m) => byMod.has(m) && m !== n.mod);
  if (ms.length) mathUp.set(n.mod, ms);
  for (const e of es) for (const m of e.mods) mathVia.set(`${n.mod}→${m}`, e.via);
}

/** 逆辺: mod → それを直接 import しているノードの mod 一覧(数学的な辺を含む)。 */
const rdeps = new Map();
const addR = (from, to) => {                 // to が from に依存している
  if (!rdeps.has(from)) rdeps.set(from, []);
  if (!rdeps.get(from).includes(to)) rdeps.get(from).push(to);
};
for (const n of nodes) {
  for (const im of n.imports) {
    if (!byMod.has(im)) continue;            // Mathlib 等、木の外は辿らない
    addR(im, n.mod);
  }
  for (const m of mathUp.get(n.mod) ?? []) addR(m, n.mod);
}

/** v から辺を辿って到達できる集合（v 自身は含めない）。 */
function reach(start, edgesOf) {
  const seen = new Set();
  const stack = [...(edgesOf(start) ?? [])];
  while (stack.length) {
    const m = stack.pop();
    if (seen.has(m)) continue;
    seen.add(m);
    for (const w of edgesOf(m) ?? []) if (!seen.has(w)) stack.push(w);
  }
  return seen;
}

const downOf = (m) => rdeps.get(m) ?? [];
const upOf = (m) => [
  ...(byMod.get(m)?.imports ?? []).filter((x) => byMod.has(x)),
  ...(mathUp.get(m) ?? []),
];

const sorryMods = new Set(nodes.filter((n) => n.hasSorry).map((n) => n.mod));

const rows = [];
for (const n of nodes) {
  if (!n.hasSorry) continue;
  const down = reach(n.mod, downOf);
  const up = reach(n.mod, upOf);
  const blockers = [...up].filter((m) => sorryMods.has(m)).sort();
  /** ★import だけでは見えない blocker（`.needs` から拾った辺で初めて見えたもの）。 */
  const upImport = reach(n.mod, (m) => (byMod.get(m)?.imports ?? []).filter((x) => byMod.has(x)));
  const mathOnly = blockers.filter((m) => !upImport.has(m));
  let dsItems = 0;
  for (const m of down) dsItems += (byMod.get(m)?.items.length ?? 0);
  rows.push({
    mod: n.mod,
    rel: n.rel,
    owner: n.owner,
    ownerKind: n.ownerKind,
    bucket: n.bucket,
    items: n.items,
    downstream: down.size,
    dsItems,
    blockers,
    mathOnly,
    startable: blockers.length === 0,
  });
}

/* ------------------------------------------------------------------ 印
 *
 * ★ここから下は **`startable` の判定に触らない**。表示だけを足す（backlog M3 と同じ設計）。 */

/** `decisions-pending.md` の**保留節**を読む。★`決定` / `採用` / `解決` だけの節は数えない。
 *
 * 当て方は 2 通りで、**両方要る**（実測 2026-09-06）:
 *   (a) 節の本文が Lean のファイル名を名指ししている … 前線 19 のうち 9 件が当たる
 *   (b) 節の**見出し**の `[論文] 項目` が、ノードの `.src` 項目と一致する … 9 件が当たる
 *   → **どちらかで当たる 11 件**。(a) だけだと `Skeleton/GenEll/VeluSemistable.lean`
 *     （D9 は「[GenEll] Lemma 3.5」としか書いておらずファイル名が無い）を落とす。
 * ★項目の照合は**ローマ数字まで見る**。見ないと D16「[FrdI] Theorem 6.4 (i)」が
 *   `Chebotarev.lean`（Theorem 6.4 **(iv)**）に当たる（実測 1 件の誤報）。 */
const KIND = 'Theorem|Proposition|Lemma|Definition|Corollary|Example|Remark|Claim|Section|Conjecture|Notation';
const ITEM_RE = new RegExp(`\\[(\\w+)\\]\\s*(${KIND})\\s*([0-9]+(?:\\.[0-9]+)*)\\s*,?\\s*(?:\\((i{1,3}|iv|vi{0,3}|ix|x)\\))?`, 'g');
const FILE_RE = /(?:lean[/\\]ABC3[/\\])?((?:Skeleton|Found|Check|Interface|Gap|Meta)[/\\][\w\-./\\]*\.lean)/g;
const itemKeys = (s) => {
  const out = [];
  ITEM_RE.lastIndex = 0;
  let m;
  while ((m = ITEM_RE.exec(s))) out.push([`[${m[1]}] ${m[2]} ${m[3]}`, m[4] ?? null]);
  return out;
};
function loadPending() {
  const byFile = new Map();
  const byItem = new Map();
  let sections = 0;
  let text;
  try { text = readFileSync(join(ROOT, 'ResearchPaper', 'decisions-pending.md'), 'utf8'); }
  catch { return { byFile, byItem, sections: -1 }; }
  let cur = null;
  const secs = [];
  for (const l of text.split(/\r?\n/)) {
    // ★`##` と `D` の間に `★` が入る形(`## ★★★★★★★D24.`)を見ていなかった
    //   ——2026-09-06 実測: D21・D23・D24 の節が丸ごと存在せず、本文が直前の節に
    //   吸収されていた(この解析器は見出し以外で節を切らないため)。
    const h = /^##\s+[★*\s]*(D\d+)\.\s*(.*)$/.exec(l);
    if (h) { cur = { id: h[1], title: h[2], state: null, body: [] }; secs.push(cur); continue; }
    if (!cur) continue;
    if (cur.state === null) { const s = /^[-*]?\s*\*\*状態\*\*\s*[:：]\s*(.*)$/.exec(l); if (s) cur.state = s[1]; }
    cur.body.push(l);
  }
  for (const s of secs) {
    if (!(s.state ?? '').includes('保留')) continue;   // ★`決定` だけの節は落とす
    sections++;
    const body = s.body.join('\n');
    let m;
    FILE_RE.lastIndex = 0;
    while ((m = FILE_RE.exec(body))) {
      const f = m[1].replace(/\\/g, '/');
      if (!byFile.has(f)) byFile.set(f, new Set());
      byFile.get(f).add(s.id);
    }
    for (const [k, roman] of itemKeys(s.title)) {   // ★見出しだけ。本文から拾うと言及が全部当たる
      if (!byItem.has(k)) byItem.set(k, []);
      byItem.get(k).push({ id: s.id, roman });
    }
  }
  return { byFile, byItem, sections };
}

/** `unwired.mjs` を 1 回だけ呼んで「配線候補」と「空撃ち」を得る。
 *  ★無ければ静かに諦める（`unwired.mjs` は D14 で取り込んだ新しい道具なので、
 *    古い木でも `frontier.mjs` が落ちないようにしておく）。 */
function loadUnwired() {
  try {
    const out = execFileSync('node', [join(ROOT, 'tools', 'unwired.mjs'), '--marks'], {
      encoding: 'utf8', maxBuffer: 1 << 28, stdio: ['ignore', 'pipe', 'ignore'],
    });
    const j = JSON.parse(out);
    const strongBy = new Map();   // rel → 「裸 かつ 80% 以上 かつ 情報量 15 以上」の件数
    const anyBy = new Map();      // rel → 候補が 1 件でもある sorry 宣言の数
    for (const r of j.report ?? []) {
      const strong = r.bare && (r.cands ?? []).some((c) => c.cov >= 0.8 && c.mass >= 15);
      if (strong) strongBy.set(r.rel, (strongBy.get(r.rel) ?? 0) + 1);
      if ((r.cands ?? []).length) anyBy.set(r.rel, (anyBy.get(r.rel) ?? 0) + 1);
    }
    const deadBy = new Map();
    for (const d of j.dead ?? []) deadBy.set(d.rel, d.kind);
    return { strongBy, anyBy, deadBy, ok: true };
  } catch {
    return { strongBy: new Map(), anyBy: new Map(), deadBy: new Map(), ok: false };
  }
}

const wantMarks = !flag('--no-marks');
const pending = wantMarks ? loadPending() : null;
const uw = wantMarks ? loadUnwired() : null;
if (wantMarks) {
  for (const r of rows) {
    const ds = new Set(pending.byFile.get(r.rel) ?? []);
    for (const it of r.items) {
      for (const [k, roman] of itemKeys(it)) {
        for (const e of pending.byItem.get(k) ?? []) {
          if (e.roman && roman && e.roman !== roman) continue;  // ★(i) と (iv) を混ぜない
          ds.add(e.id);
        }
      }
    }
    r.marks = {
      unwired: uw.strongBy.get(r.rel) ?? 0,
      cands: uw.anyBy.get(r.rel) ?? 0,
      dead: uw.deadBy.get(r.rel) ?? null,
      pending: [...ds].sort((a, b) => Number(a.slice(1)) - Number(b.slice(1))),
    };
  }
}

const ownerFilter = opt('--owner');
let sel = ownerFilter ? rows.filter((r) => r.owner === ownerFilter) : rows;
if (!flag('--all')) sel = sel.filter((r) => r.startable);

sel.sort((a, b) =>
  (b.startable - a.startable) || (b.downstream - a.downstream) ||
  (b.dsItems - a.dsItems) || a.rel.localeCompare(b.rel));

/** ★既定 5 件。1 波で配る持ち場の上限（`ResearchPaper/orchestration.md` §0）。
 *  上限を外したいときだけ `--limit 0`。 */
const shown = sel.length;
const limit = Number(opt('--limit') ?? 5);
if (limit > 0) sel = sel.slice(0, limit);

if (flag('--json')) {
  console.log(JSON.stringify({ generated: new Date().toISOString(), frontier: sel }, null, 1));
  process.exit(0);
}

const total = rows.length;
const startable = rows.filter((r) => r.startable).length;
const mathBlocked = rows.filter((r) => !r.startable && r.mathOnly.length && r.blockers.length === r.mathOnly.length);
console.log('★前線（sorry を持つノード）');
console.log(`  sorry ノード ${total} / うち着手可能 ${startable}`);
if (useMath && mathBlocked.length) {
  console.log(`  ★うち ${mathBlocked.length} 件は **import には現れない依存**だけで塞がっている`);
  console.log(`    （\`.needs\` から拾った。--no-math で切ると着手可能に見えてしまう）`);
} else if (!useMath) {
  console.log('  ★--no-math: `.needs` の辺を無視している（2026-09-05 以前の挙動）');
}
if (wantMarks) {
  const np = rows.filter((r) => r.marks.pending.length).length;
  const nu = rows.filter((r) => r.marks.unwired).length;
  const nd = rows.filter((r) => r.marks.dead).length;
  console.log(`  ★印: 保留 ${np} 件 / 配線候補 ${nu} 件 / 空撃ち・消費者なし ${nd} 件` +
    `（${total} ノード中）${uw.ok ? '' : '  ※unwired.mjs を呼べなかったので配線/空撃ちは 0'}`);
  console.log('    保留 = decisions-pending.md に**保留として**載っている（人の判断待ち）。--no-marks で消せる');
} else {
  console.log('  ★--no-marks: 保留・配線候補・空撃ちの印を出していない');
}
if (ownerFilter) console.log(`  所属で絞り込み: ${ownerFilter}`);
console.log(`  ${flag('--all') ? '★--all: 着手不可も表示' : '（着手不可を隠している。--all で全部）'}`);
console.log();
console.log('  下流  項目  ノード');
for (const r of sel) {
  const mark = r.startable ? ' ' : '×';
  console.log(`${mark} ${String(r.downstream).padStart(5)} ${String(r.dsItems).padStart(5)}  ${r.rel}  (${r.owner})`);
  if (wantMarks) {
    const tags = [];
    if (r.marks.pending.length) tags.push(`★保留 ${r.marks.pending.join(',')}（人の判断待ち）`);
    if (r.marks.unwired) tags.push(`★配線候補 ${r.marks.unwired} 件（unwired.mjs --node ${r.rel}）`);
    else if (r.marks.cands) tags.push(`配線候補 なし（弱い候補 ${r.marks.cands} 件）`);
    if (r.marks.dead) tags.push(`★${r.marks.dead}（この語を誰も使っていない）`);
    if (tags.length) console.log(`                  印: ${tags.join(' / ')}`);
  }
  if (r.items.length) console.log(`                  項目: ${r.items.slice(0, 4).join(' / ')}${r.items.length > 4 ? ' …' : ''}`);
  if (!r.startable) {
    console.log(`                  ★上流の sorry ${r.blockers.length} 件で止まる:`);
    for (const b of r.blockers.slice(0, 3)) {
      const via = mathVia.get(`${r.mod}→${b}`);
      const tag = r.mathOnly.includes(b) ? `  ★import に無い依存${via ? `（${via}）` : ''}` : '';
      console.log(`                    ${byMod.get(b)?.rel ?? b}${tag}`);
    }
    if (r.blockers.length > 3) console.log(`                    … 他 ${r.blockers.length - 3} 件`);
  }
}
console.log();
if (limit > 0 && shown > limit) {
  console.log(`  … 他 ${shown - limit} 件（★既定は 5 件まで。全部見るなら --limit 0）`);
  console.log();
}
console.log('☆次の一手は「下流」が大きい startable なノード。');
console.log('  そのノードの持ち場を sub-agent に渡すには node tools/brief.mjs --node <rel>');
console.log('★同時に起動する agent は 5 個まで。足りなければ次の波にする。');
console.log('  前線が 5 件も無いときに agent を増やしても、上流の sorry に当たって止まる。');
