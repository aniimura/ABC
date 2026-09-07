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

/* ══════════════════════════════════════════════════════════════════════════
 * ★★★M167 —— **この一覧の単位はファイルであって、配る持ち場ではない**
 *
 * ★背景（メタ第 33 回 M166 が「説明できない」と書いた食い違い）:
 *   同じ日に「lean-prover は走行時間の 54.2% で同時 2 本（天井）」なのに、
 *   ここは「保留でも空撃ちでもない着手可能は **1 件**」と出していた。
 *
 * ★メタ第 34 回が測った答え（★どちらの数字も正しい。★**単位が違うだけ**）:
 *   2026-09-07 に配った lean-prover **67 本**のうち、
 *     ・**65 本(97.0%)** は、その日の始まりに `frontier` が出していた
 *       着手可能にも着手不可の sorry ノードにも **1 つも当たらない**本を書いていた。
 *     ・**52 本**は、辿ると **`Skeleton/PGC/Section1` ただ 1 つ**の下流だった。
 *   ⇒ ★★**1 つの startable ノードから、その日のうちに 52 本の持ち場が切り出されている。**
 *   （2026-09-06 も同じ形: 39 本中 33 本(84.6%)が一覧の外、当たったのは 5 本）
 *
 * ⇒ ★★**「着手可能 N 件」を供給量として読んではいけない。**
 *   ★配る単位は「ノードの中に残っている `sorry` 1 つ」であり、その 1 つが
 *     `Found/…` の補題 1 本（= agent 1 本）になる。
 *   ★だからここは **`sorry` の残り本数**を併記する（★`startable` の判定は 1 ビットも変えない）。
 * ══════════════════════════════════════════════════════════════════════════ */

/** コメントと文字列リテラルを潰す。★`graph.mjs` の同名関数と**同じ規則**にしてある
 *  （潰さないと `.needs` の説明文「… sorry 0 …」を本物と取り違える）。 */
export const stripComments = (s) => s
  .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
  .replace(/--[^\n]*/g, ' ')
  .replace(/"(?:[^"\\\n]|\\.)*"/g, '""');

/** ★ソースに残っている `sorry` の**個数**。★`graph.mjs` の `hasSorry` は真偽しか返さない。
 *  ★実測(2026-09-08、startable 6 ノード): 潰さないと **50**、潰すと **28**。★22 件が注釈の中だった。 */
export function countSorry(src) {
  return (stripComments(src).match(/(?<![\w.])sorry(?![\w])/g) ?? []).length;
}

/** ★配線を試験できるようにするための薄い層（`read` を注入する）。
 *  ★読めなければ **null**（ゲートを落とさない。表示は `?`）。 */
export function openSorryOf(rel, read) {
  try { return countSorry(read(rel)); } catch { return null; }
}

/** ★「残り」欄の字面。★null は `?`。 */
export const fmtOpen = (v) => (v === null || v === undefined ? '?' : String(v));

/* ══════════════════════════════════════════════════════════════════════════
 * ★★★M187 —— `--jobs`「配れる持ち場の**在庫**」（メタ第 36 回）
 *
 * ★上の一覧（既定）は **ファイル**を数える。M167/M173 が「それは供給量ではない」と
 *   書いたので、ここは **宣言粒度**で「いま名指しできる持ち場」を並べる。
 *
 * ★★ただし実測の答えは「**在庫では供給量に届かない**」である（下の J0 の断り）。
 *   ★測り方は `ResearchPaper/meta-backlog.md` M187 に式ごと書いてある。
 * ══════════════════════════════════════════════════════════════════════════ */

/** コメントだけ潰す（文字列は残す）。★`.needs` の説明文を読むために要る。 */
export const stripCommentsKeepStrings = (s) => s
  .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
  .replace(/--[^\n]*/g, ' ');

/** 宣言の頭（列 0 から始まるキーワード）で切る。★行番号を保つ。 */
export const DECL_HEAD =
  /^(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|scoped\s+|noncomputable\s+|nonrec\s+|partial\s+|unsafe\s+)*(theorem|lemma|def|abbrev|instance|structure|class|inductive)\s+([^\s({[:]+)/;

export function splitDecls(src) {
  const lines = stripCommentsKeepStrings(src).split('\n');
  const out = [];
  let cur = null;
  for (let i = 0; i < lines.length; i++) {
    const m = DECL_HEAD.exec(lines[i]);
    if (m) { if (cur) out.push(cur); cur = { kind: m[1], name: m[2], line: i + 1, text: lines[i] }; }
    else if (cur) cur.text += '\n' + lines[i];
  }
  if (cur) out.push(cur);
  return out;
}

/** ★J1 —— `sorry` を持つ**宣言**。★`.needs` の説明文の中の `sorry` は数えない。 */
export function sorryDecls(src) {
  return splitDecls(src)
    .map((d) => ({ ...d, sorry: countSorry(d.text) }))
    .filter((d) => d.sorry > 0);
}

export const OBLIGATION_KINDS = ['citation', 'folklore', 'implicitStep', 'otherPaper', 'derivation'];

/** ★J2 —— `def X.needs : List ProofObligation := [ … ]` の中身。
 *  ★戻り値は `{ owner, kind, what }`。`what` は最初の 2 つの文字列リテラルを繋いだもの。 */
export function needsEntries(src) {
  const out = [];
  for (const d of splitDecls(src)) {
    if (!/\.needs$/.test(d.name)) continue;
    const owner = d.name.replace(/\.needs$/, '');
    for (const m of d.text.matchAll(
      /\.(citation|folklore|implicitStep|otherPaper|derivation)\s+((?:"(?:[^"\\]|\\.)*"\s*){0,2})/g)) {
      const what = (m[2].match(/"(?:[^"\\]|\\.)*"/g) ?? []).map((s) => s.slice(1, -1)).join(' — ');
      out.push({ owner, kind: m[1], what });
    }
  }
  return out;
}

/** `.cache/decl-index.txt` の 1 行を割る。 */
export function parseIndexLine(line) {
  const f = line.split('\t');
  if (f.length < 3) return null;
  const loc = f[2];
  return { kind: f[0].trim(), full: f[1], short: f[1].split('.').pop(),
           loc, rel: loc.slice(0, loc.lastIndexOf(':')), stmt: f[3] ?? '' };
}

/** 宣言の型を「束縛子」と「結論」に割る（トップレベルの最初の `:` で切る）。 */
export function splitSignature(stmt) {
  const m = stmt.match(/^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|scoped\s+|noncomputable\s+|nonrec\s+|partial\s+|unsafe\s+)*(?:theorem|lemma|def|abbrev|instance|structure|class|inductive)\s+\S*\s*/);
  const rest = m ? stmt.slice(m[0].length) : stmt;
  let depth = 0;
  for (let i = 0; i < rest.length; i++) {
    const c = rest[i];
    if ('([{⦃'.includes(c)) depth++;
    else if (')]}⦄'.includes(c)) depth--;
    else if (depth === 0 && c === ':' && rest[i + 1] !== '=' && rest[i - 1] !== ':') {
      return { binders: rest.slice(0, i), concl: rest.slice(i + 1).trim() };
    }
  }
  return { binders: rest, concl: '' };
}

/** トップレベルの束縛子を 1 つずつ。 */
export function eachBinder(binders) {
  const out = [];
  let depth = 0, start = -1, open = '';
  for (let i = 0; i < binders.length; i++) {
    const c = binders[i];
    if ('([{⦃'.includes(c)) { if (depth === 0) { start = i; open = c; } depth++; }
    else if (')]}⦄'.includes(c)) {
      depth--;
      if (depth === 0 && start >= 0) { out.push({ open, text: binders.slice(start + 1, i) }); start = -1; }
    }
  }
  return out;
}

/** ★J3 —— 「開いた仮説」= `def X … : Prop` であって
 *   (a) どこかで**仮説として使われ**ており（需要 ≥ 1）、かつ
 *   (b) **無条件の**証明（束縛子に仮説を 1 つも持たない theorem で結論が `X …`）が無いもの。
 *  ★(b) を「theorem が 1 本も無い」にすると `HasCoherentFunctional` を取り落とす
 *    （条件つきの `_of_smoothModelCarrier` が在るため）。★実測でそれを確かめてある。 */
export function openPropObligations(recs) {
  const propNames = new Set();
  for (const r of recs) {
    if ((r.kind === 'def' || r.kind === 'abbrev') && splitSignature(r.stmt).concl === 'Prop') propNames.add(r.short);
    else if (r.kind === 'structure' && /:\s*Prop\s*$/.test(r.stmt)) propNames.add(r.short);
  }
  const LOGIC = /[∀∃→↔=≠≤≥∈∉¬∣∧∨]/;
  const isHyp = (b) => {
    if (b.open === '[') return false;                    // インスタンス束縛子は仮説にしない
    const i = b.text.indexOf(':');
    if (i < 0) return false;
    const ty = b.text.slice(i + 1).trim();
    const head = (ty.match(/([A-Za-z_][A-Za-z0-9_'.]*)/) ?? [])[1];
    return (head && propNames.has(head)) || LOGIC.test(ty);
  };
  const defs = recs.filter((r) => (r.kind === 'def' || r.kind === 'abbrev')
    && splitSignature(r.stmt).concl === 'Prop');
  const stat = new Map();
  for (const d of defs) stat.set(d.short, { rec: d, uncond: 0, cond: 0, demand: 0 });
  for (const r of recs) {
    const { binders, concl } = splitSignature(r.stmt);
    const bs = eachBinder(binders);
    const head = (concl.match(/^([A-Za-z_][A-Za-z0-9_'.]*)/) ?? [])[1];
    if (head && stat.has(head) && (r.kind === 'theorem' || r.kind === 'lemma')) {
      if (bs.some(isHyp)) stat.get(head).cond++; else stat.get(head).uncond++;
    }
    for (const b of bs) {
      if (!isHyp(b)) continue;
      const i = b.text.indexOf(':');
      const ty = b.text.slice(i + 1);
      for (const mm of ty.matchAll(/([A-Za-z_][A-Za-z0-9_'.]*)/g)) {
        if (stat.has(mm[1])) { stat.get(mm[1]).demand++; break; }
      }
    }
  }
  return [...stat.values()]
    .filter((s) => s.uncond === 0 && s.demand > 0)
    .map((s) => ({ name: s.rec.short, loc: s.rec.loc, rel: s.rec.rel, demand: s.demand, cond: s.cond }))
    .sort((a, b) => b.demand - a.demand || a.name.localeCompare(b.name));
}

if (flag('--selftest')) {
  let pass = 0, fail = 0;
  const t = (name, got, want) => {
    if (JSON.stringify(got) === JSON.stringify(want)) { pass++; }
    else { fail++; console.log(`  ✗ ${name}: got ${JSON.stringify(got)} want ${JSON.stringify(want)}`); }
  };
  t('素の sorry 1 個', countSorry('theorem a : True := sorry'), 1);
  t('同じ行に 2 個', countSorry('  sorry  sorry'), 2);
  t('行コメントの中は数えない', countSorry('-- sorry\ntheorem a := sorry'), 1);
  t('ブロックコメントの中は数えない', countSorry('/- sorry sorry -/\nsorry'), 1);
  t('文字列の中は数えない', countSorry('def s := "sorry 0"'), 0);
  t('`.sorry` は数えない（射影）', countSorry('h.sorry'), 0);
  t('`sorryAx` は数えない', countSorry('sorryAx Nat'), 0);
  t('`resorry` は数えない', countSorry('resorry'), 0);
  t('無ければ 0', countSorry('theorem a : True := trivial'), 0);
  t('空文字列', countSorry(''), 0);
  t('入れ子のブロックコメントの後ろは数える', countSorry('/- x -/ sorry'), 1);
  t('コメント潰しは行数を変えない', stripComments('/- a\nb -/\nx').split('\n').length, 3);
  t('openSorryOf は read の中身を数える', openSorryOf('X', () => 'sorry sorry'), 2);
  t('openSorryOf は rel を read に渡す', openSorryOf('A/B.lean', (r) => (r === 'A/B.lean' ? 'sorry' : '')), 1);
  t('openSorryOf は読めなければ null', openSorryOf('X', () => { throw new Error('nope'); }), null);
  t('fmtOpen(null) は ?', fmtOpen(null), '?');
  t('fmtOpen(0) は 0（? ではない）', fmtOpen(0), '0');

  /* ── ★M187 `--jobs` ── */
  const SRC1 = 'theorem a : True := sorry\n\ndef b : Nat := 1\n\ntheorem c : True := trivial\n';
  t('J1: 宣言に割る', splitDecls(SRC1).map((d) => d.name), ['a', 'b', 'c']);
  t('J1: 行番号を保つ', splitDecls(SRC1).map((d) => d.line), [1, 3, 5]);
  t('J1: sorry を持つ宣言だけ', sorryDecls(SRC1).map((d) => d.name), ['a']);
  t('J1: sorry の個数', sorryDecls('theorem a : True := by\n  sorry\n  sorry\n')[0].sorry, 2);
  t('J1: 注釈の中の sorry は宣言にしない', sorryDecls('/- sorry -/\ndef a : Nat := 1\n').length, 0);
  t('J1: 列 0 でない行は宣言の頭にしない',
    splitDecls('def a : Nat :=\n  def b := 1\n').map((d) => d.name), ['a']);
  t('J1: noncomputable/@[simp] を剥がす',
    splitDecls('@[simp] noncomputable def a : Nat := 1\n').map((d) => d.name), ['a']);
  const SRC2 = 'theorem p : True := sorry\n\ndef p.needs : List ProofObligation :=\n' +
    '  [ .citation "[6] Serre" "上付き↔下付き" (.absent "0 件") 5,\n' +
    '    .folklore "well-known" 3,\n    .otherPaper "pGC" "Corollary 1.3" 3 ]\n';
  t('J2: .needs を数える', needsEntries(SRC2).length, 3);
  t('J2: 種別', needsEntries(SRC2).map((e) => e.kind), ['citation', 'folklore', 'otherPaper']);
  t('J2: 持ち主', needsEntries(SRC2)[0].owner, 'p');
  t('J2: 説明文を文字列から取る', needsEntries(SRC2)[0].what, '[6] Serre — 上付き↔下付き');
  t('J2: .absent の中の文字列を拾わない', needsEntries(SRC2)[1].what, 'well-known');
  t('J2: needs でない宣言は無視', needsEntries('def q : Nat := 1\n').length, 0);
  /* ★突然変異「.needs 以外も拾う」を殺す —— .needs でない宣言の中に .citation を置く */
  t('J2: .needs でない宣言の中の .citation は拾わない',
    needsEntries('def q : List X :=\n  [ .citation "a" "b" 1 ]\n').length, 0);
  /* ★突然変異「文字列を 3 つまで取る」を殺す */
  t('J2: 説明文に使う文字列は 2 つまで',
    needsEntries('def p.needs : List ProofObligation :=\n  [ .citation "a" "b" "c" 1 ]\n')[0].what, 'a — b');
  /* ★突然変異「種別の一覧を削る」を殺す（この一覧は集計欄の見出しに使う） */
  t('J2: 種別の一覧は 5 つ', OBLIGATION_KINDS,
    ['citation', 'folklore', 'implicitStep', 'otherPaper', 'derivation']);
  t('J2: 注釈の中の .citation は数えない',
    needsEntries('def p.needs : List ProofObligation :=\n  /- .citation "x" "y" 1 -/ []\n').length, 0);
  t('sig: 束縛子と結論に割る', splitSignature('theorem f (a : Nat) : True').concl, 'True');
  t('sig: 括弧の中の : では切らない',
    splitSignature('def g (h : A = B) (k : C) : Prop').binders.trim(), '(h : A = B) (k : C)');
  t('sig: := では切らない', splitSignature('def g : Nat := 1').concl, 'Nat := 1');
  /* ★突然変異「:= でも切る」を殺す —— 型注釈が無い定義では `:=` の `:` が最初に来る */
  t('sig: 型注釈が無ければ結論は空', splitSignature('def g := 1').concl, '');
  t('sig: 束縛子を 1 つずつ',
    eachBinder(' (a : Nat) {b : M} [Fact p] ').map((b) => b.open), ['(', '{', '[']);
  const IDX = [
    'def      \tN.HasX\tFound/A.lean:10\tdef HasX (K : T) : Prop',
    'def      \tN.IsY\tFound/A.lean:20\tdef IsY (K : T) : Prop',
    'def      \tN.HasZ\tFound/A.lean:30\tdef HasZ (K : T) : Prop',
    'theorem  \tN.hasX_of_isY\tFound/B.lean:5\ttheorem hasX_of_isY (K : T) (h : IsY K) : HasX K',
    'theorem  \tN.isY_all\tFound/B.lean:9\ttheorem isY_all (K : T) : IsY K',
    'theorem  \tN.use\tFound/C.lean:3\ttheorem use (K : T) (h : HasX K) (h2 : HasZ K) : True',
  ].map(parseIndexLine);
  const OP = openPropObligations(IDX);
  t('J3: 条件つきの証明しか無い Prop は「開いている」', OP.map((o) => o.name).includes('HasX'), true);
  t('J3: 無条件の証明がある Prop は載せない', OP.map((o) => o.name).includes('IsY'), false);
  t('J3: 証明が 1 本も無い Prop も載せる', OP.map((o) => o.name).includes('HasZ'), true);
  t('J3: 需要を数える', OP.find((o) => o.name === 'HasX').demand, 1);
  t('J3: 条件つきの本数を数える', OP.find((o) => o.name === 'HasX').cond, 1);
  t('J3: 需要 0 は載せない',
    openPropObligations([parseIndexLine('def      \tN.Lonely\tF/A.lean:1\tdef Lonely : Prop')]).length, 0);
  /* ★★設計の判断: インスタンス束縛子は **Prop であっても**仮説と数えない
   *   （`[Fact p.Prime]` `[DecidableEq α]` を仮説にすると木のほぼ全部が「条件つき」になる）。
   *   ★試験は `[_i : IsY K]`（IsY は Prop の def）で書く —— `[Fact p]` だと `:` が無く、
   *   ★規則を消しても素通りする（メタ第 36 回の突然変異で実際に素通りした）。 */
  t('J3: インスタンス束縛子は Prop でも仮説にしない',
    openPropObligations([...IDX,
      parseIndexLine('theorem  \tN.inst\tF/D.lean:1\ttheorem inst (K : T) [_i : IsY K] : HasZ K')])
      .map((o) => o.name).includes('HasZ'), false);
  /* ★突然変異「論理記号を見ない」を殺す —— 型が `1 = 1` の束縛子は仮説である */
  t('J3: 論理記号を含む束縛子は仮説',
    openPropObligations([...IDX,
      parseIndexLine('theorem  \tN.cnd\tF/D.lean:2\ttheorem cnd (K : T) (h : 1 = 1) : HasZ K')])
      .map((o) => o.name).includes('HasZ'), true);
  t('parseIndexLine: rel は行番号を落とす', parseIndexLine('def\tN.a\tFound/A.lean:10\tdef a : Prop').rel, 'Found/A.lean');
  t('parseIndexLine: 壊れた行は null', parseIndexLine('x'), null);
  /* ★突然変異「壊れた行を通す」を殺す —— 3 欄無い行は使えない */
  t('parseIndexLine: 2 欄しか無い行も null', parseIndexLine('def\tN.a'), null);

  console.log(`frontier.mjs selftest: ${pass}/${pass + fail} PASS`);
  process.exit(fail ? 1 : 0);
}

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
  /** ★M167: そのノードに残っている `sorry` の個数 = **そこから切り出せる持ち場の数**。
   *  ★読めなければ null（ゲートを落とさない。表示は `?` になる）。 */
  const openSorry = openSorryOf(n.rel, (r) => readFileSync(join(ROOT, 'lean', 'ABC3', r), 'utf8'));
  rows.push({
    mod: n.mod,
    rel: n.rel,
    openSorry,
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

/* ══════════════════════════════════════════════════════════════════════════
 * ★★★M187 `--jobs` —— 「配れる持ち場の**在庫**」を宣言粒度で並べる
 * ══════════════════════════════════════════════════════════════════════════ */
if (flag('--jobs')) {
  const startableRel = new Set(rows.filter((r) => r.startable).map((r) => r.rel));
  const ownerOf = new Map(rows.map((r) => [r.rel, r.owner]));
  const relOwner = new Map(nodes.map((n) => [n.rel, n.owner]));
  const leanRoot = join(ROOT, 'lean', 'ABC3');
  const wantOwner = opt('--owner');
  const lim = Number(opt('--limit') ?? 10);

  // ── J1 / J2: 木を 1 回だけ舐める ──
  const j1 = [], j2 = [];
  for (const n of nodes) {
    if (!n.hasSorry) continue;
    let src; try { src = readFileSync(join(leanRoot, n.rel), 'utf8'); } catch { continue; }
    const sd = sorryDecls(src);
    const withSorry = new Set(sd.map((d) => d.name));
    for (const d of sd) {
      j1.push({ rel: n.rel, owner: n.owner, bucket: n.bucket, name: d.name, line: d.line,
                sorry: d.sorry, startable: startableRel.has(n.rel) });
    }
    for (const e of needsEntries(src)) {
      if (!withSorry.has(e.owner)) continue;   // ★もう証明済みの宣言の .needs は持ち場ではない
      j2.push({ rel: n.rel, owner: n.owner, decl: e.owner, kind: e.kind, what: e.what,
                startable: startableRel.has(n.rel) });
    }
  }

  // ── J3: `.cache/decl-index.txt`（在庫の口。CLAUDE.md 記載）を作り直して読む ──
  let j3 = [], j3note = '';
  try {
    execFileSync('node', [join(ROOT, 'tools', 'decl-index.mjs')], { cwd: ROOT, stdio: 'ignore' });
    const idx = readFileSync(join(ROOT, '.cache', 'decl-index.txt'), 'utf8')
      .split('\n').filter(Boolean).map(parseIndexLine).filter(Boolean);
    j3 = openPropObligations(idx).map((o) => ({ ...o, owner: relOwner.get(o.rel) ?? '?' }));
  } catch (e) { j3note = `★decl-index を作れなかった（${String(e.message).slice(0, 60)}）`; }

  const pick = (xs) => (wantOwner ? xs.filter((x) => x.owner === wantOwner) : xs);
  const J1 = pick(j1), J2 = pick(j2), J3 = pick(j3);
  const kindTally = OBLIGATION_KINDS
    .map((k) => `${k} ${J2.filter((e) => e.kind === k).length}`).join(' / ');

  if (flag('--json')) {
    console.log(JSON.stringify({ generated: new Date().toISOString(),
      note: '★これは在庫であって供給量ではない（M187）', j1: J1, j2: J2, j3: J3 }, null, 1));
    process.exit(0);
  }

  console.log('★配れる持ち場の**在庫**（--jobs、宣言粒度）');
  console.log('  ★★これは**在庫**であって**供給量ではない**。★下限である。');
  console.log('  ★実測(M187。2026-09-08 の波 = lean-prover **14 本**。波の開始 = commit 3beec898):');
  console.log('    ・**持ち場そのもの**がこの在庫に載っていた …… **3 本 / 14**');
  console.log('      （prop_2_1・prop_2_2 は着手可能／cor_3_1&cor_3_3 は**着手不可**の J1 だった）');
  console.log('    ・**その持ち場が仕えるゴール**を brief が名指し … **8 本 / 14 = 57.1%**');
  console.log('      （★比較: 既定のファイル一覧では 0 本 / 73 本だった —— M183）');
  console.log('    ・**6 本 / 14 = 43%** は**名前そのものが波の開始時点に無かった** ——');
  console.log('      前の持ち場が次の仮説を作る（SmoothModelCarrier → HasCoherentFunctional → …）');
  console.log('    ・残る 5 本は sorry でも .needs でも開いた仮説でもない');
  console.log('      （外部引用の鋭い形 / 無名の束縛子の仮説 / 定義の移設）');
  console.log('  ⇒ ★★供給は**在庫ではなく流量**である。「いま何本配れるか」にこの数で答えてはいけない。');
  console.log('    ★使い方は「ゴールを選ぶ」。★1 波分の本数はここからは出ない。');
  console.log('');
  console.log(`  J1 sorry を持つ**宣言**            ${String(J1.length).padStart(4)} 件` +
    `（うち着手可能なファイル ${J1.filter((d) => d.startable).length} 件 / sorry 合計 ${J1.reduce((a, d) => a + d.sorry, 0)}）`);
  console.log(`  J2 その宣言が書いた .needs         ${String(J2.length).padStart(4)} 件（${kindTally}）`);
  console.log(`  J3 開いた仮説（def _ : Prop）      ${String(J3.length).padStart(4)} 件` +
    `${j3note ? '  ' + j3note : '（需要 ≥ 1 かつ無条件の証明が無い）'}`);
  if (wantOwner) console.log(`  所属で絞り込み: ${wantOwner}`);
  console.log('');

  const show = (title, xs, fmt) => {
    console.log(`  ── ${title}`);
    for (const x of xs.slice(0, lim > 0 ? lim : xs.length)) console.log('   ' + fmt(x));
    if (lim > 0 && xs.length > lim) console.log(`   … 他 ${xs.length - lim} 件（--limit 0 で全部）`);
    console.log('');
  };
  show('J1 sorry を持つ宣言（着手可能なファイルのものだけ）',
    J1.filter((d) => d.startable).sort((a, b) => b.sorry - a.sorry || a.rel.localeCompare(b.rel)),
    (d) => `${String(d.sorry).padStart(2)} sorry  ${d.name}  ${d.rel}:${d.line}  (${d.owner})`);
  show('J2 木自身が「これが要る」と書いたもの（着手可能なファイルのものだけ）',
    J2.filter((e) => e.startable),
    (e) => `${e.kind.padEnd(13)} ${e.decl}  ${e.what.slice(0, 70)}  (${e.owner})`);
  show('J3 開いた仮説（需要の大きい順。★条件つきの証明しか無いものを含む）',
    J3, (o) => `需要${String(o.demand).padStart(4)} 条件つき${String(o.cond).padStart(2)}  ${o.name}  ${o.loc}  (${o.owner})`);

  console.log('  ★★測れないこと（正直に）:');
  console.log('    ・J3 の「仮説」は**名前が付いているもの**だけ。★束縛子に直接書かれた無名の仮説');
  console.log('      （例 `absGalStageFiltration (hcompat : ∀ ⦃M⦄ …)`）は数えていない。実測 1,002 件あり、');
  console.log('      ★その大半は正当な引数なので、いまの規則では持ち場と区別できない。');
  console.log('    ・J3 は「定義の述語」（`IsCartierDiv` 等、全称では成り立たない分類述語）を');
  console.log('      ★**除けていない**。所属で絞って読むこと（--owner pGC）。');
  console.log('    ・定義の移設・配線は 3 つのどれにも載らない（実測 14 本中 1 本）。');
  console.log('    ・★**この一覧は「何本配れるか」に答えていない。**「どのゴールが開いているか」に答えている。');
  console.log('      ★再現手順は `ResearchPaper/meta-backlog.md` M187 に式ごと書いてある。');
  process.exit(0);
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
const startRows = rows.filter((r) => r.startable);
const openTotal = startRows.reduce((a, r) => a + (r.openSorry ?? 0), 0);
console.log('★前線（sorry を持つノード）');
console.log(`  sorry ノード ${total} / うち着手可能 ${startable}`);
/* ★★M167: この 3 行が無いと「着手可能 N」を供給量と読んでしまう（M166 が実際に読みかけた）。 */
console.log(`  ★★単位に注意 —— ここは**ファイル**を数えている。★これは供給量ではない`);
console.log(`    着手可能なノードに残る sorry は ${openTotal} 件（ノード ${startable} 件ではない）。` +
  `★★そしてこれも **下限**である`);
console.log(`    実測(2026-09-07): その日の始まりに \`Skeleton/PGC/Section1\` の残りは **1 件**。` +
  `★そこから **52 本**の持ち場が切り出された（配った 67 本中 65 本 = 97% はこの一覧に載らない本を書いた）`);
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
console.log('  下流  項目  残り  ノード      （残り = そのノードに残る sorry ＝ 切り出せる持ち場の数）');
for (const r of sel) {
  const mark = r.startable ? ' ' : '×';
  const os = fmtOpen(r.openSorry);
  console.log(`${mark} ${String(r.downstream).padStart(5)} ${String(r.dsItems).padStart(5)} ${os.padStart(5)}  ${r.rel}  (${r.owner})`);
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
console.log(`★★「同時に何本走らせられるか」を **この一覧の件数で決めてはいけない**（M167）。`);
console.log(`  ★実測: 残り 1 件のノードから 1 日に 52 本の持ち場が出た。★上限の議論にこの数を使わないこと。`);
