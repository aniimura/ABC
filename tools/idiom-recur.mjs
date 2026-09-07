#!/usr/bin/env node
// idiom-recur.mjs —— 「idiom を書いた後に、同じ失敗形の字面がログに再出現したか」を数える。
//
// 測っているのは 1 つだけ:
//   idiom S の登録時刻 T(S) より後に出た Lean の `error:` 診断のうち、
//   S の「字面」(= S 自身が名指ししている Mathlib 宣言名 / S が引用しているエラー文)を
//   本文に含むものの件数。
//
// ★これは「S が防いだ / 防がなかった」ではない。字面が一致した件数だけである。
//   ・0 でも「当たらなかっただけ」かもしれない(上界 0 とだけ言える)
//   ・逆に 0 でないからといって S の失敗形とは限らない(名指しして目で見ること)
//   ・★既知の偽陰性がある: 再発時の Lean のエラー文は「idiom が名指しした補題」ではなく
//     「詰まった現場の項」を印字するため、字面が繋がらないことがある(--calibrate 参照)。
//
// ★★メタ第 26 回で digest を v2 にした。★変わったのは 3 つ:
//   (1) ★**MCP REPL(`lean_check`)の診断を読むようになった** —— 書式が `error:` ではなく
//       `error <行>:<桁>` なので、v1 は **1 件も**読めていなかった。
//       ★実測 4,835 → **13,784 件**(lake 4,838 / repl 8,946)。★**65% を取り落としていた。**
//       ★これは「遅い経路(lake build)だけを測っていた」ということでもある
//       —— CLAUDE.md が勧める速い経路(lean_check)が丸ごと見えていなかった。
//   (2) ★各診断に**出所**を付けた(`doc` / `scratch` / `tree`。`provenance` 参照)。
//       ★族の表から「Lean のエラーですらないもの」を機械で外すため。
//   (3) ★既定の `--sig` を `both` から **`lit`** にした(M109 の結論に道具を合わせた)。
//
// ★★メタ第 27 回で digest を v3 にした。★変わったのは 2 つ:
//   (1) ★出所に **`meta`** を足した(M120) —— 改善係(`meta-optimizer`)の probe が印字した
//       診断の字面が、次の `--rescan` で「Lean のエラー」として数えられていた。
//       ★判定は **agent の種類**で行う(ファイル名で当てない)。probe の名前が変わっても漏れない。
//   (2) ★各診断に **それを出した agent の `toolUseId`** を付けた。
//       ★これは `agent-timing.mjs` が完了通知(`duration_ms`)を引く鍵と同じもの。
//       ⇒ ★`node tools/agent-timing.mjs --diag` で「族 × 所要時間」が突き合わせられる。
//
// 使い方:
//   node tools/idiom-recur.mjs --rescan        会話ログ(~/.claude/projects/D--Math-ABC3)を舐めて digest を作る
//   node tools/idiom-recur.mjs                 digest から数える(要 --rescan 済)
//   node tools/idiom-recur.mjs --hits          再出現があった idiom を名指しする
//   node tools/idiom-recur.mjs --family        失敗形の族の表(出所ごとに割る)
//   node tools/idiom-recur.mjs --family --sub "Invalid field"   ★族の中を下位に割る
//   node tools/idiom-recur.mjs --with-doc      文書を印字しただけの診断も入れる(既定は外す)
//   node tools/idiom-recur.mjs --with-meta     ★改善係自身の probe の出力も入れる(既定は外す。M120)
//   node tools/idiom-recur.mjs --sig lit|name|both              既定 lit
//   node tools/idiom-recur.mjs --calibrate     既知の実例で検出できるか見る(2 通り)
//   node tools/idiom-recur.mjs --selftest      内蔵 fixture
import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import readline from 'node:readline';
import { execFileSync } from 'node:child_process';

const ARGV = process.argv.slice(2);
const has = f => ARGV.includes(f);
const REPO = process.cwd();
const IDIOMS = path.join(REPO, 'tools', 'lean-idioms.md');
const CACHE = path.join(REPO, '.cache', 'idiom-recur-digest.json');
// ★digest の版。v1 は `error:` 書式だけ・出所なし。v2 で MCP REPL の診断と出所を足した。
//   ★v3(メタ第 27 回)で **診断を出した agent の `toolUseId`** と 出所 `meta` を足した。
//   ★版が違う digest は黙って使わない(古い数字で判定してしまうため)。
const DIGEST_V = 3;
const LOGROOT = path.join(os.homedir(), '.claude', 'projects', 'D--Math-ABC3');

// ---------- 純関数群(selftest はここだけを叩く) ----------

export const norm = s => String(s).replace(/\s+/g, ' ').trim();

export function headingKey(line) {
  let s = String(line).replace(/^#{1,6}\s*/, '').trim();
  s = s.replace(/（[^（）]*\d{4}-\d{2}-\d{2}[^（）]*）/g, '').replace(/\([^()]*\d{4}-\d{2}-\d{2}[^()]*\)/g, '');
  return s.replace(/[★☆]/g, '').replace(/\s+/g, ' ').trim();
}

// 見出しで節に割る。```fence``` の中の # は見出しにしない。
export function sections(md) {
  const out = []; let fence = false, cur = null;
  String(md).split(/\r?\n/).forEach((ln, i) => {
    if (/^\s*```/.test(ln)) fence = !fence;
    if (!fence && /^#{2,4}\s+\S/.test(ln)) { if (cur) out.push(cur); cur = { line: i + 1, head: ln.trim(), key: headingKey(ln), body: [] }; }
    else if (cur) cur.body.push(ln);
  });
  if (cur) out.push(cur);
  return out;
}

// Lean 診断の見出し行: "Path/F.lean:12:34: error: msg" / "error: msg"
export const DIAG = /^(?:(\S*\.lean):(\d+):(\d+):\s*)?(error|warning|info):\s?(.*)$/;

// tool_result の本文から error 診断だけを塊で取り出す
export function errorBlocks(text) {
  const out = []; let cur = null;
  const push = () => { if (cur) { out.push({ file: cur.file, msg: cur.buf.join('\n').replace(/\s+$/, '') }); cur = null; } };
  for (const ln of String(text).split(/\r?\n/)) {
    const m = DIAG.exec(ln);
    if (m) { push(); if (m[4] === 'error') cur = { file: m[1] || '', buf: [m[5]] }; }
    else if (cur && cur.buf.length < 40) cur.buf.push(ln);
  }
  push();
  return out;
}

// ★★MCP REPL(`mcp__abc3-lean__lean_check`)の診断は書式が違う —— `error:` ではなく
//   行頭の `error <行>:<桁>` で、本文は次の行から始まる。DIAG は 1 件も拾えない。
//   ⇒ メタ第 26 回の実測: この書式の error 診断が **8,939 件**あり、
//     `error:` 書式の 4,835 件と合わせて 13,774 件。**65% を取り落としていた。**
export const MCPDIAG = /^(error|warning|info)\s+(\d+):(\d+)\s*$/;
export function mcpErrorBlocks(text) {
  const out = []; let cur = null;
  const push = () => { if (cur) { out.push({ file: '', msg: cur.buf.join('\n').replace(/\s+$/, '') }); cur = null; } };
  for (const ln of String(text).split(/\r?\n/)) {
    const m = MCPDIAG.exec(ln.trim());
    if (m) { push(); if (m[1] === 'error') cur = { buf: [] }; }
    else if (cur && cur.buf.length < 40) cur.buf.push(ln);
  }
  push();
  return out.filter(b => b.msg.trim());
}

// ★★出所(provenance) —— 「その診断が何を検査して出たか」。族の表から
//   「失敗ではないもの」を機械で外すために要る。手で列挙すると次に増えたとき腐る。
//     doc     : 我々の文書(lean-idioms.md 等)を印字しただけ。★Lean は動いていない。
//     scratch : 検査したのが木の本ではない(REPL の断片 / `lean/ABC3/` の外の .lean)。
//               ★捨てられる試作。#158「同じ名前で書いて出させる」の合図は必ずここに落ちる。
//     tree    : `lean/ABC3/**` の本を検査した。
//   ★逆は言えない —— scratch にも本物の失敗はある。ここで言えるのは
//     「tree でない＝木は 1 行も損なわれていない」までである。
export const DOCFILES = /lean-idioms\.md|meta-backlog\.md|workflow-speedup\.md|decisions-pending\.md|orchestration\.md|autonomy-policy\.md|idiom-recur\.mjs/;
// ★`(lean[/\\])?` は冗長(どちらも `/ABC3/` で当たる)。
//   ★わざと壊しても selftest が黏ることを確かめたので落とした。
export const TREEFILE = /(^|[\/\\])ABC3[\/\\]/;
// ★★出所 `meta` —— **改善係(この道具を回す人)自身の出力**(M120。メタ第 27 回)。
//   ★メタ係は probe を回して診断の字面をそのまま印字する。その出力は会話ログに載り、
//   次の `--rescan` で「Lean のエラー」として数えられる。★第 26 回は `DOCFILES` に
//   `idiom-recur.mjs` を足して凌いだが、★**probe のファイル名が変われば漏れる。**
//   ⇒ ★ファイル名ではなく **agent の種類**で決める。`meta-optimizer` の会話から出た診断は
//     すべて `meta`。★これは字面に依存しないので、probe を何と名付けても漏れない。
//   ★`doc` より先に見る: メタ係の probe は文書名を含まないことがある(数だけ印字するなど)。
export const METATYPES = new Set(['meta-optimizer']);
export function provenance(toolName, cmd, file, agentType) {
  if (agentType && METATYPES.has(String(agentType))) return 'meta';
  if (DOCFILES.test(String(cmd || ''))) return 'doc';
  if (toolName === 'mcp__abc3-lean__lean_check') return 'scratch';
  if (file && TREEFILE.test(file)) return 'tree';
  if (file) return 'scratch';                       // 木の外の .lean(scratchpad 等)
  // ★file が空の診断(裸の `error: ...`)は命令で決める。`lake build` は木の検査。
  const c = String(cmd || '');
  return (/lake\s+build/.test(c) || TREEFILE.test(c)) ? 'tree' : 'scratch';
}

// ★族の下位割り —— `Invalid field` のように Lean の汎用エラーは
//   1 つの族に別々の原因が混ざる。族の字面の**直後の句**で機械的に割る。
export function subKey(msg, fam) {
  const i = String(msg).indexOf(fam);
  if (i < 0) return '';
  const tail = norm(String(msg).slice(i + fam.length));
  if (/^notation:/.test(tail)) return 'notation';
  const m = /does not contain `([^`]+)`/.exec(tail);
  if (m) { const f = m[1]; const j = f.lastIndexOf('.'); return 'env:' + (j > 0 ? f.slice(0, j) : f); }
  return tail.split(' ').slice(0, 3).join(' ');
}

// 「字面」その 1 —— idiom が名指ししている Mathlib 宣言名。
// backtick の中の識別子だけを見る。名前空間が 3 文字以上の大文字始まりで、全体 12 文字以上。
export const NAME_MIN = 12;
export function declNames(text) {
  const out = new Set();
  for (const m of String(text).matchAll(/`([^`\n]{2,160})`/g)) {
    for (const t of m[1].matchAll(/[A-Za-z_][A-Za-z0-9_'.]*[A-Za-z0-9_']/g)) {
      const id = t[0];
      if (id.length < NAME_MIN) continue;
      if (/\.(lean|mjs|md|json|py|txt|html|toml|yml)$/i.test(id)) continue;
      const i = id.indexOf('.');
      if (i < 0) continue;                       // 名前空間なしは捨てる(語が広すぎる)
      const ns = id.slice(0, i);
      if (ns.length < 3 || !/^[A-Z]/.test(ns)) continue;
      out.add(id);
    }
  }
  return out;
}

// 「字面」その 2 —— idiom が引用している Lean のエラー文そのもの。
export const LIT_MIN = 30;
// ★★v1 —— メタ第 37 回まで唯一の版。★**凍結する**(消さない)。
//   ★理由: M182(メタ第 35 回)は `MCP_INFRA_RE_V1` を残してあったので「なぜ動かないか」が割れた。
//   ★同じ作法で残す。`--errwords v1` で変更前の数字がいつでも再現できる。
export const ERRWORDS_V1 = /failed to synthesize|could not synthesize|type mismatch|unknown identifier|unknown constant|function expected|motive is not type correct|no goals|unsolved goals|maximum recursion depth|timeout|invalid field|invalid projection|ambiguous, possible interpretations|made no progress|has already been declared|does not contain/i;
// ★★v2 —— v1 + **Lean 4 が書式を変えた 6 語**(M189 が名指し / M192 で事前登録)。
//   ★M189 の実測: ログの実エラー 12,456 件のうち 4,096 件(32.9%)が v1 のどの語にも当たらず、
//     上位は ``Tactic `rewrite` failed``(809)/ `Lean exited with code 1`(744)/
//     ``` `exact?` could not close the goal ```(272)。★旧書式 `rewrite tactic failed` は v1 にも無い。
//   ★★ただし `ERRWORDS` の当て先は **`tools/lean-idioms.md` の節**であって
//     ログのエラー本文ではない(M189 が呼び出し元を機械で辿った)。★4,096 件が直接動くのではない。
//   ★`tactic [^ ]{1,40} failed` は M189 の提案 `Tactic .* failed` を絞った形
//     (`.*` は貪欲で 240 字の引用を跨ぐ。`` `rewrite` `` は空白を含まないのでこれで届く)。
export const ERRWORDS_V2 = /failed to synthesize|could not synthesize|type mismatch|unknown identifier|unknown constant|function expected|motive is not type correct|no goals|unsolved goals|maximum recursion depth|timeout|invalid field|invalid projection|ambiguous, possible interpretations|made no progress|has already been declared|does not contain|tactic [^ ]{1,40} failed|exited with code|could not close the goal|failed to compile definition|unexpected token|unknown namespace/i;
export const ERRWORDS_VER = (() => {
  const i = ARGV.indexOf('--errwords');
  const v = i >= 0 ? ARGV[i + 1] : 'v2';
  if (v !== 'v1' && v !== 'v2') { console.error(`--errwords は v1 か v2(受け取った: ${JSON.stringify(v)})`); process.exit(2); }
  return v;
})();
export const ERRWORDS = ERRWORDS_VER === 'v1' ? ERRWORDS_V1 : ERRWORDS_V2;
export function errLits(text, re = ERRWORDS) {
  const out = new Set();
  for (const m of String(text).matchAll(/`([^`\n]{12,240})`/g)) {
    const n = norm(m[1]);
    if (n.length >= LIT_MIN && re.test(n)) out.add(n);
  }
  let fence = false;
  for (const ln of String(text).split(/\r?\n/)) {
    if (/^\s*```/.test(ln)) { fence = !fence; continue; }
    if (!fence) continue;
    const n = norm(ln);
    if (n.length >= LIT_MIN && re.test(n)) out.add(n);
  }
  return out;
}
// ★★「エラーらしい引用」—— ★M189 の「23 節」を**再現可能にする**ための規則(M188 の規約)。
//   ★M189 が使った規則は台帳に残っていないので、★これは**別の規則**である。23 と一致する保証はない。
//   ★規則: LIT_MIN 以上の引用(backtick の内側 / fence の中の行)であって、下の語のどれかを含むもの。
export const ERRLIKE = /error|failed|unknown|invalid|unexpected|cannot|could not|no goals|mismatch|ambiguous|expected|does not|not type correct|timeout|deprecated/i;
/** ★節の本文から「エラーらしい引用」を全部返す(字面として採れたかは問わない)。 */
export function errLikeQuotes(text) {
  const out = new Set();
  for (const m of String(text).matchAll(/`([^`\n]{12,240})`/g)) {
    const n = norm(m[1]);
    if (n.length >= LIT_MIN && ERRLIKE.test(n)) out.add(n);
  }
  let fence = false;
  for (const ln of String(text).split(/\r?\n/)) {
    if (/^\s*```/.test(ln)) { fence = !fence; continue; }
    if (!fence) continue;
    const n = norm(ln);
    if (n.length >= LIT_MIN && ERRLIKE.test(n)) out.add(n);
  }
  return out;
}

//   `lean-idioms.md` を読んだ / 探した記録があるか」。★あるなら仮説 (2)(引いても防げない)、
//   ★無いなら仮説 (1)(引かれていない)。★ログを歩くので 20 秒ほどかかる(digest を使わない)。
//   ★★書きの規則は rescan より**厳しい**(`>` / `>>` / `tee` を要求)。理由: `cat tools/lean-idioms.md`
//     は rescan では「書き」に落ちるが、★これは**読み**である。読みを取り落とすと仮説 (1) に
//     有利な方へ偏るので、★こちらは読みとして数える(自分に不利な側へ倒す)。
export const IDFILE = /lean-idioms\.md/i;
export const IDWRITE = /(>>?\s*["']?[^"'\s]*lean-idioms\.md)|(\btee\b[^\n]{0,80}lean-idioms\.md)/;
export const IDLOOK = /\b(grep|rg|sed|head|tail|awk|less|cat|wc|nl)\b/;
// ★★「読んだ」だけでは仮説を分けられない —— ★**何を探したか**で分ける。
//   ★`^## #25[012]` は「次の節番号を採りに行った」であって「似た節を探した」ではない。
//   ★見出しと番号だけからなる pattern を **numbering** とし、内容語を含むものだけ **content** とする。
export const STRUCTPAT = /^[\^$#\s0-9\[\]\-|()?*+.\\/]*$/;
/** ★引用符の外の `;` `|` `&` だけで命令を切る(`"motive\|carrier"` の中で切らない)。 */
export function splitShell(cmd) {
  const out = []; let cur = '', q = '';
  for (const ch of String(cmd)) {
    if (q) { cur += ch; if (ch === q) q = ''; continue; }
    if (ch === '"' || ch === "'") { q = ch; cur += ch; continue; }
    if (ch === ';' || ch === '|' || ch === '&') { out.push(cur); cur = ''; continue; }
    cur += ch;
  }
  out.push(cur);
  return out.filter(s => s.trim());
}
// ★★命令を `;` `&&` `||` `|` で切り、★**lean-idioms.md を含む区間だけ**を見る。
//   ★そうしないと `git add … tools/lean-idioms.md && git commit … | tail -1` の `tail` を
//   「末尾を見た」と数えてしまう(★実際に 1 度そう数えて、selftest に落とされた)。
export function lookKind(kind, pat) {
  const p = String(pat || '');
  if (kind === 'Read') return 'range';
  if (kind === 'Grep') return STRUCTPAT.test(p) ? 'numbering' : 'content';
  const segs = splitShell(p).filter(s => IDFILE.test(s));
  let best = 'other';
  const rank = { other: 0, range: 1, numbering: 2, content: 3 };
  for (const s of segs) {
    const gs = [...s.matchAll(/(?:grep|rg)\b[^\n"']*"([^"]*)"/g)].map(m => m[1])
      .concat([...s.matchAll(/(?:grep|rg)\b[^\n"']*'([^']*)'/g)].map(m => m[1]));
    let k = 'other';
    if (gs.length) k = gs.some(g => !STRUCTPAT.test(g)) ? 'content' : 'numbering';
    else if (/\b(tail|head|sed|cat|less|nl|wc)\b/.test(s)) k = 'range';
    if (rank[k] > rank[best]) best = k;
  }
  return best;
}
// 語が広すぎないか —— コーパス中の出現率で切る
export const GENERIC_RATE = 0.01;
export function isDistinctive(sig, corpusHits, corpusSize) {
  return corpusSize > 0 && corpusHits <= Math.max(1, Math.floor(corpusSize * GENERIC_RATE));
}

// ---------- ★★重複して書かれた idiom(M195 / 規則は M199 に事前登録) ----------
// ★問題: 同じ失敗形が別々の言い方で何度も節になる。★M193 が `expected 'lemma'` で見つけた。
// ★「重複」を字面の**完全一致**で測ると取り落とす —— `unexpected token 'set_option'; expected 'lemma'`
//   と `unexpected token 'omit'; expected 'lemma'` は別の文字列だが同じ罠である。
// ⇒ ★**字面の長さ GRAM_N の n-gram を共有するか**で繋ぐ。上の 2 つは `; expected 'lemma'` で繋がる。
// ★汎用句(`failed to synthesize` 等)で全部が繋がるのを防ぐため、★**既存の `isDistinctive` を通した
//   n-gram だけ**を辺に使う(新しい閾値を持ち込まない)。
export const GRAM_N = 16;
export function gramsOf(lits, n = GRAM_N) {
  const out = new Set();
  for (const lit of lits) { const s = String(lit); for (let i = 0; i + n <= s.length; i++) out.add(s.slice(i, i + n)); }
  return out;
}
/** n-gram を共有する節を連結成分にまとめる。keep(g) が偽の n-gram は辺にしない。 */
export function dupClusters(gramSets, keep = () => true) {
  const g2s = new Map();
  gramSets.forEach((gs, i) => { for (const g of gs) { let a = g2s.get(g); if (!a) g2s.set(g, a = new Set()); a.add(i); } });
  const shared = [...g2s].filter(([g, set]) => set.size >= 2 && keep(g)).map(([g, set]) => ({ g, set }));
  const par = gramSets.map((_, i) => i);
  const find = x => { while (par[x] !== x) { par[x] = par[par[x]]; x = par[x]; } return x; };
  for (const { set } of shared) { const a = [...set]; for (let k = 1; k < a.length; k++) { const r = find(a[0]), q = find(a[k]); if (r !== q) par[q] = r; } }
  const comp = new Map();
  gramSets.forEach((_, i) => { const r = find(i); if (!comp.has(r)) comp.set(r, []); comp.get(r).push(i); });
  const clusters = [...comp.values()].filter(c => c.length >= 2).sort((a, b) => b.length - a.length);
  return { clusters, shared };
}
// ★見出しの中の backtick の中身 —— 「人が節を探すときに打つであろう語」の代理。
//   ★これが先行節の本文に**literal で出るか**を見ることで、「grep で見つけられたか」を機械で測る。
export function headSpans(head) {
  const out = new Set();
  for (const m of String(head).matchAll(/`([^`\n]{4,80})`/g)) { const n = norm(m[1]); if (n.length >= 4) out.add(n); }
  return out;
}

// ---------- ログ走査 ----------

function walk(d, out = []) {
  let ents; try { ents = fs.readdirSync(d, { withFileTypes: true }); } catch { return out; }
  for (const e of ents) {
    const p = path.join(d, e.name);
    if (e.isDirectory()) walk(p, out); else if (e.name.endsWith('.jsonl')) out.push(p);
  }
  return out;
}

// ★★診断を **どの agent が出したか**に結び付ける(メタ第 27 回)。
//   子 agent の会話は `<session>/subagents/agent-XXXX.jsonl` に、その素性は
//   同じ名前の `.meta.json`(`agentType` / `toolUseId`)に在る。
//   ★`toolUseId` は `agent-timing.mjs` が完了通知(duration_ms)を引く鍵と**同じもの**なので、
//     これを digest に入れておくと「診断の族」と「所要時間」が突き合わせられる。
//   ★親セッション直下の `<session>.jsonl` には meta が無い ⇒ null(= 本体セッション)。
export function metaPathOf(logPath) {
  const p = String(logPath).replace(/\\/g, '/');
  const m = /\/subagents\/([^/]+)\.jsonl$/.exec(p);
  if (!m) return null;
  return p.slice(0, p.length - 6) + '.meta.json';
}
const metaCache = new Map();
function agentOf(logPath) {
  if (metaCache.has(logPath)) return metaCache.get(logPath);
  const mp = metaPathOf(logPath);
  let v = { id: '', type: '' };
  if (mp) { try { const j = JSON.parse(fs.readFileSync(mp, 'utf8')); v = { id: j.toolUseId || '', type: j.agentType || '' }; } catch { /* meta が無い子もある */ } }
  metaCache.set(logPath, v);
  return v;
}

// ★出所を追うため tool_use_id も返す(以前は本文だけ返していた)。
function resultTexts(o) {
  const out = []; const c = o.message?.content;
  if (Array.isArray(c)) {
    for (const x of c) if (x.type === 'tool_result') {
      const id = x.tool_use_id || '';
      if (typeof x.content === 'string') out.push({ id, text: x.content });
      else if (Array.isArray(x.content)) for (const y of x.content) if (y.type === 'text') out.push({ id, text: y.text || '' });
    }
  } else if (o.type === 'user' && typeof c === 'string') out.push({ id: '', text: c });
  return out;
}

// tool_use の入力から「何を検査したか」の 1 行を作る(provenance の材料)
function useCmd(input) {
  const i = input || {};
  return String(i.command ?? i.file_path ?? i.path ?? i.snippet ?? i.code ?? '').slice(0, 2000);
}

async function rescan() {
  const errors = [], writes = [];
  let nLines = 0, nFiles = 0, nMcp = 0, nLake = 0;
  for (const f of walk(LOGROOT)) {
    nFiles++;
    const ag = agentOf(f);                        // {id: toolUseId, type: agentType}
    const uses = new Map();                       // tool_use_id -> {name, cmd}
    const rl = readline.createInterface({ input: fs.createReadStream(f), crlfDelay: Infinity });
    for await (const line of rl) {
      nLines++;
      const hasUse = line.includes('"tool_use"');
      // ★行の絞りは `error` の語だけで見る。JSONL は本文を \n で逃がすので
      //   `/\berror \d+:\d+/` は **必ず外れる**(直前の文字が `n` なので語境界が立たない)。
      //   ★これで MCP の診断を 8,939 件中 454 件しか拾えていなかった。
      const hasErr = line.includes('error');
      const hasId = line.includes('lean-idioms');
      if (!hasUse && !hasErr && !hasId) continue;
      let o; try { o = JSON.parse(line); } catch { continue; }
      const ts = o.timestamp ? Math.floor(Date.parse(o.timestamp) / 1000) : 0;
      if (hasUse) { const c = o.message?.content; if (Array.isArray(c)) for (const x of c) if (x.type === 'tool_use') uses.set(x.id, { name: x.name, cmd: useCmd(x.input) }); }
      if (!ts) continue;
      if (hasId) {
        const c = o.message?.content;
        if (Array.isArray(c)) for (const x of c) {
          if (x.type !== 'tool_use') continue;
          const fp = String(x.input?.file_path || '');
          if (/lean-idioms\.md$/i.test(fp)) {
            const b = String(x.input?.new_string ?? x.input?.content ?? '');
            if (b) writes.push({ ts, how: x.name, text: b });
          } else if (x.name === 'Bash') {
            const cmd = String(x.input?.command || '');
            if (/>>?\s*["']?[^"'\s]*lean-idioms\.md/.test(cmd) || /(cat|tee|printf|echo)[^\n]{0,80}lean-idioms\.md/.test(cmd)) writes.push({ ts, how: 'Bash', text: cmd });
          }
        }
      }
      if (hasErr) for (const { id, text: t } of resultTexts(o)) {
        const u = uses.get(id) || { name: '', cmd: '' };
        // ★MCP REPL は書式が違う。tool の名前で分ける(本文の顔で当てない)。
        const bs = u.name === 'mcp__abc3-lean__lean_check' ? mcpErrorBlocks(t) : (t.includes('error:') ? errorBlocks(t) : []);
        if (u.name === 'mcp__abc3-lean__lean_check') nMcp += bs.length; else nLake += bs.length;
        for (const b of bs) errors.push({ ts, file: b.file, msg: b.msg, src: provenance(u.name, u.cmd, b.file, ag.type), ag: ag.id });
      }
    }
  }
  fs.mkdirSync(path.dirname(CACHE), { recursive: true });
  fs.writeFileSync(CACHE, JSON.stringify({ v: DIGEST_V, builtAt: Math.floor(Date.now() / 1000), nFiles, nLines, nMcp, nLake, errors, writes }));
  return { v: DIGEST_V, nFiles, nLines, nMcp, nLake, errors, writes };
}

// ★237 commit ぶんの `git show` で 16 秒かかる。★鍵は「その本を触った commit の並び」なので
//   新しい commit が来たときだけ作り直す(0.28 秒に戻る)。
const GCACHE = path.join(REPO, '.cache', 'idiom-recur-gitfirst.json');
function gitFirstSeen() {
  const g = a => execFileSync('git', a, { encoding: 'utf8', maxBuffer: 1 << 28 });
  const log = g(['log', '--reverse', '--format=%H %ct', '--', 'tools/lean-idioms.md']).trim().split('\n').filter(Boolean).map(l => l.split(' '));
  const key = log.length + ':' + (log[log.length - 1]?.[0] ?? '');
  if (fs.existsSync(GCACHE)) {
    try { const c = JSON.parse(fs.readFileSync(GCACHE, 'utf8')); if (c.key === key) return new Map(c.rows); } catch { /* 作り直す */ }
  }
  const first = new Map();
  for (const [sha, ts] of log) {
    let txt; try { txt = g(['show', sha + ':tools/lean-idioms.md']); } catch { continue; }
    for (const s of sections(txt)) if (!first.has(s.key)) first.set(s.key, +ts);
  }
  fs.mkdirSync(path.dirname(GCACHE), { recursive: true });
  fs.writeFileSync(GCACHE, JSON.stringify({ key, rows: [...first] }));
  return first;
}

// ---------- selftest ----------
function selftest() {
  const T = []; const eq = (name, a, b) => T.push({ name, ok: JSON.stringify(a) === JSON.stringify(b), a, b });
  eq('S1 heading key strips date+star', headingKey('## #97 `K.carrier` の罠（2026-09-07、Λ7）'), '#97 `K.carrier` の罠');
  eq('S2 sections ignore # inside fence', sections('## a\n```\n## not a heading\n```\n## b').map(s => s.key), ['a', 'b']);
  eq('S3 section line numbers', sections('x\n## a\ny\n## b').map(s => s.line), [2, 4]);
  eq('S4 errorBlocks picks only error', errorBlocks('F.lean:1:1: warning: w\nF.lean:2:2: error: boom\n  detail\nF.lean:3:3: info: i').map(b => b.msg), ['boom\n  detail']);
  eq('S5 errorBlocks bare error:', errorBlocks('error: plain').map(b => b.msg), ['plain']);
  eq('S6 errorBlocks file captured', errorBlocks('A/B.lean:9:1: error: z')[0].file, 'A/B.lean');
  eq('S7 errorBlocks none', errorBlocks('all fine\nok'), []);
  eq('S8 declNames keeps dotted long', [...declNames('`IntermediateField.mem_fixingSubgroup_iff F g`')], ['IntermediateField.mem_fixingSubgroup_iff']);
  eq('S9 declNames drops short', [...declNames('`Nat.card` `K.carrier`')], []);
  eq('S10 declNames drops lowercase ns', [...declNames('`toAffine.PointStuff`')], []);
  eq('S11 declNames drops undotted', [...declNames('`maxHeartbeatsLong`')], []);
  eq('S12 declNames drops file paths', [...declNames('`tools/lean-idioms.md`')], []);
  eq('S13 declNames needs backticks', [...declNames('IntermediateField.mem_fixingSubgroup_iff')], []);
  eq('S14 errLits keeps long error quote', [...errLits('`failed to synthesize instance IsIso (P.Base a)`')], ['failed to synthesize instance IsIso (P.Base a)']);
  eq('S15 errLits drops short', [...errLits('`unsolved goals`')], []);
  eq('S16 errLits reads fenced lines', [...errLits('```\nfailed to synthesize instance of type class Foo Bar\n```')], ['failed to synthesize instance of type class Foo Bar']);
  eq('S17 errLits needs an error word', [...errLits('`this sentence is long enough but has no error word`')], []);
  eq('S18 distinctive cut', [isDistinctive('x', 10, 4831), isDistinctive('x', 60, 4831)], [true, false]);
  eq('S19 norm collapses', norm('a \n  b\t c '), 'a b c');
  eq('S20 DIAG on windows path', DIAG.exec('C:/t/x.lean:1:2: error: e')[4], 'error');
  // ★★MCP REPL の書式(メタ第 26 回。これを取り落として 8,939 件 = 全体の 65% が見えていなかった)
  eq('S21 mcp picks only error', mcpErrorBlocks('エラー 2 件 (0.4 秒)\n\nerror 15:53\nunsolved goals\n  x\n\ninfo 1:0\nnope\n\nerror 3:1\nboom').map(b => b.msg), ['unsolved goals\n  x', 'boom']);
  eq('S22 mcp ignores `error:` form', mcpErrorBlocks('F.lean:1:1: error: not this shape'), []);
  eq('S23 DIAG ignores mcp form', errorBlocks('error 15:53\nunsolved goals'), []);
  eq('S24 mcp drops empty body', mcpErrorBlocks('error 1:1\n\nerror 2:2\nreal').map(b => b.msg), ['real']);
  // ★出所 —— 族から「Lean のエラーでないもの」を機械で外す根拠
  eq('S25 prov doc beats all', provenance('Bash', 'tail -30 tools/lean-idioms.md', 'ABC3/F.lean'), 'doc');
  eq('S26 prov repl is scratch', provenance('mcp__abc3-lean__lean_check', 'theorem foo : True := trivial', ''), 'scratch');
  eq('S27 prov tree by file', provenance('Bash', 'lake build', 'ABC3/Found/PGC/X.lean'), 'tree');
  eq('S28 prov tree by lean path', provenance('Bash', 'lake build', 'lean/ABC3/Found/X.lean'), 'tree');
  eq('S29 prov scratch by outside file', provenance('Bash', 'node tools/leanfile.mjs C:/tmp/t1.lean', 'C:/tmp/t1.lean'), 'scratch');
  eq('S30 prov no file falls back to cmd', [provenance('Bash', 'cd lean && lake build ABC3', ''), provenance('Bash', 'echo hi', '')], ['tree', 'scratch']);
  eq('S36 prov ABC3 must be a directory', provenance('Bash', '', 'myABC3lib/F.lean'), 'scratch');
  // ★出所 `meta` —— 改善係自身の probe(M120)。★ファイル名ではなく **agent の種類**で決める。
  eq('S37 prov meta wins over tree', provenance('Bash', 'lake build', 'ABC3/Found/PGC/X.lean', 'meta-optimizer'), 'meta');
  eq('S38 prov meta wins over doc', provenance('Bash', 'tail tools/lean-idioms.md', '', 'meta-optimizer'), 'meta');
  eq('S39 prov meta only for that type', provenance('Bash', 'lake build', 'ABC3/F.lean', 'lean-prover'), 'tree');
  eq('S40 prov no agentType is unchanged', [provenance('Bash', 'lake build', 'ABC3/F.lean', ''), provenance('Bash', 'lake build', 'ABC3/F.lean')], ['tree', 'tree']);
  // ★診断 → agent の突き合わせ(`agent-timing --diag` の鍵)
  eq('S41 metaPath from subagent log', metaPathOf('C:/u/.claude/projects/P/S/subagents/agent-a1.jsonl'), 'C:/u/.claude/projects/P/S/subagents/agent-a1.meta.json');
  eq('S42 metaPath backslashes too', metaPathOf('C:\\u\\P\\S\\subagents\\agent-a1.jsonl'), 'C:/u/P/S/subagents/agent-a1.meta.json');
  eq('S43 metaPath none for parent session', metaPathOf('C:/u/.claude/projects/P/S.jsonl'), null);
  eq('S44 metaPath none for other dirs', metaPathOf('C:/u/P/other/agent-a1.jsonl'), null);
  // ★下位の族 —— `Invalid field` は 1 つの罠ではない(メタ第 26 回で 3 つに割れた)
  eq('S31 sub notation', subKey('Invalid field notation: Field projection operates on types of the form `C ...`', 'Invalid field'), 'notation');
  eq('S32 sub env namespace', subKey('Invalid field `mp`: The environment does not contain `Function.mp`, so it is', 'Invalid field'), 'env:Function');
  eq('S33 sub env our own tree', subKey('Invalid field `filt`: The environment does not contain `ABC3.Interface.PGC.RamificationFiltration.filt`', 'Invalid field'), 'env:ABC3.Interface.PGC.RamificationFiltration');
  eq('S34 sub falls back to 3 words', subKey('has already been declared here and there ok', 'has already been declared'), 'here and there');
  eq('S35 sub absent family', subKey('nothing', 'Invalid field'), '');
  // ★★ERRWORDS v2 —— Lean 4 の書式変更(M189 / 事前登録 M192)。★v1 を凍結したことも試す。
  const q4 = '```\nTactic `rewrite` failed: Did not find an occurrence of the pattern\n```';
  eq('S45 v2 catches new tactic format', [...errLits(q4, ERRWORDS_V2)], ['Tactic `rewrite` failed: Did not find an occurrence of the pattern']);
  eq('S46 v1 is frozen and misses it', [...errLits(q4, ERRWORDS_V1)], []);
  eq('S47 v2 catches the other five', [
    errLits('`Lean exited with code 1 Some required targets failed`', ERRWORDS_V2).size,
    errLits('```\n`exact?` could not close the goal. Try `apply?` first\n```', ERRWORDS_V2).size,
    errLits('`failed to compile definition, consider marking it noncomputable`', ERRWORDS_V2).size,
    errLits("`unexpected token 'set_option'; expected 'lemma'`", ERRWORDS_V2).size,
    errLits('`unknown namespace AlgebraicGeometry`', ERRWORDS_V2).size,
  ], [1, 1, 1, 1, 1]);
  eq('S48 v1 misses the other five', [
    errLits('`Lean exited with code 1 Some required targets failed`', ERRWORDS_V1).size,
    errLits('```\n`exact?` could not close the goal. Try `apply?` first\n```', ERRWORDS_V1).size,
    errLits('`failed to compile definition, consider marking it noncomputable`', ERRWORDS_V1).size,
    errLits("`unexpected token 'set_option'; expected 'lemma'`", ERRWORDS_V1).size,
    errLits('`unknown namespace AlgebraicGeometry`', ERRWORDS_V1).size,
  ], [0, 0, 0, 0, 0]);
  // ★`tactic … failed` が貪欲でないこと(M189 の提案 `Tactic .* failed` を絞った理由)
  eq('S49 tactic word does not span', ERRWORDS_V2.test('tactic here is a very long sentence with many words that eventually failed'), false);
  eq('S50 tactic word still reaches backticked name', ERRWORDS_V2.test('Tactic `simp` failed'), true);
  // ★「エラーらしい引用」の規則(M189 の 23 節を再現可能にするため。★M189 の規則とは別物)
  eq('S51 errLike finds what errLits misses', [
    [...errLits('`failed to read file X.olean: incompatible header`', ERRWORDS_V1)].length,
    [...errLikeQuotes('`failed to read file X.olean: incompatible header`')].length,
  ], [0, 1]);
  eq('S52 errLike ignores plain prose', [...errLikeQuotes('`this sentence is long enough but says nothing bad`')], []);
  // ★★重複クラスタ(M195 / 規則は M199 に事前登録)
  eq('S53 gramsOf window count', gramsOf(['abcdefgh'], 4).size, 5);
  eq('S54 gramsOf drops too short', gramsOf(['abc'], 4).size, 0);
  eq('S55 gramsOf unions literals', [...gramsOf(['abcd', 'abcd', 'zbcd'], 4)].sort(), ['abcd', 'zbcd']);
  // ★要: 字面が**違っても**共有 n-gram で繋がること(`'set_option'` と `'omit'` の実例の縮小版)
  eq('S56 dup links different literals via shared gram', dupClusters([
    gramsOf(["unexpected token 'set_option'; expected 'lemma'"], 16),
    gramsOf(["unexpected token 'omit'; expected 'lemma'"], 16),
    gramsOf(['something else entirely and long'], 16),
  ]).clusters, [[0, 1]]);
  // ★共有する 16-gram が **ちょうど 1 つ**になるように組んである(接尾辞 16 字だけが共通)
  eq('S57a the pair shares exactly one gram', dupClusters([gramsOf(['aaaaGENERICGENERICGE'], 16), gramsOf(['bbbbGENERICGENERICGE'], 16)]).shared.map(x => x.g), ['GENERICGENERICGE']);
  eq('S57b dup keep() drops that gram', dupClusters([
    gramsOf(['aaaaGENERICGENERICGE'], 16), gramsOf(['bbbbGENERICGENERICGE'], 16),
  ], g => g !== 'GENERICGENERICGE').clusters, []);
  eq('S58 dup is transitive', dupClusters([
    gramsOf(['AAAAAAAAAAAAAAAAxxxx'], 16), gramsOf(['AAAAAAAAAAAAAAAAxxxxBBBBBBBBBBBBBBBB'], 16), gramsOf(['BBBBBBBBBBBBBBBByyyy'], 16),
  ]).clusters.map(c => c.slice().sort()), [[0, 1, 2]]);
  eq('S59 dup ignores singletons', dupClusters([gramsOf(['AAAAAAAAAAAAAAAAAAAA'], 16), gramsOf(['BBBBBBBBBBBBBBBBBBBB'], 16)]).clusters, []);
  eq('S60 dup empty sets are not a cluster', dupClusters([new Set(), new Set(), new Set()]).clusters, []);
  // ★見出しの語(「人が打つであろう語」)—— 省略記号の綴りが違えば literal では当たらない
  eq('S61 headSpans picks backticked', [...headSpans('## `omit [Inst] in` は docstring の**前**に置く')], ['omit [Inst] in']);
  eq('S62 headSpans drops short', [...headSpans('## `in` と `x`')], []);
  // ★「読んだ」ではなく「何を探したか」で仮説を分ける(M195)
  eq('S64 lookKind numbering vs content', [
    lookKind('Bash', 'grep -nE "^## +#?25[012]" tools/lean-idioms.md; tail -20 tools/lean-idioms.md'),
    lookKind('Bash', 'grep -n "motive\\|carrier" tools/lean-idioms.md | head -20'),
    lookKind('Bash', 'tail -30 tools/lean-idioms.md'),
    lookKind('Bash', 'git add tools/lean-idioms.md && git commit -m x | tail -1'),
    lookKind('Read', 'offset=4736 limit=10'),
    lookKind('Grep', '^## #[0-9]+'),
    lookKind('Grep', 'docstring'),
  ], ['numbering', 'content', 'range', 'other', 'range', 'numbering', 'content']);
  eq('S65 lookKind single quotes too', lookKind('Bash', "grep -n 'set_option' tools/lean-idioms.md"), 'content');
  eq('S66 splitShell keeps quoted pipes', splitShell('grep -n "a\\|b" f.md | head -20'), ['grep -n "a\\|b" f.md ', ' head -20']);
  // ★D9(突然変異)が素通りした穴 —— **単引用符**の中の `|` で切ってはいけない
  eq('S67 splitShell keeps single-quoted pipes', splitShell("grep -n 'a\\|b' f.md | head"), ["grep -n 'a\\|b' f.md ", ' head']);
  eq('S68 lookKind single-quoted pattern with pipe', lookKind('Bash', "grep -n 'motive\\|carrier' tools/lean-idioms.md | head -20"), 'content');
  // ★D2(突然変異)が素通りした穴 —— **既定の GRAM_N** を通る道が 1 本も無かった
  // ★★最初の書き方(`gramsOf(['x'.repeat(GRAM_N)])`)は **GRAM_N を両辺に使っていた**ので
  //   定数を 4 に書き換えても素通りした(突然変異 D2)。★幅そのものを外から固定する。
  eq('S69 GRAM_N is the value M199 registered', GRAM_N, 16);
  const shortShare = [['AAAAAAAAAAAAAAAA' + 'SHARED12CHRS'], ['BBBBBBBBBBBBBBBB' + 'SHARED12CHRS']];
  eq('S70 default width ignores a 12-char overlap', dupClusters(shortShare.map(x => gramsOf(x))).clusters, []);
  eq('S71 width 4 would join it (so the width is what decides)', dupClusters(shortShare.map(x => gramsOf(x, 4))).clusters, [[0, 1]]);
  eq('S63 ellipsis spelling breaks literal grep', [
    '## `set_option ... in` を置けない'.includes('set_option … in'),
    '## `set_option … in` は前'.includes('set_option ... in'),
  ], [false, false]);
  const bad = T.filter(t => !t.ok);
  for (const t of bad) console.log('  NG', t.name, '\n     got', JSON.stringify(t.a), '\n     want', JSON.stringify(t.b));
  console.log(`selftest ${T.length - bad.length}/${T.length}`);
  return bad.length === 0;
}

// ---------- main ----------
if (has('--selftest')) { process.exit(selftest() ? 0 : 1); }

let dig;
if (has('--rescan')) { const t0 = Date.now(); dig = await rescan(); console.log(`rescan ${((Date.now() - t0) / 1000).toFixed(1)}s  files ${dig.nFiles} lines ${dig.nLines} errors ${dig.errors.length} (lake ${dig.nLake} / repl ${dig.nMcp}) idiom-writes ${dig.writes.length}`); }
else if (fs.existsSync(CACHE)) dig = JSON.parse(fs.readFileSync(CACHE, 'utf8'));
else { console.error('digest がない。まず `node tools/idiom-recur.mjs --rescan`'); process.exit(2); }
if ((dig.v || 1) !== DIGEST_V) { console.error(`digest が古い(v${dig.v || 1}、要 v${DIGEST_V})。\`node tools/idiom-recur.mjs --rescan\` で作り直すこと。`); process.exit(2); }

const md = fs.readFileSync(IDIOMS, 'utf8');
const secs = sections(md);
const writes = dig.writes.slice().sort((a, b) => a.ts - b.ts);
const gitFirst = has('--nogit') ? new Map() : gitFirstSeen();

for (const s of secs) {
  s.regTs = null; s.regSrc = 'none';
  for (const w of writes) if (w.text.includes(s.head)) { s.regTs = w.ts; s.regSrc = 'log:' + w.how; break; }
  if (s.regTs == null && gitFirst.has(s.key)) { s.regTs = gitFirst.get(s.key); s.regSrc = 'git'; }
  const txt = s.head + '\n' + s.body.join('\n');
  s.names = declNames(txt); s.lits = errLits(txt);
}

// ★★出所 `doc`(= 我々の文書を印字しただけ)は **Lean のエラーではない**ので既定で外す。
//   ★メタ第 26 回の実測: doc は 4,835 件中 23 件(0.5%)。★小さいが、
//   本体が「idiom はエラー文を逐語で引用する」と決めた以上、★これは**必ず増える**。
//   `--with-doc` で入れられる(前後を並べたいとき用)。
// ★★出所 `meta`(= 改善係自身の probe の出力)も同じ理由で既定で外す(M120)。
//   ★`--with-meta` で入れられる。★これは**ファイル名ではなく agent の種類**で決まるので、
//   probe を何と名付けても漏れない。
const ALLERRS = dig.errors.map(e => ({ ts: e.ts, file: e.file, src: e.src || 'tree', ag: e.ag || '', n: norm(e.msg) })).sort((a, b) => a.ts - b.ts);
const errs = ALLERRS.filter(e => (has('--with-doc') || e.src !== 'doc') && (has('--with-meta') || e.src !== 'meta'));
const N = errs.length;
const bySrcErr = {}; for (const e of ALLERRS) bySrcErr[e.src] = (bySrcErr[e.src] || 0) + 1;
const corpus = new Map();
const hitsOf = sig => { if (!corpus.has(sig)) corpus.set(sig, errs.filter(e => e.n.includes(sig)).length); return corpus.get(sig); };

// --sig lit  : idiom が引用しているエラー文そのものだけを字面とする(狭い・これが「同じ失敗形」に一番近い)
// --sig name : idiom が名指しした宣言名だけ(広い・上界)
// --sig both : 両方
// ★★既定は `lit`(メタ第 26 回で変更。それまでは `both`)。
//   理由: M109 が「(b) 補題名は上界としてすら意味が薄い。★報告には (a) を使う」と結論したのに
//   道具の既定が `both` のままで、★叩いた人が黙って `name` 込みの数字(87 件 / 625 事象)を
//   受け取っていた。★本体も「`--sig name` は使わない」と判断済み。
const SIG = ARGV.includes('--sig') ? ARGV[ARGV.indexOf('--sig') + 1] : 'lit';
let dropped = new Set();
for (const s of secs) {
  s.sigs = [];
  const pool = SIG === 'lit' ? [...s.lits] : SIG === 'name' ? [...s.names] : [...s.names, ...s.lits];
  for (const g of pool) { if (isDistinctive(g, hitsOf(g), N)) s.sigs.push(g); else dropped.add(g); }
  s.ev = [];
  if (s.regTs) for (const e of errs) { if (e.ts <= s.regTs) continue; for (const g of s.sigs) if (e.n.includes(g)) { s.ev.push({ ts: e.ts, sig: g, file: e.file, msg: e.n.slice(0, 120) }); break; } }
}

const noReg = secs.filter(s => !s.regTs).length;
const noSig = secs.filter(s => s.sigs.length === 0).length;
const measurable = secs.filter(s => s.regTs && s.sigs.length);
const withEv = secs.filter(s => s.ev.length);
const bySrc = {}; for (const s of secs) bySrc[s.regSrc] = (bySrc[s.regSrc] || 0) + 1;

console.log('== idiom-recur ==');
console.log(`  コーパス      : error 診断 ${N} 件 / idiom への書き込み ${dig.writes.length} 回  (--sig ${SIG})`);
console.log(`  出所          : ${JSON.stringify(bySrcErr)}  ${has('--with-doc') ? '(doc も入れている)' : '★doc は外している'}`);
console.log(`  idiom(節)     : ${secs.length}`);
console.log(`  登録時刻の出所: ${JSON.stringify(bySrc)}`);
console.log(`  ★照合できない : ${noReg + noSig - secs.filter(s => !s.regTs && !s.sigs.length).length} 件  (登録時刻なし ${noReg} / 字面なし ${noSig})`);
console.log(`  ★測れる       : ${measurable.length} 件`);
console.log(`  ★登録後の再出現: idiom ${withEv.length} 件 / 事象 ${secs.reduce((a, s) => a + s.ev.length, 0)} 件`);
console.log(`  広すぎて捨てた字面: ${dropped.size} 個 (コーパスの ${(GENERIC_RATE * 100).toFixed(0)}% 超)`);

// ★★`--lits-audit` —— 「エラーらしい引用はあるのに字面が 0 の節」を**名指しで**数える。
//   ★M189(メタ第 36 回)が「23 節」と書いたが規則を残さなかったので、★規則ごと道具に入れる(M188)。
//   ★規則は `ERRLIKE`(上)。★M189 の規則とは別物なので、23 と一致する保証は無い。
//   ★同時に v1 → v2 で**新たに字面が取れるようになった節**を全部出す(これが効果の主指標)。
if (has('--lits-audit')) {
  const rows = secs.map(s => {
    const txt = s.head + '\n' + s.body.join('\n');
    return { s, v1: errLits(txt, ERRWORDS_V1), v2: errLits(txt, ERRWORDS_V2), like: errLikeQuotes(txt) };
  });
  const cnt = f => rows.filter(f).length;
  console.log('\n-- lits-audit(規則: ERRLIKE / LIT_MIN = ' + LIT_MIN + ') --');
  console.log(`  節 ${rows.length}`);
  console.log(`  字面が取れた節          : v1 ${cnt(r => r.v1.size)} → v2 ${cnt(r => r.v2.size)}`);
  console.log(`  ★エラーらしい引用はあるのに字面 0 : v1 ${cnt(r => !r.v1.size && r.like.size)} → v2 ${cnt(r => !r.v2.size && r.like.size)}`);
  console.log(`  ★v1 で 0 → v2 で非 0(効果の主指標): ${cnt(r => !r.v1.size && r.v2.size)} 節`);
  console.log('\n-- v1 で 0 → v2 で非 0 の節(全部) --');
  for (const r of rows.filter(r => !r.v1.size && r.v2.size)) {
    console.log(`  L${String(r.s.line).padStart(5)} ${r.s.head.slice(0, 76)}`);
    for (const g of [...r.v2].slice(0, 2)) console.log(`         + ${JSON.stringify(g.slice(0, 100))}`);
  }
  console.log('\n-- v2 でもまだ字面 0 だが「エラーらしい引用」を持つ節(全部) --');
  for (const r of rows.filter(r => !r.v2.size && r.like.size)) {
    console.log(`  L${String(r.s.line).padStart(5)} ${r.s.head.slice(0, 76)}`);
    for (const g of [...r.like].slice(0, 1)) console.log(`         ? ${JSON.stringify(g.slice(0, 100))}`);
  }
}
if (has('--hits')) {
  console.log('\n-- 再出現があった idiom(名指し) --');
  for (const s of withEv.sort((a, b) => b.ev.length - a.ev.length)) {
    console.log(`  [${String(s.ev.length).padStart(3)}] L${s.line} (${s.regSrc} ${new Date(s.regTs * 1000).toISOString().slice(0, 16)}) ${s.head.slice(0, 88)}`);
    for (const e of s.ev.slice(0, 3)) console.log(`        ${new Date(e.ts * 1000).toISOString().slice(0, 16)} <${e.sig}> ${JSON.stringify(e.msg.slice(0, 90))}`);
  }
}
if (has('--family')) {
  // 字面が本当に一致する粒度は「補題名」ではなく「Lean の診断文そのもの」である。
  // 診断文の先頭句を族とし、その族を引用した最初の idiom の後で族が何回出たかを数える。
  console.log('\n-- 失敗形の族(= Lean の診断文の先頭句)ごとの再出現 --');
  console.log('   族 / 総数 / 族を引用した最初の idiom / その後の出現');
  const fams = [
    'Invalid projection: Projections cannot be used on functions',
    'Invalid field', 'has already been declared', 'motive is not type correct',
    'failed to synthesize', 'Application type mismatch', 'Type mismatch',
    'Did not find an occurrence of the pattern', 'unsolved goals',
    '(deterministic) timeout at `whnf`', '(kernel) deterministic timeout',
    'maximum recursion depth', 'unknown identifier', 'unknown constant',
    'made no progress', 'synthesized type class instance is not definitionally equal',
    'typeclass instance problem is stuck', 'The rfl tactic', 'linarith failed',
  ];
  const secTxt = secs.map(s => ({ s, t: s.head + '\n' + s.body.join('\n') }));
  const rows = [];
  for (const f of fams) {
    const all = ALLERRS.filter(e => e.n.includes(f)).sort((a, b) => a.ts - b.ts);
    if (!all.length) continue;
    const owners = secTxt.filter(x => x.t.includes(f) && x.s.regTs).map(x => x.s).sort((a, b) => a.regTs - b.regTs);
    const first = owners[0];
    const cnt = k => all.filter(e => e.src === k).length;
    // ★除外後 = doc(文書を印字しただけ)と scratch(捨てられる試作)を外し、
    //   ★**木の本を検査して出た失敗だけ**にしたもの。
    const kept = all.filter(e => e.src === 'tree');
    const after = first ? all.filter(e => e.ts > first.regTs).length : null;
    const afterK = first ? kept.filter(e => e.ts > first.regTs).length : null;
    rows.push({ f, all: all.length, doc: cnt('doc'), meta: cnt('meta'), scratch: cnt('scratch'), tree: kept.length, first, after, afterK });
  }
  const P = (n, w) => String(n).padStart(w);
  console.log('   総数 = doc + meta + scratch + tree。★doc は Lean のエラーではない(文書を印字しただけ)。');
  console.log('   ★meta = 改善係(meta-optimizer)自身の probe の出力。★Lean は動いていない(M120)。');
  console.log('   ★scratch = 検査したのが木の本でない(REPL の断片 / `lean/ABC3/` 外の .lean)= 捨てられる試作。');
  console.log('   族                                                     総数   doc  meta scratch  tree | 最初の idiom            → 後(総) → 後(tree)');
  for (const r of rows.sort((a, b) => b.all - a.all)) {
    console.log(`  ${r.f.slice(0, 52).padEnd(54)}${P(r.all, 5)}${P(r.doc, 6)}${P(r.meta, 6)}${P(r.scratch, 8)}${P(r.tree, 6)} | ` +
      (r.first ? `L${String(r.first.line).padEnd(5)} ${new Date(r.first.regTs * 1000).toISOString().slice(0, 16)} ${P(r.after, 6)} ${P(r.afterK, 9)}` : '(idiom に引用なし)'));
  }
  console.log('  ★「後 N」は「idiom が効かなかった回数」ではない。字面が後にも出た回数でしかない。');
  console.log('  ★★scratch を外すのは「意図した合図だから」ではない。★言えるのは「木は 1 行も損なわれていない」まで。');
  console.log('    ★#158(同じ名前で書いて `already been declared` を出させる)の合図は**必ず** scratch に落ちる。');
  console.log('    ★逆は言えない —— scratch にも本物の失敗はある。');
  // ★下位の族 —— 汎用エラーは 1 つの族に別々の原因が混ざる。次の句で機械的に割る。
  if (has('--sub')) {
    const want = ARGV[ARGV.indexOf('--sub') + 1];
    for (const f of fams) {
      if (want && !f.includes(want)) continue;
      const all = ALLERRS.filter(e => e.n.includes(f));
      if (all.length < 20) continue;
      const m = new Map();
      for (const e of all) { const k = subKey(e.n, f); m.set(k, (m.get(k) || 0) + 1); }
      const ent = [...m.entries()].sort((a, b) => b[1] - a[1]);
      console.log(`\n  -- 下位の族: ${f} (総数 ${all.length} / 下位 ${ent.length} 通り) --`);
      for (const [k, v] of ent.slice(0, 14)) console.log(`     ${P(v, 5)}  ${k}`);
      if (ent.length > 14) console.log(`     … 残り ${ent.length - 14} 通り`);
    }
  }
}
if (has('--calibrate')) {
  console.log('\n-- 校正: 本体が名指しした既知の 1 件 --');
  const key = 'mem_fixingSubgroup_iff';
  const owners = secs.filter(s => (s.head + s.body.join('\n')).includes(key));
  console.log(`  この字面を書いた idiom: ${owners.length} 件 (L${owners.map(s => s.line).join(', L')})`);
  const es = errs.filter(e => e.n.includes(key));
  console.log(`  この字面を含む error  : ${es.length} 件 ${es.map(e => new Date(e.ts * 1000).toISOString().slice(0, 16)).join(' ')}`);
  const t = Math.min(...owners.filter(s => s.regTs).map(s => s.regTs));
  console.log(`  最も早い登録時刻      : ${isFinite(t) ? new Date(t * 1000).toISOString().slice(0, 16) : '不明'}`);
  console.log(`  ⇒ 登録後の再出現(字面): ${es.filter(e => e.ts > t).length} 件`);
  console.log('  ★★メタ第 25 回(M109)はここを **0 件**と読み、「Lean は補題名ではなく現場の項を印字するから');
  console.log('    構造的な偽陰性がある」と結論した。★★メタ第 26 回でその結論は**覆った**。');
  console.log('    ★原因は「補題名で照合する設計」ではなく、★**MCP REPL の診断を 1 件も読めていなかった**こと。');
  console.log('    ★本体が「2026-09-07 に踏んだ」と言った事象は 2026-09-07T06:56 に**そのまま在る**:');
  console.log('      "Invalid projection: Projections cannot be used on functions, and');
  console.log('       IntermediateField.mem_fixingSubgroup_iff ?m.216 has function type"');
  console.log('    ⇒ ★補題名は**含まれていた**。読めていなかっただけである。');
  console.log('  ★残る偽陰性は別にある —— 上の 6 件のうち 4 件は `scratch`(REPL)で、');
  console.log('    ★書式が変われば同じことが起きる。★**書式は tool の名前で分けること**(本文の顔で当てない)。');

  // ★★校正その 2(メタ第 26 回) —— 「意図した合図」を数え切れているか。
  //   idiom #158 / #165 は「同じ名前で書いて `already been declared` を出させた」結果
  //   **8 本**が在庫だったと名指ししている。★その 8 本をログで引けるかを見る。
  console.log('\n-- 校正 2: #158 / #165 が「この手で見つけた」と名指しした 8 本 --');
  const found158 = ['fixingSubgroup_sup', 'isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot',
    'quotientMulSemiringActionOfTrivial', 'quotientSMul_mk', 'quotientSMul_mk_fixedRing',
    'ramIndex_quotient_mk', 'mem_ramificationGroupReal_quotient_mk_iff', 'herbrandPhiGroup_quotient_eq'];
  const decl = ALLERRS.filter(e => e.n.includes('has already been declared'));
  let seen = 0;
  for (const n of found158) {
    const es = decl.filter(e => e.n.includes('.' + n + '`') || e.n.includes('`' + n + '`'));
    if (es.length) seen++;
    console.log(`  ${n.padEnd(56)} ${String(es.length).padStart(2)} 件  ${[...new Set(es.map(e => e.src))].join(',')}`);
  }
  console.log(`  ⇒ ★ログで引けるのは ${seen}/8。★残り ${8 - seen} 本は `);
  console.log('    `already been declared` を **1 度も出していない**(別の手で見つけたか、合図が別の顔で出た)。');
  console.log('  ★★だから「scratch の `already been declared` = 意図した合図」と**逆向きには使えない**。');
  console.log('    ★言えるのは「合図は必ず scratch に落ちる」までで、★合図の総数は**下から**しか押さえられない。');
}
// ★★`--dupes` —— 「同じ失敗形が別々の節に何度書かれたか」と「書いた後も出続けているか」。
//   ★規則は M199 に**測る前に**事前登録した。★`--gram N` で感度を見られる(既定 16)。
if (has('--dupes') || has('--retro')) {
  const GN = ARGV.includes('--gram') ? Number(ARGV[ARGV.indexOf('--gram') + 1]) : GRAM_N;
  const build = n => {
    const sets = secs.map(s => gramsOf(s.lits, n));
    const r = dupClusters(sets, g => isDistinctive(g, hitsOf(g), N));
    return { sets, ...r };
  };
  const { sets, clusters, shared } = build(GN);
  const inClu = new Set(); for (const c of clusters) for (const i of c) inClu.add(i);
  // ★★M199 に「代表 = 繋いだ n-gram のうち**最長**」と書いたが、★共有 n-gram はどれも長さ GN で
  //   **同着**である。★事前登録した規則が ill-defined だった。★勝手に選び直さず、両方出す:
  //   (a) 事前登録の規則 + 決定的な同着処理(節数が最多、同数なら辞書順の最初)
  //   (b) 直した規則(クラスタの共有 n-gram の**どれか**に当たる診断を数える)
  const repOf = c => {
    const own = new Set(c); let best = null;
    for (const { g, set } of shared) { if (![...set].every(i => own.has(i))) continue; if (!best || set.size > best.n || (set.size === best.n && g < best.g)) best = { g, n: set.size }; }
    return best ? best.g : '';
  };
  const rows = clusters.map((c, ci) => {
    const ss = c.map(i => secs[i]).sort((a, b) => (a.regTs || 0) - (b.regTs || 0));
    const own = new Set(c);
    const gs = shared.filter(x => [...x.set].every(i => own.has(i))).sort((a, b) => b.set.size - a.set.size || (a.g < b.g ? -1 : 1));
    const t0 = Math.min(...ss.filter(s => s.regTs).map(s => s.regTs));
    return { ci, c, ss, gs, rep: repOf(c), t0, totalAny: 0, afterAny: 0, totalRep: 0, afterRep: 0 };
  });
  // ★診断を 1 度だけ舐めて、どのクラスタの共有 n-gram に当たるかを数える(gram ごとに 14k 件を舐めない)
  const g2c = new Map(); rows.forEach(r => { for (const { g } of r.gs) if (!g2c.has(g)) g2c.set(g, r); });
  for (const e of errs) {
    const touched = new Set();
    for (let i = 0; i + GN <= e.n.length; i++) { const r = g2c.get(e.n.slice(i, i + GN)); if (r) touched.add(r); }
    for (const r of touched) { r.totalAny++; if (isFinite(r.t0) && e.ts > r.t0) r.afterAny++; }
    for (const r of rows) { if (r.rep && e.n.includes(r.rep)) { r.totalRep++; if (isFinite(r.t0) && e.ts > r.t0) r.afterRep++; } }
  }
  console.log(`\n-- dupes(規則 M199: 字面の ${GN}-gram を共有 + isDistinctive / --errwords ${ERRWORDS_VER}) --`);
  console.log(`  字面を持つ節            : ${secs.filter(s => s.lits.size).length} / ${secs.length}`);
  console.log(`  ★重複クラスタ           : ${clusters.length} 個`);
  console.log(`  ★クラスタに入る節       : ${inClu.size} 節`);
  console.log(`  ★最大クラスタ           : ${clusters[0] ? clusters[0].length : 0} 節`);
  console.log(`  ★書いた後も再発(代表 gram / 事前登録の規則): ${rows.filter(r => r.afterRep > 0).length} クラスタ / 事象 ${rows.reduce((a, r) => a + r.afterRep, 0)} 件`);
  console.log(`  ★書いた後も再発(共有 gram のどれか / 直した規則): ${rows.filter(r => r.afterAny > 0).length} クラスタ / 事象 ${rows.reduce((a, r) => a + r.afterAny, 0)} 件`);
  console.log('\n  -- クラスタ(全部。★標本ではない) --');
  for (const r of rows) {
    console.log(`  [${String(r.ss.length).padStart(2)} 節] 代表 ${JSON.stringify(r.rep)} (共有 gram ${r.gs.length} 個)  総数 rep ${r.totalRep} / any ${r.totalAny}  ―  T0 後 rep ${r.afterRep} / any ${r.afterAny}`);
    // ★「後から書いた節の見出しの語で、先行節を grep できたか」を機械で見る(M195 の仮説 1 の検査)
    for (let k = 0; k < r.ss.length; k++) {
      const s = r.ss[k];
      const prior = r.ss.slice(0, k).map(p => p.head + '\n' + p.body.join('\n')).join('\n');
      const sp = [...headSpans(s.head)];
      const hit = k === 0 ? '' : (sp.some(x => prior.includes(x)) ? '見つかる' : '★見つからない');
      console.log(`      L${String(s.line).padStart(5)} ${s.regTs ? new Date(s.regTs * 1000).toISOString().slice(0, 10) : '   ?      '} ${hit.padEnd(12)} ${s.head.slice(0, 62)}`);
      if (k > 0 && sp.length) console.log(`               見出しの語: ${sp.slice(0, 4).map(x => JSON.stringify(x)).join(' ')}`);
    }
  }
  const later = rows.flatMap(r => r.ss.slice(1).map((s, k) => ({ s, prior: r.ss.slice(0, k + 1) })));
  const grepOk = later.filter(x => { const prior = x.prior.map(p => p.head + '\n' + p.body.join('\n')).join('\n'); return [...headSpans(x.s.head)].some(g => prior.includes(g)); }).length;
  console.log(`\n  ★後から書かれた節 ${later.length} 件のうち、見出しの語で先行節を grep できたのは ${grepOk} 件 (${later.length ? (grepOk / later.length * 100).toFixed(0) : 0}%)`);
  if (has('--gram-sweep')) {
    console.log('\n  -- 感度(★選び直しではなく全部見せる) --');
    console.log('     N   クラスタ   節   最大');
    for (const n of [12, 16, 24, 32]) { const b = build(n); const s2 = new Set(); for (const c of b.clusters) for (const i of c) s2.add(i); console.log(`   ${String(n).padStart(3)}${String(b.clusters.length).padStart(9)}${String(s2.size).padStart(6)}${String(b.clusters[0] ? b.clusters[0].length : 0).padStart(6)}`); }
  }
  // ★★`--retro` —— 「節を書く**前**に『似た節がある』と言う口」が、実際に先行節を出せたかを
  //   **その時点の `lean-idioms.md`(git の版)**に対して確かめる。★口の価値の主指標。
  if (has('--retro')) {
    const g = a => { try { return execFileSync('git', a, { encoding: 'utf8', maxBuffer: 1 << 28 }); } catch { return ''; } };
    console.log('\n  -- retro: その時点の lean-idioms.md に対して口を当てたら先行節を出せたか --');
    let ok = 0, no = 0, nohist = 0, self = 0;
    for (const r of rows) for (let k = 1; k < r.ss.length; k++) {
      const s = r.ss[k]; if (!s.regTs) { nohist++; continue; }
      const sha = g(['log', '--before=' + new Date((s.regTs - 1) * 1000).toISOString(), '-1', '--format=%H', '--', 'tools/lean-idioms.md']).trim();
      if (!sha) { nohist++; continue; }
      const md0 = g(['show', sha + ':tools/lean-idioms.md']); if (!md0) { nohist++; continue; }
      const hs = sections(md0).map(x => ({ x, lits: errLits(x.head + '\n' + x.body.join('\n')) })).filter(x => headingKey(x.x.head) !== headingKey(s.head));
      const mine = gramsOf(s.lits, GN);
      const found = hs.filter(x => { for (const q of gramsOf(x.lits, GN)) if (mine.has(q) && isDistinctive(q, hitsOf(q), N)) return true; return false; });
      if (sections(md0).some(x => headingKey(x.head) === headingKey(s.head))) self++;
      if (found.length) ok++; else no++;
      console.log(`      L${String(s.line).padStart(5)} ${sha.slice(0, 8)} ${found.length ? '★出せた ' + found.length + ' 節: L' + found.slice(0, 3).map(x => x.x.line).join(', L') : '出せない'}   ${s.head.slice(0, 46)}`);
    }
    console.log(`  ⇒ ★出せた ${ok} / 出せない ${no} / 版が取れない ${nohist}   (★書いた版に自分が既に居た: ${self} 件)`);
  }
}
// ★★`--reads` —— M195 の 2 仮説を分ける。「重複した節を書く前に、同じセッションで
if (has('--reads')) {
  const GN = ARGV.includes('--gram') ? Number(ARGV[ARGV.indexOf('--gram') + 1]) : GRAM_N;
  const sets = secs.map(s => gramsOf(s.lits, GN));
  const { clusters } = dupClusters(sets, g => isDistinctive(g, hitsOf(g), N));
  const target = new Map(); for (const c of clusters) for (const i of c) target.set(secs[i].head, secs[i]);
  const found = new Map();
  for (const f of walk(LOGROOT)) {
    const ev = [];
    const rl = readline.createInterface({ input: fs.createReadStream(f), crlfDelay: Infinity });
    for await (const line of rl) {
      if (!line.includes('lean-idioms')) continue;
      let o; try { o = JSON.parse(line); } catch { continue; }
      const ts = o.timestamp ? Math.floor(Date.parse(o.timestamp) / 1000) : 0;
      const c = o.message?.content; if (!ts || !Array.isArray(c)) continue;
      for (const x of c) {
        if (x.type !== 'tool_use') continue;
        const inp = x.input || {}; const fp = String(inp.file_path || inp.path || ''); const cmd = String(inp.command || '');
        const body = String(inp.new_string ?? inp.content ?? '');
        if (IDFILE.test(fp) && body) ev.push({ ts, kind: 'write', text: body });
        else if (x.name === 'Bash' && IDWRITE.test(cmd)) ev.push({ ts, kind: 'write', text: cmd });
        else if (x.name === 'Read' && IDFILE.test(fp)) ev.push({ ts, kind: 'Read', pat: `offset=${inp.offset ?? '-'} limit=${inp.limit ?? '-'}` });
        else if (x.name === 'Grep' && (IDFILE.test(fp) || IDFILE.test(String(inp.glob || '')))) ev.push({ ts, kind: 'Grep', pat: String(inp.pattern || '') });
        else if (x.name === 'Bash' && IDFILE.test(cmd) && IDLOOK.test(cmd)) ev.push({ ts, kind: 'Bash', pat: cmd.replace(/\s+/g, ' ') });
      }
    }
    ev.sort((a, b) => a.ts - b.ts);
    for (let i = 0; i < ev.length; i++) {
      if (ev[i].kind !== 'write') continue;
      for (const [head] of target) {
        if (!ev[i].text.includes(head)) continue;
        let prev = null; for (let j = i - 1; j >= 0; j--) if (ev[j].kind !== 'write') { prev = ev[j]; break; }
        const rec = { ts: ev[i].ts, prev, nBefore: ev.slice(0, i).filter(e => e.kind !== 'write').length };
        if (!found.has(head) || found.get(head).ts > rec.ts) found.set(head, rec);
      }
    }
  }
  const rows = [...target.values()].sort((a, b) => a.line - b.line);
  const wrote = rows.filter(s => found.has(s.head));
  const withLook = wrote.filter(s => found.get(s.head).prev);
  console.log(`\n-- reads(重複クラスタの ${rows.length} 節。★書きの規則は rescan より厳しい) --`);
  console.log(`  書きの記録が取れた節              : ${wrote.length} / ${rows.length}`);
  console.log(`  ★書く前に同じセッションで読み/探しがあった: ${withLook.length} / ${wrote.length}`);
  const byKind = {}; for (const s of withLook) { const k = found.get(s.head).prev.kind; byKind[k] = (byKind[k] || 0) + 1; }
  console.log(`  直前の手段                        : ${JSON.stringify(byKind)}`);
  const byWhat = {}; for (const s of withLook) { const p = found.get(s.head).prev; const k = lookKind(p.kind, p.pat); byWhat[k] = (byWhat[k] || 0) + 1; }
  console.log(`  ★★直前に**何を**探したか           : ${JSON.stringify(byWhat)}`);
  console.log(`     content=内容語で探した / numbering=見出し番号を採りに行った / range=末尾か行範囲を見た / other=探していない(git 等)`);
  console.log('\n  -- ★content と判定されたもの(★全部。目で見るため) --');
  for (const s of withLook) { const p = found.get(s.head).prev; if (lookKind(p.kind, p.pat) !== 'content') continue; console.log(`  L${String(s.line).padStart(5)} ${JSON.stringify(p.pat).slice(0, 150)}`); }
  const only = ARGV.includes('--only') ? ARGV[ARGV.indexOf('--only') + 1] : '';
  console.log('\n  -- 節ごと(★--only <語> で見出しを絞れる) --');
  for (const s of rows) {
    if (only && !s.head.includes(only)) continue;
    const r = found.get(s.head);
    if (!r) { console.log(`  L${String(s.line).padStart(5)} 書きの記録なし        ${s.head.slice(0, 60)}`); continue; }
    const p = r.prev;
    console.log(`  L${String(s.line).padStart(5)} ${new Date(r.ts * 1000).toISOString().slice(0, 16)} ${p ? '★読/探 ' + Math.round((r.ts - p.ts) / 60) + ' 分前 ' + p.kind : '読み無し           '}  ${s.head.slice(0, 52)}`);
    if (p) console.log(`             ${JSON.stringify(p.pat).slice(0, 130)}`);
  }
}
// ★★`--similar <file>` —— 提案の口そのもの。新しく書こうとしている節の下書きを渡すと、
//   同じ字面の n-gram を持つ既存の節を名指しする。★`--dupes` と**同じ規則**を使う。
if (has('--similar')) {
  const p = ARGV[ARGV.indexOf('--similar') + 1];
  const GN = ARGV.includes('--gram') ? Number(ARGV[ARGV.indexOf('--gram') + 1]) : GRAM_N;
  if (!p) { console.error('--similar <下書きのファイル>'); process.exit(2); }
  const txt = p === '-' ? fs.readFileSync(0, 'utf8') : fs.readFileSync(p, 'utf8');
  const mine = gramsOf(errLits(txt), GN);
  const hit = [];
  for (const s of secs) { const sh = []; for (const q of gramsOf(s.lits, GN)) if (mine.has(q) && isDistinctive(q, hitsOf(q), N)) sh.push(q); if (sh.length) hit.push({ s, g: sh.sort((a, b) => b.length - a.length)[0] }); }
  console.log(`-- similar(${GN}-gram / --errwords ${ERRWORDS_VER}) 下書きの字面 ${errLits(txt).size} 個 --`);
  if (!hit.length) console.log('  似た節は無い。');
  for (const h of hit) console.log(`  L${String(h.s.line).padStart(5)} <${h.g}> ${h.s.head.slice(0, 76)}`);
  console.log(`  ⇒ ${hit.length} 節`);
}
if (has('--json')) fs.writeFileSync(ARGV[ARGV.indexOf('--json') + 1], JSON.stringify(secs.map(s => ({ line: s.line, head: s.head, regTs: s.regTs, regSrc: s.regSrc, sigs: s.sigs, ev: s.ev }))));
