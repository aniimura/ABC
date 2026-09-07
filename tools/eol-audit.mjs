#!/usr/bin/env node
// eol-audit.mjs —— 木全体の改行を数える口(メタ第 31 回。M152)
//
// ★なぜ要るか
//   改善係は「触る前に改行を数えろ」と毎回言われているのに、数える口が無かった。
//   その結果 **CRLF の本を LF で書き戻して 1 ファイル丸ごと差分になる**事故が
//   第 29・30 回で起きている(`tools/check.mjs` / `tools/decl-index.mjs` /
//   `lean/ABC3/Found.lean` / `lean/ABC3/Found/PGC/FilteredGroup.lean`)。
//
// ★★`grep -c $'\r'` は Git Bash では嘘をつく(パス変換と `\r` の扱いで 0 を返す)。
//   ⇒ **バイトを直接数える**。この本は node で `\r\n` / 裸の `\n` / 裸の `\r` を数える。
//
// 使い方:
//   node tools/eol-audit.mjs                 … CRLF を含む本だけを出す(既定)
//   node tools/eol-audit.mjs --all           … 全部の本を出す
//   node tools/eol-audit.mjs --mixed         … ★CRLF と LF が混ざった本だけ(いちばん危ない)
//   node tools/eol-audit.mjs --ext .mjs,.md  … 拡張子を絞る
//   node tools/eol-audit.mjs --dir tools     … 部分木を絞る
//   node tools/eol-audit.mjs --root D:/Math_ABC3  … ★別のチェックアウトを見る(本体を読むだけ)
//   node tools/eol-audit.mjs --json
//   node tools/eol-audit.mjs --selftest
//
// ★見落とさないための約束:
//   - `.git` / `.lake` / `node_modules` / `.cache` / `ResearchPaper/0_Source` は**降りない**
//     (0_Source は 212MB の junction。降りると本体を舐める)。
//   - 2MB を超える本は「大きすぎる」として**行数だけ**数える(中身は読むが保持しない)。

import fs from 'node:fs';
import path from 'node:path';
import { spawnSync } from 'node:child_process';

const SKIP_DIR = new Set(['.git', '.lake', 'node_modules', '.cache', '0_Source', '__pycache__', '.claude']);
const DEF_EXT = ['.mjs', '.js', '.md', '.lean', '.json', '.py', '.yml', '.yaml', '.toml', '.txt', '.html'];

/**
 * ★純関数。バッファから改行を数える。
 * 返り値: { crlf, lf, cr, total, kind }
 *   crlf … `\r\n` の数 / lf … `\r` を伴わない `\n` の数 / cr … `\n` を伴わない `\r` の数
 *   kind … 'crlf' | 'lf' | 'mixed' | 'cr' | 'none'
 */
export function countEol(buf) {
  let crlf = 0, lf = 0, cr = 0;
  for (let i = 0; i < buf.length; i++) {
    const b = buf[i];
    if (b === 13) {                       // \r
      if (i + 1 < buf.length && buf[i + 1] === 10) { crlf++; i++; }
      else cr++;
    } else if (b === 10) {                // 裸の \n
      lf++;
    }
  }
  const total = crlf + lf + cr;
  let kind = 'none';
  if (total > 0) {
    const nz = [crlf > 0, lf > 0, cr > 0].filter(Boolean).length;
    if (nz > 1) kind = 'mixed';
    else if (crlf > 0) kind = 'crlf';
    else if (lf > 0) kind = 'lf';
    else kind = 'cr';
  }
  return { crlf, lf, cr, total, kind };
}

/** ★純関数。バイト列が binary っぽいか(NUL を含む)。 */
export function looksBinary(buf) {
  const n = Math.min(buf.length, 8192);
  for (let i = 0; i < n; i++) if (buf[i] === 0) return true;
  return false;
}

export function* walk(root, rel = '') {
  let ents;
  try { ents = fs.readdirSync(path.join(root, rel), { withFileTypes: true }); } catch { return; }
  for (const e of ents) {
    const r = rel ? path.join(rel, e.name) : e.name;
    if (e.isDirectory()) {
      if (SKIP_DIR.has(e.name)) continue;
      yield* walk(root, r);
    } else if (e.isFile()) {
      yield r;
    }
  }
}

/**
 * ★★★M164(メタ第 33 回)—— **git が追跡している本の集合**を返す。
 * ★なぜ要るか(実測): 既定の `audit()` は D:\Math_ABC3 で **15,662 本**を舐めるが、
 *   ★**git が追跡しているのは 2,638 本(16.8%)だけ**である。残り 13,024 本の内訳は
 *   `external/` 10,285 / `scratch/` 2,380 / `tools/`(生成物) 355。
 * ★★M155 の公表値 `crlf 9,534` は、★**その 9,249(97%)が git の知らない本**だった。
 *   ⇒ ★「プロジェクトの姿」ではなく「他人の repo の姿」を数えていた。
 * ★`git ls-files` は実測 **0.1 秒**。★木が git でないときは null を返す(呼び手が既定に落ちる)。
 */
export function trackedSet(root) {
  const r = spawnSync('git', ['ls-files'], {
    cwd: root, encoding: 'utf8', maxBuffer: 256 * 1024 * 1024, shell: false,
  });
  if (r.status !== 0) return null;
  const s = new Set();
  for (const line of String(r.stdout || '').split('\n')) {
    const p = line.trim();
    if (p) s.add(p);
  }
  return s.size ? s : null;
}

export function audit(root, opts = {}) {
  const exts = opts.exts ?? DEF_EXT;
  const sub = opts.dir ?? '';
  // ★M164: `tracked` を渡すと git が追跡している本だけを見る。★既定は従来どおり(公表値を動かさない)。
  const only = opts.tracked ? (opts.trackedSet ?? trackedSet(root)) : null;
  const rows = [];
  for (const rel of walk(root, sub)) {
    if (exts.length && !exts.includes(path.extname(rel).toLowerCase())) continue;
    if (only && !only.has(rel.split(path.sep).join('/'))) continue;
    let st;
    try { st = fs.statSync(path.join(root, rel)); } catch { continue; }
    let buf;
    try { buf = fs.readFileSync(path.join(root, rel)); } catch { continue; }
    if (looksBinary(buf)) continue;
    const c = countEol(buf);
    rows.push({ file: rel.split(path.sep).join('/'), size: st.size, ...c });
  }
  rows.sort((a, b) => (a.file < b.file ? -1 : 1));
  return rows;
}

// ══════════════════════════════════════════════════════════════════════════
// ★★★M161 —— `.gitattributes` を **入れずに** その費用を測る(メタ第 32 回。持ち場 2)
// ══════════════════════════════════════════════════════════════════════════
/**
 * ★問い(本体が判断するために要る数字): ★**`* text=auto eol=lf` を入れたら差分はどれだけ出るか。**
 *   ★第 31 回は「入れると全ファイルが 1 度だけ巨大な差分になる恐れがある」と書いた。★測る。
 *
 * ★★**入れずに測れる**理由: `git ls-files --eol` は 1 行につき
 *     `i/<index の改行>  w/<作業ツリーの改行>  attr/<いま効いている属性>  <path>`
 *   を返す。★入れたときに起きることは、この 2 列だけで決まる:
 *
 *   | いまの姿 | `eol=lf` を入れると |
 *   |---|---|
 *   | `i/crlf` / `i/mixed` / `i/cr` | ★**index の blob に `\r` がある** ⇒ `add --renormalize` で**巨大 commit** |
 *   | `i/lf` かつ `w/crlf` | index は既に LF ⇒ ★**commit の差分は 0**。作業ツリーは次の checkout で LF になる |
 *   | `i/-text`(バイナリ) | 触られない |
 *
 * ★★**この関数は 1 バイトも書かない。**★`git ls-files` を読むだけ。
 */
export function gitattrCost(root) {
  const r = spawnSync('git', ['ls-files', '--eol'], {
    cwd: root, encoding: 'utf8', maxBuffer: 256 * 1024 * 1024, shell: false,
  });
  if (r.status !== 0) return { ok: false, note: ((r.stderr || '') + '').trim().split('\n')[0] };
  return { ok: true, ...parseEolLines(r.stdout || '') };
}

/** ★純関数(git を呼ばない ⇒ selftest で較正できる)。`git ls-files --eol` の出力を読む。 */
export function parseEolLines(text) {
  const rows = [];
  let unparsed = 0;
  for (const line of String(text || '').split('\n')) {
    if (!line.trim()) continue;
    const m = /^i\/(\S+)\s+w\/(\S+)\s+attr\/(\S*)\s*\t(.*)$/.exec(line);
    if (!m) { unparsed++; continue; }
    rows.push({ i: m[1], w: m[2], attr: m[3], path: m[4] });
  }
  const tally = (k) => {
    const t = {};
    for (const x of rows) t[x[k]] = (t[x[k]] || 0) + 1;
    return Object.entries(t).sort((a, b) => b[1] - a[1]);
  };
  // ★1 度きりの巨大 commit になるのは index に \r がある本だけ
  const indexDirty = rows.filter((x) => x.i === 'crlf' || x.i === 'mixed' || x.i === 'cr');
  // ★作業ツリーだけが書き換わる本(commit の差分は 0)
  const workOnly = rows.filter((x) => x.i === 'lf' && (x.w === 'crlf' || x.w === 'mixed' || x.w === 'cr'));
  const byExt = (list) => {
    const t = {};
    for (const x of list) {
      const e = (x.path.match(/(\.[A-Za-z0-9]+)$/) || ['(なし)'])[0];
      t[e] = (t[e] || 0) + 1;
    }
    return Object.entries(t).sort((a, b) => b[1] - a[1]).slice(0, 12);
  };
  return {
    tracked: rows.length, unparsed, i: tally('i'), w: tally('w'),
    attrs: [...new Set(rows.map((x) => x.attr))].slice(0, 8),
    indexDirty: { n: indexDirty.length, byExt: byExt(indexDirty), sample: indexDirty.slice(0, 15).map((x) => x.path) },
    workOnly: { n: workOnly.length, byExt: byExt(workOnly) },
    inert: rows.length - indexDirty.length - workOnly.length,
  };
}

// ══════════════════════════════════════════════════════════════════════════
// ★★★M162 —— 「混在だけが危ないのか」を**道具の側**で数える(メタ第 32 回)
// ══════════════════════════════════════════════════════════════════════════
/**
 * ★本体の見立て:「★混在は道具を黙って壊すが、★一貫している限り CRLF でも LF でも実害が無い」。
 * ★★**測ったら、この切り分けは成り立たない。**
 *
 *   素の `split(改行)` に 3 つの姿を食わせると(`--fragile` が実演する):
 *     LF 一貫   … `\r` が残る行 **0**
 *     CRLF 一貫 … `\r` が残る行 **全部**   ← ★★混在より**悪い**
 *     混在      … `\r` が残る行 **混ざった分だけ**
 *   ⇒ ★**危ないのは「混在か否か」ではなく「LF か否か」である。**
 *     ★混在が目立つのは**壊れ方がまだらで気づきにくい**からであって、被害が大きいからではない。
 *
 * ★だから `--mixed` をゲートに入れるのは**必要だが十分ではない**。
 * ★ここは「道具の字面」を数えるだけの粗い測り方である(★どの `split` が
 *   ファイルを食っていて、どれが `git` の出力を食っているかは字面では割れない)。⇒ ★**上限**である。
 */
export function splitFragility(dir) {
  // ★正規表現はここに直接書く(★Bash の heredoc に書くとフックが CR を潰す。M30)
  const RE_NAIVE = [new RegExp("split\\(\\s*['\"]\\\\n['\"]\\s*\\)", 'g'), new RegExp('split\\(\\s*/\\\\n/', 'g')];
  const RE_SAFE = [new RegExp('split\\(\\s*/\\\\r\\?\\\\n/', 'g'), new RegExp('split\\(\\s*/\\\\r\\\\n\\|\\\\n/', 'g')];
  const rows = [];
  let files = 0, naive = 0, safe = 0;
  let ents = [];
  try { ents = fs.readdirSync(dir); } catch { return { files: 0, naive: 0, safe: 0, rows: [] }; }
  for (const f of ents) {
    if (!f.endsWith('.mjs') || f.startsWith('_')) continue;   // 使い捨ては数えない(M1)
    let src; try { src = fs.readFileSync(path.join(dir, f), 'utf8'); } catch { continue; }
    files++;
    const n = RE_NAIVE.reduce((a, re) => a + (src.match(re) || []).length, 0);
    const s = RE_SAFE.reduce((a, re) => a + (src.match(re) || []).length, 0);
    naive += n; safe += s;
    if (n || s) rows.push({ f, n, s });
  }
  rows.sort((a, b) => b.n - a.n || (a.f < b.f ? -1 : 1));
  return { files, naive, safe, rows };
}

/** ★純関数。素の `split(改行)` が 3 つの姿をどう壊すかを実演する。 */
export function fragilityDemo() {
  const NL = String.fromCharCode(10), CR = String.fromCharCode(13);
  const body = ['alpha', 'beta', 'gamma'];
  const cases = {
    'LF 一貫': body.join(NL),
    'CRLF 一貫': body.join(CR + NL),
    '混在': ['alpha' + NL, 'beta' + CR + NL, 'gamma'].join(''),
  };
  return Object.entries(cases).map(([name, text]) => {
    const lines = text.split(NL);
    return {
      name, n: lines.length,
      dirty: lines.filter((l) => l.endsWith(CR)).length,
      eq: lines.filter((l, i) => l === body[i]).length,
    };
  });
}

function selftest() {
  let ok = 0, ng = 0;
  const t = (n, c) => { if (c) ok++; else { ng++; console.log('  NG ' + n); } };
  const B = (s) => Buffer.from(s, 'utf8');

  t('countEol: 純 LF', (() => { const c = countEol(B('a\nb\nc\n')); return c.lf === 3 && c.crlf === 0 && c.cr === 0 && c.kind === 'lf'; })());
  t('countEol: 純 CRLF', (() => { const c = countEol(B('a\r\nb\r\n')); return c.crlf === 2 && c.lf === 0 && c.cr === 0 && c.kind === 'crlf'; })());
  t('countEol: 混在', (() => { const c = countEol(B('a\r\nb\nc\r\n')); return c.crlf === 2 && c.lf === 1 && c.kind === 'mixed'; })());
  t('countEol: 裸の CR(旧 Mac)', (() => { const c = countEol(B('a\rb\rc')); return c.cr === 2 && c.crlf === 0 && c.lf === 0 && c.kind === 'cr'; })());
  t('countEol: 改行なし', (() => { const c = countEol(B('abc')); return c.total === 0 && c.kind === 'none'; })());
  t('countEol: 末尾が裸の CR', (() => { const c = countEol(B('a\r')); return c.cr === 1 && c.crlf === 0; })());
  t('countEol: CRLF を \\r + \\n に二重計上しない', (() => { const c = countEol(B('\r\n')); return c.crlf === 1 && c.lf === 0 && c.cr === 0; })());
  t('countEol: \\n\\r は CRLF ではない', (() => { const c = countEol(B('\n\r')); return c.crlf === 0 && c.lf === 1 && c.cr === 1; })());
  // ★CRLF が 1 行だけ混ざった LF ファイル(第 30 回が踏んだ形)
  t('countEol: LF の中に CRLF 1 行', (() => { const c = countEol(B('a\nb\r\nc\n')); return c.kind === 'mixed' && c.crlf === 1 && c.lf === 2; })());
  // ★walk が降りない場所(第 31 回の突然変異 E5 が素通りしたので足した)
  {
    const os2 = process.env.TEMP || process.env.TMP || '.';
    const base = path.join(os2, 'eol-audit-selftest-' + process.pid);
    fs.rmSync(base, { recursive: true, force: true });
    for (const d of ['.git', '.lake', 'node_modules', '.cache', '0_Source', 'ok']) fs.mkdirSync(path.join(base, d), { recursive: true });
    for (const d of ['.git', '.lake', 'node_modules', '.cache', '0_Source', 'ok']) fs.writeFileSync(path.join(base, d, 'a.md'), 'x\r\n');
    fs.writeFileSync(path.join(base, 'top.md'), 'y\n');
    const got = [...walk(base)].map(s => s.split(path.sep).join('/')).sort();
    t('walk: ★降りないディレクトリを本当に飛ばす', got.length === 2 && got[0] === 'ok/a.md' && got[1] === 'top.md');
    const rows = audit(base, { exts: ['.md'] });
    t('audit: ★飛ばした先は表に出ない', rows.length === 2 && !rows.some(r => /\.git|\.lake|node_modules|\.cache|0_Source/.test(r.file)));
    t('audit: ★CRLF の本を crlf と判定する', rows.find(r => r.file === 'ok/a.md').kind === 'crlf');
    t('audit: ★LF の本を lf と判定する', rows.find(r => r.file === 'top.md').kind === 'lf');
    // ★★M164(メタ第 33 回)—— `tracked` で母集団を絞る
    t('M164: ★tracked の集合に無い本は落ちる',
      audit(base, { exts: ['.md'], tracked: true, trackedSet: new Set(['top.md']) })
        .map(r => r.file).join(',') === 'top.md');
    t('M164: ★tracked を渡さなければ従来どおり(既定を動かしていない)',
      audit(base, { exts: ['.md'] }).length === 2);
    t('M164: ★tracked が空集合なら 0 本(黙って全部通さない)',
      audit(base, { exts: ['.md'], tracked: true, trackedSet: new Set() }).length === 0);
    t('M164: ★path の区切りは / で突き合わせる(Windows の \\ でも当たる)',
      audit(base, { exts: ['.md'], tracked: true, trackedSet: new Set(['ok/a.md']) })
        .map(r => r.file).join(',') === 'ok/a.md');
    fs.rmSync(base, { recursive: true, force: true });
  }
  // ★★M161 —— `git ls-files --eol` の読み(★git を呼ばない純関数として較正する)
  {
    const TAB = String.fromCharCode(9);
    const L = (i, w, p, attr) => `i/${i}\tw/${w}\tattr/${attr || ''}${TAB}${p}`.replace(/\t/g, '  ').replace(/  ([^ ]*)$/, TAB + '$1');
    const mk = (rows) => rows.map(([i, w, p]) => `i/${i}  w/${w}  attr/  ${TAB}${p}`).join(String.fromCharCode(10));
    void L;
    const c1 = parseEolLines(mk([['lf', 'crlf', 'a.mjs'], ['lf', 'lf', 'b.md'], ['-text', '-text', 'c.png']]));
    t('gitattr: 3 本読める', c1.tracked === 3 && c1.unparsed === 0);
    t('gitattr: ★index に CR が無ければ巨大 commit は 0', c1.indexDirty.n === 0);
    t('gitattr: ★作業ツリーだけ書き換わる本を数える', c1.workOnly.n === 1);
    t('gitattr: 何も起きない本を数える', c1.inert === 2);
    const c2 = parseEolLines(mk([['crlf', 'crlf', 'x.md'], ['mixed', 'crlf', 'y.md'], ['lf', 'crlf', 'z.md']]));
    t('gitattr: ★★index が CRLF/mixed なら巨大 commit に数える', c2.indexDirty.n === 2);
    t('gitattr: ★その内訳(拡張子)を出す', c2.indexDirty.byExt[0][0] === '.md' && c2.indexDirty.byExt[0][1] === 2);
    t('gitattr: ★見本を出す', c2.indexDirty.sample.includes('x.md') && c2.indexDirty.sample.includes('y.md'));
    t('gitattr: 空でも落ちない', parseEolLines('').tracked === 0);
    t('gitattr: ★読めない行を黙って数に入れない',
      (() => { const c = parseEolLines('ゴミ' + String.fromCharCode(10) + mk([['lf', 'lf', 'a.md']])); return c.tracked === 1 && c.unparsed === 1; })());
    t('gitattr: ★空白を含むパスも 1 本として読む',
      parseEolLines(mk([['lf', 'crlf', 'ResearchPaper/1_Structured/A Version of X/index.html']])).workOnly.n === 1);
    t('gitattr: i/ と w/ の内訳を出す',
      c1.i.find(([k]) => k === 'lf')[1] === 2 && c1.w.find(([k]) => k === 'crlf')[1] === 1);
  }
  // ★★M162 —— 「混在だけが危ない」を否定する実演そのものを試験にする
  {
    const d = Object.fromEntries(fragilityDemo().map((x) => [x.name, x]));
    t('fragile: LF 一貫は 1 行も汚れない', d['LF 一貫'].dirty === 0 && d['LF 一貫'].eq === 3);
    t('fragile: ★★CRLF 一貫は混在より多くの行を汚す', d['CRLF 一貫'].dirty > d['混在'].dirty);
    t('fragile: ★CRLF 一貫は素の比較を通さない', d['CRLF 一貫'].eq < 3);
    t('fragile: 混在も素の比較を通さない', d['混在'].eq < 3);
    t('fragile: 3 つの姿とも行数は同じ(壊れ方は行数に出ない)',
      d['LF 一貫'].n === 3 && d['CRLF 一貫'].n === 3 && d['混在'].n === 3);
    const fr = splitFragility(path.join(new URL('.', import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, '$1')));
    t('fragile: tools/ を数えられる(0 本ではない)', fr.files > 10);
    t('fragile: 素と安全を別に数える', fr.naive > 0 && fr.safe > 0);
    t('fragile: ★使い捨て `_*.mjs` を数えない', !fr.rows.some((r) => r.f.startsWith('_')));
    t('fragile: 無いディレクトリでも落ちない', splitFragility('/no/such/dir').files === 0);
  }
  t('looksBinary: NUL あり', looksBinary(Buffer.from([65, 0, 66])));
  t('looksBinary: NUL なし', !looksBinary(B('hello')));
  // ★UTF-8 の多バイト(日本語)を \r/\n と取り違えない
  t('countEol: 日本語を誤検出しない', (() => { const c = countEol(B('抽象核\n改行\n')); return c.lf === 2 && c.cr === 0 && c.crlf === 0; })());

  console.log(`eol-audit selftest: ok ${ok} / ng ${ng}`);
  return ng === 0;
}

function main() {
  const argv = process.argv.slice(2);
  const has = (f) => argv.includes(f);
  const val = (f, d) => { const i = argv.indexOf(f); return i >= 0 && argv[i + 1] ? argv[i + 1] : d; };
  if (has('--selftest')) process.exit(selftest() ? 0 : 1);

  const root = path.resolve(val('--root', null)
    ?? path.join(new URL('.', import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, '$1'), '..'));
  // ★★M161 —— `.gitattributes` を入れたときの費用(★入れずに測る)
  if (has('--gitattr')) {
    const c = gitattrCost(root);
    if (!c.ok) { console.error('★git ls-files --eol が読めない: ' + c.note); process.exit(2); }
    if (has('--json')) { console.log(JSON.stringify({ root, ...c }, null, 1)); return; }
    console.log('== .gitattributes(`* text=auto eol=lf`)の費用 —— ★測っただけ。入れていない ==');
    console.log(`  ★数えた木 : ${root}`);
    console.log('  ★★「どの木で数えたか」を必ず添えること(第 31 回の教訓。木ごとに改行が真逆になる)');
    console.log(`  git の追跡下 : ${c.tracked} 本(解釈できなかった行 ${c.unparsed})`);
    console.log(`  いま効いている attr : ${JSON.stringify(c.attrs)}`);
    console.log('');
    console.log(`  index 側(i/) : ${c.i.map(([k, v]) => `${k} ${v}`).join(' / ')}`);
    console.log(`  作業ツリー(w/): ${c.w.map(([k, v]) => `${k} ${v}`).join(' / ')}`);
    console.log('');
    console.log(`  ★★1 度きりの巨大 commit になる本(index に CR がある) : ${c.indexDirty.n} 本`);
    if (c.indexDirty.n) {
      console.log(`      拡張子 : ${c.indexDirty.byExt.map(([k, v]) => `${k} ${v}`).join(' / ')}`);
      for (const p of c.indexDirty.sample) console.log(`      · ${p}`);
      if (c.indexDirty.n > 15) console.log(`      · …他 ${c.indexDirty.n - 15} 本`);
    } else {
      console.log('      ⇒ ★**0 本**。`core.autocrlf=true` が入口で正規化し続けてきたので、');
      console.log('        ★`.gitattributes` を足しても **commit の差分は出ない**。');
    }
    console.log(`  ★作業ツリーだけ書き換わる本(commit の差分 0) : ${c.workOnly.n} 本`);
    if (c.workOnly.n) console.log(`      拡張子 : ${c.workOnly.byExt.map(([k, v]) => `${k} ${v}`).join(' / ')}`);
    console.log(`  何も起きない本 : ${c.inert} 本`);
    console.log('');
    console.log('  ★注意: 足しただけでは作業ツリーは**書き換わらない**(次の checkout から LF になる)。');
    console.log('    ⇒ 既存の木を今すぐ揃えたいなら再 checkout が要る。★未 commit の作業がある木では危険。');
    return;
  }

  // ★★M162 —— 「混在だけが危ないのか」を道具の側で数える
  if (has('--fragile')) {
    const fr = splitFragility(path.join(root, 'tools'));
    const demo = fragilityDemo();
    if (has('--json')) { console.log(JSON.stringify({ root, ...fr, demo }, null, 1)); return; }
    console.log('== ★道具の側の「改行の食べ方」(★M162。粗い上限) ==');
    console.log(`  ★数えた木 : ${root}`);
    console.log(`  tools/*.mjs(使い捨て _*.mjs を除く) : ${fr.files} 本`);
    console.log(`  ★素の split(改行) —— CRLF を食うと **全行に CR が残る** : ${fr.naive} 箇所`);
    console.log(`  ★安全な split(CR? + 改行)                              : ${fr.safe} 箇所`);
    console.log(`  ⇒ 素の割合 : ${(100 * fr.naive / Math.max(1, fr.naive + fr.safe)).toFixed(1)}%`);
    console.log('');
    console.log('    素   安全  ファイル');
    for (const r of fr.rows.slice(0, 14)) {
      console.log(`   ${String(r.n).padStart(4)} ${String(r.s).padStart(5)}  ${r.f}`);
    }
    console.log('');
    console.log('== ★★実演: 素の split(改行) に 3 つの姿を食わせる ==');
    for (const d of demo) {
      console.log(`  ${d.name.padEnd(10)} 行 ${d.n} / ★CR が残った行 ${d.dirty} / 素の比較が通った行 ${d.eq}/3`);
    }
    console.log('');
    console.log('  ★★結論: ★**「一貫していれば安全」は成り立たない。**');
    console.log('    ★CRLF 一貫は混在より**多くの行**を汚す。混在が怖いのは「まだらで気づけない」から。');
    console.log('    ⇒ ★`--mixed` をゲートに入れるのは**必要だが十分ではない**。');
    console.log('  ★★ただしこれは字面の上限である。★どの split がファイルを食っているかは字面では割れない。');
    return;
  }

  const extArg = val('--ext', null);
  const exts = extArg ? extArg.split(',').map(s => (s.startsWith('.') ? s : '.' + s).toLowerCase()) : DEF_EXT;
  // ★★M164(メタ第 33 回)—— `--tracked` で「git が追跡している本」だけに絞る。
  //   ★既定は**変えない**。M155 / M161 / M162 の公表値がこの母集団で出ているため。
  const wantTracked = has('--tracked');
  const only = wantTracked ? trackedSet(root) : null;
  if (wantTracked && !only) {
    console.error('★--tracked: `git ls-files` が読めない(この木は git ではない?)。既定の母集団に落とす。');
  }
  const rows = audit(root, { exts, dir: val('--dir', ''), tracked: !!only, trackedSet: only });

  let show = rows;
  if (has('--mixed')) show = rows.filter(r => r.kind === 'mixed');
  else if (!has('--all')) show = rows.filter(r => r.crlf > 0 || r.cr > 0);

  if (has('--json')) { console.log(JSON.stringify({ root, tracked: !!only, n: rows.length, rows: show }, null, 1)); return; }

  const by = {};
  for (const r of rows) by[r.kind] = (by[r.kind] || 0) + 1;
  console.log(`== eol-audit —— ${root}`);
  console.log(`   走査 ${rows.length} 本 (${exts.join(' ')})`
    + (only ? `  ★--tracked: git の追跡下 ${only.size} 本に絞った` : ''));
  if (!only) {
    console.log('   ★★母集団は**木の中の全ファイル**である(M164)。'
      + '★`external/` や `scratch/` を含むので');
    console.log('     ★「プロジェクトの姿」を見たいときは **`--tracked`** を付けること。');
  }
  console.log(`   内訳: ` + Object.entries(by).sort().map(([k, v]) => `${k} ${v}`).join(' / '));
  console.log('');
  if (!show.length) { console.log('   ★CRLF / 裸の CR を含む本は無い。'); return; }
  console.log('   種別   CRLF    LF    CR  ファイル');
  for (const r of show) {
    console.log(`   ${r.kind.padEnd(6)}${String(r.crlf).padStart(5)}${String(r.lf).padStart(6)}${String(r.cr).padStart(6)}  ${r.file}`);
  }
  console.log('');
  console.log('   ★★この一覧に載っている本を編集するときは、書き戻す改行を合わせること。');
  console.log('     LF で書き戻すと 1 ファイル丸ごと差分になり、採用の単位が読めなくなる。');
}

if (process.argv[1] && process.argv[1].endsWith('eol-audit.mjs')) main();
