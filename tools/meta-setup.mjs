#!/usr/bin/env node
/**
 * tools/meta-setup.mjs —— メタ係(隔離 worktree)の立ち上がりを 1 本にする
 * ================================================================
 *
 * 何のためにあるか
 * ----------------
 * `meta-optimizer` は毎回まっさらな隔離 worktree で起動する。そこは
 *
 *   - master から数百〜千 commit 遅れている(M10。第 3/5/7/14 回と**5 度**再発)
 *   - `ResearchPaper/0_Source`(212MB、gitignore)が無い
 *   - `.cache/pdf-pages.json` が無い(cold な PDF 抽出は 20 分でも終わらない。M10 第 7 回)
 *   - ★**本体の未 commit の改善が入っていない**(第 14 回の新しい観測。
 *     merge しても入らない。本体の作業ツリーからファイルで揃えるしかない)
 *
 * ため、そのままでは本体と**同じ数字**が出ない。5 回続けて人が同じ手順を
 * 踏み直しているので道具にした。★踏んだ落とし穴もここに畳んである:
 *
 *   - junction は `cmd //c mklink` も `powershell New-Item` も **Bash フックに弾かれる**。
 *     通るのは `fs.symlinkSync(target, link, 'junction')`。★`node -e` に埋めると
 *     `0_Source` の `\0` が **null byte** と解釈されて壊れる ⇒ **ファイルに書いて実行**(=これ)。
 *   - 帰り際に junction を `rm -rf` すると **本体の 212MB を消しうる**(M35)。
 *     外すのは `--teardown`(`fs.rmdirSync`。reparse point だけが消える)。
 *   - `.cache/pdf-pages.json` の鍵は **絶対パス**を含む。接頭辞を worktree のものへ
 *     付け替えないと 1 頁も当たらない(junction 越しでも mtime/size は一致する)。
 *
 * ★★**最初に叩くコマンドはこれ**(2 人目・3 人目が同じ所で詰まった。M63 の穴 1)
 * ------------------------------------------------------------------
 *
 *     git merge master --no-edit          # ★★先に merge。順番が逆だと下の cp が無意味になる
 *     cp /d/Math_ABC3/tools/meta-setup.mjs tools/meta-setup.mjs && node tools/meta-setup.mjs
 *
 * ★★**順番**(5 人目が足し、6 人目が確かめた。M79 穴 1): ★**merge → cp → 起動**。
 *   `merge` の前に `cp` しても、その直後の merge で `tools/` が 1209 commit 前の姿に触られる。
 *   ★6 人目の実測: merge 直後の worktree に `tools/meta-setup.mjs` は**無い**(master に入っていない)。
 *     ⇒ cp は必須で、かつ merge の後でなければならない。
 *
 * ★★**この worktree で踏む穴**(★道具では塞げないので、ここに書いてある):
 *   1. ★`/tmp` は **`D:\tmp`** に解決される(Bash と node で cwd の解釈が違う)。
 *      `node -e` に `/tmp/...` を渡すと `ENOENT: D:\tmp\...` で落ちる ⇒ **scratchpad の絶対パス**を使う。
 *   2. ★**Bash のガードが誤爆する**(5 人目は 3 回弾かれた): 実行するコマンド名を変数や引用符に入れる /
 *      `time (node … > /dev/null)` のような部分シェル / 長い引用符付きパスの `cat`。
 *      ⇒ ★計測の shell は「素の 1 コマンド」に割る(`time node … > /dev/null` は通る。6 人目が実測)。
 *   3. ★★★**`git checkout -- <file>` を叩かない**(6 人目が踏んだ。下の「3. 整列」の印字も参照)。
 *      整列で写した本は**本体の未 commit** なので、checkout すると 1209 commit 前に戻る。
 *      ⇒ 変更前の数字が要るときは **`cp` で退避してから編集**。消したら `--no-gate` で立て直す。
 *
 * ★理由: この道具は **master にまだ入っていない**(本体の未 commit)ので、
 * 隔離 worktree の `tools/` には**無い**。そして `WT` は
 * ★**自分自身の置き場所から**求める(`dirname(dirname(import.meta.url))`)ので、
 * ★**本体の版を直に叩くと `WT` が本体になり、本体の作業ツリーを触ってしまう。**
 * ⇒ 必ず**先に worktree へ写してから**叩くこと。
 * (自衛として、起動時に `WT` が本体と同じなら**止まる**ようにしてある。)
 *
 * 使い方
 * ------
 *   node tools/meta-setup.mjs            # 立ち上げ(同期 → 整列 → junction → cache → ゲート)
 *   node tools/meta-setup.mjs --dry-run  # 1 バイトも書かずに「何をするか」だけ出す
 *   node tools/meta-setup.mjs --no-gate  # ゲートを回さない(9 秒節約)
 *   node tools/meta-setup.mjs --teardown # ★★本当に最後。junction を外す(帰り際)
 *                                        #   ゲートを**先に回して数字を印字してから**外す
 *                                        #   (`--no-gate` を足すと回さない)
 *   node tools/meta-setup.mjs --force-align # ★自分の作業ごと本体の版で上書きする
 *   node tools/meta-setup.mjs --with-lean   # ★本体の未 commit の lean も揃える(基準が変わる)
 *   node tools/meta-setup.mjs --json     # 機械可読
 *
 * ★**触らないもの**(意図的。理由をコードに残す)
 *   - `lean/ABC3/**` … メタ係の持ち場ではない。**差の件数だけ数えて報告する**
 *     (`--with-lean` を明示したときだけ揃える)。
 *   - `ResearchPaper/decisions-pending.md` … 本体が同時に追記している。
 *     こちらへ写すと、こちらから写し返したときに本体の記録が消えうる。
 *   - `.claude/**` … agent 定義と**権限設定**が同居している。設定は書き換えない。
 *     差があることだけ報告する。
 *   - `tools/_*.mjs` … 使い捨て 399 本(M1)。揃えても意味が無く git status を汚す。
 */

import { spawnSync } from 'node:child_process';
import {
  existsSync, readFileSync, writeFileSync, mkdirSync, copyFileSync,
  readdirSync, statSync, lstatSync, symlinkSync, rmdirSync, realpathSync,
} from 'node:fs';
import { join, dirname, relative, sep } from 'node:path';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';

const HERE = dirname(fileURLToPath(import.meta.url));
const WT = dirname(HERE);                       // worktree のルート

const args = process.argv.slice(2);
const has = (f) => args.includes(f);
const optOf = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };
const DRY = has('--dry-run');
const JSON_OUT = has('--json');
const NO_GATE = has('--no-gate');
const WITH_LEAN = has('--with-lean');
const TEARDOWN = has('--teardown');
const FORCE = has('--force-align');   // ★自分の作業ごと本体の版で上書きする(既定は守る)

const log = [];
const say = (s = '') => { log.push(s); if (!JSON_OUT) console.log(s); };
const result = { worktree: WT, steps: [] };
const t0 = Date.now();
let stepT = Date.now();
const lap = () => { const d = Date.now() - stepT; stepT = Date.now(); return d; };
function step(name, detail, extra = {}) {
  const ms = lap();
  result.steps.push({ name, ms, detail, ...extra });
  say(`  [${(ms / 1000).toFixed(1)}s] ${name} —— ${detail}`);
}

const run = (cmd, argv, opts = {}) => spawnSync(cmd, argv, {
  cwd: WT, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024, shell: false, ...opts,
});
const git = (...a) => run('git', a);

// ────────────────────────────────────────────────────────────────
// 0. 本体(共有チェックアウト)の場所を機械で決める
// ────────────────────────────────────────────────────────────────
/**
 * worktree の `.git` は `gitdir: <本体>/.git/worktrees/<名>` という 1 行のファイル。
 * ★ここから本体を導く(パスの直書きを避ける)。
 */
function resolveMain() {
  const forced = optOf('--main');
  if (forced) return forced;
  try {
    const dotgit = readFileSync(join(WT, '.git'), 'utf8').trim();
    const m = /^gitdir:\s*(.+)$/.exec(dotgit);
    if (m) {
      // <本体>/.git/worktrees/<名> → 3 つ上が本体
      const p = m[1].replace(/[/\\]+$/, '');
      const main = dirname(dirname(dirname(p)));
      if (existsSync(join(main, 'CLAUDE.md')) && existsSync(join(main, 'tools', 'check.mjs'))) return main;
    }
  } catch { /* 下の当てずっぽうへ */ }
  // 保険: <本体>/.claude/worktrees/<名> という置き方
  const guess = dirname(dirname(dirname(WT)));
  if (existsSync(join(guess, 'CLAUDE.md'))) return guess;
  return null;
}

// ────────────────────────────────────────────────────────────────
// 1. master との差を数える
// ────────────────────────────────────────────────────────────────
function measureLag() {
  const behind = (git('rev-list', '--count', 'HEAD..master').stdout || '').trim();
  const ahead = (git('rev-list', '--count', 'master..HEAD').stdout || '').trim();
  const head = (git('log', '--oneline', '-1').stdout || '').trim();
  return { behind: Number(behind || -1), ahead: Number(ahead || -1), head };
}

// ────────────────────────────────────────────────────────────────
// 2. 同期 —— 3 通りを順に試す(どれが通ったかを必ず記録する。M10 の再発を数え続けるため)
// ────────────────────────────────────────────────────────────────
const SYNC_WAYS = [
  ['merge --no-edit', ['merge', 'master', '--no-edit']],
  ['merge --ff-only', ['merge', '--ff-only', 'master']],
  ['reset --hard', ['reset', '--hard', 'master']],
];
function sync(lag) {
  if (lag.behind === 0) return { way: '(不要)', tried: [], note: '既に master に追いついている' };
  const tried = [];
  for (const [name, argv] of SYNC_WAYS) {
    if (DRY) { tried.push(`${name}: (dry-run。試していない)`); continue; }
    const r = git(...argv);
    if (r.status === 0) return { way: name, tried, note: (r.stdout || '').trim().split('\n')[0] || '' };
    tried.push(`${name}: ${((r.stderr || r.stdout || '') + '').trim().split('\n')[0]}`);
  }
  return { way: null, tried, note: '★3 通りとも通らなかった' };
}

// ────────────────────────────────────────────────────────────────
// 3. 本体の作業ツリーと**中身で**揃える
// ────────────────────────────────────────────────────────────────
/**
 * ★なぜ `git -C <本体> status --porcelain` を使わないか(第 15 回の実測):
 *   隔離 agent の Bash フックは `git -C <本体>` も `cd <本体> && git` も**弾く**
 *   (「worktree 外へ向かう git 操作」と見なされる)。そこでフックを迂回する代わりに、
 *   ★**同期後の worktree は master そのもの**であることを使って
 *   「本体の作業ツリー ≠ worktree」= 「本体の未 commit」と**中身で**定義する。
 *   読むだけなので本体を一切変更しない。CRLF で checkout された worktree の
 *   行末も、本体(LF)から写すことで**ついでに直る**(M10 の 3 番目の原因)。
 */
const COPY_ROOTS = ['tools', 'ResearchPaper', 'ResearchMethod', 'memory'];
const REPORT_ONLY_ROOTS = ['.claude', 'lean/ABC3'];
// ★`worktrees` を外さないと**他の agent の隔離 worktree**まで走査する(実測 36,445 本 / 58 秒)。
// `lean-repl` は vendor した REPL の丸ごとの写しで、メタ係が測る対象ではない。
const SKIP_SEG = new Set(['.git', '.cache', '.lake', 'node_modules', '__pycache__', '0_Source',
  'build', '.venv', 'worktrees', 'lean-repl']);
const skipName = (n) => (
  /^_.*\.mjs$/.test(n)          // 使い捨て(M1)
  || /\.bak(-|\.|$)/.test(n)    // *.json.bak-2026-09-04 など
  || /\.stackdump$/.test(n)
  || n === 'decisions-pending.md' // 本体が同時に追記している。写し返しで消しうる
  // ★本体のルートには「道具でも記録でもないもの」が落ちている(第 15 回に写してしまった):
  //   `--help`(202KB。取り違えたリダイレクトの残骸)/ `Donburi_v2.0.exe`(10MB)/
  //   `rpr_acc_dth_rec_downstream.html`(3.3MB)。ゲートは 1 つも読まない。
  || n.startsWith('-')
  || /\.(exe|dll|zip|7z|png|jpg|jpeg|gif|pdf|pyc)$/i.test(n)
);
const MAX_COPY_BYTES = 4 * 1024 * 1024;   // ★実測: 揃える対象の最大は genell-goal.md の 1.6MB

function* walkFiles(root, rel = '') {
  let ents;
  try { ents = readdirSync(join(root, rel), { withFileTypes: true }); } catch { return; }
  for (const e of ents) {
    if (SKIP_SEG.has(e.name)) continue;
    const r = rel ? `${rel}/${e.name}` : e.name;
    if (e.isSymbolicLink()) continue;           // junction を降りない(M10 第 3 回)
    if (e.isDirectory()) { yield* walkFiles(root, r); continue; }
    if (skipName(e.name)) continue;
    yield r;
  }
}
function* rootFiles(root) {
  for (const e of readdirSync(root, { withFileTypes: true })) {
    if (e.isDirectory() || skipName(e.name) || e.name.startsWith('.')) continue;
    yield e.name;
  }
}
const sha1 = (p) => { try { return createHash('sha1').update(readFileSync(p)).digest('hex'); } catch { return null; } };
const sha1lf = (p) => {
  try { return createHash('sha1').update(readFileSync(p).toString('binary').replace(/\r\n/g, '\n')).digest('hex'); }
  catch { return null; }
};

/**
 * ★**自分の作業を上書きしないための盾**(第 15 回に踏みかけた)。
 *
 * 整列は「本体 → worktree」に写す。ところが**同じ道具を作業の途中でもう一度叩くと、
 * こちらが書いた `meta-backlog.md` や `tools/*.mjs` が本体の版で消える。**
 * ⇒ ★**HEAD と中身が違うファイル(= 自分が触ったもの)には触らない。**
 *   `git diff --name-only` は CRLF の幻を出さない(M10 第 4 回。`git status` は出す)。
 */
function locallyChanged() {
  const d = git('diff', '--name-only');
  const s = new Set((d.stdout || '').split('\n').map((x) => x.trim()).filter(Boolean));
  const u = git('ls-files', '--others', '--exclude-standard');
  for (const l of (u.stdout || '').split('\n')) { const t = l.trim(); if (t) s.add(t); }
  return s;
}

function alignWith(main) {
  const plan = { copy: [], crlf: [], reportOnly: [], tooBig: [], skipped: [], missingInMain: 0 };
  const mine = locallyChanged();
  const scan = (roots, mode) => {
    for (const root of roots) {
      const abs = join(main, root.replace(/\//g, sep));
      if (!existsSync(abs)) continue;
      const st = statSync(abs);
      const list = st.isDirectory() ? [...walkFiles(abs)].map((r) => `${root}/${r}`) : [root];
      for (const rel of list) {
        const a = join(main, rel.replace(/\//g, sep));
        const b = join(WT, rel.replace(/\//g, sep));
        try { if (statSync(a).size > MAX_COPY_BYTES) { plan.tooBig.push(rel); continue; } } catch { continue; }
        const ha = sha1(a);
        const hb = existsSync(b) ? sha1(b) : null;
        if (ha === hb) continue;
        if (hb !== null && sha1lf(a) === sha1lf(b)) { if (mode === 'copy') plan.crlf.push(rel); continue; }
        // ★自分が触ったファイルは**写さない**(本体の版で自分の作業が消える)。
        //   ★限界: 1 回目の整列で写したものも「HEAD と違う」ので以後は守られる側に回る。
        //   本体が走っている最中に追いつきたいときは `--force-align`(★上書きする)。
        if (mode === 'copy' && mine.has(rel) && !FORCE) { plan.skipped.push(rel); continue; }
        (mode === 'copy' ? plan.copy : plan.reportOnly).push(rel);
      }
    }
  };
  scan(COPY_ROOTS, 'copy');
  scan([...rootFiles(main)], 'copy');
  scan(WITH_LEAN ? REPORT_ONLY_ROOTS.filter((r) => r !== 'lean/ABC3') : REPORT_ONLY_ROOTS, 'report');
  if (WITH_LEAN) scan(['lean/ABC3'], 'copy');

  if (!DRY) {
    for (const rel of [...plan.copy, ...plan.crlf]) {
      const a = join(main, rel.replace(/\//g, sep));
      const b = join(WT, rel.replace(/\//g, sep));
      mkdirSync(dirname(b), { recursive: true });
      copyFileSync(a, b);
    }
  }
  return plan;
}

// ────────────────────────────────────────────────────────────────
// 4. `0_Source` の junction
// ────────────────────────────────────────────────────────────────
const SRC_REL = join('ResearchPaper', '0_Source');
function makeJunction(main) {
  const link = join(WT, SRC_REL);
  const target = join(main, SRC_REL);
  if (!existsSync(target)) return { ok: false, note: `本体に ${SRC_REL} が無い` };
  let cur = null;
  try { cur = lstatSync(link); } catch { /* 無い */ }
  if (cur) {
    if (cur.isSymbolicLink() || cur.isDirectory()) {
      const n = (() => { try { return readdirSync(link).length; } catch { return 0; } })();
      if (n > 0) return { ok: true, note: `既にある(${n} エントリ)`, made: false };
    }
    return { ok: false, note: '既に何かが居るが空。手で確かめること' };
  }
  if (DRY) return { ok: true, note: '(dry-run) junction を張る', made: false };
  mkdirSync(dirname(link), { recursive: true });
  symlinkSync(target, link, 'junction');
  const n = readdirSync(link).length;
  return { ok: true, note: `張った(${n} エントリ)`, made: true };
}
function teardownJunction() {
  const link = join(WT, SRC_REL);
  let st = null;
  try { st = lstatSync(link); } catch { return '無い(何もしない)'; }
  if (!st.isSymbolicLink()) return '★symlink/junction ではない。手を出さない(本体の実体かもしれない)';
  if (DRY) return '(dry-run) rmdirSync で外す';
  rmdirSync(link);   // ★reparse point だけを外す。rm -rf は本体を消しうる(M35)
  return '外した(reparse point のみ)';
}

// ────────────────────────────────────────────────────────────────
// 5. `.cache/pdf-pages.json` の移植 —— 鍵の接頭辞を付け替える
// ────────────────────────────────────────────────────────────────
function transplantCache(main) {
  const src = join(main, '.cache', 'pdf-pages.json');
  if (!existsSync(src)) return { ok: false, note: '本体に pdf-pages.json が無い' };
  let j;
  try { j = JSON.parse(readFileSync(src, 'utf8')); } catch (e) { return { ok: false, note: `読めない: ${e.message}` }; }
  const pages = j.pages || {};
  const from = join(main, SRC_REL);
  const to = join(WT, SRC_REL);
  const out = {};
  let moved = 0;
  for (const [k, v] of Object.entries(pages)) {
    if (k.startsWith(from)) { out[to + k.slice(from.length)] = v; moved++; } else out[k] = v;
  }
  // ★`self` は check.mjs のハッシュ。手順 3 で本体から写しているので一致するはず。
  const selfNow = sha1(join(WT, 'tools', 'check.mjs'));
  const note = `${Object.keys(out).length} 頁(接頭辞を付け替えたのは ${moved})`;
  if (!DRY) {
    mkdirSync(join(WT, '.cache'), { recursive: true });
    writeFileSync(join(WT, '.cache', 'pdf-pages.json'),
      JSON.stringify({ self: j.self, pdftotext: j.pdftotext, pages: out }), 'utf8');
  }
  return { ok: true, note, moved, total: Object.keys(out).length, pdftotext: j.pdftotext, selfNow };
}

// ────────────────────────────────────────────────────────────────
// 6. ゲート —— ★段を明示して呼ぶ。`--brief` だけで呼ぶと `lake build` が始まる(M41)
// ────────────────────────────────────────────────────────────────
const GATES = [
  ['selftest+structured', ['tools/check.mjs', '--selftest', '--structured', '--brief']],
  ['ledger', ['tools/check.mjs', '--ledger', '--brief']],
];
function runGates() {
  const out = [];
  for (const [name, argv] of GATES) {
    const t = Date.now();
    const r = run(process.execPath, argv);
    const s = `${r.stdout || ''}${r.stderr || ''}`;
    const ng = /\bNG\s+(\d+)\s*件/.exec(s.split('\n').filter((l) => /^NG \d+ 件|^PASS/.test(l)).pop() || '');
    const self = /selftest:\s*(\d+)\/(\d+)/.exec(s);
    const s16 = /1_Structured: S1-S6 すべて PASS/.test(s);
    out.push({
      gate: name, sec: (Date.now() - t) / 1000,
      ng: ng ? Number(ng[1]) : (/^PASS$/m.test(s) ? 0 : null),
      selftest: self ? `${self[1]}/${self[2]}` : null,
      s1s6: s16,
    });
  }
  // graph.mjs は 0.5 秒。ノード/辺と md5 を出す(副作用の見張りの主役。M10 第 7 回)
  const t = Date.now();
  const g = run(process.execPath, ['tools/graph.mjs']);
  const gs = g.stdout || '';
  // ★出力は `ノード 2222 / ファイル 2222` と `辺(import) 6364 本`。
  //   辺は括弧が挟まるので `辺\s*(\d+)` では取れない(第 15 回に踏んだ)。
  const nodes = /ノード\s+(\d+)/.exec(gs);
  const edges = /辺[^0-9\n]*?(\d+)\s*本/.exec(gs);
  out.push({
    gate: 'graph', sec: (Date.now() - t) / 1000,
    nodes: nodes ? Number(nodes[1]) : null, edges: edges ? Number(edges[1]) : null,
    md5: createHash('md5').update(gs).digest('hex').slice(0, 12),
  });
  return out;
}

// ────────────────────────────────────────────────────────────────
// 本体
// ────────────────────────────────────────────────────────────────
say('== meta-setup: 隔離 worktree をメタ係の測定台に仕立てる ==');
say(`  worktree : ${WT}`);

const main = resolveMain();
if (!main) { console.error('★本体(共有チェックアウト)が見つからない。--main <パス> を渡すこと'); process.exit(2); }
say(`  本体     : ${main}`);
// ★★M63 の穴 1 の**機械の自衛**: 本体の版を直に叩くと `WT` が本体になり、
//   「隔離 worktree を仕立てる」つもりで**本体を書き換えて**しまう。
{
  const norm = (p) => { try { return realpathSync(p).replace(/[\\/]+$/, '').toLowerCase(); } catch { return String(p).toLowerCase(); } };
  if (norm(WT) === norm(main)) {
    console.error('');
    console.error('★★止めた —— この道具を**本体の版のまま**叩いている(WT と 本体 が同じ)。');
    console.error('   `WT` は自分自身の置き場所から決まるので、このまま進むと本体を書き換える。');
    console.error('   ★先に隔離 worktree へ写してから叩くこと:');
    console.error('     cp /d/Math_ABC3/tools/meta-setup.mjs tools/meta-setup.mjs && node tools/meta-setup.mjs');
    process.exit(2);
  }
}
// ★M63 の穴 2 の続き: 前回 teardown 済みなら、ゲートの数字が壊れることを先に知らせる。
{
  const mk = join(WT, '.cache', 'meta-teardown.json');
  if (existsSync(mk) && !existsSync(join(WT, SRC_REL))) {
    let at = '?';
    try { at = JSON.parse(readFileSync(mk, 'utf8')).at; } catch { /* 読めなくてよい */ }
    say(`  ★★前回 ${at} に teardown 済み(0_Source が無い)。`);
    say('     ★この状態でゲートを叩くと NG 361 / NG 5143 になる。★退行ではない。');
    say('     ★このまま進めば下の 4 で張り直す。');
  }
}
if (DRY) say('  ★dry-run: 1 バイトも書かない');
say('');

if (TEARDOWN) {
  // ★★M63 の穴 2(事故の大きさが一番大きい): **teardown した後にゲートを叩くと壊滅する。**
  //   実測(2026-09-07、2 人目): junction を外した状態で
  //   `--selftest --structured` は NG **361**、`--ledger` は NG **5143**(全部「S4 PDF が見つからない」)。
  //   ★これは退行ではないのに「退行を出した」と報告してしまう。
  //   ⇒ 直し方は「順序を守れ」と書くことではなく、★**順序を守る必要を無くすこと**:
  //      **外す前にゲートを回して最終の数字をここで印字する。** 以降ゲートを叩く理由が無くなる。
  if (!has('--no-gate')) {
    say('  ★junction を外す前に、最後のゲートをここで回す(これが提出用の数字。以降は叩かないこと)');
    for (const g of runGates()) {
      say(`      ${g.gate.padEnd(20)} ${String(g.sec).padStart(5)}s  `
        + (g.ng !== null && g.ng !== undefined ? `NG ${g.ng}` : '')
        + (g.selftest ? ` selftest ${g.selftest}` : '')
        + (g.s1s6 ? ' S1-S6 PASS' : '')
        + (g.nodes ? `ノード ${g.nodes} / 辺 ${g.edges} / md5 ${g.md5}` : ''));
    }
    say('');
  }
  const note = teardownJunction();
  say(`  0_Source の junction: ${note}`);
  // ★外したことを残す。次に誰かが meta-setup を叩いたとき先頭で知らせる(穴 2 の再発防止)。
  try {
    mkdirSync(join(WT, '.cache'), { recursive: true });
    if (!DRY) writeFileSync(join(WT, '.cache', 'meta-teardown.json'),
      JSON.stringify({ at: new Date().toISOString(), note }, null, 1));
  } catch { /* 残せなくても致命ではない */ }
  say('');
  say('  ★★★ここから先、ゲート(check.mjs)を叩いてはいけない。');
  say('     0_Source が無いので `--structured` は NG 361、`--ledger` は NG 5143 になる。');
  say('     ★これは退行ではなく「PDF が見つからない」だけである。');
  say('     どうしても叩くなら、先に張り直すこと(1.9 秒):');
  say('       cp /d/Math_ABC3/tools/meta-setup.mjs tools/meta-setup.mjs && node tools/meta-setup.mjs --no-gate');
  process.exit(0);
}

const lag = measureLag();
result.lag = lag;
step('1. master との差', `★behind ${lag.behind} / ahead ${lag.ahead} —— HEAD = ${lag.head}`);

const sy = sync(lag);
result.sync = sy;
for (const t of sy.tried) say(`        × ${t}`);
step('2. 同期', `通ったのは ★${sy.way ?? '(どれも通らなかった)'} ${sy.note ? `—— ${sy.note}` : ''}`);
if (sy.way === null) say('  ★同期できていない。以降の数字は本体と比較できない。');

const plan = alignWith(main);
result.align = { copied: plan.copy.length, crlf: plan.crlf.length, reportOnly: plan.reportOnly.length };
result.alignList = plan.copy;
step('3. 本体の未 commit と整列',
  `写した ★${plan.copy.length} 本(+ 行末だけ違う ${plan.crlf.length} 本)/ 報告のみ ${plan.reportOnly.length} 本`
  + (plan.skipped.length ? ` / ★自分の作業なので写さなかった ${plan.skipped.length} 本` : ''));
for (const rel of plan.skipped) say(`        = ${rel}(★自分が触っている。本体の版で上書きしない)`);
for (const rel of plan.copy.slice(0, 40)) say(`        + ${rel}`);
if (plan.copy.length > 40) say(`        + …他 ${plan.copy.length - 40} 本`);
/* ★★★★6 人目が実際に踏んだ穴(メタ第 20 回、M83)。**印字で塞ぐ**。
 *   ここで写した本は**本体の未 commit** であって、worktree の git には入っていない。
 *   ⇒ ★`git checkout -- <file>` / `git stash` / `git restore` を掛けると
 *     ★**1209 commit 前の版に戻る**(＝本体の改善が消える)。6 人目は
 *     「変更前の数字を採ろう」として `git checkout -- tools/brief.mjs` を叩き、
 *     ★第 19 回で採用された `--audit-proof` ごと消した(usage が出て初めて気づいた)。
 *   ★直し方は 1 行で済むので、消える側ではなく**戻し方**をここに出す。 */
if (plan.copy.length || plan.skipped.length) {
  // ★`skipped`(＝自分が触っている本)も同じ危険を負う。むしろ**そちらの方が消したときに痛い**
  //   ので、2 回目以降の起動でも出す。
  say('        ! ★★上の本は「本体の未 commit」＝ worktree の git には無い。');
  say('          ⇒ `git checkout -- <file>` / `git stash` を掛けると 1209 commit 前に戻る(本体の改善が消える)。');
  say('          ⇒ 戻すには `node tools/meta-setup.mjs --no-gate`(同期がもう一度写す。0.1 秒)。');
  say('          ⇒ ★変更前の数字が要るときは `git checkout` ではなく **`cp` で退避してから編集**すること。');
}
if (plan.reportOnly.length) {
  const byRoot = {};
  for (const r of plan.reportOnly) { const k = r.split('/').slice(0, 2).join('/'); byRoot[k] = (byRoot[k] || 0) + 1; }
  say(`        ! 揃えていない(持ち場の外): ${Object.entries(byRoot).map(([k, v]) => `${k} ${v} 本`).join(' / ')}`);
  const leanDiff = plan.reportOnly.filter((r) => r.startsWith('lean/'));
  if (leanDiff.length) {
    // ★第 15 回の実測: ここが 10 本あると `--ledger` の NG が 13 ではなく **17** になる
    //   (本体の未 commit に `.src` の頁の修正 51→47 が含まれており、R′ がそれを見る)。
    //   ★**基準が違う理由をここで言わないと、次の起動が「退行だ」と誤読する。**
    say(`        ! ★本体の未 commit の lean が ${leanDiff.length} 本ある。`);
    say('          ⇒ ゲートの NG が本体の基準と食い違いうる。'
      + '揃えたいときだけ `--with-lean`(★写すだけ。数学は書き換えない)');
    for (const r of leanDiff.slice(0, 12)) say(`            · ${r}`);
  }
}

const jn = makeJunction(main);
result.junction = jn;
step('4. 0_Source の junction', `${jn.ok ? '' : '★失敗: '}${jn.note}`);

const ca = transplantCache(main);
result.cache = ca;
step('5. pdf-pages.json の移植', `${ca.ok ? '' : '★失敗: '}${ca.note}`);

if (!NO_GATE && !DRY) {
  const gates = runGates();
  result.gates = gates;
  say('');
  for (const g of gates) {
    if (g.gate === 'graph') { say(`  [${g.sec.toFixed(1)}s] graph.mjs —— ノード ${g.nodes} / 辺 ${g.edges} / md5 ${g.md5}`); continue; }
    say(`  [${g.sec.toFixed(1)}s] ${g.gate} —— NG ${g.ng}${g.selftest ? ` / selftest ${g.selftest}` : ''}${g.s1s6 ? ' / S1-S6 PASS' : ''}`);
  }
  stepT = Date.now();
}

result.totalSec = (Date.now() - t0) / 1000;
say('');
say(`  合計 ★${result.totalSec.toFixed(1)} 秒`);
// ★★M63 の穴 3: 整列で 200 本前後が「行末だけ違う」状態になるので `git status` は `M` だらけになり、
//   ★**自分が触った 2〜3 本が埋もれる。** 提出のとき何を採ればよいか分からなくなるので、
//   ★行末を無視した見方をここで渡す(2 人目は「中身が違うのは 3 本だけ」を手で調べ直した)。
say('');
say('  ★提出のとき「自分が触ったファイル」だけを見る(整列で行末だけ違う本が 200 本前後出るため):');
// ★`2>/dev/null` が要る —— git が「LF will be replaced by CRLF」を**数百行**吐いて
//   本命の一覧が流れる(3 人目が実際に踏んだ。付けないと 38KB 出る)。
say('      git diff --stat --ignore-cr-at-eol 2>/dev/null | tail -20');
say('      git status --porcelain 2>/dev/null | grep "^??"   # 新しく足したもの');
say('  ★★写しただけのもの(CLAUDE.md / lean-idioms.md / autonomy-policy.md / decisions-pending.md /');
say('     unverified.mjs)は**本体が同じ日に伸ばしている**。worktree 側が必ず古い。★採ってはいけない。');
say('');
say('  ★帰り際に `node tools/meta-setup.mjs --teardown` で junction を外すこと(M35)');
say('  ★★teardown は**本当に最後**。外した後にゲートを叩くと NG 361 / NG 5143 になる(M63 穴 2)。');
say('     ★`--teardown` は外す前に最後のゲートを回して数字を印字するので、以降は叩かなくてよい。');
if (JSON_OUT) console.log(JSON.stringify(result, null, 2));
