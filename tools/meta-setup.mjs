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
 *   node tools/meta-setup.mjs --selftest # ★M10 の原因判定だけを較正する(git を呼ばない)
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
  // ★★★M151(メタ第 31 回)—— ★**4 人続けてここで止まった**(第 27・28・29・30 回)。
  //   ★ここが**本体そのもの**なら `.git` は**ファイルではなくディレクトリ**である。
  //   従来は下の `readFileSync` が EISDIR で落ち、保険の当てずっぽうも外れて `null` を返し、
  //   ★「本体が見つからない。`--main <パス>` を渡すこと」という**見当違いの**文で終わっていた
  //   (★本体は目の前にある。足りないのは worktree のほうである)。
  //   ⇒ ★**WT 自身を返す**。すると下の「WT === 本体」の番人が発火し、
  //     `noWorktreeMsg()` が **worktree の作り方**を名指しで出す(それが要る文である)。
  try { if (statSync(join(WT, '.git')).isDirectory()) return WT; } catch { /* .git が無ければ下へ */ }
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

/**
 * ★★M10(worktree が毎回 1209〜1210 commit 遅れる)の**原因**を測る。
 *
 * ★第 21 回までの見立て:「push すれば直る」→ ★**外れ**(push したら 1209 → 1210 に増えた)。
 * ★第 22 回の測定: worktree の枝は **リポジトリの既定枝(origin/main)** から切られている。
 *   ところが本体が押しているのは `master` なので `origin/main` は動かない。
 *   ⇒ ★**push では永久に解消しない。** 本体側で 1 度 `master` を `main` に流すしかない。
 * ★ここは**測って印字するだけ**。worktree の中からは直せない(押す先が worktree の外)。
 */
function diagnoseLag(lag, g = git) {
  if (!(lag.behind > 0)) return null;
  const out = [];
  for (const ref of ['origin/main', 'origin/master', 'main', 'master']) {
    if (g('rev-parse', '--verify', '--quiet', ref).status !== 0) continue;
    out.push({
      ref,
      cutHere: g('merge-base', '--is-ancestor', 'HEAD', ref).status === 0,
      // ★HEAD からその枝の先端までの距離。★**切り元なら 0**(HEAD がその枝の先端そのもの)。
      dist: Number((g('rev-list', '--count', `HEAD..${ref}`).stdout || '').trim() || '-1'),
      gap: Number((g('rev-list', '--count', `${ref}..master`).stdout || '').trim() || '-1'),
      tip: (g('log', '-1', '--format=%h %ad %s', '--date=short', ref).stdout || '').trim().slice(0, 60),
    });
  }
  // ★★M98 の宿題: 「HEAD はこの枝の中にある」は **HEAD を含む枝すべて**に付くので識別に効かない
  //   (使い捨て repo で 4/4 に付いた)。★HEAD を含む枝のうち **距離が最小のもの** だけを候補にする。
  //   ★同点なら同点のまま出す(同じ commit を指す別名のことがある)。★1 本に絞らない = 断定しない。
  const contains = out.filter((r) => r.cutHere && r.dist >= 0);
  const near = contains.length ? Math.min(...contains.map((r) => r.dist)) : -1;
  for (const r of out) r.nearest = r.cutHere && r.dist === near;
  const cause = out.find((r) => r.cutHere && Math.abs(r.gap - lag.behind) <= 2);
  return { refs: out, cause, near };
}

/**
 * ★★M135(第 28 回の宿題)—— **遅れているなら、まず大きく警告して逃げ道を出す**。
 *
 * ★なぜ要るか: **3 回続けて立ち上がりで落ちている**。
 *   第 27 回 = worktree が渡されていない / 第 28 回 = `main` から切られ master より 1,210 commit 前
 *   (`tools/meta-setup.mjs` すら無く MODULE_NOT_FOUND) / 第 29 回 = また worktree が渡されていない。
 * ★逃げ道は **`git reset --hard master`(実測 0.4 秒)** だが、これは台帳の奥にしか書いていなかった。
 * ★ここは**純関数**(行の配列を返すだけ)。★git を呼ばない ⇒ selftest で較正できる。
 */
export function lagBanner(lag) {
  if (!(lag && lag.behind > 0)) return [];
  const n = lag.behind;
  const loud = n >= 100 ? '★★★★' : n >= 10 ? '★★★' : '★★';
  return [
    '',
    `  ${loud} 警告 —— この worktree は master より ${n} commit 前である。`,
    '     ⇒ このままでは **本体に無い道具・無い台帳**を見ることになり、数字が本体と食い違う。',
    `     ⇒ ★**逃げ道は 1 つ**(自分の枝を動かすだけ。実測 0.4 秒):`,
    '',
    '           git reset --hard master && node tools/meta-setup.mjs',
    '',
    n >= 1000
      ? '     ★遅れが 1,000 を超える = **`main` から切られている**(M92 / M135)。この場合 `tools/` に'
        + '\n       `meta-setup.mjs` すら無く、最初の起動が MODULE_NOT_FOUND で落ちる。★先に reset すること。'
      : '     ★下の「2. 同期」も merge で埋めにいくが、遅れが大きいときは reset の方が速い。',
  ];
}

/**
 * ★★隔離 worktree が**そもそも渡されていない**ときの案内。★純関数(行の配列)。
 * ★第 27 回と第 29 回がこれ。★従来の文は「worktree へ写してから叩け」だけで、
 *   ★**worktree の作り方**が書いていなかった(無いものへは写せない)。
 */
export function noWorktreeMsg(wt, nth = 'NN') {
  return [
    '',
    '★★止めた —— cwd が **本体そのもの**である(WT と 本体 が同じ)。',
    `   いま : ${wt}`,
    '   `WT` は自分自身の置き場所から決まるので、このまま進むと**本体を書き換える**。',
    '',
    '   ★隔離 worktree が渡されていないなら、**自分で 1 本作る**(実測 1 秒):',
    '',
    `       git worktree add -b meta${nth} /d/Math_ABC3/.claude/worktrees/meta${nth} master`,
    `       cd /d/Math_ABC3/.claude/worktrees/meta${nth}`,
    '       node tools/meta-setup.mjs',
    '',
    '   ★`-b` で**新しい枝**を切るので master は動かない。帰り際は `--teardown` のあと',
    '     `git worktree remove` は**しない**(本体が採否を読む)。',
  ];
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
 * ★**自分の作業を上書きしないための盾**(第 15 回に踏みかけた)。★純関数(git を呼ばない)。
 *
 * 整列は「本体 → worktree」に写す。ところが**同じ道具を作業の途中でもう一度叩くと、
 * こちらが書いた `meta-backlog.md` や `tools/*.mjs` が本体の版で消える。**
 * ⇒ ★**HEAD と中身が違うファイル(= 自分が触ったもの)には触らない。**
 *   `git diff --name-only` は CRLF の幻を出さない(M10 第 4 回。`git status` は出す)。
 *
 * ★★**M157(第 31 回が踏んだ)—— `git diff` だけ見ると commit で盾が外れる。**
 *   事前登録は「データを見る前に書いた」ことを示すため ★**commit で時刻を刻むのが作法**である
 *   (M152 はそうした)。ところが commit した瞬間にその本は `git diff --name-only` から消え、
 *   ★**「自分の作業」でなくなって次の整列で本体の版に上書きされる。**
 *   ★実測(第 32 回が再現): commit 前は `= tools/meta-setup.mjs`、commit 直後は `+ tools/meta-setup.mjs`。
 *   ⇒ ★**枝側の差分 `master...HEAD` を足す**。3 点(merge-base からの差)なので、
 *     master が進んでも**自分の commit だけ**が入る(本体の前進を「自分の作業」と誤認しない)。
 *
 * @param workOut      `git diff --name-only` の stdout
 * @param untrackedOut `git ls-files --others --exclude-standard` の stdout
 * @param branchOut    `git diff --name-only master...HEAD` の stdout(読めなければ '')
 */
export function changedSets(workOut, untrackedOut, branchOut) {
  const toSet = (s) => new Set(String(s ?? '').split('\n').map((x) => x.trim()).filter(Boolean));
  const work = toSet(workOut);
  for (const x of toSet(untrackedOut)) work.add(x);
  const branch = toSet(branchOut);
  const all = new Set(work);
  for (const x of branch) all.add(x);
  // ★commit したから盾が外れかけていた本(= 枝にしか居ない本)。印字して見せる。
  const onlyCommitted = [...branch].filter((x) => !work.has(x)).sort();
  return { all, work, branch, onlyCommitted };
}

/** ★上を本物の git で呼ぶ薄い皮。★`g` を差し替えると selftest で較正できる。 */
function locallyChanged(g = git) {
  const out = (r) => (r && r.status === 0 ? (r.stdout || '') : '');
  // ★`master` が無い / HEAD が master そのもの、でも落ちないこと(status≠0 は空として扱う)。
  return changedSets(
    out(g('diff', '--name-only')),
    out(g('ls-files', '--others', '--exclude-standard')),
    out(g('diff', '--name-only', 'master...HEAD')),
  );
}

function alignWith(main) {
  const plan = {
    copy: [], crlf: [], reportOnly: [], tooBig: [], skipped: [], missingInMain: 0, savedByCommit: [],
  };
  const changed = locallyChanged();
  const mine = changed.all;
  plan.onlyCommitted = changed.onlyCommitted;
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
        if (mode === 'copy' && mine.has(rel) && !FORCE) {
          plan.skipped.push(rel);
          // ★M157: 「commit したから `git diff` から消えていた」本を名指しで数える。
          //   ここが 0 のままなら盾の新しい半分は一度も効いていない ⇒ 効果を測れる。
          if (!changed.work.has(rel)) plan.savedByCommit.push(rel);
          continue;
        }
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
    for (const l of noWorktreeMsg(WT)) console.error(l);
    console.error('   ★既に worktree があるのにここへ来たなら、先に写してから叩くこと:');
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

/**
 * ★★M157 の盾の較正。★本物の git を呼ばない(スタブを注入する)。
 * ★試験の穴を作らないための約束(第 30 回・第 31 回の反省):
 *   - ★**純関数 `changedSets` だけでなく、`locallyChanged` の配線も試す**
 *     (第 31 回は `--since` を CLI に埋めていて selftest が届かなかった)。
 *   - ★**逆(守らない側)も試す**。守りすぎると worktree が本体に追いつけなくなる。
 */
function selftestShield() {
  // git のスタブ。呼ばれた引数を記録する。
  const calls = [];
  const mk = (work, untracked, branch, branchStatus = 0) => (...a) => {
    calls.push(a.join(' '));
    if (a[0] === 'ls-files') return { status: 0, stdout: untracked };
    if (a[2] === 'master...HEAD') return { status: branchStatus, stdout: branch };
    return { status: 0, stdout: work };
  };
  const S = (o) => [...o.all].sort().join(',');
  const both = changedSets('a.md\nb.md', 'u.md', 'c.mjs');
  return [
    // —— 純関数として ——
    ['★M157 commit 済みの変更も「自分の作業」に入る', both.all.has('c.mjs')],
    ['★M157 未 commit の変更は今までどおり入る', both.all.has('a.md') && both.all.has('b.md')],
    ['★M157 untracked も今までどおり入る', both.all.has('u.md')],
    ['★M157 全部で 4 本(重複なし)', both.all.size === 4],
    ['★M157 枝にしか居ない本を名指しできる', S({ all: new Set(both.onlyCommitted) }) === 'c.mjs'],
    ['★M157 work と枝で重なる本は onlyCommitted に出さない',
      changedSets('x.md', '', 'x.md').onlyCommitted.length === 0
      && changedSets('x.md', '', 'x.md').all.size === 1],
    ['★M157 空でも落ちない', changedSets('', '', '').all.size === 0],
    ['★M157 null/undefined でも落ちない',
      changedSets(null, undefined, null).all.size === 0],
    ['★M157 空行と余白を落とす', changedSets('  a.md  \n\n\n', '', '').all.size === 1],
    // —— ★配線として(locallyChanged 経由。ここが第 31 回の穴) ——
    ['★★M157 locallyChanged が 3 点 master...HEAD を引く',
      (() => { calls.length = 0; locallyChanged(mk('', '', '')); return calls.includes('diff --name-only master...HEAD'); })()],
    ['★★M157 locallyChanged が commit 済みを盾に入れる',
      locallyChanged(mk('', '', 'k.mjs')).all.has('k.mjs')],
    ['★★M157 master が無くて(status≠0)も落ちず、work だけで動く',
      (() => { const r = locallyChanged(mk('w.md', '', 'ゴミ', 128)); return r.all.has('w.md') && !r.all.has('ゴミ') && r.all.size === 1; })()],
    ['★M157 従来の 2 本(diff / ls-files)を引き続き引く',
      (() => { calls.length = 0; locallyChanged(mk('', '', '')); return calls.includes('diff --name-only') && calls.includes('ls-files --others --exclude-standard'); })()],
    // —— ★逆(守りすぎないこと)——
    ['★★M157 触っていない本は盾に入らない(整列は今までどおり進む)',
      !locallyChanged(mk('a.md', 'u.md', 'c.mjs')).all.has('CLAUDE.md')],
    ['★M157 何も触っていなければ盾は空(全部写す)',
      locallyChanged(mk('', '', '')).all.size === 0],
  ];
}

/** ★M10 の原因判定の較正。★本物の git を呼ばない(スタブを注入する)。 */
function selftestLag() {
  const stub = (refs, cutRef) => (...a) => {
    const [cmd, x, y, z] = a;
    if (cmd === 'rev-parse') return { status: refs[z] === undefined ? 1 : 0, stdout: '' };
    if (cmd === 'merge-base') return { status: z === cutRef ? 0 : 1, stdout: '' };
    if (cmd === 'rev-list') {
      const arg = String(y);
      if (arg.startsWith('HEAD..')) return { status: 0, stdout: arg.slice(6) === cutRef ? '0' : '99' };
      return { status: 0, stdout: String(refs[arg.replace('..master', '')] ?? 0) };
    }
    return { status: 0, stdout: 'deadbeef 2026-09-03 tip' };
  };
  const R = { 'origin/main': 1210, 'origin/master': 0, master: 0 };
  // ★★M98 が使い捨て repo で観測した形(HEAD が **全枝の祖先**)を作る。
  //   dist: origin/main と main は 0(切り元)、origin/master と master は 3。
  const allContain = (dists) => (...a) => {
    const [cmd, x, y, z] = a;
    if (cmd === 'rev-parse') return { status: dists[z] === undefined ? 1 : 0, stdout: '' };
    if (cmd === 'merge-base') return { status: 0, stdout: '' };          // ★どの枝も HEAD を含む
    if (cmd === 'rev-list') {
      const arg = String(y);
      if (arg.startsWith('HEAD..')) return { status: 0, stdout: String(dists[arg.slice(6)] ?? 0) };
      return { status: 0, stdout: '3' };
    }
    void x; void cmd;
    return { status: 0, stdout: 'deadbeef 2026-09-03 tip' };
  };
  const D = { 'origin/main': 0, 'origin/master': 3, main: 0, master: 3 };
  const marked = diagnoseLag({ behind: 3 }, allContain(D)).refs.filter((r) => r.nearest);
  const checks = [
    // ★★M98 の宿題: 全枝が HEAD を含む形でも、印は**最も近い枝だけ**に付くこと。
    //   ★これが壊れる(印を cutHere に戻す)と 4 本になって鳴る。
    ['★全枝が HEAD を含んでも印は最短の枝だけに付く(4 → 2)', marked.length === 2],
    ['★印が付くのは距離 0 の枝(origin/main と main)',
      marked.map((r) => r.ref).sort().join(',') === 'main,origin/main'],
    ['★距離を測って持っている', diagnoseLag({ behind: 3 }, allContain(D)).refs
      .every((r) => Number.isFinite(r.dist))],
    ['★同点はどちらも残す(1 本に絞って断定しない)', marked.length > 1],
    ['★最短の距離を返す', diagnoseLag({ behind: 3 }, allContain(D)).near === 0],
    ['遅れが 0 なら何も言わない', diagnoseLag({ behind: 0 }, stub(R, 'origin/main')) === null],
    ['origin/main から切られていれば原因と断定する',
      diagnoseLag({ behind: 1210 }, stub(R, 'origin/main'))?.cause?.ref === 'origin/main'],
    ['遅れの数と枝の古さが食い違えば断定しない',
      !diagnoseLag({ behind: 7 }, stub(R, 'origin/main'))?.cause],
    ['どの枝にも入っていなければ断定しない',
      !diagnoseLag({ behind: 1210 }, stub(R, 'nowhere'))?.cause],
    ['無い枝は表に出さない',
      diagnoseLag({ behind: 1210 }, stub(R, 'origin/main')).refs.every((r) => r.ref !== 'main')],
    // ★★M135(第 29 回で実装)—— 立ち上がりで 3 回落ちているので、警告そのものを試験する。
    ['★遅れ 0 なら警告を出さない', lagBanner({ behind: 0 }).length === 0],
    ['★遅れ 0 は ahead があっても黙る', lagBanner({ behind: 0, ahead: 5 }).length === 0],
    ['★lag が無くても落ちない', lagBanner(null).length === 0 && lagBanner(undefined).length === 0],
    ['★遅れ 1 でも警告する', lagBanner({ behind: 1 }).length > 0],
    ['★警告に遅れの数が入る', lagBanner({ behind: 1210 }).join('\n').includes('1210 commit 前')],
    ['★★警告に逃げ道(reset --hard master)が入る',
      lagBanner({ behind: 1210 }).join('\n').includes('git reset --hard master')],
    ['★遅れが大きいほど星が増える',
      lagBanner({ behind: 1210 }).join('\n').includes('★★★★')
      && !lagBanner({ behind: 3 }).join('\n').includes('★★★★')],
    ['★1000 超なら main から切られた形だと言う',
      lagBanner({ behind: 1210 }).join('\n').includes('MODULE_NOT_FOUND')
      && !lagBanner({ behind: 30 }).join('\n').includes('MODULE_NOT_FOUND')],
    // ★worktree が渡されていない場合(第 27 回・第 29 回)
    ['★★worktree の作り方を出す', noWorktreeMsg('D:/Math_ABC3').join('\n').includes('git worktree add -b')],
    ['★作る場所に master を指定している', noWorktreeMsg('D:/Math_ABC3').join('\n').includes('worktrees/metaNN master')],
    ['★いまの cwd を出す', noWorktreeMsg('D:/Math_ABC3').join('\n').includes('D:/Math_ABC3')],
    ['★番号を差し替えられる', noWorktreeMsg('X', 29).join('\n').includes('meta29')],
    // ★★★M157(第 31 回が踏み、第 32 回が再現した)—— commit すると盾が外れる。
    //   ★実測の再現: commit 前 `= tools/meta-setup.mjs` / commit 直後 `+ tools/meta-setup.mjs`。
    ...selftestShield(),
  ];
  let ok = 0;
  for (const [label, pass] of checks) { say(`  ${pass ? 'ok ' : 'NG '} ${label}`); if (pass) ok++; }
  say(`\n  selftest(M10 の原因判定): ${ok}/${checks.length} PASS`);
  return ok === checks.length;
}
if (has('--selftest')) process.exit(selftestLag() ? 0 : 1);

const lag = measureLag();
result.lag = lag;
step('1. master との差', `★behind ${lag.behind} / ahead ${lag.ahead} —— HEAD = ${lag.head}`);
// ★★M135: 遅れているなら **原因の診断より先に** 大きく警告し、逃げ道を出す。
for (const l of lagBanner(lag)) say(l);
const lagWhy = diagnoseLag(lag);
result.lagWhy = lagWhy;
if (lagWhy) {
  for (const r of lagWhy.refs) {
    const mark = r.nearest ? '  ★HEAD に最も近い(= 切り元の候補)'
      : r.cutHere ? `  (この枝にも含まれる。HEAD から ${r.dist} commit 先)`
      : '  (HEAD を含まない)';
    say(`        ${r.ref.padEnd(14)} master より ${String(r.gap).padStart(5)} commit 古い` +
      `${mark}  ${r.tip}`);
  }
  if (lagWhy.cause) {
    say(`  ! ★★M10 の原因はこれ: worktree は **${lagWhy.cause.ref}** から切られており、`);
    say(`        その枝が master より ${lagWhy.cause.gap} commit 古い。`);
    say('        ⇒ ★**本体が master へ push しても解消しない**(第 21 回に実測。1209 → 1210 に増えただけ)。');
    say("        ⇒ 直すのは本体側で 1 度だけ: `git push origin master:main`(または PR を 1 本 merge)。");
    say('        ⇒ それまでは下の「2. 同期」が毎回 merge して埋める(★7 回連続、競合 0)。');
  } else {
    say('  ! ★M10 の原因は特定できなかった(切り元の枝が見つからない)。★そう書いて次へ進むこと。');
  }
}

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
for (const rel of plan.skipped) {
  // ★M157: commit で `git diff` から消えていた本は、そう名指しする(盾のどちら半分が効いたか)。
  const why = plan.savedByCommit.includes(rel) ? '★commit 済みだが自分の作業' : '★自分が触っている';
  say(`        = ${rel}(${why}。本体の版で上書きしない)`);
}
if (plan.savedByCommit.length) {
  say(`        ! ★★M157 の盾が効いた —— commit 済みの ${plan.savedByCommit.length} 本を守った`
    + '(旧版はここで本体の版に戻していた)。');
}
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
