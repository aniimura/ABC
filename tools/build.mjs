#!/usr/bin/env node
/**
 * tools/build.mjs —— `lake build` を **1 回だけ**回してログをファイルに落とし、以後はそこを見る
 * =========================================================================
 *
 * 何のためにあるか（M130 / M134 の実測。メタ第 28・29 回）
 * -----------------------------------------------------
 * 本体の道具待ち 86.2 時間のうち `lake build` が **41.6 時間**（実測。M130 の 44.1 時間は
 * heredoc の本文に書かれた "lake build" を 1,112 回ぶん数え込んでいた。第 29 回が補正）。
 * そのうち機械で「無駄」と言い切れるのは 2 つだけである:
 *
 *   ★1. **同じ命令の中で 2 度以上建てている** … 313 命令 / うち同じ対象の繰り返し 217 回。
 *        典型は `lake build ABC3 2>&1 | grep error | head -3; …; lake build ABC3 2>&1 | tail -2`
 *        —— ★**出力を保存していないので、別の切り口を見るためにもう一度建てている。**
 *   ★2. **間に 1 文字も書かずにもう一度建てている** … 609 回 / 合計 **6.0 時間**（中央 10.4 秒）。
 *        —— ★これも「さっきの出力をもう一度見たい」だけのことが多い。
 *
 * ⇒ ★**建てるのは 1 回、出力はファイル、以後は切り出す。** それだけの道具。
 *    実測: 出力をファイルに落としていたのは **84 / 4,739 命令**（1.8%）しかない。
 *
 * ★★この道具が「安全側を下げない」理由
 * -----------------------------------
 * ★対象を絞る（`lake build ABC3.Found.X`）のは**速いが見落とす**（下の `--why` を読むこと）。
 * ★この道具の既定は **対象を絞らない**。減らすのは「**同じものを 2 度建てる**」だけなので、
 * ★**見落としは増えない。** 対象を絞りたいときは自分で引数に書く（判断は人がする）。
 *
 * 使い方
 * ------
 *   node tools/build.mjs                     # 木全部を 1 回建てて .cache/build-ABC3.log に落とし、要約を出す
 *   node tools/build.mjs ABC3.Found.PGC.X    # 対象を指定（★見落としうる。--why を読むこと）
 *   node tools/build.mjs --if-stale          # ★.lean が 1 本も新しくなければ**建てずに**前のログを読む
 *   node tools/build.mjs --errors            # ★建てない。直近のログから error だけ
 *   node tools/build.mjs --sorry             # ★建てない。`declaration uses 'sorry'` を数えて列挙
 *                                            #   （★`grep sorry` は当てにならない —— docstring を拾う）
 *   node tools/build.mjs --jobs              # ★建てない。成否と jobs 数
 *   node tools/build.mjs --grep RE [--tail N] --context N
 *   node tools/build.mjs --files             # ★建てない。error を出しているファイルの一覧
 *   node tools/build.mjs --why               # ★木全部と対象指定の交換（実測の数字）を出す
 *   node tools/build.mjs --dry-run           # 何を叩くかだけ出す
 *   node tools/build.mjs --selftest          # 純関数の試験（lake を呼ばない）
 *   node tools/build.mjs --json              # 機械可読
 *
 * ★**読み出しの口（--errors / --sorry / --jobs / --grep / --files）は lake を呼ばない。**
 *   ⇒ 2 度目の切り口は **建て直しではなくファイル読み**になる（実測 10.4 秒 → 0.0 秒）。
 *
 * ★置き場: `.cache/build-<対象>.log`（`.cache/` は gitignore 済み）。
 *   ★`<対象>` は `ABC3` / `ABC3.Found.PGC.X` などをそのまま使う（複数なら `+` で繋ぐ）。
 */

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
export const ROOT = path.dirname(HERE);
const LEAN = path.join(ROOT, 'lean');
const CACHE = path.join(ROOT, '.cache');

// ────────────────────────────────────────────────────────────
// 純関数（★selftest はここだけを見る。lake は呼ばない）
// ────────────────────────────────────────────────────────────

/** 対象の並びからログの鍵を作る。★引数なしは `ABC3`（= 木全部）。★純関数。 */
export function logKey(targets) {
  const t = (targets || []).filter(x => x && !x.startsWith('-'));
  const k = (t.length ? t : ['ABC3']).join('+');
  return k.replace(/[^A-Za-z0-9_.+-]/g, '_');
}
export function logPath(targets, cache = CACHE) {
  return path.join(cache, `build-${logKey(targets)}.log`);
}

/**
 * `lake build` の出力を読む。★純関数（文字列 → 要約）。
 * ★`declaration uses 'sorry'` は **warning** で出る。★`grep sorry` では docstring を拾ってしまう。
 */
export function parseLog(text) {
  const s = String(text || '');
  const lines = s.split(/\r?\n/);
  const errors = [];      // { file, line, col, msg }
  const sorries = [];     // { file, line, col, decl }
  const warns = [];
  let jobs = null, ok = null, built = 0, failed = [];
  const LOC = /^(error|warning):\s+([A-Za-z0-9_.\/\\-]+\.lean):(\d+):(\d+):\s*(.*)$/;
  const SORRY = /declaration uses ['`‘]sorry['`’]/;
  for (let i = 0; i < lines.length; i++) {
    const ln = lines[i];
    const m = ln.match(LOC);
    if (m) {
      const rec = { file: m[2].replace(/\\/g, '/'), line: +m[3], col: +m[4], msg: m[5] };
      if (m[1] === 'error') errors.push(rec);
      else {
        warns.push(rec);
        // ★`declaration uses 'sorry'` は **その warning 自身の本文** に出る。
        // ★★次の行を覗いてはいけない —— 直後が別の warning のとき、手前の
        //   `unused variable` まで sorry に数えてしまう(第 29 回が偽の lake で踏んだ)。
        //   次行を見るのは「本文が空で、次行が新しい診断の見出しでない」ときだけ。
        const nxt = lines[i + 1] || '';
        const body = m[5].trim() !== '' ? m[5] : (LOC.test(nxt) ? '' : nxt);
        if (SORRY.test(body)) sorries.push({ ...rec, decl: body.trim() });
      }
      continue;
    }
    let j;
    if ((j = ln.match(/Build completed successfully\s*\((\d+)\s*jobs?\)/))) { jobs = +j[1]; ok = true; }
    else if (/^error:.*build failed/i.test(ln) || /\bBuild failed\b/.test(ln)) ok = false;
    else if ((j = ln.match(/^✖\s*\[(\d+)\/(\d+)\]\s*Built\s+(\S+)/))) { failed.push(j[3]); jobs = +j[2]; ok = false; }
    else if ((j = ln.match(/^✔\s*\[(\d+)\/(\d+)\]/))) { built++; jobs = jobs ?? +j[2]; }
  }
  const files = [...new Set(errors.map(e => e.file))];
  return { ok, jobs, built, errors, warns, sorries, failedTargets: [...new Set(failed)], files, lines: lines.length };
}

/** ログより新しい `.lean` があるか。★無ければ建て直す必要がない（`--if-stale`）。★純関数。 */
export function isStale(logMtimeMs, fileMtimes) {
  if (!(logMtimeMs > 0)) return { stale: true, why: 'ログが無い', newest: null };
  let newest = null, nm = 0;
  for (const [f, t] of fileMtimes) if (t > nm) { nm = t; newest = f; }
  if (nm > logMtimeMs) return { stale: true, why: `${newest} がログより新しい`, newest };
  return { stale: false, why: 'ログより新しい .lean が無い', newest };
}

/** ログの中身を切り出す。★純関数。 */
export function slice(text, { grep = null, tail = 0, context = 0 } = {}) {
  let lines = String(text || '').split(/\r?\n/);
  if (grep) {
    const re = new RegExp(grep);
    const keep = new Set();
    for (let i = 0; i < lines.length; i++) if (re.test(lines[i]))
      for (let d = -context; d <= context; d++) if (i + d >= 0 && i + d < lines.length) keep.add(i + d);
    lines = [...keep].sort((a, b) => a - b).map(i => lines[i]);
  }
  if (tail > 0) lines = lines.slice(-tail);
  return lines.join('\n');
}

// ────────────────────────────────────────────────────────────
// 実行
// ────────────────────────────────────────────────────────────

function lakeBin() {
  if (process.env.ABC3_LAKE) return process.env.ABC3_LAKE;
  const elan = path.join(os.homedir(), '.elan', 'bin', 'lake.exe');
  if (fs.existsSync(elan)) return elan;
  return 'lake';
}

function leanMtimes(root = path.join(LEAN, 'ABC3')) {
  const out = [];
  if (!fs.existsSync(root)) return out;
  (function walk(d) {
    for (const e of fs.readdirSync(d, { withFileTypes: true })) {
      const p = path.join(d, e.name);
      if (e.isDirectory()) walk(p);
      else if (e.name.endsWith('.lean')) { try { out.push([p, fs.statSync(p).mtimeMs]); } catch { /* 消えた */ } }
    }
  })(root);
  const extra = path.join(LEAN, 'ABC3.lean');
  if (fs.existsSync(extra)) out.push([extra, fs.statSync(extra).mtimeMs]);
  return out;
}

/**
 * 叩くもの。★`ABC3_LAKE` が `.mjs`/`.js` を指していれば node で走らせる
 * （★試験で lake を差し替えるための口。shell を挟まないので Windows でも同じ経路を通る）。
 */
export function lakeCmd(bin = lakeBin(), exe = process.execPath) {
  return /\.(mjs|cjs|js)$/i.test(bin) ? { cmd: exe, pre: [bin] } : { cmd: bin, pre: [] };
}

function runBuild(targets, { dryRun = false, timeoutS = 1800 } = {}) {
  const { cmd: bin, pre } = lakeCmd();
  const args = [...pre, 'build', ...targets];
  if (dryRun) return { cmd: `${bin} ${args.join(' ')}  (cwd=${LEAN})`, dryRun: true };
  const t0 = Date.now();
  const r = spawnSync(bin, args, { cwd: LEAN, encoding: 'utf8', maxBuffer: 512 * 1024 * 1024, timeout: timeoutS * 1000 });
  const text = (r.stdout || '') + (r.stderr || '');
  return { text, code: r.status, secs: (Date.now() - t0) / 1000, cmd: `${bin} ${args.join(' ')}` };
}

const WHY = `
★木全部 (\`lake build ABC3\`) と 対象指定 (\`lake build ABC3.Found.X\`) の交換 —— 実測（メタ第 29 回、13.1 日）
  ・直前に .lean を 1 本だけ触った後の値段: ★木全部 中央 37.7 秒 (n=1,098) / 対象指定 中央 14.8 秒 (n=1,229)
  ・全部を対象指定にすると 13.2 時間 → 4.5 時間（★節約 8.7 時間）
  ・★その代わり: 木全部が error を報告した 140 件のうち ★25 件 (18%) は
    「触った本 + その import」の外の error だった ⇒ ★対象指定なら**その場では見落とす**。
    行き先の上位は Skeleton.GenEll.EllModuliWitness (5) / ABC3.Found (3) / Found.GenEll.QuotClassExists (3)。
  ・★見落としても「消える」わけではない —— 次の木全部 build で必ず出る。失うのは**気づく早さ**。
  ⇒ ★この道具は既定で対象を絞らない。絞るかどうかは人が決める。
`.trim();

function main() {
  const argv = process.argv.slice(2);
  const has = f => argv.includes(f);
  const val = f => { const i = argv.indexOf(f); return i >= 0 ? argv[i + 1] : null; };
  const json = has('--json');
  if (has('--selftest')) return selftest();
  if (has('--why')) { console.log(WHY); return; }

  const FLAGS = new Set(['--json', '--errors', '--sorry', '--jobs', '--files', '--if-stale', '--dry-run',
                         '--why', '--selftest', '--grep', '--tail', '--context', '--quiet']);
  const targets = [];
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (FLAGS.has(a)) { if (a === '--grep' || a === '--tail' || a === '--context') i++; continue; }
    if (a.startsWith('-')) continue;
    targets.push(a);
  }
  const tg = targets.length ? targets : ['ABC3'];
  const lp = logPath(tg);
  const readOnly = has('--errors') || has('--sorry') || has('--jobs') || has('--files') || has('--grep');

  fs.mkdirSync(CACHE, { recursive: true });

  // ---- 読み出しだけ（★lake を呼ばない） ----
  if (readOnly) {
    if (!fs.existsSync(lp)) {
      console.log(`  ★ログが無い: ${path.relative(ROOT, lp)}\n  先に \`node tools/build.mjs ${tg.join(' ')}\` を 1 回だけ回すこと。`);
      process.exitCode = 2; return;
    }
    const text = fs.readFileSync(lp, 'utf8');
    const p = parseLog(text);
    const age = ((Date.now() - fs.statSync(lp).mtimeMs) / 60000).toFixed(1);
    if (json) { console.log(JSON.stringify({ log: lp, ageMin: +age, ...p, errors: p.errors, sorries: p.sorries })); return; }
    console.log(`  ログ ${path.relative(ROOT, lp)}（${age} 分前 / ${p.lines} 行）`);
    if (has('--jobs')) console.log(`  ${p.ok === true ? 'ok  Build completed' : p.ok === false ? '★失敗' : '?'}  jobs ${p.jobs ?? '?'}  失敗した対象 ${p.failedTargets.join(' ') || '(なし)'}`);
    if (has('--errors')) {
      console.log(`  error ${p.errors.length} 件`);
      for (const e of p.errors.slice(0, +(val('--tail') || 40))) console.log(`    ${e.file}:${e.line}:${e.col}  ${e.msg.slice(0, 140)}`);
    }
    if (has('--files')) { console.log(`  error を出しているファイル ${p.files.length} 本`); for (const f of p.files) console.log(`    ${f}`); }
    if (has('--sorry')) {
      console.log(`  ★declaration uses \`sorry\` ${p.sorries.length} 件（★grep sorry ではなく build の warning を数えている）`);
      for (const s of p.sorries) console.log(`    ${s.file}:${s.line}  ${s.decl.slice(0, 80)}`);
    }
    if (has('--grep')) console.log(slice(text, { grep: val('--grep'), tail: +(val('--tail') || 0), context: +(val('--context') || 0) }));
    return;
  }

  // ---- 建てる ----
  if (has('--if-stale')) {
    const mt = fs.existsSync(lp) ? fs.statSync(lp).mtimeMs : 0;
    const st = isStale(mt, leanMtimes());
    if (!st.stale) {
      const p = parseLog(fs.readFileSync(lp, 'utf8'));
      const age = ((Date.now() - mt) / 60000).toFixed(1);
      if (json) { console.log(JSON.stringify({ rebuilt: false, why: st.why, log: lp, ageMin: +age, ok: p.ok, jobs: p.jobs, errors: p.errors.length, sorries: p.sorries.length })); return; }
      console.log(`  ★建てなかった（${st.why}）—— ${path.relative(ROOT, lp)} は ${age} 分前`);
      console.log(`  ${p.ok === true ? 'ok  Build completed' : p.ok === false ? '★失敗' : '?'}  jobs ${p.jobs ?? '?'}  error ${p.errors.length}  sorry ${p.sorries.length}`);
      return;
    }
  }
  const r = runBuild(tg, { dryRun: has('--dry-run') });
  if (r.dryRun) { console.log(`  叩くもの: ${r.cmd}\n  落とす先: ${path.relative(ROOT, lp)}`); return; }
  fs.writeFileSync(lp, r.text);
  const p = parseLog(r.text);
  if (json) { console.log(JSON.stringify({ rebuilt: true, secs: r.secs, log: lp, ok: p.ok, jobs: p.jobs, errors: p.errors.length, sorries: p.sorries.length, files: p.files, failedTargets: p.failedTargets })); return; }
  console.log(`  ${r.cmd}  —— ${r.secs.toFixed(1)} 秒`);
  console.log(`  ${p.ok === true ? 'ok  Build completed' : p.ok === false ? '★失敗' : '?'}  jobs ${p.jobs ?? '?'}  error ${p.errors.length}  warning ${p.warns.length}  sorry ${p.sorries.length}`);
  for (const e of p.errors.slice(0, 10)) console.log(`    ${e.file}:${e.line}:${e.col}  ${e.msg.slice(0, 140)}`);
  if (p.errors.length > 10) console.log(`    … 残り ${p.errors.length - 10} 件は \`node tools/build.mjs ${tg.join(' ')} --errors\`（★建て直さない）`);
  console.log(`  → ${path.relative(ROOT, lp)}  以後は --errors / --sorry / --jobs / --grep で切り出す（★lake を呼ばない）`);
  if (p.ok === false) process.exitCode = 1;
}

// ────────────────────────────────────────────────────────────
// selftest
// ────────────────────────────────────────────────────────────
function selftest() {
  let ok = 0, ng = 0;
  const T = (name, cond) => { if (cond) ok++; else { ng++; console.log(`  NG ${name}`); } };

  // logKey
  T('logKey 引数なしは ABC3', logKey([]) === 'ABC3');
  T('logKey 単一', logKey(['ABC3.Found.PGC.X']) === 'ABC3.Found.PGC.X');
  T('logKey 複数は + で繋ぐ', logKey(['ABC3', 'ABC3.Found']) === 'ABC3+ABC3.Found');
  T('logKey は旗を無視', logKey(['--json', 'ABC3.A']) === 'ABC3.A');
  T('logKey は危ない字を潰す', logKey(['a/b:c']) === 'a_b_c');

  const LOG = [
    'info: stdout',
    '✔ [10/6905] Built ABC3.Found.A',
    'warning: ABC3/Found/A.lean:12:9: declaration uses \'sorry\'',
    'warning: ABC3/Found/A.lean:30:1: unused variable `x`',
    'error: ABC3/Found/B.lean:5:2: unknown identifier `foo`',
    'error: ABC3/Found/B.lean:9:4: unsolved goals',
    '✖ [6904/6905] Built ABC3.Found.B',
    'error: Lean exited with code 1',
  ].join('\n');
  const p = parseLog(LOG);
  T('parseLog error 2 件', p.errors.length === 2);
  T('parseLog error の在処', p.errors[0].file === 'ABC3/Found/B.lean' && p.errors[0].line === 5);
  T('parseLog sorry 1 件', p.sorries.length === 1);
  T('parseLog sorry の在処', p.sorries[0].file === 'ABC3/Found/A.lean' && p.sorries[0].line === 12);
  T('parseLog ★sorry は warning から取る（error に混ぜない）', !p.errors.some(e => /sorry/.test(e.msg)));
  T('parseLog warning 2 件', p.warns.length === 2);
  T('parseLog 失敗した対象', p.failedTargets.join() === 'ABC3.Found.B');
  T('parseLog jobs', p.jobs === 6905);
  T('parseLog ok=false', p.ok === false);
  T('parseLog error のファイル一覧は重複しない', p.files.length === 1);

  const OKLOG = ['✔ [6905/6905] Built ABC3', 'Build completed successfully (6905 jobs).'].join('\n');
  const q = parseLog(OKLOG);
  T('parseLog 成功', q.ok === true && q.jobs === 6905 && q.errors.length === 0);
  // ★成功行だけでも jobs が取れること（★✔ 行に隠されて素通りしていた。第 29 回の突然変異 #13）
  T('parseLog ★成功行だけで jobs が取れる', parseLog('Build completed successfully (42 jobs).').jobs === 42);
  // ★行の途中の "error:" を診断と間違えないこと（★ゴール表示に混ざる。第 29 回の突然変異 #8）
  const MID = ['error: ABC3/A.lean:1:1: unsolved goals',
               '  h : foo error: ABC3/B.lean:9:9: これは本文であって診断ではない'].join('\n');
  T('parseLog ★行の途中の error: は診断でない', parseLog(MID).errors.length === 1);
  T('parseLog ★行の途中の在処を拾わない', parseLog(MID).files.join() === 'ABC3/A.lean');
  T('parseLog 空でも落ちない', parseLog('').ok === null);
  T('parseLog CRLF でも読める', parseLog(LOG.replace(/\n/g, '\r\n')).errors.length === 2);
  // ★docstring の中の "sorry" を拾わないこと（本体の 134 件はこれ）
  const DOC = 'info: -- sorry と書いてあるだけの docstring\nBuild completed successfully (3 jobs).';
  T('parseLog ★docstring の sorry を数えない', parseLog(DOC).sorries.length === 0);
  // ★★隣の warning に引きずられないこと（偽の lake で実際に踏んだ。14 件 → 7 件）
  const ADJ = [
    "warning: ABC3/Check/X.lean:99:10: unused variable `x`",
    "warning: ABC3/Found/Y.lean:12:2: declaration uses 'sorry'",
    "warning: ABC3/Check/Z.lean:5:1: unused variable `y`",
  ].join('\n');
  const a = parseLog(ADJ);
  T('parseLog ★隣の warning を sorry に数えない', a.sorries.length === 1);
  T('parseLog ★sorry の在処は自分の行', a.sorries[0].file === 'ABC3/Found/Y.lean' && a.sorries[0].line === 12);
  // ★本文が空で次行に続く形（Lean が折り返した場合）は拾う
  const WRAP = "warning: ABC3/Found/W.lean:3:1:\n  declaration uses 'sorry'";
  T('parseLog 折り返した sorry は拾う', parseLog(WRAP).sorries.length === 1);
  T('parseLog ★バッククォート版も拾う', parseLog('warning: ABC3/A.lean:1:1: declaration uses `sorry`').sorries.length === 1);

  // isStale
  T('isStale ログが無い', isStale(0, []).stale === true);
  T('isStale 新しい .lean がある', isStale(100, [['a.lean', 200]]).stale === true);
  T('isStale 新しい .lean が無い', isStale(300, [['a.lean', 200], ['b.lean', 100]]).stale === false);
  T('isStale 空の木', isStale(300, []).stale === false);
  T('isStale ★同時刻は建て直さない', isStale(200, [['a.lean', 200]]).stale === false);

  // slice
  const S = 'a1\nb2\nc3\nd4\ne5';
  T('slice grep', slice(S, { grep: 'b|d' }) === 'b2\nd4');
  T('slice tail', slice(S, { tail: 2 }) === 'd4\ne5');
  T('slice context', slice(S, { grep: 'c3', context: 1 }) === 'b2\nc3\nd4');
  T('slice 何もしなければそのまま', slice(S) === S);
  T('slice ★context が重なっても行が重複しない', slice('x\nh\nh\ny', { grep: 'h', context: 1 }) === 'x\nh\nh\ny');

  // lakeCmd（★試験の差し替えが本番の経路を変えないこと）
  T('lakeCmd 実体はそのまま', lakeCmd('C:/x/lake.exe', 'node').cmd === 'C:/x/lake.exe');
  T('lakeCmd 実体は pre が空', lakeCmd('C:/x/lake.exe', 'node').pre.length === 0);
  T('lakeCmd .mjs は node で走る', lakeCmd('C:/x/fake.mjs', 'node').cmd === 'node');
  T('lakeCmd .mjs は pre に入る', lakeCmd('C:/x/fake.mjs', 'node').pre[0] === 'C:/x/fake.mjs');

  console.log(`  build.mjs selftest ${ok}/${ok + ng}`);
  if (ng) process.exitCode = 1;
}

if (process.argv[1] && process.argv[1].replace(/\\/g, '/').endsWith('build.mjs')) main();
