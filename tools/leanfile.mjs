#!/usr/bin/env node
// 単一 Lean ファイルの検査。**`lake build` と同じオプション**で走らせる。
//
// ★★★なぜ要るか(2026-08-17 実測)
//
// `lake env lean FILE` は `LEAN_PATH` を張るだけで、
// **`lakefile.toml` の `[leanOptions]` を引き継がない**。
// そのため
//
//     structure Foo where
//       Point : Type u        -- `universe u` を書いていない
//
// が `lake env lean` では**通り**、`lake build` では
// `unknown universe level u` で**落ちた**。
//
// ★食い違いの実例はこれで 3 種目である
// (セクション変数の自動包含 / auto-bound universe / `cd` の欠落)。
// 共通点は「`lake env lean` の方が緩い」。
//
// ★★本スクリプトは `lakefile.toml` を**読んで**オプションを渡すので、
// lakefile を変えても自動で追従する(写経による drift が起きない)。
//
// ★★★★2026-09-08: **ログをキャッシュし、既定では要点だけ出す**ようにした。
//
// 理由(実測): 本ファイルは 86 行でキャッシュを持たず、**出力が丸ごと呼び出し側の文脈に入っていた**。
// 実装 1 体あたり 7〜14 往復あり、エラー 1 件が長い(型不一致は `inst✝` 込みで 8 行)。
// ★圧縮が起きると「どのタクティクをなぜ失敗したか」が消え、**同じタクティクを再試行する**。
// ⇒ 全文は `.cache/leanfile-<モジュール>.log` に置き、既定は
//   **エラー/警告の行とその直後の文脈だけ**を上限つきで出す。
//   後から `--errors` / `--grep` / `--full` で切り出せる(★`lean` を呼び直さない)。
//
// 使い方(`lean/` から):
//     node ../tools/leanfile.mjs ABC3/Found/GenEll/Foo.lean
// あるいはリポジトリ根から:
//     node tools/leanfile.mjs lean/ABC3/Found/GenEll/Foo.lean
//
// 切り出し(★lean を呼ばない。直前の実行のログを読むだけ):
//     node tools/leanfile.mjs --errors lean/ABC3/Found/GenEll/Foo.lean
//     node tools/leanfile.mjs --grep 'unknown identifier' lean/ABC3/Found/GenEll/Foo.lean
//     node tools/leanfile.mjs --full lean/ABC3/Found/GenEll/Foo.lean
//
// ★`lake build` と違い **olean を書かない**ので、
// 並行セッションと同じワークツリーを共有していても安全である。

import { readFileSync, writeFileSync, existsSync, mkdirSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { dirname, join, relative, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const LEAN_DIR = join(ROOT, 'lean');
const LAKEFILE = join(LEAN_DIR, 'lakefile.toml');
const CACHE = join(ROOT, '.cache');

/** 既定で出す最大行数。これを超えたぶんはログにだけ残る。 */
export const HEAD_LINES = 60;

/** `lakefile.toml` の `[leanOptions]` を読む。 */
export function readLeanOptions(path = LAKEFILE) {
  if (!existsSync(path)) return [];
  const src = readFileSync(path, 'utf8');
  const out = [];
  let inSection = false;
  for (const raw of src.split('\n')) {
    const line = raw.trim();
    if (line.startsWith('[')) { inSection = line === '[leanOptions]'; continue; }
    if (!inSection) continue;
    if (!line || line.startsWith('#')) continue;
    const m = /^([\w.]+)\s*=\s*(.+?)\s*$/.exec(line);
    if (m) out.push(`${m[1]}=${m[2].replace(/^"(.*)"$/, '$1')}`);
  }
  return out;
}

/** `lean/` からの相対パスを、ログ名に使える安全な名前にする。 */
export function logNameFor(rel) {
  return `leanfile-${rel.replace(/\.lean$/, '').replace(/[\\/]/g, '.')}.log`;
}

/**
 * Lean の診断の「見出し行」か。
 * ★Lean 4 は書式が 2 通りある(`FILE:L:C: error: ...` と `error: FILE:L:C: ...`)。
 * ★どちらも拾わないと取り落とす(2026-09-07 に実測、メタ第 37 回)。
 */
export function isDiagHead(line) {
  return /^\S.*?:\d+:\d+:\s*(error|warning)\b/.test(line) ||
         /^(error|warning)\b.*?:\s*\S.*?:\d+:\d+:/.test(line) ||
         /^(error|warning)(\(|:)/.test(line);
}

/**
 * 全文から「診断の行とその直後の文脈」だけを抜く。
 * 見出しから次の見出し(または空行 2 連)までを 1 ブロックとし、
 * 1 ブロックあたり `ctx` 行までに切る。
 */
export function extractDiagnostics(text, { ctx = 12 } = {}) {
  const lines = text.split('\n');
  const out = [];
  let i = 0;
  while (i < lines.length) {
    if (!isDiagHead(lines[i])) { i++; continue; }
    const block = [lines[i]];
    let j = i + 1;
    let blank = 0;
    while (j < lines.length && block.length < ctx) {
      if (isDiagHead(lines[j])) break;
      if (lines[j].trim() === '') { blank++; if (blank >= 2) break; } else blank = 0;
      block.push(lines[j]);
      j++;
    }
    while (block.length && block[block.length - 1].trim() === '') block.pop();
    out.push(block.join('\n'));
    i = j;
  }
  return out;
}

/** 出力を上限つきで出し、切ったぶんの案内を返す。 */
export function printCapped(blocks, cap = HEAD_LINES) {
  let used = 0;
  let shown = 0;
  for (const b of blocks) {
    const n = b.split('\n').length;
    if (used + n > cap) break;
    console.log(b);
    used += n;
    shown++;
  }
  return { shown, total: blocks.length, hidden: blocks.length - shown };
}

// ---------------------------------------------------------------- selftest

function selftest() {
  let pass = 0, fail = 0;
  const eq = (name, got, want) => {
    const ok = JSON.stringify(got) === JSON.stringify(want);
    if (ok) pass++; else { fail++; console.log(`  NG ${name}\n     got  ${JSON.stringify(got)}\n     want ${JSON.stringify(want)}`); }
  };
  const t = (name, cond) => { if (cond) pass++; else { fail++; console.log(`  NG ${name}`); } };

  // logNameFor
  eq('logNameFor 1', logNameFor('ABC3/Found/PGC/Foo.lean'), 'leanfile-ABC3.Found.PGC.Foo.log');
  eq('logNameFor 2', logNameFor('ABC3\\Found\\Bar.lean'), 'leanfile-ABC3.Found.Bar.log');
  eq('logNameFor 3', logNameFor('Foo.lean'), 'leanfile-Foo.log');

  // isDiagHead —— ★2 書式 + 括弧つき
  t('diag 書式A', isDiagHead('ABC3/Found/Foo.lean:12:4: error: unknown identifier'));
  t('diag 書式B', isDiagHead('error: ABC3/Found/Foo.lean:12:4: import failed'));
  t('diag 括弧', isDiagHead('error(lean.unknownIdentifier): Unknown constant `Foo`'));
  t('diag warning', isDiagHead('ABC3/Foo.lean:3:0: warning: declaration uses sorry'));
  t('diag 非診断', !isDiagHead('  have h : x = y := by'));
  t('diag 空行', !isDiagHead(''));
  t('diag ok 行', !isDiagHead('ok  ABC3/Found/Foo.lean  (オプション: pp.unicode.fun=true)'));

  // extractDiagnostics
  const sample = [
    'building ABC3.Found',
    'ABC3/Foo.lean:1:0: error: unknown identifier',
    '  context line 1',
    '  context line 2',
    '',
    '',
    'unrelated tail',
    'ABC3/Foo.lean:9:2: warning: unused variable',
    '  ctx',
  ].join('\n');
  const blocks = extractDiagnostics(sample);
  eq('extract 個数', blocks.length, 2);
  t('extract 1 本目に見出し', blocks[0].startsWith('ABC3/Foo.lean:1:0: error'));
  t('extract 1 本目に文脈', blocks[0].includes('context line 2'));
  t('extract 1 本目に無関係な行が入らない', !blocks[0].includes('unrelated tail'));
  t('extract 2 本目', blocks[1].includes('unused variable'));
  eq('extract 空入力', extractDiagnostics(''), []);
  eq('extract 診断なし', extractDiagnostics('all good\nok  X'), []);

  // ★見出しが連続する場合、ブロックが融合しないこと
  const back2back = 'A.lean:1:1: error: a\nB.lean:2:2: error: b';
  eq('extract 連続見出し', extractDiagnostics(back2back).length, 2);

  // ctx 上限
  const long = 'A.lean:1:1: error: a\n' + Array.from({ length: 40 }, (_, i) => `  l${i}`).join('\n');
  t('extract ctx 上限', extractDiagnostics(long, { ctx: 5 })[0].split('\n').length <= 5);

  // readLeanOptions —— 実ファイルで
  const opts = readLeanOptions();
  t('lakefile の leanOptions が読める', Array.isArray(opts) && opts.length > 0);

  console.log(`selftest ${pass}/${pass + fail}`);
  return fail === 0 ? 0 : 1;
}

// ---------------------------------------------------------------- main

const argv = process.argv.slice(2);
if (argv.includes('--selftest')) process.exit(selftest());

const wantFull = argv.includes('--full');
const wantErrors = argv.includes('--errors');
const grepIdx = argv.indexOf('--grep');
const grepRe = grepIdx >= 0 && argv[grepIdx + 1] ? new RegExp(argv[grepIdx + 1], 'i') : null;
const files = argv.filter((a, i) =>
  !a.startsWith('--') && !(grepIdx >= 0 && i === grepIdx + 1));

if (files.length === 0) {
  console.error('使い方: node tools/leanfile.mjs <FILE.lean> [...]');
  console.error('  ★lakefile.toml の [leanOptions] を渡すので `lake build` と食い違わない。');
  console.error('  --errors / --grep <re> / --full  直前のログから切り出す(★lean を呼ばない)');
  console.error('  --selftest');
  process.exit(2);
}

const opts = readLeanOptions();
if (opts.length === 0) {
  console.error(`★lakefile.toml の [leanOptions] を読めなかった: ${LAKEFILE}`);
  process.exit(2);
}
if (!existsSync(CACHE)) mkdirSync(CACHE, { recursive: true });

let bad = 0;
for (const a of files) {
  const abs = resolve(a);
  if (!existsSync(abs)) { console.error(`NG  ファイルが無い: ${a}`); bad++; continue; }
  const rel = relative(LEAN_DIR, abs).replace(/\\/g, '/');
  const logPath = join(CACHE, logNameFor(rel));

  // --- 切り出しモード: lean を呼ばずにログを読む ---
  if (wantErrors || grepRe) {
    if (!existsSync(logPath)) {
      console.error(`★ログが無い: ${logPath}  —— 先にオプション無しで 1 度走らせること`);
      bad++; continue;
    }
    const text = readFileSync(logPath, 'utf8');
    if (grepRe) {
      const hits = text.split('\n').filter((l) => grepRe.test(l));
      console.log(hits.length ? hits.join('\n') : '(該当なし)');
    } else {
      const blocks = extractDiagnostics(text);
      const { hidden, total } = printCapped(blocks, Number.MAX_SAFE_INTEGER);
      console.log(`  診断 ${total} 件  ${hidden ? `(${hidden} 件は省略)` : ''}`);
    }
    continue;
  }

  // --- 実行 ---
  const args2 = ['env', 'lean', ...opts.map((o) => `-D${o}`), rel];
  const t0 = Date.now();
  const r = spawnSync('lake', args2, { cwd: LEAN_DIR, encoding: 'utf8', shell: true });
  const out = `${r.stdout ?? ''}${r.stderr ?? ''}`.trim();
  const secs = ((Date.now() - t0) / 1000).toFixed(1);

  writeFileSync(logPath, out + '\n', 'utf8');

  const isNG = r.status !== 0 || /^\S+:\d+:\d+: error/m.test(out) || /^error(\(|:)/m.test(out);

  if (wantFull) {
    if (out) console.log(out);
  } else {
    const blocks = extractDiagnostics(out);
    if (blocks.length) {
      const { hidden, total } = printCapped(blocks, HEAD_LINES);
      if (hidden) console.log(`  … 診断 ${total} 件中 ${hidden} 件は省略(全文はログ)`);
    }
  }

  if (isNG) { console.log(`NG  ${rel}  —— ${secs} 秒`); bad++; }
  else console.log(`ok  ${rel}  —— ${secs} 秒  (オプション: ${opts.join(' ')})`);
  console.log(`  → ${relative(ROOT, logPath).replace(/\\/g, '/')}  以後は --errors / --grep <re> / --full で切り出す（★lean を呼ばない）`);
}
process.exit(bad ? 1 : 0);
