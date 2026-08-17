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
// 使い方(`lean/` から):
//     node ../tools/leanfile.mjs ABC3/Found/GenEll/Foo.lean
// あるいはリポジトリ根から:
//     node tools/leanfile.mjs lean/ABC3/Found/GenEll/Foo.lean
//
// ★`lake build` と違い **olean を書かない**ので、
// 並行セッションと同じワークツリーを共有していても安全である。

import { readFileSync, existsSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { dirname, join, relative, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const LEAN_DIR = join(ROOT, 'lean');
const LAKEFILE = join(LEAN_DIR, 'lakefile.toml');

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

const args = process.argv.slice(2);
if (args.length === 0) {
  console.error('使い方: node tools/leanfile.mjs <FILE.lean> [...]');
  console.error('  ★lakefile.toml の [leanOptions] を渡すので `lake build` と食い違わない。');
  process.exit(2);
}

const opts = readLeanOptions();
if (opts.length === 0) {
  console.error(`★lakefile.toml の [leanOptions] を読めなかった: ${LAKEFILE}`);
  process.exit(2);
}

let bad = 0;
for (const a of args) {
  const abs = resolve(a);
  if (!existsSync(abs)) { console.error(`NG  ファイルが無い: ${a}`); bad++; continue; }
  const rel = relative(LEAN_DIR, abs).replace(/\\/g, '/');
  const argv = ['env', 'lean', ...opts.map((o) => `-D${o}`), rel];
  const r = spawnSync('lake', argv, { cwd: LEAN_DIR, encoding: 'utf8', shell: true });
  const out = `${r.stdout ?? ''}${r.stderr ?? ''}`.trim();
  if (out) console.log(out);
  if (r.status !== 0 || /^\S+:\d+:\d+: error/m.test(out)) {
    console.log(`NG  ${rel}`);
    bad++;
  } else {
    console.log(`ok  ${rel}  (オプション: ${opts.join(' ')})`);
  }
}
process.exit(bad ? 1 : 0);
