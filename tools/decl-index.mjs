#!/usr/bin/env node
/**
 * 宣言インデックス —— 在庫検索の対策(2026-08-21)
 *
 * ★動機(実測): 宣言が 9,937 個になり、「使える補題がすでにあるか」を探すのに
 * 木全体(198,000 行)を `grep` する往復が増えていた。
 * ★1 本のテキストにしておけば **1 回の `grep` で済む**。
 *
 * 出力(既定 `.cache/`、gitignore 済み):
 *   `decl-index.txt`  1 行 1 宣言 —— `種別 <TAB> 名前 <TAB> ファイル:行`
 *   `src-index.txt`   1 行 1 locator —— `論文 pPDF頁 <TAB> item <TAB> ファイル:行 <TAB> 宣言`
 *
 * 使い方:
 *   node tools/decl-index.mjs
 *   grep -i "rootMap" .cache/decl-index.txt
 *   grep "Proposition 5.5" .cache/src-index.txt
 */

import { readFileSync, readdirSync, statSync, mkdirSync, writeFileSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const LEAN_SRC = join(ROOT, 'lean', 'ABC3');
const OUT_DIR = join(ROOT, '.cache');

function walk(dir, acc = []) {
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) walk(p, acc);
    else if (p.endsWith('.lean')) acc.push(p);
  }
  return acc;
}

const DECL_RE =
  /^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*(theorem|lemma|def|abbrev|instance|structure|inductive|class)\s+([A-Za-z_][\w'!?₀-₉]*(?:\.[\w'!?₀-₉]+)*)/;
const NS_RE = /^\s*namespace\s+([A-Za-z_][\w'.₀-₉]*)/;
const END_RE = /^\s*end\s+([A-Za-z_][\w'.₀-₉]*)\s*$/;

const decls = [];
const srcs = [];

for (const file of walk(LEAN_SRC)) {
  const rel = relative(LEAN_SRC, file).replace(/\\/g, '/');
  const lines = readFileSync(file, 'utf8').split(/\r?\n/);
  const nsStack = [];
  let pendingSrcDecl = null;
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    const ns = NS_RE.exec(line);
    if (ns) { nsStack.push(ns[1]); continue; }
    const en = END_RE.exec(line);
    if (en && nsStack.length && nsStack[nsStack.length - 1] === en[1]) { nsStack.pop(); continue; }
    const m = DECL_RE.exec(line);
    if (!m) continue;
    const full = [...nsStack, m[2]].join('.');
    decls.push(`${m[1].padEnd(9)}\t${full}\t${rel}:${i + 1}`);
    if (m[2].endsWith('.src')) pendingSrcDecl = { name: full, rel, line: i + 1, at: i };
    // `.src` の中身(paper / pdfPage / item)を後続 6 行から拾う
    if (pendingSrcDecl && i - pendingSrcDecl.at <= 6) {
      const blob = lines.slice(pendingSrcDecl.at, Math.min(pendingSrcDecl.at + 7, lines.length)).join(' ');
      const paper = /paper\s*:=\s*"([^"]*)"/.exec(blob);
      const page = /pdfPage\s*:=\s*(\d+)/.exec(blob);
      const item = /item\s*:=\s*"([^"]*)"/.exec(blob);
      if (paper && item) {
        srcs.push(`${paper[1]} p${page ? page[1] : '?'}\t${item[1]}\t${pendingSrcDecl.rel}:${pendingSrcDecl.line}\t${pendingSrcDecl.name}`);
        pendingSrcDecl = null;
      }
    }
  }
}

decls.sort();
srcs.sort();
mkdirSync(OUT_DIR, { recursive: true });
writeFileSync(join(OUT_DIR, 'decl-index.txt'), decls.join('\n') + '\n', 'utf8');
writeFileSync(join(OUT_DIR, 'src-index.txt'), srcs.join('\n') + '\n', 'utf8');
console.log(`decl-index.txt: ${decls.length} 宣言`);
console.log(`src-index.txt:  ${srcs.length} locator`);
