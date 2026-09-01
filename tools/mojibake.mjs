#!/usr/bin/env node
// UTF-8 のバイト列を latin1 として読んでしまった文字化け(mojibake)を検出・修復する。
//
//   node tools/mojibake.mjs           走査のみ(終了コード 1 = 検出あり)
//   node tools/mojibake.mjs --fix     復元して書き戻す
//
// 機構: 化けた行は「latin1 で符号化 → UTF-8 で復号」が成功して別の文字列になる。
// 正常な日本語行は latin1 で符号化できない(U+0080 以上の非 latin1 文字を含む)ので素通りする。
import fs from 'fs';
import path from 'path';

const ROOT = 'D:/Math_ABC3';
const TARGETS = ['lean', 'tools', 'ResearchPaper'];
const EXT = new Set(['.lean', '.mjs', '.md', '.json']);

function* walk(dir) {
  for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
    if (e.name === '.lake' || e.name === 'node_modules' || e.name === '.cache') continue;
    const p = path.join(dir, e.name);
    if (e.isDirectory()) yield* walk(p);
    else if (EXT.has(path.extname(e.name))) yield p;
  }
}

// 行が化けているなら復元後の文字列を、そうでなければ null を返す。
function recover(line) {
  if (![...line].some((c) => c.codePointAt(0) > 127)) return null;
  let bytes;
  try {
    bytes = Buffer.from(line, 'latin1');
    // latin1 は 0-255 しか表せない。範囲外があれば化けていない。
    if ([...line].some((c) => c.codePointAt(0) > 255)) return null;
  } catch {
    return null;
  }
  const dec = new TextDecoder('utf-8', { fatal: true });
  let out;
  try {
    out = dec.decode(bytes);
  } catch {
    return null;
  }
  return out === line ? null : out;
}

const fix = process.argv.includes('--fix');
let files = 0;
let lines = 0;

for (const t of TARGETS) {
  const dir = path.join(ROOT, t);
  if (!fs.existsSync(dir)) continue;
  for (const f of walk(dir)) {
    const src = fs.readFileSync(f, 'utf8');
    const arr = src.split('\n');
    let hit = [];
    for (let i = 0; i < arr.length; i++) {
      const r = recover(arr[i]);
      if (r !== null) {
        hit.push(i + 1);
        if (fix) arr[i] = r;
      }
    }
    if (hit.length) {
      files++;
      lines += hit.length;
      const rel = path.relative(ROOT, f).split(path.sep).join('/');
      console.log(`${rel}  行 ${hit.join(', ')}`);
      if (fix) fs.writeFileSync(f, arr.join('\n'), 'utf8');
    }
  }
}

if (files === 0) {
  console.log('ok  文字化けなし');
  process.exit(0);
}
console.log(`${fix ? '復元した' : '検出した'}: ${files} ファイル / ${lines} 行`);
process.exit(fix ? 0 : 1);
