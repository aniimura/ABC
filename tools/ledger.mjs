#!/usr/bin/env node
// 台帳を 1 本の CLI で更新する。
//
// これまで台帳の 1 行を書くたびに tools/_ledger-NNNN.mjs のような
// 使い捨てスクリプトを作っていた(2026-09-03 の実測で 397 本 / 10,528 行)。
// 同じことを毎回書き直す理由は無いので、ここに集約する。
//
//   node tools/ledger.mjs blocked --key "<key>" --field <名前> --body-file <file>
//   node tools/ledger.mjs gap     --field <名前> --body-file <file>
//   node tools/ledger.mjs goal    --body-file <file>
//   node tools/ledger.mjs show    blocked --key "<key>"
//   node tools/ledger.mjs keys    blocked
//
// 本文は --body-file(UTF-8)か --body で渡す。★シェルの引用を通さずに済むので
// 「`」や「★」や改行を含む日本語がそのまま書ける。
// --dry-run を付けると書かずに差分だけ出す。

import fs from 'fs';
import path from 'path';

// ★LEDGER_RP で差し替えられる(試験のときに複製へ向けるため)。
const RP = process.env.LEDGER_RP || 'D:/Math_ABC3/ResearchPaper';
const TARGETS = {
  blocked: { file: `${RP}/blocked-leaves.json`, kind: 'keyed', array: 'blocked', idField: 'key' },
  gap:     { file: `${RP}/mathlib-gap.json`,    kind: 'flat' },
  tree:    { file: `${RP}/obligation-tree.json`, kind: 'keyed', array: 'obligations', idField: 'name' },
  goal:    { file: `${RP}/genell-goal.md`,      kind: 'markdown' },
};

function usage(msg) {
  if (msg) console.error(`ledger: ${msg}\n`);
  console.error(`使い方:
  node tools/ledger.mjs blocked --key "<key>" --field <名前> (--body-file <file> | --body "<text>")
  node tools/ledger.mjs gap                   --field <名前> (--body-file <file> | --body "<text>")
  node tools/ledger.mjs tree    --key "<key>" --field <名前> (--body-file <file> | --body "<text>")
  node tools/ledger.mjs goal                                 (--body-file <file> | --body "<text>")
  node tools/ledger.mjs show <blocked|gap|tree> [--key "<key>"] [--field <名前>]
  node tools/ledger.mjs keys <blocked|gap|tree>

  --dry-run   書かずに何をするかだけ出す
  --replace   既にその欄がある場合、既定は拒否。上書きするならこれを付ける`);
  process.exit(msg ? 2 : 0);
}

function parseArgs(argv) {
  const out = { _: [] };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a === '--dry-run') out.dryRun = true;
    else if (a === '--replace') out.replace = true;
    else if (a === '--reformat') out.reformat = true;
    else if (a.startsWith('--')) out[a.slice(2)] = argv[++i];
    else out._.push(a);
  }
  return out;
}

function readBody(args) {
  if (args['body-file']) {
    const p = path.resolve(args['body-file']);
    if (!fs.existsSync(p)) usage(`--body-file が無い: ${p}`);
    return fs.readFileSync(p, 'utf8').replace(/\r\n/g, '\n').replace(/\n+$/, '');
  }
  if (typeof args.body === 'string') return args.body;
  usage('--body-file か --body が要る');
}

function loadJson(t) {
  const raw = fs.readFileSync(t.file, 'utf8');
  // ★書き戻したときに**関係ない行まで整形が変わる**台帳がある
  //   (2026-09-03 実測: obligation-tree.json は 6,482 → 6,982 文字に変わる)。
  //   差分が読めなくなるので、既定では拒否する。
  t._roundTrips = raw === JSON.stringify(JSON.parse(raw), null, 2) + '\n';
  return JSON.parse(raw);
}

function saveJson(t, j, args) {
  if (!t._roundTrips && !args.reformat) {
    console.error(
      `${path.basename(t.file)} は現在 JSON.stringify(...,null,2) の整形と一致していない。\n` +
      `そのまま書くと**触っていない行まで差分に出る**。承知のうえなら --reformat を付ける`);
    process.exit(1);
  }
  const text = JSON.stringify(j, null, 2) + '\n';
  if (args.dryRun) { console.log(`(dry-run) ${t.file} へ ${text.length} 文字書く`); return; }
  fs.writeFileSync(t.file, text, 'utf8');
  console.log(`ok  ${path.basename(t.file)} を更新した (${text.length} 文字)`);
}

const args = parseArgs(process.argv.slice(2));
const cmd = args._[0];
if (!cmd) usage();

// ── show / keys ────────────────────────────────────────────────────────────
if (cmd === 'show' || cmd === 'keys') {
  const name = args._[1];
  const t = TARGETS[name];
  if (!t || t.kind === 'markdown') usage(`show/keys は blocked|gap|tree のみ`);
  const j = loadJson(t);
  if (cmd === 'keys') {
    const ks = t.kind === 'keyed' ? j[t.array].map((e) => e[t.idField]) : Object.keys(j);
    ks.forEach((k) => console.log(k));
    process.exit(0);
  }
  let obj;
  if (t.kind === 'keyed') {
    if (!args.key) usage('show <keyed> には --key が要る');
    obj = j[t.array].find((e) => e[t.idField] === args.key);
    if (!obj) { console.error(`該当なし: ${args.key}`); process.exit(1); }
  } else {
    obj = j;
  }
  if (args.field) console.log(obj[args.field] ?? '(その欄は無い)');
  else Object.keys(obj).forEach((k) => console.log(`${k}\t(${String(obj[k]).length} 文字)`));
  process.exit(0);
}

// ── 追記 ───────────────────────────────────────────────────────────────────
const t = TARGETS[cmd];
if (!t) usage(`知らない台帳: ${cmd}`);
const body = readBody(args);
if (!body.trim()) usage('本文が空である');

if (t.kind === 'markdown') {
  if (args.dryRun) { console.log(`(dry-run) ${t.file} へ ${body.length} 文字を追記する`); process.exit(0); }
  const cur = fs.readFileSync(t.file, 'utf8');
  const sep = cur.endsWith('\n') ? '\n' : '\n\n';
  fs.writeFileSync(t.file, cur + sep + body + '\n', 'utf8');
  console.log(`ok  ${path.basename(t.file)} へ ${body.length} 文字を追記した`);
  process.exit(0);
}

if (!args.field) usage('--field が要る');
const j = loadJson(t);

let target;
if (t.kind === 'keyed') {
  if (!args.key) usage('--key が要る');
  target = j[t.array].find((e) => e[t.idField] === args.key);
  if (!target) { target = { [t.idField]: args.key }; j[t.array].push(target); console.log(`新しい葉を作った: ${args.key}`); }
} else {
  target = j;
}

if (target[args.field] !== undefined && !args.replace) {
  console.error(`欄 ${args.field} は既にある (${String(target[args.field]).length} 文字)。上書きするなら --replace`);
  process.exit(1);
}
target[args.field] = body;
saveJson(t, j, args);
