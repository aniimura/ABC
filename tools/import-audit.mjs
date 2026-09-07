#!/usr/bin/env node
// import-audit.mjs —— 「ビルドが通っている」を「依存が正しい」と取り違える失敗を止めるための観測点。
//                     ★測るだけ。1 行も書き換えない(メタ第 28 回、持ち場 3)。
//
// ★何を見つけたくて作ったか(実際に 2 度起きた失敗形)
//   2026-09-07 に本体は同じ失敗を 2 度踏んだ ——
//   「`Skeleton` の**定義**が定理ファイルに置かれていて、`Found` がそれを引くために
//    その定理ファイルを import し、配線した瞬間に **import 循環** になる」。
//   ★1 度目 `Skeleton/PGC/Section1` → `Section1Defs` / ★2 度目 `Section2` → `Section2Defs`。
//   ⇒ ★同じ形の本が他にもあるなら、**踏む前に**名指しできるはずである。
//
// ★★循環するかどうかは「引く側ごと」に決まる(ここが肝)
//   F が S を import して循環するのは **F ∈ reach(S)** のときだけ。
//   ⇒ 「この Skeleton は危ない」ではなく「**この Skeleton を引けない Found が N 本ある**」と数える。
//
// 使い方:
//   node tools/import-audit.mjs                 … 罠(定義+定理が同居し Found へ推移 import)を並べる
//   node tools/import-audit.mjs --edge <F> <S>  … その 1 本を張ったら循環するかを判定する
//   node tools/import-audit.mjs --fragile       … 「使っているのに直接 import していない」を経路の本数で割る
//   node tools/import-audit.mjs --selftest      … 自己検査(★純関数だけを叩く)
//   共通: --root <lean/ABC3 の場所>  --ns Skeleton  --limit N
//
// ★★測れないこと(先に書く)
//   - ★「その定義を Found が**引きたがっている**か」は測れない。ここで数えるのは
//     **名前がコードに現れるか**だけ。★docstring の言及は剥がすが、同名の別宣言は残る
//     (`WeilDiv` は Skeleton と Found の両方にある)。⇒ ★同名がある場合は印字して区別する。
//   - ★Lean の名前解決(`open` / 部分修飾 / 記法)は再現していない。★語で当てているだけ。
//   - ★「直接 import していない」は**誤りではない**(Lean では推移 import が普通)。
//     ★危ないのは「到達経路が 1 本しかない」場合で、その 1 本を切ると落ちる(ArtinMap の型)。

import fs from 'node:fs';
import path from 'node:path';

// ────────────────────────────────────────────────────────────
// 純関数(selftest はここだけを叩く)
// ────────────────────────────────────────────────────────────
export const DECL = /^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+|partial\s+|unsafe\s+|scoped\s+)*(theorem|lemma|def|abbrev|structure|class|instance|inductive|axiom|opaque)\b\s+([A-Za-z_][A-Za-z0-9_'!?.]*)/;
export const THMK = new Set(['theorem', 'lemma']);
export const DEFK = new Set(['def', 'abbrev', 'structure', 'class', 'instance', 'inductive', 'axiom', 'opaque']);
/** ★台帳の付随宣言。`def foo.src` は「定義」ではない —— 数えると罠の判定が甘くなる。 */
export const LEDGER = /\.(src|needs|deps|obligation|obligations)$/;

/** `/- … -/`(ネスト)と行末 `--` を落とす。★行数は保つ。 */
export function stripComments(src) {
  let out = '', depth = 0;
  const s = String(src);
  for (let i = 0; i < s.length; i++) {
    if (!depth && s[i] === '/' && s[i + 1] === '-') { depth = 1; i++; continue; }
    if (depth) {
      if (s[i] === '/' && s[i + 1] === '-') { depth++; i++; continue; }
      if (s[i] === '-' && s[i + 1] === '/') { depth--; i++; continue; }
      if (s[i] === '\n') out += '\n';
      continue;
    }
    if (s[i] === '-' && s[i + 1] === '-') { while (i < s.length && s[i] !== '\n') i++; out += '\n'; continue; }
    out += s[i];
  }
  return out;
}

/** 1 本の .lean から import と宣言を読む。★名前は点まで読む(`def A.B.c` の名は `A.B.c`)。 */
export function parseLean(src) {
  const code = stripComments(src);
  const imports = [], decls = [], ns = [];
  code.split(/\r?\n/).forEach((ln, i) => {
    const t = ln.trim();
    let m;
    if ((m = /^import\s+([A-Za-z_][A-Za-z0-9_.]*)/.exec(t))) { imports.push(m[1]); return; }
    if ((m = /^namespace\s+([A-Za-z_][A-Za-z0-9_'.]*)/.exec(t))) { ns.push(m[1]); return; }
    if ((m = /^end\s+([A-Za-z_][A-Za-z0-9_'.]*)/.exec(t))) { if (ns.length && ns[ns.length - 1] === m[1]) ns.pop(); return; }
    const d = DECL.exec(ln);
    if (d) {
      const name = d[2].replace(/\.$/, '');
      const abs = name.includes('.') && /^[A-Z]/.test(name);
      decls.push({ kind: d[1], name, full: abs ? name : ((ns.join('.') ? ns.join('.') + '.' : '') + name),
                   ledger: LEDGER.test(name), line: i + 1 });
    }
  });
  return { imports, decls, code };
}

/** 推移閉包。★循環があっても止まらない(訪問済みで打ち切る)。 */
export function closure(graph) {
  const memo = new Map();
  const walk = (m, stack) => {
    if (memo.has(m)) return memo.get(m);
    if (stack.has(m)) return new Set();
    stack.add(m);
    const out = new Set();
    for (const im of graph.get(m) || []) { out.add(im); for (const x of walk(im, stack)) out.add(x); }
    stack.delete(m);
    memo.set(m, out);
    return out;
  };
  const all = new Map();
  for (const m of graph.keys()) all.set(m, walk(m, new Set()));
  return all;
}

/** F → S を張ったら循環するか。★循環するのは F ∈ reach(S) のときだけ。 */
export function wouldCycle(reach, F, S) {
  if (F === S) return true;
  const r = reach.get(S);
  return !!(r && r.has(F));
}

/** その名前へ到達できる**直接 import** の本数。★1 なら、その 1 本を切った瞬間に落ちる。 */
export function pathCount(reach, direct, target) {
  return direct.filter(h => h === target || (reach.get(h) && reach.get(h).has(target))).length;
}

// ────────────────────────────────────────────────────────────
// 読み出し
// ────────────────────────────────────────────────────────────
export function loadTree(root) {
  const files = [];
  (function walk(d) {
    for (const e of fs.readdirSync(d, { withFileTypes: true })) {
      const p = path.join(d, e.name);
      if (e.isDirectory()) walk(p); else if (e.name.endsWith('.lean')) files.push(p.replace(/\\/g, '/'));
    }
  })(root);
  const base = path.basename(root);
  const info = new Map(), graph = new Map();
  for (const f of files) {
    const rel = f.slice(root.length + 1);
    const mod = base + '.' + rel.replace(/\.lean$/, '').replace(/\//g, '.');
    const p = parseLean(fs.readFileSync(f, 'utf8'));
    info.set(mod, { file: f, rel, ...p });
    graph.set(mod, p.imports);
  }
  return { info, reach: closure(graph), base };
}

const realDefs = v => v.decls.filter(d => DEFK.has(d.kind) && !d.ledger);
const theorems = v => v.decls.filter(d => THMK.has(d.kind));
const pad = (s, n) => String(s).padStart(n);
const padr = (s, n) => { s = String(s); return s + ' '.repeat(Math.max(0, n - s.length)); };

function cmdTraps({ info, reach, base }, opts) {
  const NS = opts.ns || 'Skeleton', TO = opts.to || 'Found';
  const pre = `${base}.${NS}.`, tpre = `${base}.${TO}.`;
  const consumers = [...info.keys()].filter(m => m.startsWith(tpre));
  // 短い名前 → それを宣言している本(同名の見分けに使う)
  const owners = new Map();
  for (const [mod, v] of info) for (const d of v.decls) {
    const last = d.name.split('.').pop();
    if (!owners.has(last)) owners.set(last, new Set());
    owners.get(last).add(mod);
  }
  const rows = [];
  for (const [mod, v] of info) {
    if (!mod.startsWith(pre)) continue;
    const df = realDefs(v), th = theorems(v);
    if (!df.length || !th.length) continue;
    const blocked = consumers.filter(m => reach.get(mod).has(m));
    if (!blocked.length) continue;
    const wanted = new Map();                       // ★同じ名前が 2 度宣言されている本があるので name で畳む
    for (const d of df) {
      const last = d.name.split('.').pop();
      if (last.length < 4 || wanted.has(d.name)) continue;
      const re = new RegExp('(^|[^A-Za-z0-9_.\'])' + last.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '($|[^A-Za-z0-9_\'])');
      const users = consumers.filter(m => re.test(info.get(m).code));
      if (users.length) wanted.set(d.name, { name: d.name, n: users.length, dup: [...(owners.get(last) || [])].filter(x => x !== mod) });
    }
    rows.push({ mod, rel: v.rel, df: df.length, th: th.length, blocked: blocked.length, wanted: [...wanted.values()] });
  }
  rows.sort((a, b) => b.wanted.length - a.wanted.length || b.blocked - a.blocked);
  console.log(`== 罠: ${NS} の中で「定義と定理が同居し、${TO} を推移 import している」本 ==`);
  console.log(`   ⇒ ★その本を import できない ${TO} が居る。★定義だけ別の本に割れば引けるようになる。\n`);
  console.log(`  定義 定理  引けない${TO}  名前がコードに出る定義  ファイル`);
  for (const r of rows.slice(0, opts.limit ?? 25)) {
    console.log(`  ${pad(r.df, 4)}${pad(r.th, 5)}${pad(r.blocked, 12)}${pad(r.wanted.length, 12)}      ${r.rel}`);
    for (const w of r.wanted.slice(0, 4)) {
      console.log(`        ${padr(w.name, 30)} ${TO} ${pad(w.n, 4)} 本`
        + (w.dup.length ? `  ★同名が別の本にもある(${w.dup.length} 本: ${w.dup.slice(0, 2).map(x => x.replace(base + '.', '')).join(' ')})`
                        : '  ★同名なし ⇒ 引くならこの本'));
    }
  }
  console.log(`\n  ★該当 ${rows.length} 本 / ${NS} 全 ${[...info.keys()].filter(m => m.startsWith(pre)).length} 本`);
  console.log(`  ★★「名前がコードに出る」は**引きたがっている証拠ではない**(同名の別宣言を引いている場合がある)。`);
  let cyc = 0;
  for (const m of info.keys()) if (reach.get(m).has(m)) { cyc++; console.log(`  ★循環がある: ${m}`); }
  console.log(`  ★いまの木にある循環: ${cyc} 件(0 でなければ木は建たない)`);
}

function cmdEdge({ info, reach, base }, F, S) {
  const norm = (x) => x.replace(/\.lean$/, '').replace(/[\\/]/g, '.').replace(new RegExp('^.*?(?=' + base + '\\.)'), '');
  const f = info.has(F) ? F : norm(F), s = info.has(S) ? S : norm(S);
  if (!info.has(f)) { console.error(`引く側が見つからない: ${F} → ${f}`); process.exit(2); }
  if (!info.has(s)) { console.error(`引き先が見つからない: ${S} → ${s}`); process.exit(2); }
  const cyc = wouldCycle(reach, f, s);
  console.log(`  ${f}\n    → ${s}`);
  console.log(`  ★${cyc ? '循環する(引き先が引く側を推移的に import している)' : '循環しない'}`);
  if (cyc) {
    const via = (info.get(s).imports || []).filter(h => h === f || (reach.get(h) && reach.get(h).has(f)));
    console.log(`  ★戻ってくる経路の入口: ${via.join(' / ') || '(直接)'}`);
    const df = realDefs(info.get(s));
    console.log(`  ★引き先の定義 ${df.length} 件: ${df.map(d => d.name).join(' ').slice(0, 160)}`);
    console.log(`  ⇒ ★定義だけを ${info.get(s).rel.replace(/\.lean$/, 'Defs.lean')} に割れば引ける(Section1Defs / Section2Defs と同じ手)。`);
  }
}

/**
 * ★★★先回り —— 「この Skeleton を `XDefs.lean` に割ったら、本当に引けるようになるか」を**機械で確かめる**。
 *
 * ★M132 は「割れば直る」と書いたが、**割って直るとは限らない**:
 *   定義そのものが `Found` 側の名前を使っていたら、`XDefs.lean` へ移しても同じ import が要る。
 * ⇒ ここでは 3 つを分けて出す。★純関数(git も fs も触らない)。
 *
 *   1. **引けない `Found`**  … F ∈ reach(S) の F(= `F → S` を張ると循環する本)
 *   2. **毒の import**      … S の import のうち、引けない F へ到達するもの
 *                             (★`XDefs.lean` がこれを持ち越したら**循環は消えない**)
 *   3. **持ち越さないと失う名前** … 毒側にしか無い短名。定義の本体がこれに触れていたら
 *                             ★**そのままでは割れない**(先にその名前ごと動かす必要がある)
 *
 * 最後に **割った後の姿を組み立てて `wouldCycle` をやり直し**、
 * ★「N 本中 M 本が引けるようになる」を**数えて**返す(見込みではなく計算)。
 */
export function planSplit({ info, reach }, S, toPre) {
  const v = info.get(S);
  if (!v) return null;
  const defs = realDefs(v), thms = theorems(v);
  const consumers = [...info.keys()].filter(m => m.startsWith(toPre));
  const blocked = consumers.filter(m => reach.get(S).has(m));
  const hits = (h, f) => h === f || !!(reach.get(h) && reach.get(h).has(f));
  const imports = v.imports || [];
  const poison = imports.filter(h => blocked.some(f => hits(h, f)));
  const safe = imports.filter(h => !poison.includes(h));

  // 割った後に見える短名 / 失う短名
  const namesOf = (mods) => {
    const s = new Set();
    for (const m of mods) { const w = info.get(m); if (!w) continue; for (const d of w.decls) s.add(d.name.split('.').pop()); }
    return s;
  };
  const expand = (mods) => { const s = new Set(); for (const m of mods) { s.add(m); for (const x of (reach.get(m) || [])) s.add(x); } return s; };
  const keep = namesOf(expand(safe));
  const lost = new Set([...namesOf(expand(poison))].filter(n => !keep.has(n) && n.length >= 4));

  // 定義の本体(次の宣言の手前まで)が失う名前に触れているか
  const lines = v.code.split(/\r?\n/);
  const marks = [...v.decls].sort((a, b) => a.line - b.line);
  const bodyOf = (d) => {
    const nx = marks.find(x => x.line > d.line);
    return lines.slice(d.line - 1, nx ? nx.line - 1 : lines.length).join('\n');
  };
  // 失う名前を**どの本が宣言しているか**(★「では何を import すればよいか」に答えるため)。
  const owner = new Map();
  for (const m of expand(poison)) {
    const w = info.get(m); if (!w) continue;
    for (const d of w.decls) { const n = d.name.split('.').pop(); if (lost.has(n) && !owner.has(n)) owner.set(n, m); }
  }
  const isPoison = (m) => blocked.some(f => hits(m, f));

  const movable = [], stuck = [];
  for (const d of defs) {
    const used = new Set();
    for (const id of bodyOf(d).match(/[A-Za-z_][A-Za-z0-9_']*/g) || []) if (lost.has(id)) used.add(id);
    const needs = [...used].slice(0, 6).map(n => {
      const m = owner.get(n) || null;
      return { name: n, from: m, poison: m ? isPoison(m) : null };
    });
    (used.size ? stuck : movable).push({ name: d.name, line: d.line, needs });
  }

  // ★★割った後を組み立て直して数え直す。
  // ★★★ここは「safe だけ持ち越す」で数えては**意味が無い**(定義上 safe は blocked に届かないので
  //   必ず「全部 引ける」になる ―― 第 29 回が突然変異 #5 で気づいた **同語反復**)。
  //   ⇒ ★**移せない定義が要求している名前の宣言元も import する**ものとして数える。
  //   これなら「要る名前が毒の本にしか無い」場合に **freed < blocked** になる(= 割っても直らない)。
  const needImports = new Set();
  for (const d of stuck) for (const n of d.needs) if (n.from) needImports.add(n.from);
  const proposed = [...new Set([...safe, ...needImports])];
  const reachDefs = expand(proposed);
  const freed = blocked.filter(f => !reachDefs.has(f) && f !== S);
  const blockers = [...needImports].filter(m => isPoison(m));
  return { S, rel: v.rel, defs: defs.length, thms: thms.length, blocked, poison, safe,
           movable, stuck, freed, blockers, proposed, lost: lost.size };
}

function cmdPlan({ info, reach, base }, targets, opts) {
  const TO = opts.to || 'Found', toPre = `${base}.${TO}.`;
  const norm = (x) => x.replace(/\.lean$/, '').replace(/[\\/]/g, '.').replace(new RegExp('^.*?(?=' + base + '\\.)'), '');
  let list = targets.map(t => (info.has(t) ? t : norm(t)));
  if (opts.all) {
    const pre = `${base}.${opts.ns || 'Skeleton'}.`;
    list = [...info.keys()].filter(m => m.startsWith(pre)
      && realDefs(info.get(m)).length && theorems(info.get(m)).length
      && [...info.keys()].some(c => c.startsWith(toPre) && reach.get(m).has(c)));
  }
  console.log(`== 先回り: 「${TO} が引けるようにするには、どの定義を割るか」 ==`);
  console.log('   ★「引けない」= その本を import したら循環する。★割った後を組み直して数え直している。\n');
  let totBlocked = 0, totFreed = 0;
  for (const S of list) {
    const p = planSplit({ info, reach }, S, toPre);
    if (!p) { console.log(`  ★見つからない: ${S}`); continue; }
    totBlocked += p.blocked.length; totFreed += p.freed.length;
    console.log(`  ── ${p.rel}   定義 ${p.defs} / 定理 ${p.thms}`);
    console.log(`     引けない ${TO}: ★${p.blocked.length} 本`);
    console.log(`     割り先   : ${p.rel.replace(/\.lean$/, 'Defs.lean')}`);
    console.log(`     ★その import に持ち越してよい: ${p.safe.length ? p.safe.join(' ') : '(なし ⇒ import 無しで立つ)'}`);
    if (p.poison.length) console.log(`     ★★持ち越すと循環が残る(毒): ${p.poison.join(' ')}`);
    console.log(`     ★そのまま移せる定義 ${p.movable.length} 件: ${p.movable.map(d => d.name).join(' ').slice(0, 150) || '(なし)'}`);
    for (const d of p.stuck) {
      console.log(`     ★★そのままでは移せない: ${d.name}(${d.line} 行目) —— 毒側の名前を使っている`);
      for (const n of d.needs)
        console.log(`          ${padr(n.name, 26)} ${n.from ? n.from.replace(base + '.', '') : '(宣言元が木に無い = Mathlib か namespace)'}`
          + (n.from ? (n.poison ? '  ★★これも毒 ⇒ この本も割る必要がある' : '  ★毒ではない ⇒ Defs から import してよい') : ''));
    }
    console.log(`     ★Defs が実際に import することになる: ${p.proposed.length ? p.proposed.join(' ') : '(なし)'}`);
    const ok = p.freed.length === p.blocked.length;
    console.log(`     ⇒ ★★それで引けるようになる: ${p.freed.length} / ${p.blocked.length} 本${ok ? '  ★全部' : '  ★★全部ではない'}`);
    if (!ok) {
      console.log(`        ★★★直らない理由: 要る名前が **毒の本にしか無い** —— ${p.blockers.map(x => x.replace(base + '.', '')).join(' ')}`);
      console.log(`        ⇒ ★先にその本を割る(あるいは名前を Interface / Setup 側へ移す)必要がある。`);
    }
    console.log('');
  }
  console.log(`  ★合計: 引けない ${totBlocked} 本 → 割ると ${totFreed} 本が引けるようになる(${list.length} 本の Skeleton について)`);
}

function cmdFragile({ info, reach, base }, opts) {
  const byFull = new Map(), byShort = new Map();
  for (const [mod, v] of info) for (const d of v.decls) {
    byFull.set(d.full, { mod, ...d });
    const last = d.name.split('.').pop();
    if (!byShort.has(last)) byShort.set(last, []);
    byShort.get(last).push({ mod, ...d });
  }
  const IDENT = /[A-Za-z_][A-Za-z0-9_'.]*/g;
  let pairs = 0, one = 0;
  const rows = [];
  for (const [mod, v] of info) {
    const direct = v.imports, dset = new Set(direct), trans = reach.get(mod);
    const need = new Map();
    for (const m of v.code.replace(/^import .*$/gm, '').matchAll(IDENT)) {
      const id = m[0];
      let d = byFull.get(id);
      if (!d) { const c = byShort.get(id); if (!c || c.length !== 1) continue; d = c[0]; }
      if (d.mod === mod || dset.has(d.mod) || !trans.has(d.mod)) continue;
      if (!need.has(d.mod)) need.set(d.mod, new Set());
      need.get(d.mod).add(d.name);
    }
    for (const [g, names] of need) {
      pairs++;
      const k = pathCount(reach, direct, g);
      // ★`reach.get(h)` は木の外(Mathlib など)の import で undefined になる。★必ず守る。
      if (k === 1) { one++; rows.push({ rel: v.rel, g, names: [...names], via: direct.find(h => h === g || (reach.get(h) && reach.get(h).has(g))) }); }
    }
  }
  console.log('== 使っているのに直接 import していない —— ★到達経路の本数で割る ==\n');
  console.log(`  (使う本, 宣言している本)の組 ${pairs} 組 / ★経路が 1 本だけ ${one} 組(${(100 * one / Math.max(1, pairs)).toFixed(1)}%)`);
  console.log('  ★経路が 1 本 = その 1 本の import を切り替えた瞬間に `unknown identifier` で落ちる。');
  const byFile = new Map();
  for (const r of rows) byFile.set(r.rel, (byFile.get(r.rel) || 0) + 1);
  console.log(`  ★1 本経路を持つ本: ${byFile.size} 本\n`);
  console.log('  件数  ファイル');
  for (const [f, n] of [...byFile.entries()].sort((a, b) => b[1] - a[1]).slice(0, opts.limit ?? 15)) console.log(`  ${pad(n, 4)}  ${f}`);
}

// ────────────────────────────────────────────────────────────
function selftest() {
  let ok = 0, ng = 0;
  const t = (n, c) => { if (c) ok++; else { ng++; console.log('  NG ' + n); } };
  t('strip: 行末 --', stripComments('def a := 1 -- コメント\ndef b := 2').includes('def b'));
  t('strip: -- の中身は消える', !stripComments('def a := 1 -- filteredGroupOf').includes('filteredGroupOf'));
  t('strip: /- … -/', !stripComments('/- filteredGroupOf -/\ndef b := 2').includes('filteredGroupOf'));
  t('strip: 入れ子の /- -/', !stripComments('/- a /- b -/ c -/\ndef b := 2').includes('c'));
  t('strip: 行数を保つ', stripComments('a\n/- x -/\nb').split('\n').length === 3);
  const P = parseLean([
    'import ABC3.Found.X', 'import ABC3.Skeleton.Y', 'namespace ABC3.Skeleton.PGC',
    'def foo : Nat := 1', 'theorem bar : True := trivial', 'def bar.src : Source := {}',
    'end ABC3.Skeleton.PGC', 'noncomputable def ABC3.Interface.PGC.RamificationFiltration.filt : Nat := 0',
  ].join('\n'));
  t('parse: import を 2 本', P.imports.length === 2);
  t('parse: 名前空間を前に付ける', P.decls.some(d => d.full === 'ABC3.Skeleton.PGC.foo'));
  t('parse: 点を含む名前を最後まで読む', P.decls.some(d => d.name === 'ABC3.Interface.PGC.RamificationFiltration.filt'));
  t('parse: 完全修飾の宣言に名前空間を足さない', P.decls.some(d => d.full === 'ABC3.Interface.PGC.RamificationFiltration.filt'));
  t('parse: .src は台帳印', P.decls.find(d => d.name === 'bar.src').ledger === true);
  t('parse: 台帳印は定義に数えない', realDefs({ decls: P.decls }).length === 2);
  t('parse: 定理を拾う', theorems({ decls: P.decls }).length === 1);
  const g = new Map([['A', ['B']], ['B', ['C']], ['C', []], ['D', ['A']]]);
  const R = closure(g);
  t('closure: 推移で C まで', R.get('A').has('C'));
  t('closure: 逆向きは入らない', !R.get('A').has('D'));
  t('cycle: F ∈ reach(S) なら循環', wouldCycle(R, 'C', 'A') === true);
  t('cycle: そうでなければ循環しない', wouldCycle(R, 'D', 'A') === false);
  t('cycle: 自分自身は循環', wouldCycle(R, 'A', 'A') === true);
  const gc = new Map([['A', ['B']], ['B', ['A']]]);
  t('closure: 循環があっても止まらない', closure(gc).get('A').has('A'));
  t('path: 経路 1 本', pathCount(R, ['B'], 'C') === 1);
  t('path: 経路 0 本', pathCount(R, ['C'], 'A') === 0);   // ★C は何も import していない
  t('path: D は A 経由で C に届く(1 本)', pathCount(R, ['D'], 'C') === 1);
  t('path: 直接も 1 本と数える', pathCount(R, ['C'], 'C') === 1);
  t('path: 2 経路', pathCount(new Map([['X', new Set(['Z'])], ['Y', new Set(['Z'])]]), ['X', 'Y'], 'Z') === 2);

  // ── planSplit(★第 29 回)。★形は実機と同じ: Skeleton.S が Found.F を推移 import している ──
  const mk = (imports, src) => ({ file: '', rel: 'x.lean', ...parseLean(src), imports });
  const I = new Map([
    // ★`safeName` は **両側にある**(同名の別宣言)。★safe 側との差を取らないと
    //   cleanDef まで「移せない」に落ちる ―― 実機の 96 短名がこの形(M132)。
    ['B.Found.F',        mk([], 'def onlyInF : Nat := 0\ndef safeName : Nat := 9')],
    ['B.Interface.Safe', mk([], 'def safeName : Nat := 0')],
    ['B.Skeleton.Mid',   mk(['B.Found.F'], 'def midName : Nat := 0')],
    ['B.Skeleton.S',     mk(['B.Skeleton.Mid', 'B.Interface.Safe'],
      'def cleanDef : Nat := safeName\ndef dirtyDef : Nat := onlyInF\ntheorem thm : True := trivial')],
  ]);
  const G = new Map([...I].map(([k, v]) => [k, v.imports]));
  const RR = closure(G);
  const pl = planSplit({ info: I, reach: RR }, 'B.Skeleton.S', 'B.Found.');
  t('plan: 引けない Found を数える', pl.blocked.join() === 'B.Found.F');
  t('plan: 毒の import を名指す', pl.poison.join() === 'B.Skeleton.Mid');
  t('plan: 毒でない import は持ち越してよい', pl.safe.join() === 'B.Interface.Safe');
  t('plan: 毒側の名前を使わない定義は移せる', pl.movable.map(d => d.name).join() === 'cleanDef');
  t('plan: ★毒側の名前を使う定義は移せない', pl.stuck.map(d => d.name).join() === 'dirtyDef');
  t('plan: ★何が要るかを名指す', pl.stuck[0].needs.some(n => n.name === 'onlyInF'));
  t('plan: ★その名前の宣言元を出す', pl.stuck[0].needs.find(n => n.name === 'onlyInF').from === 'B.Found.F');
  t('plan: ★宣言元が毒かを判定する', pl.stuck[0].needs.find(n => n.name === 'onlyInF').poison === true);
  // ★★ここが要 —— dirtyDef が要る `onlyInF` は **B.Found.F(毒)にしか無い**ので、
  //   Defs はそれを import せざるを得ず、★**割っても引けるようにならない**。
  t('plan: ★★割っても直らない場合を捕まえる', pl.freed.length === 0 && pl.blocked.length === 1);
  t('plan: ★直らない理由(毒の本)を名指す', pl.blockers.join() === 'B.Found.F');
  t('plan: ★Defs が持つことになる import を出す', pl.proposed.sort().join() === 'B.Found.F,B.Interface.Safe');
  t('plan: 定理は移す対象に数えない', pl.defs === 2 && pl.thms === 1);
  t('plan: 無い本には null', planSplit({ info: I, reach: RR }, 'B.Nope', 'B.Found.') === null);
  // ★毒の import を持ち越したら循環は消えない(= safe に毒が混ざっていないこと)
  t('plan: ★safe に毒は入らない', !pl.safe.some(h => pl.poison.includes(h)));
  // ★引けない Found が 0 本なら、割る必要が無いこと
  const pl2 = planSplit({ info: I, reach: RR }, 'B.Interface.Safe', 'B.Found.');
  t('plan: 引けない Found が無ければ毒も無い', pl2.blocked.length === 0 && pl2.poison.length === 0);
  t('plan: そのときは全部の定義が移せる', pl2.stuck.length === 0);
  // ★★割ると本当に直る形 —— 要る名前が **毒でない本** にある場合
  const I3 = new Map([
    ['B.Found.F',        mk([], 'def onlyInF : Nat := 0')],
    ['B.Interface.Safe', mk([], 'def safeName : Nat := 0')],
    ['B.Setup',          mk([], 'def setupName : Nat := 0')],
    ['B.Skeleton.Mid',   mk(['B.Found.F'], 'def midName : Nat := 0')],
    ['B.Skeleton.T',     mk(['B.Skeleton.Mid', 'B.Setup'], 'def okDef : Nat := setupName\ntheorem th : True := trivial')],
  ]);
  const R3 = closure(new Map([...I3].map(([k, v]) => [k, v.imports])));
  const pl3 = planSplit({ info: I3, reach: R3 }, 'B.Skeleton.T', 'B.Found.');
  t('plan: ★★直る場合は freed が blocked に一致する', pl3.freed.join() === 'B.Found.F' && pl3.blocked.length === 1);
  t('plan: ★直る場合は詰まる本が無い', pl3.blockers.length === 0);
  t('plan: ★毒でない本の名前は失わない(移せる側に居る)', pl3.movable.map(d => d.name).join() === 'okDef');

  console.log(`\nselftest: ${ok}/${ok + ng}`);
  return ng === 0;
}

// ────────────────────────────────────────────────────────────
function main() {
  const argv = process.argv.slice(2);
  const has = f => argv.includes(f);
  const val = (f, d) => { const i = argv.indexOf(f); return i >= 0 && argv[i + 1] ? argv[i + 1] : d; };
  if (has('--selftest')) process.exit(selftest() ? 0 : 1);
  const here = path.resolve(path.join(new URL('.', import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, '$1'), '..'));
  const root = val('--root', path.join(here, 'lean', 'ABC3')).replace(/\\/g, '/');
  if (!fs.existsSync(root)) { console.error(`木が無い: ${root}(--root で指す)`); process.exit(2); }
  const t0 = Date.now();
  const tree = loadTree(root);
  const opts = { ns: val('--ns', 'Skeleton'), to: val('--to', 'Found'), limit: Number(val('--limit', '25')) };
  if (has('--edge')) {
    const i = argv.indexOf('--edge');
    cmdEdge(tree, argv[i + 1], argv[i + 2]);
  } else if (has('--plan') || has('--plan-all')) {
    const i = argv.indexOf('--plan');
    const t = [];
    for (let j = i + 1; i >= 0 && j < argv.length && !argv[j].startsWith('-'); j++) t.push(argv[j]);
    cmdPlan(tree, t, { ...opts, all: has('--plan-all') });
  } else if (has('--fragile')) cmdFragile(tree, opts);
  else cmdTraps(tree, opts);
  console.error(`  (${root} の .lean ${tree.info.size} 本 / ${((Date.now() - t0) / 1000).toFixed(1)} 秒)`);
}
if (process.argv[1] && process.argv[1].replace(/\\/g, '/').endsWith('import-audit.mjs')) main();
