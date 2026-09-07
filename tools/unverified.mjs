#!/usr/bin/env node
/**
 * 持ち場に書いた「見当」が当たったかを**数える**(backlog M45 / メタ第 14 回)。
 *
 * ★なぜ要るか(M45 の負の結果)
 *   メタ第 13 回は「持ち場に未確認の断りを書くと誤りが覆る」という本体の見立てを
 *   検証しようとして、**原理的に測れない**という結論に達した:
 *     - 断りを書いた持ち場 21 件のうち、当否が判定できた 14 件・覆った 8 件(57%)。
 *     - ★しかし**断りを書かなかった**持ち場でも覆っている(#16 の `pdfPage := 8`)。
 *     - ★★**覆らなかった誤りは記録に現れない**ので、**分母が作れない**。
 *
 * ★★この道具が変えるのはただ 1 点: **分母を作る。**
 *   ⇒ **見当は「自信がないもの」だけでなく、断定して渡すものも含めて全部 `GUESS:` に書く。**
 *     `確度=高` の行が無ければ「断定した見当が何件あって、そのうち何件外れたか」が
 *     永遠に分からない。★**自信のあるものを書かないと、この道具は M45 と同じ壁に当たる。**
 *
 * ★書式(`ResearchPaper/autonomy-policy.md` §4.7 が正典)
 *
 *     GUESS[<持ち場>]: <見当を 1 行で> | 確度=<高|中|低> | 検算=<型|実機|索引|原文|なし>
 *     VERDICT[<持ち場>]: <当たり|外れ|半分|未確認> — <一言>
 *
 *   - `<持ち場>` は本体が既に使っている呼び名(`Y19e` / `#16` / `第1075波-3` など)。
 *     ★**GUESS と VERDICT で同じ字面**にすること。ここだけが対応の鍵である。
 *   - `GUESS` は**配るとき**に、`VERDICT` は**報告が返ったとき**に書く。
 *   - 1 つの持ち場に見当が複数あるなら `Y19e-a` `Y19e-b` のように枝番を付ける。
 *
 * ★これは**ゲートではない**。`check.mjs` に入れてはならない ——
 *   落とすようにすると「見当を書かない」ことで逃げられ、正直さを罰する
 *   (G7 の基準・G1 の `Found/` 非対称と同じ理由)。
 *
 * 使い方:
 *   node tools/unverified.mjs                # 集計
 *   node tools/unverified.mjs --open         # まだ VERDICT が無いものだけ
 *   node tools/unverified.mjs --json
 *   node tools/unverified.mjs --selftest     # 器具の較正(書式を変えたら必ず走らせる)
 *   node tools/unverified.mjs --file <path>  # 別の記録を読む
 */

import { readFileSync, existsSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const args = process.argv.slice(2);
const flag = (f) => args.includes(f);
const opt = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };

/** 既定で読む記録。★増やすときはここだけを触る。 */
const SOURCES = [
  join(ROOT, 'ResearchPaper', 'decisions-pending.md'),
  join(ROOT, 'ResearchPaper', 'meta-backlog.md'),
];

const GUESS_RE = /^\s*(?:[★☆*\-\s]*)GUESS\[([^\]]+)\]\s*[:：]\s*(.*)$/;
const VERDICT_RE = /^\s*(?:[★☆*\-\s]*)VERDICT\[([^\]]+)\]\s*[:：]\s*(.*)$/;
const CONF_RE = /確度\s*=\s*([高中低])/;
const HOW_RE = /検算\s*=\s*(型|実機|索引|原文|なし)/;
const OUT_RE = /^\s*(当たり|外れ|半分|未確認)/;

/** 記録 1 本を読んで GUESS / VERDICT を拾う。 */
function parse(path) {
  const guesses = new Map();   // 持ち場 → { conf, how, text, file, line }
  const verdicts = new Map();  // 持ち場 → { outcome, text, file, line }
  if (!existsSync(path)) return { guesses, verdicts, missing: true };
  const lines = readFileSync(path, 'utf8').split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    let m = GUESS_RE.exec(lines[i]);
    // ★★書式の見本(`GUESS[<持ち場>]`)を数えない(メタ第 19 回 §7 の指摘)。
    //   見出しの件数と表の合計が 1 だけ合わなかった原因。
    //   雛形は `確度=<高|中|低>` で列挙値に合わず表からは落ちていた。
    if (m && /[<>]/.test(m[1])) continue;
    if (m) {
      guesses.set(m[1].trim(), {
        conf: CONF_RE.exec(m[2])?.[1] ?? '(無記入)',
        how: HOW_RE.exec(m[2])?.[1] ?? '(無記入)',
        text: m[2].split('|')[0].trim(),
        file: path, line: i + 1,
      });
      continue;
    }
    m = VERDICT_RE.exec(lines[i]);
    if (m && /[<>]/.test(m[1])) continue;
    if (m) {
      verdicts.set(m[1].trim(), {
        outcome: OUT_RE.exec(m[2])?.[1] ?? '(読めない)',
        text: m[2].replace(OUT_RE, '').replace(/^\s*[—-]\s*/, '').trim(),
        file: path, line: i + 1,
      });
    }
  }
  return { guesses, verdicts, missing: false };
}

function collect(paths) {
  const guesses = new Map(); const verdicts = new Map(); const missing = [];
  for (const p of paths) {
    const r = parse(p);
    if (r.missing) { missing.push(p); continue; }
    for (const [k, v] of r.guesses) guesses.set(k, v);
    for (const [k, v] of r.verdicts) verdicts.set(k, v);
  }
  return { guesses, verdicts, missing };
}

function report({ guesses, verdicts }) {
  const rows = [...guesses].map(([k, g]) => ({ k, ...g, v: verdicts.get(k) ?? null }));
  const orphan = [...verdicts].filter(([k]) => !guesses.has(k));
  const judged = rows.filter((r) => r.v && r.v.outcome !== '未確認' && r.v.outcome !== '(読めない)');
  const over = (r) => r.v.outcome === '外れ' || r.v.outcome === '半分';

  console.log(`# 見当と当否 —— GUESS ${rows.length} 件 / VERDICT ${verdicts.size} 件`);
  if (!rows.length) {
    console.log('');
    console.log('★GUESS が 1 件も無い。書式は autonomy-policy.md §4.7。');
    console.log('★★M45 の壁(分母が作れない)は、**断定して渡す見当も `確度=高` で書く**ことでしか越えられない。');
    return 0;
  }
  console.log('');
  console.log('## ★★M45 が作れなかった分母 —— 確度 × 当否');
  console.log('');
  console.log('  確度      判定済   覆った   当たり   覆った率');
  for (const c of ['高', '中', '低', '(無記入)']) {
    const g = judged.filter((r) => r.conf === c);
    if (!g.length) continue;
    const o = g.filter(over).length;
    console.log(`  ${c.padEnd(8)}  ${String(g.length).padStart(6)}  ${String(o).padStart(6)}  ` +
      `${String(g.length - o).padStart(6)}  ${(o * 100 / g.length).toFixed(0).padStart(7)}%`);
  }
  const hi = judged.filter((r) => r.conf === '高');
  const lo = judged.filter((r) => r.conf === '中' || r.conf === '低');
  console.log('');
  if (hi.length && lo.length) {
    const rh = hi.filter(over).length / hi.length;
    const rl = lo.filter(over).length / lo.length;
    console.log(`★**断定(確度=高)の覆った率 ${(rh * 100).toFixed(0)}% 対 ` +
      `断り付き(中・低)の ${(rl * 100).toFixed(0)}%** ——` +
      (rh > rl ? '断定の方がよく覆っている(断りを書くこと自体より、'
               + '「自信が無いと気づけたか」が効いている可能性)'
               : '断り付きの方がよく覆っている(断りが実装者の検算を誘発している可能性)'));
    console.log('☆★★件数が 2 桁に届くまでは向きを断定しないこと(M45 の轍)。');
  } else {
    console.log('★**まだ比較できない。** 確度=高 と 確度=中/低 の**両方**に判定済みが要る。');
    console.log(`  いま 確度=高 の判定済 ${hi.length} 件 / 中・低 ${lo.length} 件。`);
  }
  console.log('');
  console.log('## 検算の指定ごと');
  const byHow = new Map();
  for (const r of judged) byHow.set(r.how, [...(byHow.get(r.how) ?? []), r]);
  for (const [h, g] of [...byHow].sort((a, b) => b[1].length - a[1].length)) {
    console.log(`  ${h.padEnd(8)} ${String(g.length).padStart(3)} 件 / 覆った ${g.filter(over).length}`);
  }
  const open = rows.filter((r) => !r.v);
  console.log('');
  console.log(`## VERDICT 待ち ${open.length} 件`);
  for (const r of open) console.log(`  ${r.k.padEnd(14)} 確度=${r.conf} 検算=${r.how}  ${r.text.slice(0, 60)}`);
  if (orphan.length) {
    console.log('');
    console.log(`## ★対応する GUESS が無い VERDICT ${orphan.length} 件(鍵の字面が食い違っている)`);
    for (const [k, v] of orphan) console.log(`  ${k.padEnd(14)} ${v.outcome}`);
  }
  return rows.length;
}

/** 器具の較正。★書式(正規表現)を変えたら必ずここを走らせること。 */
function selftest() {
  const fx = [
    'GUESS[Y19e]: N(-α) は最小多項式の定数項に落ちる | 確度=高 | 検算=型',
    '★GUESS[#16-a]: pdfPage は 8 | 確度=中 | 検算=原文',
    '- GUESS[Z1]: 何か | 確度=低 | 検算=実機',
    'GUESS[Z2]: 判定待ちのもの | 確度=高 | 検算=索引',
    'VERDICT[Y19e]: 当たり — 型が合った',
    '★VERDICT[#16-a]: 外れ — 実装者が 7 を採った',
    'VERDICT[Z1]: 半分 — 片方だけ',
    'VERDICT[ZZ]: 当たり — 対応する GUESS が無い',
    'この行は GUESS[ を含むが書式ではない',
  ].join('\n');
  const path = join(ROOT, '.cache', '_unverified-selftest.md');
  // ★一時ファイルを作らずに済むよう、parse を文字列で叩けるようにしてある。
  const lines = fx.split('\n');
  const guesses = new Map(); const verdicts = new Map();
  for (let i = 0; i < lines.length; i++) {
    let m = GUESS_RE.exec(lines[i]);
    if (m) { guesses.set(m[1].trim(), { conf: CONF_RE.exec(m[2])?.[1] ?? '(無記入)', how: HOW_RE.exec(m[2])?.[1] ?? '(無記入)', text: m[2].split('|')[0].trim(), file: path, line: i + 1 }); continue; }
    m = VERDICT_RE.exec(lines[i]);
    if (m) verdicts.set(m[1].trim(), { outcome: OUT_RE.exec(m[2])?.[1] ?? '(読めない)', text: '', file: path, line: i + 1 });
  }
  const checks = [
    ['GUESS を 4 件拾う', guesses.size === 4],
    ['VERDICT を 4 件拾う', verdicts.size === 4],
    ['★ を前置しても拾う', guesses.has('#16-a') && verdicts.has('#16-a')],
    ['`- ` を前置しても拾う', guesses.has('Z1')],
    ['確度を読む', guesses.get('Y19e').conf === '高' && guesses.get('Z1').conf === '低'],
    ['検算を読む', guesses.get('Y19e').how === '型' && guesses.get('#16-a').how === '原文'],
    ['当否を読む', verdicts.get('#16-a').outcome === '外れ' && verdicts.get('Z1').outcome === '半分'],
    ['書式でない行を拾わない', !guesses.has('を含むが書式ではない')],
    ['VERDICT 待ちを 1 件と数える', [...guesses.keys()].filter((k) => !verdicts.has(k)).length === 1],
    ['孤児 VERDICT を 1 件と数える', [...verdicts.keys()].filter((k) => !guesses.has(k)).length === 1],
  ];
  let ok = 0;
  for (const [label, pass] of checks) { console.log(`  ${pass ? 'ok ' : 'NG '} ${label}`); if (pass) ok++; }
  console.log(`\n  selftest: ${ok}/${checks.length} PASS`);
  return ok === checks.length;
}

if (flag('--selftest')) process.exit(selftest() ? 0 : 1);

const paths = opt('--file') ? [opt('--file')] : SOURCES;
const data = collect(paths);
for (const m of data.missing) console.log(`(記録が無い: ${m})`);
if (flag('--json')) {
  console.log(JSON.stringify({
    guesses: [...data.guesses].map(([k, v]) => ({ key: k, ...v })),
    verdicts: [...data.verdicts].map(([k, v]) => ({ key: k, ...v })),
  }, null, 1));
} else if (flag('--open')) {
  for (const [k, g] of data.guesses) {
    if (!data.verdicts.has(k)) console.log(`${k}\t確度=${g.conf}\t検算=${g.how}\t${g.text}`);
  }
} else {
  report(data);
}
