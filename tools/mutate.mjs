#!/usr/bin/env node
// mutate.mjs —— 「わざと壊して、ゲートが鳴るか」を**再現可能に**回す台。
//
// ★★これを作った理由(メタ第 36 回が自分で見つけた失敗形):
//   第 36 回の突然変異 19 通りのうち **3 通りが「当たらなかった」**。
//   原因は「その木の `tools/*.mjs` が CRLF なので、探索文字列の `\n` が当たらなかった」ことで、
//   ★**「壊せなかった(当たらない)」と「壊したのに鳴らなかった(素通り)」が見分けられなかった。**
//   ⇒ ★この台は **置換が実際に起きたことを md5 で検算し、起きていなければ試験を失敗させる**。
//
// ★★この台が保証すること(4 つ):
//   (1) 置換前後で対象ファイルの md5 が**変わったこと**を確かめる。変わらなければ `NOTAPPLIED`。
//       ★`NOTAPPLIED` は「素通り」ではなく **★試験の穴**として、終了コードを 1 にする。
//   (2) 探索は **CRLF を LF に正規化してから**行い、書き戻すときに**元の改行に戻す**(M155 / M191)。
//   (3) 走り終わったら**必ず元のバイト列に戻し**、md5 が元に戻ったことを確かめる
//       (戻らなければ `RESTORE-FAILED` として大声で落ちる)。
//   (4) 「鳴った」の定義は **基準線の出力との差**である。exit code か stdout+stderr が変われば鳴った。
//       ★ゲートごとに「何を見れば鳴ったと言えるか」を書かなくてよい。
//
// 使い方:
//   node tools/mutate.mjs --spec <spec.json>      仕様に従って回す
//   node tools/mutate.mjs --spec <spec.json> --json
//   node tools/mutate.mjs --selftest              内蔵 fixture(★「当たらない置換」が落ちることも試す)
//
// spec.json の形:
//   {
//     "gate": ["node", "tools/idiom-recur.mjs", "--selftest"],   // 回すゲート(cwd = リポジトリ)
//     "mutations": [
//       { "id": "X1", "file": "tools/idiom-recur.mjs",
//         "find": "探索する字面(LF 正規化後)", "repl": "置換後", "all": false,
//         "note": "この規則を消す" }
//     ]
//   }
//   ★`find` に `\n` を書いてよい。★CRLF の木でも当たる(それがこの台の仕事)。
//   ★`all` を true にすると全部置換する(既定は最初の 1 つだけ)。
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { execFileSync } from 'node:child_process';

const ARGV = process.argv.slice(2);
const has = f => ARGV.includes(f);
const arg = f => { const i = ARGV.indexOf(f); return i >= 0 ? ARGV[i + 1] : null; };
const md5 = b => crypto.createHash('md5').update(b).digest('hex');

// ---------- 純関数(ここが selftest で較正できる本体) ----------

/** 改行の種類を数える。★`grep -c $'\r'` を使わない(本体が 2 度誤った) */
export function eolOf(buf) {
  const s = buf.toString('binary');
  const crlf = (s.match(/\r\n/g) || []).length;
  const lf = (s.match(/(?<!\r)\n/g) || []).length;
  return { crlf, lf, kind: crlf && lf ? 'mixed' : crlf ? 'crlf' : 'lf' };
}

/** CRLF → LF に正規化して探索し、書き戻すときは元の改行に戻す。
 *  戻り値 { text, applied } —— `applied` は置換が起きた件数。 */
export function applyMutation(rawText, find, repl, all = false) {
  const norm = rawText.replace(/\r\n/g, '\n');
  if (!find) return { text: rawText, applied: 0 };
  let applied = 0, out;
  if (all) {
    out = norm.split(find).join(repl);
    applied = norm.split(find).length - 1;
  } else {
    const i = norm.indexOf(find);
    if (i < 0) { out = norm; applied = 0; }
    else { out = norm.slice(0, i) + repl + norm.slice(i + find.length); applied = 1; }
  }
  return { text: out, applied };
}

/** 正規化した本文を、元の改行の種類で書き戻す形にする。 */
export function restoreEol(text, kind) {
  const lf = text.replace(/\r\n/g, '\n');
  return kind === 'crlf' ? lf.replace(/\n/g, '\r\n') : lf;
}

/** 判定 —— ★これが持ち場の要。3 通りしかない。 */
export function verdict({ applied, md5Same, rang }) {
  if (applied === 0 || md5Same) return 'NOTAPPLIED';   // ★試験の穴。素通りではない。
  return rang ? 'RANG' : 'SILENT';                     // 発火 / 素通り
}

// ---------- selftest ----------
function selftest() {
  const T = []; const eq = (n, a, b) => T.push({ n, ok: JSON.stringify(a) === JSON.stringify(b), a, b });
  eq('T1 eol lf', eolOf(Buffer.from('a\nb\n')).kind, 'lf');
  eq('T2 eol crlf', eolOf(Buffer.from('a\r\nb\r\n')).kind, 'crlf');
  eq('T3 eol mixed', eolOf(Buffer.from('a\r\nb\n')).kind, 'mixed');
  eq('T4 eol counts', [eolOf(Buffer.from('a\r\nb\nc\n')).crlf, eolOf(Buffer.from('a\r\nb\nc\n')).lf], [1, 2]);
  // ★★核心: CRLF の本文でも `\n` を含む探索が当たる(第 36 回が落ちた穴)
  eq('T5 crlf body is searchable with \\n', applyMutation('x\r\ncatch { return null; }\r\ny', 'catch { return null; }\ny', 'Z').applied, 1);
  eq('T6 lf body too', applyMutation('x\ncatch { return null; }\ny', 'catch { return null; }\ny', 'Z').applied, 1);
  // ★★「当たらない置換」は applied 0 になる —— これが NOTAPPLIED の入口
  eq('T7 miss gives 0', applyMutation('abc', 'zzz', 'Z').applied, 0);
  eq('T8 miss leaves text', applyMutation('abc', 'zzz', 'Z').text, 'abc');
  eq('T9 all replaces all', applyMutation('a.a.a', 'a', 'b', true).applied, 3);
  eq('T10 single replaces one', applyMutation('a.a.a', 'a', 'b', false).text, 'b.a.a');
  eq('T11 restore crlf', restoreEol('a\nb', 'crlf'), 'a\r\nb');
  eq('T12 restore lf', restoreEol('a\r\nb', 'lf'), 'a\nb');
  // ★★判定の表 —— 3 通りしかないこと
  eq('T13 verdict notapplied by count', verdict({ applied: 0, md5Same: false, rang: true }), 'NOTAPPLIED');
  eq('T14 verdict notapplied by md5', verdict({ applied: 1, md5Same: true, rang: true }), 'NOTAPPLIED');
  eq('T15 verdict rang', verdict({ applied: 1, md5Same: false, rang: true }), 'RANG');
  eq('T16 verdict silent', verdict({ applied: 1, md5Same: false, rang: false }), 'SILENT');
  // ★★「置換したのに md5 が同じ」= 恒等置換。★これも試験の穴として落とす
  eq('T17 identity replacement is a hole', (() => {
    const r = applyMutation('abc', 'b', 'b');
    return verdict({ applied: r.applied, md5Same: md5(r.text) === md5('abc'), rang: false });
  })(), 'NOTAPPLIED');
  // ★実際に本物のファイルで往復して、バイト列が戻ることを確かめる
  eq('T18 round trip keeps bytes', (() => {
    const p = path.join(process.cwd(), '.cache', 'mutate-selftest.tmp');
    fs.mkdirSync(path.dirname(p), { recursive: true });
    const orig = Buffer.from('alpha\r\nbeta\r\ngamma\r\n', 'utf8');
    fs.writeFileSync(p, orig);
    const before = md5(fs.readFileSync(p));
    const raw = fs.readFileSync(p, 'utf8');
    const r = applyMutation(raw, 'beta\ngamma', 'BETA\nGAMMA');
    fs.writeFileSync(p, Buffer.from(restoreEol(r.text, eolOf(orig).kind), 'utf8'));
    const mid = md5(fs.readFileSync(p));
    fs.writeFileSync(p, orig);
    const after = md5(fs.readFileSync(p));
    fs.unlinkSync(p);
    return [r.applied, before !== mid, before === after];
  })(), [1, true, true]);
  const bad = T.filter(t => !t.ok);
  for (const t of bad) console.log('  NG', t.n, '\n     got', JSON.stringify(t.a), '\n     want', JSON.stringify(t.b));
  console.log(`selftest ${T.length - bad.length}/${T.length}`);
  return bad.length === 0;
}

if (has('--selftest')) process.exit(selftest() ? 0 : 1);

// ---------- 本体 ----------
const specPath = arg('--spec');
if (!specPath) {
  console.error('使い方: node tools/mutate.mjs --spec <spec.json> [--json]');
  console.error('        node tools/mutate.mjs --selftest');
  process.exit(2);
}
const spec = JSON.parse(fs.readFileSync(specPath, 'utf8'));
const REPO = process.cwd();
const runGate = () => {
  const [cmd, ...rest] = spec.gate;
  try {
    const out = execFileSync(cmd === 'node' ? process.execPath : cmd, rest,
      { cwd: REPO, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'], maxBuffer: 1 << 28 });
    return { code: 0, out };
  } catch (e) { return { code: e.status ?? -1, out: (e.stdout || '') + (e.stderr || '') }; }
};

const base = runGate();
const baseSig = base.code + '|' + md5(base.out);
const rows = [];
for (const m of spec.mutations) {
  const p = path.join(REPO, m.file);
  const orig = fs.readFileSync(p);
  const before = md5(orig);
  const kind = eolOf(orig).kind;
  let r, after = before, g = { code: 0, out: '' };
  try {
    r = applyMutation(orig.toString('utf8'), m.find, m.repl, !!m.all);
    fs.writeFileSync(p, Buffer.from(restoreEol(r.text, kind), 'utf8'));
    after = md5(fs.readFileSync(p));
    if (r.applied > 0 && after !== before) g = runGate();
  } finally {
    fs.writeFileSync(p, orig);
    const back = md5(fs.readFileSync(p));
    if (back !== before) { console.error(`★★RESTORE-FAILED ${m.file} (${m.id}) —— 元に戻せていない。手で直すこと。`); process.exit(3); }
  }
  const sig = g.code + '|' + md5(g.out);
  const v = verdict({ applied: r.applied, md5Same: after === before, rang: sig !== baseSig });
  rows.push({ id: m.id, file: m.file, note: m.note || '', applied: r.applied, eol: kind, verdict: v });
}

const n = k => rows.filter(r => r.verdict === k).length;
if (has('--json')) console.log(JSON.stringify({ baseSig, rows }, null, 1));
else {
  console.log('== mutate —— わざと壊してゲートが鳴るか ==');
  console.log(`  ゲート    : ${spec.gate.join(' ')}`);
  console.log(`  基準線    : exit ${base.code} / md5 ${md5(base.out).slice(0, 12)}`);
  console.log(`  突然変異  : ${rows.length} 通り`);
  console.log('');
  for (const r of rows) {
    const mark = r.verdict === 'RANG' ? '  発火' : r.verdict === 'SILENT' ? '★素通り' : '★★当たらない';
    console.log(`  ${String(r.id).padEnd(5)} ${mark.padEnd(14)} 置換 ${r.applied} / ${r.eol.padEnd(5)} ${r.file}  ${r.note}`);
  }
  console.log('');
  console.log(`  発火 ${n('RANG')} / ★素通り ${n('SILENT')} / ★★当たらない(=試験の穴) ${n('NOTAPPLIED')}`);
  if (n('NOTAPPLIED')) console.log('  ★★「当たらない」は**素通りではない**。★探索文字列が届いていないので、試験を書き直すこと。');
}
process.exit(n('NOTAPPLIED') ? 1 : 0);
