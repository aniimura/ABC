// Gained 族 26 件が AxLemmaGraded / AxSenTate まで繋がっているかを測る。
// 使い方: node tools/gained-family-reach.mjs
import fs from 'node:fs';
import path from 'node:path';

const DIR = 'lean/ABC3/Found/PGC';
const files = fs.readdirSync(DIR).filter(f => f.endsWith('.lean'));

// --- 1. 全宣言を切り出す（theorem / def、シグネチャは `:= ` か `:= by` まで） ---
const decls = []; // {name, file, line, sig}
const srcOf = new Map();
for (const f of files) {
  const text = fs.readFileSync(path.join(DIR, f), 'utf8');
  srcOf.set(f, text);
  const lines = text.split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    const m = /^(?:noncomputable\s+)?(theorem|lemma|def)\s+([A-Za-z_][A-Za-z0-9_'.]*)/.exec(lines[i]);
    if (!m) continue;
    let sig = lines[i];
    let j = i;
    while (j + 1 < lines.length && !/:=/.test(sig) && j - i < 40) {
      j++;
      sig += '\n' + lines[j];
    }
    decls.push({ kind: m[1], name: m[2], file: f, line: i + 1, sig });
  }
}

// --- 2. 出口（結論に AxSenTate / AxLemmaGraded / AxLemma / AxWildDescent が出る宣言） ---
const exitPat = /\b(AxSenTate|AxLemmaGraded|AxLemma|AxWildDescent)\b/;
const exits = decls.filter(d => {
  const concl = d.sig.split(/:=/)[0];
  // 結論はシグネチャの最後の `:` 以降を粗く見る
  const idx = concl.lastIndexOf(') :');
  const tail = idx >= 0 ? concl.slice(idx) : concl;
  return exitPat.test(tail);
});

// --- 3. Gained 族 26 件 ---
const FAM = /^exists_norm_sub_algebraMap_le_prod_axDecay/;
const fam = decls.filter(d => FAM.test(d.name));

// 仮説の数（シグネチャ中の先頭が `(` か `{` か `[` の塊を数える。粗い近似）
function hypCount(sig) {
  const head = sig.split(/:=/)[0];
  // 名前より後ろだけ
  const after = head.replace(/^(?:noncomputable\s+)?(theorem|lemma|def)\s+[^\s]+/, '');
  let depth = 0, groups = 0;
  for (let i = 0; i < after.length; i++) {
    const c = after[i];
    if (c === '(' || c === '{' || c === '[') { if (depth === 0) groups++; depth++; }
    else if (c === ')' || c === '}' || c === ']') depth--;
  }
  return groups;
}

// --- 4. 逆参照: 各名前が「自分の定義行以外」で何回出るか ---
function usageCount(name) {
  let n = 0;
  const re = new RegExp('\\b' + name.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '\\b', 'g');
  for (const [f, text] of srcOf) {
    const lines = text.split(/\r?\n/);
    for (let i = 0; i < lines.length; i++) {
      if (/^(?:noncomputable\s+)?(theorem|lemma|def)\s/.test(lines[i])) {
        // 定義行そのものは数えない（.src は数える）
        const m = /^(?:noncomputable\s+)?(theorem|lemma|def)\s+([A-Za-z_][A-Za-z0-9_'.]*)/.exec(lines[i]);
        if (m && m[2] === name) continue;
      }
      const hits = lines[i].match(re);
      if (hits) n += hits.length;
    }
  }
  return n;
}

console.log('=== 1. Gained 族（^exists_norm_sub_algebraMap_le_prod_axDecay*） ===');
console.log('件数 = ' + fam.length + ' / ファイル数 = ' + new Set(fam.map(d => d.file)).size);
console.log('');
console.log('仮説数\t参照数\tfile:line\tname');
const famRows = fam.map(d => ({ ...d, h: hypCount(d.sig), u: usageCount(d.name) }));
famRows.sort((a, b) => a.h - b.h);
for (const d of famRows) {
  console.log(`${d.h}\t${d.u}\t${d.file}:${d.line}\t${d.name}`);
}

console.log('');
console.log('=== 2. 結論に AxSenTate / AxLemmaGraded / AxLemma / AxWildDescent が出る宣言 ===');
for (const d of exits) {
  const concl = d.sig.split(/:=/)[0];
  const idx = concl.lastIndexOf(') :');
  const tail = (idx >= 0 ? concl.slice(idx + 3) : concl).replace(/\s+/g, ' ').trim();
  console.log(`${d.file}:${d.line}\t${d.name}\t→ ${tail.slice(0, 90)}`);
}

console.log('');
console.log('=== 3. Gained 族の名前が、出口を出す宣言の中で使われているか ===');
const exitNames = new Set(exits.map(d => d.name));
// 出口宣言の「証明本体」を取り、その中に Gained 族の名前が出るかを見る
function bodyOf(d) {
  const text = srcOf.get(d.file);
  const lines = text.split(/\r?\n/);
  let out = [];
  for (let i = d.line - 1; i < lines.length; i++) {
    if (i > d.line - 1 && /^(?:noncomputable\s+)?(theorem|lemma|def)\s/.test(lines[i])) break;
    out.push(lines[i]);
  }
  return out.join('\n');
}
let found = 0;
for (const e of exits) {
  const body = bodyOf(e);
  const used = famRows.filter(d => body.includes(d.name));
  if (used.length) {
    found++;
    console.log(`${e.file}:${e.line} ${e.name}  ← ${used.map(u => u.name).join(', ')}`);
  }
}
if (found === 0) console.log('★ 0 件 —— Gained 族は出口の証明の中で一度も使われていない。');

console.log('');
console.log('=== 4. 私の AxWildDescent 側 3 本の仮説数（比較用） ===');
for (const n of ['axSenTate_of_deep_bound', 'axSenTate_of_deep_loss_bound',
                 'axSenTate_of_cyclotomic_loss', 'axLemma_of_eventually_Icc',
                 'axSenTate_of_axLemmaGraded', 'axSenTate_of_axDecay',
                 'axLemma_of_axDecay', 'axWildDescent_prime']) {
  const d = decls.find(x => x.name === n);
  if (d) console.log(`${hypCount(d.sig)}\t${usageCount(n)}\t${d.file}:${d.line}\t${n}`);
  else console.log(`----\t----\t(見つからない)\t${n}`);
}
