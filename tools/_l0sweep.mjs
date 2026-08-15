// 使い捨て: 層0の §0 語を mathlib に当てて件数を出す(一次スクリーニング)
// ★件数は「語の一致」でしかない。概念の一致は別に確かめること。
import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join } from 'node:path';

const ROOT = 'lean/.lake/packages/mathlib/Mathlib';
const files = [];
(function walk(d) {
  for (const e of readdirSync(d)) {
    const p = join(d, e);
    if (statSync(p).isDirectory()) walk(p);
    else if (p.endsWith('.lean')) files.push(p);
  }
})(ROOT);
const blob = files.map((f) => readFileSync(f, 'utf8'));

const WORDS = [
  ['AbsTopI', 'generically scheme-like', /generically scheme-?like/i],
  ['AbsTopI', 'scheme-like', /scheme-?like/i],
  ['AbsTopI', 'slim', /\bslim\b/i],
  ['AbsTopI', 'partial coarsification', /coarsification/i],
  ['AbsTopI', 'order (§0)', /\bIsOrder\b|\bLinearOrder\b/],
  ['AbsTopII', 'sturdy', /\bsturdy\b/i],
  ['AbsTopIII', 'countably ordered set', /countably ordered/i],
  ['AbsTopIII', 'id-rigid', /id-?rigid/i],
  ['AbsTopIII', 'nexus', /\bnexus\b/i],
  ['AbsTopIII', 'saturated', /\bsaturated\b/i],
  ['AbsTopIII', 'symmetrically saturated', /symmetrically saturated/i],
  ['EtTh', 'essential image', /essential image|essImage|EssImage/],
  ['EtTh', 'group-saturated', /group-?saturated/i],
  ['EtTh', 'isomorph', /\bisomorph\b/i],
  ['EtTh', 'isomorphism-full', /isomorphism-?full/i],
  ['EtTh', 'perf-saturation', /perf-?saturation/i],
  ['FrdI', 'iso-subanchor', /subanchor/i],
  ['FrdI', 'characteristic (monoid)', /characteristic type/i],
  ['FrdI', 'fiberwise-surjective', /fiberwise-?surjective/i],
  ['FrdI', 'fsm-morphism', /fsm-?morphism/i],
  ['FrdI', 'fsmi-morphism', /fsmi-?morphism/i],
  ['FrdI', 'integral (monoid)', /IsCancelMul|IsCancelAdd|CancelCommMonoid/],
  ['FrdI', 'irreducible (monoid)', /\bIrreducible\b/],
  ['FrdI', 'primary (monoid)', /\bIsPrimary\b|\bprimary\b/i],
  ['FrdI', 'rigid', /\brigid\b/i],
  ['FrdI', 'saturated (monoid)', /IsSaturated|Saturated/],
  ['FrdI', 'sharp', /\bsharp\b/i],
  ['FrdI', 'torsion-free', /IsAddTorsionFree|IsMulTorsionFree|torsion-?free/i],
  ['FrdI', 'perfect (monoid)', /perfect monoid|IsPerfect/i],
  ['IUTchI', 'cyclotomic', /\bcyclotomic\b/i],
  ['IUTchI', 'poly-isomorphism', /poly-?isomorphism/i],
  ['IUTchI', 'pseudo-monoid', /pseudo-?monoid/i],
  ['SemiAnbd', 'indissectible', /indissectible/i],
  ['SemiAnbd', 'smooth log curve', /smooth log curve|LogScheme|log scheme/i],
  ['SemiAnbd', 'commensurably terminal', /commensurably terminal|Commensurable/i],
];

console.log(`mathlib ファイル数: ${files.length}`);
console.log('paper        word                        hits  例');
for (const [paper, w, re] of WORDS) {
  let n = 0; let ex = '';
  for (let i = 0; i < blob.length; i++) {
    const ms = blob[i].match(new RegExp(re.source, re.flags.includes('i') ? 'gi' : 'g'));
    if (ms) { n += ms.length; if (!ex) ex = files[i].replace(ROOT + '\\', '').replace(/\\/g, '/'); }
  }
  console.log(`${paper.padEnd(12)} ${w.padEnd(27)} ${String(n).padStart(5)}  ${ex}`);
}
