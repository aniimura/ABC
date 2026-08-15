// 番号を持たない**概念**を節点として拾い、番号付き項目のグラフに合流させる。
//
// ★なぜ要るか(2026-08-15 実測)
//   我々のグラフは番号付き項目しか節点にしないので、**名前で運ばれる概念**が出ない。
//   実害: 半グラフ(番号なしの地の文で定義)、anabelioid、半グラフ的アナベリオイド
//   (Definition 2.1 は存在するが、論文中で番号で参照される箇所が1つも無い)。
//   ★目的は**工数の把握**であって網羅ではない。語彙表は tools/concepts.json(人手)。
//
// ★辺の張り方(ここが設計)
//   概念 → 番号付き項目 : 概念の**定義箇所**が引用している項目
//   概念 → 概念         : 概念の定義箇所に**他の概念名が出る**とき
//   項目 → 概念         : 項目の**冒頭(= statement 相当)**に概念名が出るとき
//     ★本文全体で名前を拾うと、Frobenioid のような語がどこにでも出て
//      入次数が爆発し、循環だらけになる。**冒頭に限る**のはそのため。
//
// ★限界: 語の一致で辺を張るので、番号による引用より精度が低い。
//   「定義箇所」は「we shall refer to … as X / we define … X / X is defined」の
//   最初の出現で近似している。**目視していない。**

import { readFileSync } from 'node:fs';

// ★定義箇所の当て方。狭すぎると「定義未特定」が量産され、その概念が層 0 に落ちて
//   「依存が無い=すぐ着手できる」と**誤読**される。実測で 31 中 14 が未特定だったので広げた。
const DEF_PATS = (re) => [
  new RegExp(`(?:we shall refer to[^.]{0,240}?as (?:an? |the )?|shall be referred to as (?:an? |the )?)(?:“|")?${re}`, 'i'),
  new RegExp(`we (?:shall )?call[^.]{0,200}?(?:“|")?${re}`, 'i'),
  new RegExp(`we (?:shall )?define[^.]{0,200}?(?:“|")?${re}`, 'i'),
  new RegExp(`(?:“|")?${re}(?:”|")?[^.]{0,60}?(?:is|are) defined`, 'i'),
  new RegExp(`(?:a|an|the)\\s+(?:“|")?${re}(?:”|")?\\s+(?:is|consists of|to be)\\b`, 'i'),
  new RegExp(`(?:Definition|Notation)[^\\n]{0,40}\\n?[^.]{0,120}?${re}`, 'i'),
];

/**
 * @param {object} o
 * @param {(tag:string)=>({decls:Map,lines:string[]}|null)} o.load  論文を読む
 * @param {string[]} o.tags        探索対象の論文タグ
 * @param {RegExp}  o.NOT_A_DEP    依拠でない文脈
 * @param {string}  o.conceptsPath tools/concepts.json
 * @param {(tag:string,key:string)=>string|null} o.resolveTag
 * @param {string}  o.KIND
 */
export function conceptNodes({ load, tags, NOT_A_DEP, conceptsPath, resolveTag, KIND }) {
  const CS = JSON.parse(readFileSync(conceptsPath, 'utf8')).concepts;
  const nodes = new Map();   // "CONCEPT / term" -> {page, tag, defLine}
  const edges = [];          // [from, to]
  const key = (t) => `CONCEPT / ${t}`;

  // ── ① 各概念の定義箇所を見つける ────────────────────────
  const defBody = new Map();
  for (const c of CS) {
    const search = c.defIn ? [c.defIn] : tags;
    let found = null;
    const pats = DEF_PATS(c.re);
    const at = (P, tg, i, approx) => {
      let page = 0;
      for (const d of P.decls.values()) if (d.line <= i && i < d.end) { page = d.page; break; }
      return { tag: tg, line: i, page, approx,
               // ★概念どうしの辺は**定義文の範囲だけ**から取る。
               //   ±36 行の本体から取ると、IUT の語彙は解説文中で相互参照が濃いため
               //   23 節点の塊になり、順序が消える(2026-08-15 実測)。
               body: P.lines.slice(i, i + 5).join('\n') };
    };
    for (const tg of search) {
      const P = load(tg); if (!P) continue;
      for (let i = 0; i < P.lines.length && !found; i++) {
        const win = P.lines.slice(i, i + 3).join(' ');   // pdftotext は文を折り返すので3行つなぐ
        if (pats.some((p) => p.test(win))) found = at(P, tg, i, false);
      }
      if (found) break;
    }
    // ★定義文型に当たらない場合は、その語の**初出**で近似する。
    //   「未特定のまま層 0 に置く」より正直で、誤読も少ない(approx として印を残す)。
    if (!found) {
      const term = new RegExp(c.re, 'i');
      for (const tg of search) {
        const P = load(tg); if (!P) continue;
        for (let i = 0; i < P.lines.length && !found; i++) if (term.test(P.lines[i])) found = at(P, tg, i, true);
        if (found) break;
      }
    }
    if (!found) { nodes.set(key(c.term), { tag: '?', page: 0, missing: true }); continue; }
    nodes.set(key(c.term), { tag: found.tag, page: found.page, approx: found.approx });
    defBody.set(c.term, found);
  }

  // ── ② 概念 → 番号付き項目 / 概念 → 概念 ──────────────────
  // ★概念 → **番号付き項目** の辺は張らない(2026-08-15 の判断)。
  //   理由: 概念は番号付き項目の**中で**定義されることが多いので、その項目へ辺を張ると
  //   「項目 → 概念(項目が概念を使う)」と往復して**循環が量産される**。実測で
  //   循環ボックスが 13 → 17 に増え、45 節点の EtTh 循環が新たに生まれた。
  //   概念は語彙の層として、**概念どうしの順序だけ**を持たせる。
  //   ★これは情報を捨てている。捨てた理由は「工数把握には順序が要る」からで、
  //    網羅性より順序を優先した。
  for (const c of CS) {
    const f = defBody.get(c.term); if (!f) continue;
    for (const d of CS) {
      if (d.term === c.term) continue;
      // 自分の名前を含む語(Frobenioid ⊂ pre-Frobenioid 等)は辺にしない
      if (new RegExp(d.re, 'i').test(c.term) || new RegExp(c.re, 'i').test(d.term)) continue;
      if (new RegExp(d.re, 'i').test(f.body)) edges.push([key(c.term), key(d.term)]);
    }
  }

  // ── ③ ★向きを固定して非循環にする ────────────────────────
  //   仮定: **概念は導入順に積み上がる**。すなわち、後から定義された概念は
  //   先に定義された概念に依拠しうるが、逆は無い。
  //   順位 = (定義論文の順) × 1000 + (物理ページ)。論文の順は下の PAPER_ORDER。
  //   ★これは**モデル化の仮定であって実測ではない**。この仮定を置かないと
  //    概念どうしが1つの塊になり、工数把握に使えない。仮定を置いた事実を残す。
  const PAPER_ORDER = ['@The Geometry of Anabelioids', 'AbsAnab', 'SemiAnbd', 'pGC',
    'AbsTopI', 'AbsTopII', 'AbsTopIII', 'FrdI', 'FrdII', 'EtTh', 'GenEll',
    'IUTchI', 'IUTchII', 'IUTchIII', 'IUTchIV'];
  const rankOf = (t) => {
    const v = nodes.get(key(t));
    if (!v || v.missing) return 1e9;
    const pi = PAPER_ORDER.indexOf(v.tag);
    return (pi < 0 ? PAPER_ORDER.length : pi) * 10000 + (v.page ?? 0);
  };
  const kept = [], dropped = [];
  for (const [a, b] of edges) {
    const ta = a.slice(10), tb = b.slice(10);
    (rankOf(ta) > rankOf(tb) ? kept : dropped).push([a, b]);
  }

  return { nodes, edges: kept, droppedByOrder: dropped, concepts: CS, key, rankOf };
}

/** 番号付き項目 → 概念。項目の**冒頭**に概念名が出るときだけ張る。 */
export function itemToConceptEdges({ load, itemKeys, concepts, key, HEAD = 1200 }) {
  const out = [];
  const cache = new Map();
  for (const k of itemKeys) {
    const [tag, name] = k.split(' / ');
    const P = load(tag); if (!P) continue;
    const d = P.decls.get(name); if (!d) continue;
    const head = P.lines.slice(d.line, d.end).join(' ').slice(0, HEAD);
    for (const c of concepts) {
      let re = cache.get(c.term);
      if (!re) { re = new RegExp(c.re, 'i'); cache.set(c.term, re); }
      if (re.test(head)) out.push([k, key(c.term)]);
    }
  }
  return out;
}
