// 各論文の **§0(Notations and Conventions)** を節点として拾う。
//
// ★なぜ要るか(2026-08-15 実測)
//   依存グラフの「層 0」を読むと、そこには Frobenioid・anabelioid・semi-graph といった
//   **中層の語彙**が並ぶ。層 0 なのは「引用が無い」からで、「語彙が要らない」からではない。
//   実際に層 0 の項目を読むと:
//     [FrdI] Definition 1.1: … is pre-divisorial if it is integral [cf. §0],
//                            saturated [cf. §0], and of characteristic type [cf. §0].
//   ★**番号のない §0 を指している。** §0 は実在するが番号付き項目ではないので
//   我々のグラフに1つも現れていなかった。
//   §0 への参照回数(実測): FrdI 52 / AbsTopI 30 / IUTchI 29 / AbsTopIII 19 /
//                          SemiAnbd 18 / EtTh 17 / ★IUTchIII 0
//   ★左の層ほど §0 依存が濃い。層 0 を「着手可能」と読むには §0 を入れるしかない。
//
// ★辺の張り方(精度優先)
//   ① §0 の本文から、定義文型で**定義されている語**を抜く。
//   ② 番号付き項目の本文に `[cf. §0]` が出る箇所を探し、その**直前 70 字**に
//      ① の語があれば、項目 → §0 語 の辺を張る。
//   ★`[cf. §0]` が無い箇所では張らない。語が偶然出ただけでは辺にしない。
//     (Frobenioid のような語がどこにでも出て入次数が爆発するのを避けるため)

// ★行末ハイフン分割を戻す。pdftotext は "char- acteristically" のように割る。
//   戻さないと、抽出側でも照合側でも語が一致しない(2026-08-15 実測: 取りこぼしの主因)。
export const dehyph = (s) => s.replace(/(\w)-\s+(\w)/g, '$1$2').replace(/\s+/g, ' ');

const TERM = '([a-z][a-z0-9 ()≥-]{2,38}?)';
const END = '(?=[,.;:]|\\s+(?:and|or|if|which|whose|that|when|for|in|of|to)\\b)';
const DEF_LINES = [
  // we shall refer to … as (a|an|the)? X   ★冠詞は任意("as monoprime" の形がある)
  new RegExp(`we shall refer to[^.]{0,220}?\\bas (?:a |an |the )?${TERM}${END}`, 'gi'),
  // we shall say that … is (a|an)? X
  new RegExp(`we shall say that[^.]{0,180}?\\bis (?:a |an )?${TERM}${END}`, 'gi'),
  // … will be referred to as (a|an|the)? X
  new RegExp(`will be referred to as (?:a |an |the )?${TERM}${END}`, 'gi'),
  // we (shall )?write/denote by X
  new RegExp(`we (?:shall )?(?:write|denote by)\\s+${TERM}${END}`, 'gi'),
  // X is defined to be / X will be called
  new RegExp(`${TERM}\\s+(?:is defined to be|will be called)`, 'gi'),
];

const STOP = new Set(['set', 'element', 'elements', 'map', 'maps', 'object', 'objects',
  'category', 'morphism', 'morphisms', 'group', 'groups', 'ring', 'field', 'monoid',
  'result', 'following', 'above', 'below', 'case', 'situation', 'notation', 'usual',
  'natural', 'evident', 'same', 'given', 'such', 'unique', 'sub', 'the', 'a', 'an']);

/** §0 の行範囲を返す。無ければ null。 */
export function section0Range(lines) {
  let s = -1;
  for (let i = 0; i < lines.length; i++) {
    if (/^\s*(§0\.\s*Notations|Section 0:\s*Notations)/.test(lines[i])) s = i;
  }
  if (s < 0) return null;
  for (let i = s + 1; i < lines.length; i++) {
    if (/^\s*(§1\.|Section 1:)/.test(lines[i])) return [s, i];
  }
  return [s, Math.min(lines.length, s + 900)];
}

/** §0 で定義されている語を抜く。term -> {line, page} */
export function section0Terms(lines, pageOf) {
  const r = section0Range(lines);
  if (!r) return new Map();
  const [s, e] = r;
  const out = new Map();
  // 6行ずつの窓で走査(pdftotext は文を折り返す)
  for (let i = s; i < e; i++) {
    const win = dehyph(lines.slice(i, i + 6).join(' '));
    for (const re of DEF_LINES) {
      re.lastIndex = 0;
      let m;
      while ((m = re.exec(win))) {
        const t = m[1].trim().toLowerCase().replace(/\s+/g, ' ');
        if (t.length < 4 || t.length > 38) continue;
        if (STOP.has(t)) continue;
        if (/^\d/.test(t)) continue;
        if (!out.has(t)) out.set(t, { line: i, page: pageOf[i] ?? 0 });
      }
    }
  }
  return out;
}

/**
 * 項目 → §0 語 の辺。`[cf. §0]` の直前 WINDOW 字に §0 語があるときだけ張る。
 * @returns [{term, page}] の配列
 */
export function section0Refs(body, terms, WINDOW = 70) {
  const flat = dehyph(body);
  const out = new Set();
  const re = /\[cf\.\s*§0|\bcf\.\s*§0|§0\s*[\],]/g;
  let m, seen = 0, matched = 0;
  while ((m = re.exec(flat))) {
    seen++;
    const pre = flat.slice(Math.max(0, m.index - WINDOW), m.index).toLowerCase();
    let hit = false;
    for (const t of terms.keys()) if (pre.includes(t)) { out.add(t); hit = true; }
    if (hit) matched++;
  }
  // ★語を同定できなかった §0 参照も、依存であることには変わりない。
  //   捨てると「§0 を引いていない」と誤読され、層 0 が実際より広く見える。
  //   同定できた分と分けて、`(§0・語は未同定)` という節点に集める。
  //   実測(2026-08-15): §0 参照 165 件のうち語まで同定できたのは 64 件(約 39%)。
  if (seen > matched) out.add('(§0・語は未同定)');
  return [...out];
}
