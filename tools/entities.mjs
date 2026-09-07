/**
 * tools/entities.mjs —— 名前つき HTML 実体の表。**1 体系 1 表**にする。
 * ================================================================
 *
 * ★なぜ 1 本にしたか(メタ第 13〜15 回の実測)
 * -------------------------------------------------
 * 同じ表が `check.mjs`(**100 種**)と `brief.mjs`(**90 種**)に**別々に**あった。
 * 突き合わせると:
 *
 *   - 両方にある鍵で値が食い違うのは ★**`mu` の 1 件だけ**
 *     (`check` = U+00B5 MICRO / `brief` = U+03BC GREEK)
 *   - `check` にしか無い鍵 **23** / `brief` にしか無い鍵 **13** / 和集合 **113**
 *   - ★**生きている穴**: 生きた `1_Structured`(69 本)で使われる実体 53 種のうち
 *     `brief` が開けないものが **3 種・延べ 8 回**
 *     (`&otimes;` 4 / `&eacute;` 3 / `&ouml;` 1)。
 *     ★`&otimes;`(テンソル積)は**数学の記号**なので、これは体裁の問題ではない。
 *
 * ★`mu` の向きは**数えて**決めた(メタ第 14 回 M53):
 *
 *   | | U+00B5 MICRO | U+03BC GREEK |
 *   |---|---|---|
 *   | `0_Source/*.txt`(143 本) | 358 | ★**4,350** |
 *   | `lean/ABC3`(2,222 本) | 90 | ★**3,202** |
 *
 *   ⇒ **U+03BC を採る**(`brief` 側が正しかった)。
 *   ★`normalize()` は 2 つを畳まないので、legacy が生き返れば即座に効く。
 *
 * ★なぜ「どちらかの CLI から import する」形にしなかったか
 * -------------------------------------------------------
 * `check.mjs` も `brief.mjs` も**モジュール本体が CLI**である
 * (`check.mjs` の末尾は `if (all || only('--selftest')) selftest();` と `process.exit`)。
 * 一方から他方を import すると、★**相手の CLI が自分の argv で走り出す**。
 * よって共有するものは**副作用の無い第 3 のモジュール**に置くしかない。
 * ★ここには import も副作用も置かないこと。
 */

/**
 * 名前つき実体 → 文字。★**113 種**(`check` 100 ∪ `brief` 90)。
 *
 * ★経緯のコメントは消さないこと —— 「なぜその批が足されたか」は
 * 表そのものより価値がある(取りこぼしは**手で維持している**ことから来るため)。
 */
export const ENTITIES = {
  amp: '&', lt: '<', gt: '>', quot: '"', apos: "'", nbsp: ' ',
  ldquo: '“', rdquo: '”', lsquo: '‘', rsquo: '’',
  minus: '−', times: '×', ge: '≥', le: '≤',
  rarr: '→', larr: '←', harr: '↔',
  sube: '⊆', supe: '⊇', sub: '⊂', sup: '⊃',
  cong: '≅', or: '∨', and: '∧',
  Gamma: 'Γ', alpha: 'α', chi: 'χ', prime: '′',
  sect: '§', hellip: '…', mdash: '—', ndash: '–',
  // ★★2026-09-07 追加。★理由: Yoshida §2–§4 の構造化(68 件)で **S4 が 12 件落ちた**が、
  //   原因は逐語の誤りではなく **この表が標準の名前つき実体を取りこぼしていた**ことだった。
  //   ★実体はコーパス全体で使われている(`&middot;` 既存 55 / `&sigma;` 45 / `&isin;` 22)ので、
  //   ★**追加の前後で NG 件数を測ってから採った**(増えたら戻す、という手順を踏んだ)。
  //   ★`&ne;` は `≠` に開くが、`pdftotext` は斜線を落とすので
  //   **`data-txt="="` を併記しないと通らない**([[pdftotext-drops-negation]])。
  isin: '∈', ni: '∋', middot: '·', cap: '∩', cup: '∪',
  equiv: '≡', ne: '≠', empty: '∅', infin: '∞', bull: '•',
  rArr: '⇒', lArr: '⇐', hArr: '⇔', Prime: '″',
  prod: '∏', sum: '∑', part: '∂', radic: '√',
  pi: 'π', theta: 'θ', sigma: 'σ', beta: 'β', psi: 'ψ', phi: 'φ',
  // ★`mu` は U+03BC(GREEK SMALL LETTER MU)。U+00B5(MICRO SIGN)ではない。上の測定。
  lambda: 'λ', mu: 'μ', nu: 'ν', tau: 'τ', rho: 'ρ', delta: 'δ',
  Lambda: 'Λ', Sigma: 'Σ', Theta: 'Θ', Phi: 'Φ', Delta: 'Δ', Omega: 'Ω',
  // ★★2026-09-07 メタ第 13 回。★**手で足したのではなく `--entities` が数えたものを足した。**
  //   実測(1_Structured 87 本): 名前つき実体は異なり 76 種。この表は 67 種しか持たず、
  //   ★ゲート対象 69 本で **4 種 / 延べ 10 件**、legacy 18 本で **18 種 / 延べ 311 件**
  //   (最多は `&Pi;` 214)が**黙って `&Pi;` のまま残っていた**。
  //   ★同じ形の穴が今日 29 件の偽 NG を出したので、**穴の残りを機械で数えてから塞ぐ**。
  //   ★以後は `node tools/check.mjs --entities` が「表に無い実体」を数える。
  Pi: 'Π', gamma: 'γ', epsilon: 'ε', zeta: 'ζ', eta: 'η', iota: 'ι',
  kappa: 'κ', omega: 'ω', otimes: '⊗', darr: '↓', uarr: '↑', bot: '⊥', top: '⊤',
  plusmn: '±', Qopf: 'ℚ', Zopf: 'ℤ', Nopf: 'ℕ', Ropf: 'ℝ', Copf: 'ℂ',
  eacute: 'é', Eacute: 'É', ouml: 'ö', uuml: 'ü', auml: 'ä', egrave: 'è', agrave: 'à',
  // ★`&sup1;` `&sup2;` `&sup3;` は**旧い正規表現では 1 件も当たらなかった** ——
  //   `/&([a-zA-Z]+);/` が数字を含まないため。名前の文字類を広げて初めて効く。
  sup1: '¹', sup2: '²', sup3: '³', frac12: '½', frac13: '⅓', frac14: '¼',
  // ★★2026-09-07 メタ第 15 回。ここから下の **13 種は `brief.mjs` にしか無かった**もの。
  //   `check.mjs` 側に入れても**ゲートの NG は 1 件も動かない**ことを実測してから足した
  //   (生きた HTML に 1 回も出ないため)。★出たときに黙って素通りしないことに意味がある。
  mapsto: '↦', notin: '∉', not: '¬', forall: '∀', exist: '∃', prop: '∝',
  deg: '°', divide: '÷', sdot: '⋅', xi: 'ξ', upsilon: 'υ', Xi: 'Ξ', Psi: 'Ψ',
};

/**
 * 数値参照と名前つき実体を開く。
 *
 * ★名前は数字を含みうる(`sup2` `frac12`)。`brief.mjs` 側は長らく
 * `/&([a-zA-Z]+);/` で、★**表に持っていても開けない**状態だった
 * (legacy に `&sup1;` `&sup2;` `&sup3;` が実在する)。ここで 1 つに揃える。
 */
export function decodeEntities(s) {
  return s
    .replace(/&#x([0-9a-fA-F]+);/g, (_, h) => String.fromCodePoint(parseInt(h, 16)))
    .replace(/&#(\d+);/g, (_, d) => String.fromCodePoint(parseInt(d, 10)))
    .replace(/&([a-zA-Z][a-zA-Z0-9]*);/g, (m, name) => (name in ENTITIES ? ENTITIES[name] : m));
}
