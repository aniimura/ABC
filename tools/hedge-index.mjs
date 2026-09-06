// 原文の「省略の合図」(immediately / formally / routine / one verifies …)を数え、
// 項目ごとに地図にする道具。
//
// ★動機(2026-08-20)
//   原典は証明を「it follows immediately」「a routine verification」で畳むことが多い。
//   我々はそれを 1 つずつ開いてきたが、**どこにいくつ残っているかを見ていなかった**。
//   CLAUDE.md の「進め方」は「必要なものが出たらスケルトンを足す」だが、
//   ★**省略の合図は「これから必要になるものの予告」である**——先に数えておけば、
//   分解(chain)を書く前に規模の当たりが付く。
//
// ★2026-09-05 の一般化(メタ backlog M7)
//   それまでこの道具は**論文 5 本を名前で埋め込んでいた**。CLAUDE.md は
//   「着手前に必ず数える」と定めているのに、`--paper pGC` は FrdI の地図を出し、
//   `--item "Lemma 3.5"` は「見つからない」と「合図 0 件」を**同じ空表**で返していた。
//   ——手順が論文によって効いたり効かなかったりするのは、手順書の穴である。直した点:
//     1) 論文の表を `ResearchPaper/papers.json`(49 本)から引く
//     2) 頁の同定を `===== [page N] =====` と**改頁文字 \f** の両方に対応させる
//        (実測: 47 本中 26 本がマーカー式、**21 本が \f 式**。後者は今まで頁が出なかった)
//     3) 見出しの書式を 3 通り(行独立 `N.` / 行末コロン `N:` / 走り込み `N. 本文`)にし、
//        **走り込みは論文の書式が走り込み優勢のときだけ**有効にする
//        (FrdI で無条件に有効にすると折り返し行を見出しと誤る——実測 1 件)
//     4) `.src` の照合を **`paper` 欄で分ける**(それまで論文をまたいで混ぜていた。
//        実測: 見出し 129 件のうち **25 件(19%)を 2 本以上の論文が名乗る**——
//        `Theorem 4.2` は FrdI / CorrHyp / Falt1 / pGC の 4 本)
//     5) `--all` で全論文の被覆表を出す(「効かない論文」が黙って隠れないように)
//
// ★出力
//   1) 合図の語ごとの件数(論文全体)
//   2) 項目(Proposition N.M 等)ごとの合図の数と、我々の実装状況(条なし .src / 条つき / 未)
//   3) 未着手の項目のうち合図が多いもの = 「畳まれた量が多い」順の作業候補
//
// ★限界
//   * 合図の語が必ずしも省略とは限らない(定義の中の "clearly" など)。
//   * 「項目に属する」の判定は**直前の見出し**で行う近似である。
//   * 合図の語は英語なので、仏語の原典(Asterisque / EGA / Del)では 0 件になる。
//     これは「省略が無い」ではなく「測れていない」である——`--all` の `語` 欄で見分ける。
//   * 原文 txt は gitignore 下にあるので、無ければ静かに終わる(CI で落とさない)。
//
// 使い方: node tools/hedge-index.mjs [--paper FrdI] [--json] [--cite] [--item "Proposition 1.10"]
//         node tools/hedge-index.mjs --all          … 全論文の被覆表
//         node tools/hedge-index.mjs --papers       … 使える論文の鍵の一覧

import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(HERE, '..');
const SRC_DIR = path.join(REPO, 'ResearchPaper', '0_Source');

/** ★論文の表は `papers.json` が正典(以前はここに 5 本を埋め込んでいた)。 */
function loadPapers() {
  const j = JSON.parse(readFileSync(path.join(REPO, 'ResearchPaper', 'papers.json'), 'utf8'));
  const out = new Map();
  for (const [key, v] of Object.entries(j.papers)) out.set(key, { ...v, txt: `${v.file}.txt` });
  return out;
}
const PAPERS = loadPapers();

const args = process.argv.slice(2);
const asJson = args.includes('--json');
const paperIdx = args.indexOf('--paper');
const paper = paperIdx >= 0 ? args[paperIdx + 1] : 'FrdI';
const itemIdx = args.indexOf('--item');
const onlyItem = itemIdx >= 0 ? args[itemIdx + 1] : null;

/** 合図の語。★分類は「実際に開けてみた経験」に基づく(README を見よ)。 */
const HEDGES = [
  { key: 'routine', re: /\broutine\b/i, note: '★最大級。畳まれた量がいちばん多い' },
  // ★★2026-09-06 修正: 旧版は `/\bformal(?:ly)?\b/i` で **裸の `formal` にも当たっていた**。
  //   Lubin–Tate の章では "formal group law" が数十件あり、[MilneCFT] Chapter I では
  //   **57 件中 44 件が偽陽性**だった(合図 131 → 実質 84)。
  //   ★合図としての `formal` は「formally, this follows」の副詞形なので `formally` に絞る。
  //   `formal consequence` のような名詞修飾は稀で、取りこぼしても下界が下がるだけ
  //   (`.needs` は下界という規約に沿う)。偽陽性は**見積を膨らませる**ので害が大きい。
  { key: 'formally', re: /\bformally\b/i, note: '新しい材料は要らない。手数はある' },
  { key: 'immediately', re: /\bimmediat(?:e|ely)\b/i, note: '玉石混交。要点検' },
  { key: 'one verifies', re: /\b(?:one|One)\s+(?:verifies|checks|computes)\b/, note: '短いことが多い' },
  { key: 'easily', re: /\beasil(?:y|ier)\b/i, note: '短いことが多い' },
  { key: 'clearly', re: /\bclearl(?:y)\b/i, note: '短いことが多い' },
  { key: 'well-known', re: /\bwell-known\b/i, note: '★外部の在庫を指す。mathlib を測ること' },
  { key: 'similarly', re: /\b(?:similarly|in a similar (?:way|manner))\b/i, note: '前の議論の再演' },
];

// ───────────────── 第 2 の水路: 角括弧の傍注 ─────────────────
// ★2026-09-05(メタ第 5 回、M12)。上の HEDGES は **FrdI 較正の語彙**であり、
//   「語で畳む」論文しか測れない。実測(47 本):
//     傍注/KB ≥ 0.8 の 18 本(すべて Mochizuki、2004 年以降)… 既定 0.60/KB に対し傍注 1.87/KB
//     残り 29 本(前期 Mochizuki + 他著者すべて)      … 既定 0.26/KB に対し傍注 0.13/KB
//   すなわち**後期 Mochizuki は畳む量の 3.1 倍を角括弧に入れており、道具はそこを見ていなかった**。
//   例: [GenEll] Lemma 3.5 は既定 8 語が 0 件なので「原文が畳んでいない」と表示されていたが、
//   実際は `[Thus, (EH)F is isogenous to EF , hence has semi-stable reduction at all the
//   finite primes of F and extends to a one-dimensional semi-abelian scheme …]` が
//   **3 主張を 1 文に畳んでいる**(引用は 1 つも抱えていない)。
//
// ★なぜ HEDGES に語を足さないか(実測)。`Thus,` / `it follows` / `we conclude` を
//   HEDGES に足すと FrdI は 259 → 580 件(+124%)に膨れ、**Example 6.1 が 4 → 5 件**になる
//   (増えたのは行 6186 の `Thus,` ——書き下された結論であって省略ではない)。
//   GenEll Lemma 3.5 でも 6 件出るが、真の畳みは**そのうち 1 件だけ**(精度 1/6)。
//   傍注として数えれば同じ箇所で**精度 1/1**。よって**別の水路として数え、合図には混ぜない**。
//
// ★入れ子 1 段まで許す局所的な正規表現を使う。**大域の depth 走査は使えない**
//   ——FrdI は行 1216 に閉じ忘れの `[` が 1 個あり、深さで数えると**そこから 5798 行**が
//   depth≥1 に沈んで傍注が 698 → 196 件に化ける(実測)。M11 の「分母が実装依存」の実例。
const ASIDE_RE = /\[(?:[^\[\]]|\[[^\[\]]*\])*\]/g;
const ASIDE_CAP = 600;   // これを超える一致は抽出の壊れ。全 40,525 件中 98 件(0.24%)だけ落ちる
const ASIDE_MINW = 4;    // [12] / [Mzk1] / [FC] は傍注ではなく書誌の鍵
const BIB_RE = /^\[?[A-Z][A-Za-z]{0,7}\]?,?\s*(?:Chapter|Ch\.|Theorem|Proposition|Lemma|Corollary|Remark|Definition|§|pp?\.|vol|Introduction)/i;
const IE_RE = /^i\.e\.?\b/i;
// ★入れ子の書誌鍵は **2 文字以上** を要求する。`\[[A-Za-z]` だと数学の `V [L]` /
//   `V [M]`(実測: FrdI Example 6.1 で 3 件)を引用と誤り、
//   `[which thus induces a morphism V [L] →V [M] over V]` が推論から引用に化ける。
const CITE_IN = /\bcf\.|\b(?:Proposition|Theorem|Lemma|Corollary|Definition|Remark|Chapter|Ch\.|§)\s*[0-9IVX]|\[[A-Za-z]{2,}[0-9]*\]/;
const INFER_IN = /\b(?:[Tt]hus|hence|therefore|it follows|consequently)\b/;
/** 傍注/KB がこれ以上なら「傍注式に畳む論文」。★実測の谷は 0.29 → 1.00(3.4 倍)。 */
const ASIDE_STYLE_MIN = 0.8;

/** 原文から傍注を取り、種別を付けて返す。種別は排他:
 *  引用 = 手順書を抱えている / 推論 = 引用なしで Thus・hence 等を抱えている(★最重要) / 傍注 = その他の地の文。
 *  i.e.(言い換え)と書誌の鍵は傍注として数えない。 */
function scanAsides(text) {
  const starts = [0];
  for (let i = 0; i < text.length; i++) if (text[i] === '\n') starts.push(i + 1);
  const lineOf = (idx) => {
    let lo = 0, hi = starts.length - 1;
    while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (starts[mid] <= idx) lo = mid; else hi = mid - 1; }
    return lo + 1;
  };
  const out = [];
  let dropped = 0, ie = 0, ref = 0;
  ASIDE_RE.lastIndex = 0;
  let m;
  while ((m = ASIDE_RE.exec(text)) !== null) {
    if (m[0].length > ASIDE_CAP) { dropped += 1; continue; }
    const body = m[0].slice(1, -1).replace(/\s+/g, ' ').trim();
    const w = body ? body.split(' ').length : 0;
    if (w < ASIDE_MINW) { ref += 1; continue; }
    if (IE_RE.test(body)) { ie += 1; continue; }
    if (BIB_RE.test(body)) { ref += 1; continue; }
    const cite = CITE_IN.test(body);
    const kind = cite ? '引用' : (INFER_IN.test(body) && w >= 6 ? '推論' : '傍注');
    out.push({ line: lineOf(m.index), body, words: w, kind });
  }
  return { list: out, dropped, ie, ref };
}

const KIND = '(?:Proposition|Theorem|Definition|Corollary|Lemma|Example|Remark)';
const NUM = '[0-9]+(?:\\.[0-9]+)*';
const LEAD = '^[ \\t]{0,6}';

// ★見出しの 3 書式。番号のあとがどう閉じるかで分ける。
//   DISPLAY … 行が見出しだけで終わる(または '(' で題が続く)。Mochizuki の全論文がこれ。
//             'Proposition 5.5, (iii), below]' のような**前方参照を拾わない**ための形。
//   COLON   … 'Proposition 1.1:' と**行末で閉じる**もの。pGC がこれ(実測 11 件)。
//             ★行末を要求しないと 'Theorem 5.6: By the log purity theorem …' のような
//             折り返し行を拾う(Config / LocProP / IUTchII で実測 4 件)。
//   RUNIN   … 'Definition 2.10. Let S ⊆ Node(G) …' と本文が続くもの。NodNon / HASurII /
//             BC / ABSLean / Stacks / MilneAV がこれ。★**書式の多数決で有効にする**。
const DISPLAY_RE = new RegExp(`${LEAD}${KIND}\\s+${NUM}\\.(?:\\s*$|\\s*\\()`);
const COLON_RE = new RegExp(`${LEAD}${KIND}\\s+${NUM}\\s*:\\s*$`);
const RUNIN_RE = new RegExp(`${LEAD}${KIND}\\s+${NUM}\\.\\s+\\S`);
const KEY_RE = new RegExp(`${LEAD}(${KIND})\\s+(${NUM})`);
const PAGE_RE = /^=====\s*\[page\s+(\d+)\]\s*=====/;

function readSource(key) {
  const rec = PAPERS.get(key);
  if (!rec) return { err: 'unknown', want: null };
  const p = path.join(SRC_DIR, rec.txt);
  if (!existsSync(p)) return { err: 'missing', want: path.join('ResearchPaper', '0_Source', rec.txt) };
  return { text: readFileSync(p, 'utf8'), rec };
}

/** ★原文を走査して「行 → 頁」「行 → 属する項目」を作る。
 *  頁は `[page N]` マーカーを最優先し、無ければ改頁文字 \f を数える
 *  (実測: \f の総数は papers.json の pdfPages と 21 本すべてで一致する)。 */
function scanSource(text) {
  const lines = text.split(/\r?\n/);
  const marked = lines.some((s) => PAGE_RE.test(s));
  // ★書式の多数決: 走り込みだけで拾える見出しが、行独立の見出しより多ければ走り込み式。
  let nDisplay = 0, nRunInOnly = 0;
  for (const s of lines) {
    const t = s.replace(/\f/g, '');
    if (DISPLAY_RE.test(t)) nDisplay += 1;
    else if (RUNIN_RE.test(t)) nRunInOnly += 1;
  }
  const runIn = nRunInOnly > nDisplay;
  return { lines, marked, runIn, nDisplay, nRunInOnly };
}

function mapItems(scan) {
  const { lines, marked, runIn } = scan;
  let page = marked ? null : 1;
  const rows = []; // 行ごとの {page, item}
  let item = null;
  const seen = new Map(); // item → 最初に現れた頁
  for (const s of lines) {
    const pm = s.match(PAGE_RE);
    if (pm) { page = Number(pm[1]); rows.push({ page, item }); continue; }
    const t = s.replace(/\f/g, '');
    if (DISPLAY_RE.test(t) || COLON_RE.test(t) || (runIn && RUNIN_RE.test(t))) {
      const m = t.match(KEY_RE);
      if (m) {
        item = `${m[1]} ${m[2]}`;
        if (!seen.has(item)) seen.set(item, page);
        rows.push({ page, item, heading: true });
        const nf = (s.match(/\f/g) ?? []).length;
        if (nf && !marked) page += nf;
        continue;
      }
    }
    rows.push({ page, item });
    const nf = (s.match(/\f/g) ?? []).length;
    if (nf && !marked) page += nf;
  }
  return { rows, seen };
}

/** Lean 側の `.src` を集める: paper → (item 見出し → 条なしか条つきか)
 *  ★`Source` は 1 行の `{ paper := "…", pdfPage := N, item := "…", … }` リテラルなので、
 *  レコードごとに paper と item を**対にして**取る。以前は item だけを見ていたため、
 *  `Theorem 4.2` が FrdI / CorrHyp / Falt1 / pGC のどれのものか分からなかった。 */
function scanLean() {
  const roots = [path.join(REPO, 'lean', 'ABC3')];
  const byPaper = new Map();
  const REC = /\{[^{}]*item\s*:=\s*"([^"]*)"[^{}]*\}/g;
  const HEAD = new RegExp(`^(?:Proposition|Theorem|Definition|Corollary|Lemma|Example|Remark)\\s+${NUM}`);
  const walk = (d) => {
    for (const e of readdirSync(d)) {
      const p = path.join(d, e);
      const st = statSync(p);
      if (st.isDirectory()) walk(p);
      else if (e.endsWith('.lean')) {
        const s = readFileSync(p, 'utf8');
        let m;
        REC.lastIndex = 0;
        while ((m = REC.exec(s)) !== null) {
          const raw = m[1];
          const head = raw.match(HEAD);
          if (!head) continue;
          const pm = m[0].match(/paper\s*:=\s*"([^"]*)"/);
          const pk = pm ? pm[1] : '?';
          const key = head[0];
          if (!byPaper.has(pk)) byPaper.set(pk, new Map());
          const items = byPaper.get(pk);
          const rec = items.get(key) ?? { bare: false, qualified: 0 };
          // 条なし = item がちょうど見出しだけ(「, (ii)」や「 — …」が付かない)
          if (raw.trim() === key) rec.bare = true;
          else rec.qualified += 1;
          items.set(key, rec);
        }
      }
    }
  };
  for (const r of roots) if (existsSync(r)) walk(r);
  return byPaper;
}

/** 1 本の論文を測る。 */
function measure(key, leanByPaper) {
  const got = readSource(key);
  if (got.err) return { key, ...got };
  const scan = scanSource(got.text);
  const { rows: lineRows, seen } = mapItems(scan);
  const perItem = new Map();
  const perHedge = new Map();
  let hedgeTotal = 0, attributed = 0;
  for (let i = 0; i < scan.lines.length; i++) {
    const L = scan.lines[i];
    const info = lineRows[i];
    if (!info || info.heading) continue;
    if (PAGE_RE.test(L)) continue;
    for (const h of HEDGES) {
      if (h.re.test(L)) {
        perHedge.set(h.key, (perHedge.get(h.key) ?? 0) + 1);
        hedgeTotal += 1;
        if (info.item) {
          if (!perItem.has(info.item)) perItem.set(info.item, { page: seen.get(info.item), hits: [], total: 0 });
          const rec = perItem.get(info.item);
          rec.hits.push({ key: h.key, line: i + 1, page: info.page });
          rec.total += 1;
          attributed += 1;
        }
      }
    }
  }
  const lean = leanByPaper.get(key) ?? new Map();
  const state = (k) => {
    const li = lean.get(k);
    return li?.bare ? '済' : li?.qualified ? `条つき${li.qualified}` : '未';
  };
  // ★第 2 の水路。傍注を行から項目へ帰属させる(合図の数え方には一切触らない)。
  const A = scanAsides(got.text);
  const perAside = new Map();
  let asideTotal = 0, inferTotal = 0, asideAttr = 0;
  for (const a of A.list) {
    asideTotal += 1;
    if (a.kind === '推論') inferTotal += 1;
    const info = lineRows[a.line - 1];
    const it = info?.item;
    if (!it) continue;
    if (!perAside.has(it)) perAside.set(it, []);
    perAside.get(it).push({ ...a, page: info.page });
    asideAttr += 1;
  }
  const kb = Math.max(1, got.text.length) / 1000;
  const asidePerKb = asideTotal / kb;

  const rows = [...seen.keys()].map((k) => {
    const v = perItem.get(k) ?? { page: seen.get(k), hits: [], total: 0 };
    const as = perAside.get(k) ?? [];
    return {
      item: k, page: seen.get(k), total: v.total, state: state(k), hits: v.hits,
      asides: as, aside: as.length, infer: as.filter((a) => a.kind === '推論').length,
    };
  });
  return {
    key, scan, seen, rows, perHedge, hedgeTotal, attributed, lean,
    headings: seen.size, runIn: scan.runIn, pageMode: scan.marked ? '[page N]' : '\\f',
    asideTotal, inferTotal, asideAttr, asidePerKb, asideDropped: A.dropped,
    asideStyle: asidePerKb >= ASIDE_STYLE_MIN,
  };
}

// ───────────────────────── 口 ─────────────────────────

if (args.includes('--papers')) {
  console.log(`★使える論文の鍵(ResearchPaper/papers.json、${PAPERS.size} 本)\n`);
  for (const [k, v] of PAPERS) {
    const ok = existsSync(path.join(SRC_DIR, v.txt));
    console.log(`   ${ok ? ' ' : '×'} ${k.padEnd(11)} ${v.title ?? v.file}`);
  }
  console.log('\n   × = 0_Source に .txt が無い(pdftotext -layout で作る。M11 を見よ)');
  process.exit(0);
}

if (args.includes('--all')) {
  const leanByPaper = scanLean();
  const out = [];
  for (const k of PAPERS.keys()) out.push(measure(k, leanByPaper));
  if (asJson) {
    console.log(JSON.stringify(out.map((r) => r.err ? { key: r.key, err: r.err } : {
      key: r.key, headings: r.headings, hedges: r.hedgeTotal, attributed: r.attributed,
      pageMode: r.pageMode, runIn: r.runIn, leanItems: r.lean.size,
      asides: r.asideTotal, infer: r.inferTotal, asidePerKb: +r.asidePerKb.toFixed(2),
      style: r.asideStyle ? 'aside' : 'word',
    }), null, 1));
    process.exit(0);
  }
  console.log('★全論文の被覆(node tools/hedge-index.mjs --all)\n');
  console.log('   論文        見出し    合図   帰属  帰属率   傍注  推論 傍注/KB  畳み方  頁の出所   書式   .src');
  const done = out.filter((r) => !r.err).sort((a, b) => b.asidePerKb - a.asidePerKb);
  for (const r of done) {
    const rate = r.hedgeTotal ? `${((100 * r.attributed) / r.hedgeTotal).toFixed(0)}%` : '—';
    console.log(
      `   ${r.key.padEnd(11)}${String(r.headings).padStart(5)}${String(r.hedgeTotal).padStart(8)}` +
      `${String(r.attributed).padStart(7)}${rate.padStart(7)}${String(r.asideTotal).padStart(7)}` +
      `${String(r.inferTotal).padStart(6)}${r.asidePerKb.toFixed(2).padStart(8)}` +
      `  ${(r.asideStyle ? '傍注式' : '語式').padEnd(7)}${r.pageMode.padEnd(9)}` +
      `${(r.runIn ? '走り込み' : '行独立').padEnd(7)}${String(r.lean.size).padStart(4)}`
    );
  }
  for (const r of out.filter((r) => r.err)) {
    console.log(`   ${r.key.padEnd(11)}   —— ${r.err === 'missing' ? `.txt が無い(${r.want})` : '未登記'}`);
  }
  const T = done.reduce((a, r) => ({
    h: a.h + r.hedgeTotal, at: a.at + r.attributed, hd: a.hd + r.headings,
    as: a.as + r.asideTotal, inf: a.inf + r.inferTotal,
  }), { h: 0, at: 0, hd: 0, as: 0, inf: 0 });
  console.log(`\n   計 ${done.length} 本 / 見出し ${T.hd} / 合図 ${T.h} / 帰属 ${T.at}` +
    ` (${((100 * T.at) / T.h).toFixed(1)}%) / 傍注 ${T.as}(うち推論 ${T.inf})`);
  const st = done.filter((r) => r.asideStyle);
  console.log(`   ★畳み方が「傍注式」の論文: ${st.length} 本 —— ${st.map((r) => r.key).join(' ')}`);
  console.log('     この ${n} 本では**既定 8 語だけでは下界を取り落とす**(傍注/KB 1.87 対 既定/KB 0.60)。'
    .replace('${n}', String(st.length)));
  console.log('   ★合図 0 の論文は「省略が無い」ではなく「合図の語が英語なので測れていない」');
  console.log('     ことがある(仏語の原典 —— Asterisque / EGA / Del)。');
  process.exit(0);
}

const leanByPaper = scanLean();
const M = measure(paper, leanByPaper);
if (M.err === 'unknown') {
  console.log(`(papers.json に「${paper}」は無い。鍵の一覧は node tools/hedge-index.mjs --papers)`);
  process.exit(0);
}
if (M.err === 'missing') {
  console.log(`(原文 txt が無いので測れない: ${M.want})`);
  process.exit(0);
}

const lines = M.scan.lines;
const allRows = M.rows;
const rows = allRows.filter((r) => !onlyItem || r.item === onlyItem);

if (asJson) {
  console.log(JSON.stringify({
    paper, pageMode: M.pageMode, runIn: M.runIn, headings: M.headings,
    hedgeTotal: M.hedgeTotal, attributed: M.attributed,
    asideTotal: M.asideTotal, inferTotal: M.inferTotal,
    asidePerKb: +M.asidePerKb.toFixed(2), style: M.asideStyle ? 'aside' : 'word',
    perHedge: Object.fromEntries(M.perHedge),
    rows: rows.filter((r) => r.total > 0 || r.aside > 0 || onlyItem),
  }, null, 1));
  process.exit(0);
}

// ★`--item` は「見つからない」と「合図 0 件」を区別する。以前はどちらも空表だった
//   ——CLAUDE.md が「着手前に必ず数える」と言う道具が、黙って何も言わないのが穴だった。
if (onlyItem && rows.length === 0) {
  const near = allRows
    .filter((r) => r.item.split(' ')[1] === onlyItem.split(' ')[1] || r.item.startsWith(onlyItem.split(' ')[0]))
    .slice(0, 8).map((r) => r.item);
  console.log(`★[${paper}] 「${onlyItem}」という見出しは原文に見つからない。`);
  console.log(`   (${M.headings} 件の見出しを ${M.runIn ? '走り込み' : '行独立'}式で拾っている。頁は ${M.pageMode})`);
  if (near.length) console.log(`   近いもの: ${near.join(' / ')}`);
  console.log('   ★「合図 0 件」ではなく「見出しを取り違えている」可能性がある。');
  console.log('     --papers で論文の鍵を、引数なしで見出しの一覧を確かめること。');
  process.exit(1);
}

console.log(`★[${paper}] 原文の「省略の合図」の地図`);
console.log(`   見出し ${M.headings} 件(${M.runIn ? '走り込み' : '行独立'}式)/ 頁の出所 ${M.pageMode}` +
  ` / 合図 ${M.hedgeTotal} 件中 ${M.attributed} 件を項目に帰属`);
console.log(`   畳み方 **${M.asideStyle ? '傍注式' : '語式'}**(傍注/KB ${M.asidePerKb.toFixed(2)}、閾 ${ASIDE_STYLE_MIN})`
  + ` / 角括弧の傍注 ${M.asideTotal} 件(うち推論 ${M.inferTotal} 件)、${M.asideAttr} 件を項目に帰属`);
if (M.asideStyle) {
  console.log('   ★この論文は**語ではなく角括弧で畳む**。既定 8 語だけを数えると下界を取り落とす。');
}
console.log('');
console.log('-- 語ごとの件数(論文全体)');
for (const h of HEDGES) {
  const n = M.perHedge.get(h.key) ?? 0;
  if (n) console.log(`   ${h.key.padEnd(13)} ${String(n).padStart(4)}   ${h.note}`);
}
if (M.hedgeTotal === 0) console.log('   (0 件。合図の語は英語なので、仏語の原典では測れない)');

if (onlyItem) {
  for (const r of rows) {
    console.log(`\n-- ${r.item}(物理 p.${r.page ?? '?'}、状態 ${r.state})`);
    for (const h of r.hits) console.log(`   ${h.key.padEnd(13)} 行 ${h.line}  p.${h.page ?? '?'}`);
    // ★第 2 の水路。「既定 8 語が 0 件」を「原文が畳んでいない」と読ませないための欄。
    //   実測: [GenEll] Lemma 3.5 は既定 0 件だが推論傍注 1 件(3 主張を 1 文に畳んでいる)。
    if (r.aside) {
      console.log(`   ---- 角括弧の傍注 ${r.aside} 件(うち推論 ${r.infer} 件)`);
      for (const a of r.asides) {
        const s = a.body.length > 96 ? `${a.body.slice(0, 96)}…` : a.body;
        console.log(`   ${a.kind}  行 ${a.line}  p.${a.page ?? '?'}  [${s}]`);
      }
    }
    if (r.total === 0 && r.aside === 0) {
      console.log('   ★合図も傍注も 0 件。**この項目は原文が畳んでいない**(見つからないのではない)。');
    } else if (r.total === 0) {
      console.log(`   ★既定 8 語は 0 件だが、**傍注が ${r.aside} 件ある**。`);
      console.log(`     この論文は傍注/KB ${M.asidePerKb.toFixed(2)} の**${M.asideStyle ? '傍注式' : '語式'}**である`
        + `(閾 ${ASIDE_STYLE_MIN})。「原文が畳んでいない」と読まないこと。`);
    }
    if (r.infer) {
      console.log('   ★「推論」= 引用を抱えずに Thus / hence で結論だけ述べた傍注。');
      console.log('     **手順書が無いまま主張が畳まれている**ので、節点に起こす第一候補。');
    }
  }
} else {
  console.log('\n-- 項目ごと(合図の多い順、上位 25)');
  console.log('   状態  合図  物理p  項目                        内訳');
  for (const r of rows.filter((r) => r.total > 0).sort((a, b) => b.total - a.total).slice(0, 25)) {
    const brk = Object.entries(
      r.hits.reduce((m, h) => ({ ...m, [h.key]: (m[h.key] ?? 0) + 1 }), {})
    ).map(([k, n]) => `${k}×${n}`).join(' ');
    console.log(
      `   ${r.state.padEnd(6)}${String(r.total).padStart(3)}  ${String(r.page ?? '?').padStart(5)}  ${r.item.padEnd(26)}${brk}`
    );
  }
}

const withHedge = rows.filter((r) => r.total > 0);
const undone = withHedge.filter((r) => r.state === '未' || r.state.startsWith('条つき'));
const sum = undone.reduce((a, r) => a + r.total, 0);
console.log(`\n-- 未実装(条つき含む)の項目に残る合図: ${sum} 件 / ${undone.length} 項目`);
console.log(`   (合図を持つ項目 ${withHedge.length} / 見出し ${M.headings}。`
  + `合図 0 の見出しは「原文が畳んでいない」)`);
console.log('   ★これが「まだ開いていない省略」の下界である(合図の無い省略は写らない)。');

// ★--cite: 合図の文が抱えている「[cf. …]」と番号つき参照を、候補の依存として出す。
//   ★経験則: **合図の文に括弧の引用があれば、それが手順書である**。
//   引用が無い合図は、公理の向きが合っていない可能性を先に疑うこと。
if (args.includes('--cite')) {
  const REF = /\b(Proposition|Theorem|Definition|Corollary|Lemma|Example|Remark)s?\s+([0-9]+\.[0-9]+(?:\.[0-9]+)?)/g;
  console.log('\n-- 合図の文が抱えている引用(候補の依存)');
  for (const r of withHedge) {
    if (onlyItem === null && r.state === '済') continue;
    const seenRef = new Set();
    let noCite = 0;
    for (const h of r.hits) {
      const ctx = lines.slice(Math.max(0, h.line - 3), h.line + 2).join(' ');
      const refs = [...ctx.matchAll(REF)].map((m) => `${m[1]} ${m[2]}`);
      if (refs.length === 0) noCite += 1;
      for (const x of refs) seenRef.add(x);
    }
    console.log(`   ${r.item}(合図 ${r.total}、うち引用なし ${noCite})`);
    if (seenRef.size) console.log(`      → ${[...seenRef].join(' / ')}`);
    if (noCite) console.log('      ★引用の無い合図がある —— 公理の向きを先に確かめること');
  }
}

// ★傍注の側の残工数。合図の側とは別に出す(混ぜると過去の数字と比べられなくなる)。
if (M.asideStyle) {
  const un = rows.filter((r) => r.aside > 0 && (r.state === '未' || r.state.startsWith('条つき')));
  const sa = un.reduce((a, r) => a + r.aside, 0);
  const si = un.reduce((a, r) => a + r.infer, 0);
  console.log(`-- 未実装(条つき含む)の項目に残る傍注: ${sa} 件 / ${un.length} 項目(うち推論 ${si} 件)`);
  console.log('   ★推論傍注 = 引用を抱えない傍注。手順書が無いまま畳まれているので、まずここを開く。');
}

console.log('\n★使い方: 項目に着手する前に `--item "Proposition 1.10"` で内訳を見て、');
console.log('  合図 1 つを分解(frdi-decomposition.json)の節点 1 つに対応させると、');
console.log('  「畳まれた量」を最初に見積もれる。全論文の被覆は `--all`。');
console.log('  ★「傍注式」の論文(後期 Mochizuki 18 本)では**合図 0 件でも畳んでいる**。');
console.log('    --item の「角括弧の傍注」欄を必ず見ること。');
