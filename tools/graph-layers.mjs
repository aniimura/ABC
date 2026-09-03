// 依存グラフを **層状のボックス図** として書き出す(自己完結 HTML、外部ライブラリ無し)。
//
// ★このレイアウトが答える問い: 「どこから手を付ければよいか」
//   右端 = [IUTchIII] Corollary 3.12。左へ行くほど「先に要るもの」。
//   **左端の層は依存を持たない**ので、そこから始められる。
//
// ★循環の扱い(ここが設計の核心)
//   循環があると層に並べられないので、**強連結成分(SCC)に潰す**。
//   潰した結果が凝縮 DAG で、これは必ず層に並ぶ。
//   ★そして SCC は「まとめてしか扱えない塊」なので、**そのままグルーピングの単位**になる。
//
//   ★★訂正の記録(2026-08-15): 当初この欄には
//   「824 節点 → 最大 SCC 262 節点(IUTchI/II/III/IV にまたがる)。
//    IUT 本体は1つの循環であり、『IUTchI から順に』という順序付けは存在しない」
//   と書いていた。**誤りである。**
//   その 262 節点の循環は、辺の定義が「本文が名指しした」だったことの副作用だった。
//   辺を「証明が依拠している」に定め直した実測(ResearchPaper/cycle-analysis.md):
//     節点 824 → 653 / 最大 SCC 262 → **17** / 循環する成分 30 → 12 / 層 19 → **55**
//   ★残る最大 SCC 17 は IUTchIII + IUTchIV のみで、**IUTchI・IUTchII は循環から外れる。**
//   すなわち **IUTchI → II → III → IV の順序付けは存在する。**
//
// ★ボックスの決め方
//   ・サイズ>1 の SCC は、それ自体で1ボックス(分割できないので)
//   ・残りの単独節点は (層, 論文) でまとめて1ボックス
//     = 「依存関係が少ないものはグルーピング」
//
// ★分母の限界: 原文側は pdftotext 由来で目視していない。番号の無い依存を数え落とし、
//   「the discussion of ...」型の案内を数え過ぎる。**下界でも上界でもない。**
//
// 使い方: node tools/graph-layers.mjs [出力先.html]

import { readFileSync, writeFileSync, existsSync, readdirSync, statSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const SRC = join(ROOT, 'ResearchPaper', '0_Source');
const LEAN = join(ROOT, 'lean', 'ABC3');
const OUT = process.argv[2] ?? join(ROOT, 'dependency-graph.html');
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';
const FILE_OF = JSON.parse(readFileSync(join(ROOT, 'tools', '_fileof.json'), 'utf8'));

// ★引用キーは論文ごとに独立である(2026-08-15 実測)。
//   IUT の4本は [IUTchI] のような記号的キーだが、古い論文は [Mzk3] のような数字型で、
//   しかも SemiAnbd の [Mzk4] と EtTh の [Mzk4] は別物でありうる。
//   各論文の参考文献欄を解いて解決する。tools/bibmap.mjs を参照。
//   ★実害だったもの: anabelioid の定義元([Mzk4] = The Geometry of Anabelioids)が
//    グラフに1度も現れていなかった。
//   ★★衝突に注意: IUTchI の [pGC] は "The Local Pro-p Anabelian Geometry of Curves" で、
//    我々がタグ pGC に割り当てた "A Version of the Grothendieck Conjecture for
//    p-adic Local Fields" とは**別の論文**である。論文ごとの解決はこれも正す。
const { buildBibs } = await import('./bibmap.mjs');
const BIBS = buildBibs(SRC, FILE_OF);
const FILE_TO_TAG = new Map(Object.entries(FILE_OF).map(([t, f]) => [f, t]));
/** 引用キーを、引用している論文の参考文献欄で解決する。 */
function resolveTag(fromTag, key) {
  const b = BIBS.get(fromTag);
  const v = b?.get(key);
  if (v && !v.startsWith('FILE:')) return v;          // 我々のタグに解決
  if (v) {                                            // 0_Source にあるがタグ未登録 → 題名をタグ代わりに
    const file = v.slice(5);
    return FILE_TO_TAG.get(file) ?? `@${file}`;
  }
  return FILE_OF[key] ? key : null;                   // 記号的キーがそのまま我々のタグなら採る
}
// `@題名` 形式のタグも読めるようにする
// ★§0: 接頭辞の節点は**出辺を持たない**(真の葉)。load が null を返せばそうなる。
const fileOfTag = (tag) => (tag.startsWith('§0:') ? null : (tag.startsWith('@') ? tag.slice(1) : FILE_OF[tag]));

// ── 原文側の抽出 ──────────────────────────────────────────
const loaded = new Map();
function load(tag) {
  if (loaded.has(tag)) return loaded.get(tag);
  const f = fileOfTag(tag); const p = f ? join(SRC, `${f}.txt`) : null;
  if (!p || !existsSync(p)) { loaded.set(tag, null); return null; }
  const lines = readFileSync(p, 'utf8').split(/\r?\n/);
  const pageOf = new Array(lines.length).fill(0);
  { let pg = 0; for (let i = 0; i < lines.length; i++) {
      const m = lines[i].match(/^===== \[page (\d+)\]/); if (m) pg = Number(m[1]); pageOf[i] = pg; } }
  const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\.[ \\t]*(?:$|\\()`);
  const list = [];
  for (let i = 0; i < lines.length; i++) { const m = lines[i].match(declRe); if (m) list.push({ name: `${m[1]} ${m[2]}`, line: i, page: pageOf[i] }); }
  for (let i = 0; i < list.length; i++) list[i].end = i + 1 < list.length ? list[i + 1].line : lines.length;
  const decls = new Map(); for (const d of list) if (!decls.has(d.name)) decls.set(d.name, d);
  const o = { decls, lines }; loaded.set(tag, o); return o;
}
// ★辺の意味(2026-08-15 に定めた): 「原文の証明が実際に依拠しているもの」。
//   次は辺ではない —— 前方参照(cf. …, below)/ 解説への案内(the discussion of …)。
//   実測(ResearchPaper/cycle-analysis.md): この2つを落とすだけで
//   最大 SCC が 262 → 16 になり、**論文をまたぐ循環が消える**。
//   ★Remark を種別で落とすことはしない。[IUTchIII] Remark 3.9.5, (i) は
//    正則包を**定義**しており、我々は実際にそれに依拠している。
const NOT_A_DEP = /\bbelow\b|the discussion of|the discussion surrounding|\bsee also\b/i;

function outs(tag, name) {
  const P = load(tag); if (!P) return null;
  const d = P.decls.get(name); if (!d) return null;
  const body = P.lines.slice(d.line, d.end).join('\n');
  const ctx = (i, len) => body.slice(Math.max(0, i - 90), i + len + 60).replace(/\s+/g, ' ');
  const es = new Map(); const spans = [];
  const keep = (to, c) => {
    const soft = NOT_A_DEP.test(c);
    const cur = es.get(to);
    // 同じ先が複数回出るなら、1回でも「依拠」の文脈があれば辺として残す
    if (cur === undefined || (cur && !soft)) es.set(to, soft);
  };
  // ★タグに数字を許すこと(2026-08-15 修正)。FrdI は [Mzk15] の形で引用するので、
  //   [A-Za-z]+ だと masking し損ね、直後の同一論文内スキャンが偽の辺を作る。
  //   実測: FrdI に偽の辺 Definition 1.2 → Definition 3.1 が生まれ、偽の循環に見えていた。
  const xre = new RegExp(`\\[([A-Za-z][A-Za-z0-9]*)\\],?\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(xre)) {
    spans.push([m.index, m.index + m[0].length]);
    // ★引用キーは引用元の参考文献欄で解決する。解決できないものは辺にしない
    //   (0_Source に無い外部文献、および題名照合が 60% を満たさなかったもの)。
    //   ★誤った辺は欠けた辺より悪い、という判断による。
    const rt = resolveTag(tag, m[1]);
    if (!rt) continue;
    keep(`${rt} / ${m[2]} ${m[3]}`, ctx(m.index, m[0].length));
  }
  const ire = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(ire)) {
    if (spans.some(([a, b]) => m.index >= a && m.index < b)) continue;
    const to = `${m[1]} ${m[2]}`;
    if (to !== name && P.decls.has(to)) keep(`${tag} / ${to}`, ctx(m.index, m[0].length));
  }
  return { edges: [...es].filter(([, soft]) => !soft).map(([k]) => k), page: d.page,
           dropped: [...es].filter(([, soft]) => soft).length };
}

const ROOTK = 'IUTchIII / Corollary 3.12';
const adj = new Map(), page = new Map();
{
  const q = [ROOTK], seen = new Set();
  while (q.length) {
    const k = q.shift(); if (seen.has(k)) continue; seen.add(k);
    const [t, n] = k.split(' / ');
    const r = outs(t, n);
    adj.set(k, r ? r.edges : []); page.set(k, r ? r.page : 0);
    for (const e of (r ? r.edges : [])) if (!seen.has(e)) q.push(e);
  }
  for (const k of seen) if (!adj.has(k)) { adj.set(k, []); page.set(k, 0); }
}

// ── ★§0(Notations and Conventions)を節点として合流させる ────────────
//   層 0 に中層の語彙が並ぶのは「引用が無い」からで「語彙が要らない」からではなかった。
//   実測: [FrdI] Definition 1.1 は integral / saturated / characteristic type を
//   すべて `[cf. §0]` で指している。§0 は番号付き項目でないのでグラフに無かった。
//   詳細は tools/section0.mjs。★§0 語は出辺を持たない = **真の葉**になる。
{
  const S0 = await import('./section0.mjs');
  const terms = new Map();     // tag -> Map(term -> {page})
  for (const tag of Object.keys(FILE_OF)) {
    const P = load(tag); if (!P) continue;
    const pageOf = [];
    { let pg = 0; for (let i = 0; i < P.lines.length; i++) {
        const m = P.lines[i].match(/^===== \[page (\d+)\]/); if (m) pg = Number(m[1]); pageOf[i] = pg; } }
    const t = S0.section0Terms(P.lines, pageOf);
    if (t.size) terms.set(tag, t);
  }
  let added = 0, edged = 0;
  for (const k of [...adj.keys()]) {
    const [tag, name] = k.split(' / ');
    const t = terms.get(tag); if (!t) continue;
    const P = load(tag); if (!P) continue;
    const d = P.decls.get(name); if (!d) continue;
    const body = P.lines.slice(d.line, d.end).join('\n');
    for (const term of S0.section0Refs(body, t)) {
      const nk = `§0:${tag} / ${term}`;
      if (!adj.has(nk)) { adj.set(nk, []); page.set(nk, t.get(term)?.page ?? 0); added++; }
      if (!adj.get(k).includes(nk)) { adj.get(k).push(nk); edged++; }
    }
  }
  console.log(`  §0 を合流: ${added} 節点(${terms.size} 論文)/ 項目→§0 の辺 ${edged}`);
}

// ── ★概念(番号を持たない語彙)を節点として合流させる ─────────────
//   目的は工数の把握。番号付き項目だけでは中層が出ないので、tools/concepts.json の
//   語彙を節点にし、定義箇所から辺を張る。詳細は tools/concept-graph.mjs。
{
  const CG = await import('./concept-graph.mjs');
  const tags = Object.keys(FILE_OF);
  const { nodes: cn, edges: ce, droppedByOrder, concepts, key: ckey } =
    CG.conceptNodes({ load, tags, NOT_A_DEP, conceptsPath: join(ROOT, 'tools', 'concepts.json'),
                      resolveTag, KIND });
  // 概念を節点として登録(定義が見つからなかったものも「未特定」として置く)
  for (const [k, v] of cn) { if (!adj.has(k)) adj.set(k, []); page.set(k, v.page ?? 0); }
  const ie = CG.itemToConceptEdges({ load, itemKeys: [...adj.keys()].filter((k) => !k.startsWith('CONCEPT / ')),
                                     concepts, key: ckey });
  for (const [a, b] of [...ce, ...ie]) {
    if (!adj.has(a)) continue;
    if (!adj.has(b)) continue;                 // 到達していない項目へは張らない
    if (!adj.get(a).includes(b)) adj.get(a).push(b);
  }
  console.log(`  概念を合流: ${cn.size} 節点(定義未特定 ${[...cn.values()].filter((v) => v.missing).length}、近似 ${[...cn.values()].filter((v) => v.approx).length})`);
  console.log(`    概念どうしの辺 ${ce.length}(導入順に反するとして落とした ${droppedByOrder.length})/ 項目→概念の辺 ${ie.length}`);
}

// ── ★★基礎理論(「理論のまとまり」)を節点として合流させる ──────────────
//   ★2026-08-16 追加。動機の実測:
//   [GenEll] の必要分 24 件はグラフ上 24 節点だが、その下には **Arakelov 理論**と
//   **l 捩れ点への Galois 表現**という「mathlib に 1 件も無い古典理論 2 本」が横たわっている。
//   ★それは番号付き項目でも §0 語でも概念語でもないので、既存の 3 層のどれにも出なかった——
//   合流前の本グラフを検索すると `Arakelov` 0 件 / `torsion` 0 件だった。
//   ★概念層が拾えないのは、原文に**定義文が無い**からである。原文は
//   "the Galois representation … associated to E_L" と**使う**だけで、定義文型では導入しない。
//
//   ★出辺を持たせない = **真の葉**にする。層 0(左端)に出るのが正しい——
//   「これを作らないと右へ進めない」ものだから。
//
//   ★限界: 登記 `ResearchPaper/foundations.json` は**人手**である。機械が
//   「この項目にはこの理論が要る」と判定しているのではない。**`neededBy` は下界**。
{
  const FJ = join(ROOT, 'ResearchPaper', 'foundations.json');
  if (existsSync(FJ)) {
    const F = JSON.parse(readFileSync(FJ, 'utf8'));
    const SHORT = { inMathlib: 'mathlibにある', partial: '枠組みが違う',
                    inProject: '公開プロジェクトあり', building: '構築中', absent: '無い' };
    let added = 0, edged = 0, unreached = 0;
    for (const f of F.foundations ?? []) {
      const nk = `基礎 / ${f.name} [${SHORT[f.status] ?? f.status}]`;
      if (!adj.has(nk)) { adj.set(nk, []); page.set(nk, 0); added++; }
      for (const n of f.neededBy ?? []) for (const it of n.items ?? []) {
        const k = `${n.paper} / ${it}`;
        if (!adj.has(k)) { unreached++; continue; }   // 到達していない項目へは張らない
        if (!adj.get(k).includes(nk)) { adj.get(k).push(nk); edged++; }
      }
    }
    console.log(`  基礎理論を合流: ${added} 節点 / 項目→基礎理論の辺 ${edged}` +
      `(到達していない項目 ${unreached} 件へは張らない)`);
  }
}

// ── ★★Interface の義務(構造体)を節点として合流させる ──────────────
//   ★2026-08-17 追加。動機の実測: グラフの節点語彙は**原文の項目番号**なので、
//   `GenEll Definition 1.1` に **8 本**の obligation がぶら下がっていても
//   **1 節点に潰れて見えなかった**(検索して 0 件を確認)。
//   ★check.mjs の「Interface 実装待ち」と**同じ単位**を図にも出す。
//
//   ★基礎理論の合流と同じく**出辺を持たせない** = 真の葉にする。
//   層 0(左端)に出るのが正しい——「これを作らないと右へ進めない」ものだから。
//
//   ★限界: `.src` を持たない構造体(`.waiting` だけのもの)は**錨が無い**ので
//   孤立した葉になる。それも数には入れる。
{
  const wk = (d, a = []) => {
    for (const f of readdirSync(d)) {
      const q = join(d, f);
      if (statSync(q).isDirectory()) wk(q, a); else if (q.endsWith('.lean')) a.push(q);
    }
    return a;
  };
  const files = wk(LEAN);
  const witness = new Set();
  for (const f of files) {
    const t = readFileSync(f, 'utf8');
    for (const m of t.matchAll(/^theorem\s+([A-Za-z_][\w'.]*)\.nonvacuous/gm))
      witness.add(m[1].split('.').pop());
  }
  let added = 0, edged = 0, done = 0, unanchored = 0;
  for (const f of files) {
    const rel = f.slice(ROOT.length + 1).replace(/\\/g, '/');
    if (!/\/Interface\//.test(rel)) continue;
    const t = readFileSync(f, 'utf8');
    for (const sm of t.matchAll(/^structure\s+([A-Za-z_][A-Za-z0-9_']*)\s/gm)) {
      const nm = sm[1];
      const ok = witness.has(nm);
      const nk = `義務 / ${nm} [${ok ? '埋まった' : '実装待ち'}]`;
      if (!adj.has(nk)) { adj.set(nk, []); page.set(nk, 0); added++; if (ok) done++; }
      const sm2 = new RegExp(`${nm}\\.src\\s*:[^=]*:=\\s*\\{([\\s\\S]{0,400}?)\\}`).exec(t);
      const pm = sm2 ? /paper\s*:=\s*"([^"]*)"/.exec(sm2[1]) : null;
      const im = sm2 ? /item\s*:=\s*"([^"]*)"/.exec(sm2[1]) : null;
      const km = (pm && im) ? new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`).exec(im[1]) : null;
      if (!km || !pm) { unanchored++; continue; }
      const k = `${pm[1]} / ${km[1]} ${km[2]}`;
      if (!adj.has(k)) { unanchored++; continue; }
      if (!adj.get(k).includes(nk)) { adj.get(k).push(nk); edged++; }
    }
  }
  console.log(`  Interface の義務を合流: ${added} 節点(埋まった ${done})/ 項目→義務の辺 ${edged}` +
    `(錨の無いもの ${unanchored} 件)`);
}

// ── ★★★Arakelov / Galois の義務の**木**を展開する ──────────────
//   ★2026-08-24 追加。動機の実測: Arakelov と Galois は `foundations.json` では
//   **1 行ずつ**しか無く、`義務 /` 節点も出辺を持たないので**孤立した葉**だった。
//   ★★実際は **Arakelov 9 obligation・76 条件**、**Galois 8 obligation・39 条件**の
//   **2 本の木**である(`ResearchPaper/obligation-tree.json`、2026-08-24 実測)。
//   ★★★ここで (a) 義務どうしの依存辺 (b) 未実装の条件の節点 を足し、
//   **本来の大きさ**で出す。
//
//   ★辺の向き: 「A が B を要する」= A → B(既存の合流と同じ向き)。
//   ★限界: 条件の節点は**個数だけ**が実測で、中身は Interface の
//   フィールド名に対応する(名前は出さない——長すぎるため)。
{
  const OJ = join(ROOT, 'ResearchPaper', 'obligation-tree.json');
  if (existsSync(OJ)) {
    const O = JSON.parse(readFileSync(OJ, 'utf8'));
    const byName = new Map((O.obligations ?? []).map((o) => [o.name, o]));
    const key = (o) => `義務 / ${o.name} [${o.fieldsDone >= o.fieldsTotal ? '埋まった' : '実装待ち'}]`;
    let edged = 0, fadded = 0;
    for (const o of O.obligations ?? []) {
      const nk = key(o);
      if (!adj.has(nk)) { adj.set(nk, []); page.set(nk, 0); }
      for (const dn of o.deps ?? []) {
        const d = byName.get(dn); if (!d) continue;
        const dk = key(d);
        if (!adj.has(dk)) { adj.set(dk, []); page.set(dk, 0); }
        if (!adj.get(nk).includes(dk)) { adj.get(nk).push(dk); edged++; }
      }
      const rem = Math.max(0, (o.fieldsTotal ?? 0) - (o.fieldsDone ?? 0));
      for (let i = 0; i < rem; i++) {
        // ★追分類(`D1` 等)は**論文ではない**ので、タグはトラック名までにする。
        //   ★以前は `Arakelov:D1 …` をタグにしていたため、論文一覧が分裂して見えていた。
        const fk = `${o.track} / ${o.code} ${o.ja} 条件 ${(o.fieldsDone ?? 0) + i + 1} of ${o.fieldsTotal}`;
        if (!adj.has(fk)) { adj.set(fk, []); page.set(fk, 0); fadded++; }
        if (!adj.get(nk).includes(fk)) { adj.get(nk).push(fk); }
      }
    }
    console.log(`  義務の木を展開: 義務→義務の辺 ${edged} / 未実装の条件の節点 ${fadded}`);
  }
}


// ★**現在作業中**の節点(分解チェーンの現在の葉 = 未着手で依存が全て済)。
//   ★ボックスの左上に ☆ を出すために使う。
const WIP = new Set();

// ── ★★★分解チェーンを**木のまま**展開する ──────────────
//   ★2026-08-18 追加。動機の実測: `ResearchPaper/frdi-decomposition.json` に
//   分解を登記しても、**このグラフはそれを読んでいなかった**。
//   ★そのため `Theorem 3.4`(実際は 18 節点・9 層)や
//   `Proposition 4.4`(13 節点・5 層)が **1 箱に潰れて見えていた**。
//   ★★★ここで展開し、**本来の大きさ**で出す。
//
//   ★辺の向き: 「A が B を要する」= A → B(既存の合流と同じ向き)。
//   ★奉じる項目 → チェーンの**終端節点**(他が依存していないもの)へ張る。
//   ★限界: 登記は**人手**であり、節点の粒度も我々の設計である。
//   機械が原典から抽出したものではない。**下界**として見ること。
{
  const DJ = join(ROOT, 'ResearchPaper', 'frdi-decomposition.json');
  if (existsSync(DJ)) {
    const D = JSON.parse(readFileSync(DJ, 'utf8'));
    let added = 0, edged = 0, unreached = 0, chains = 0;
    for (const c of D.chains ?? []) {
      const st = new Map((c.nodes ?? []).map((n) => [n.id, n.status]));
      // ★分解の節点は**奉じる論文の項目に帰属する**(2026-09-03)。
      //   以前は `分解 / …` を独立した「論文」のように扱っていたので、
      //   論文一覧に 320 件の `分解` が並び、FrdI から切り離されて見えていた。
      const sp = c.serves?.paper ?? '分解';
      const si = c.serves?.item ?? c.id;
      const key = (id) => `${sp} / ${si} \u25b8 ${c.id}/${id} [${st.get(id) === 'done' ? '済' : '未'}]`;
      const pg = c.serves?.page ?? 0;
      for (const n of c.nodes ?? []) {
        const nk = key(n.id);
        if (!adj.has(nk)) { adj.set(nk, []); page.set(nk, pg); added++; }
      }
      // ★現在の葉(未着手かつ依存が全て済)を WIP に入れる
      {
        const doneIds = new Set((c.nodes ?? []).filter((n) => n.status === 'done').map((n) => n.id));
        for (const n of c.nodes ?? []) {
          if (n.status === 'done') continue;
          if ((n.deps ?? []).every((d) => !st.has(d) || doneIds.has(d))) WIP.add(key(n.id));
        }
      }
      for (const n of c.nodes ?? []) for (const dep of n.deps ?? []) {
        if (!st.has(dep)) continue;
        const a = key(n.id), b = key(dep);
        if (!adj.get(a).includes(b)) { adj.get(a).push(b); edged++; }
      }
      const depended = new Set((c.nodes ?? []).flatMap((n) => n.deps ?? []));
      const tops = (c.nodes ?? []).filter((n) => !depended.has(n.id)).map((n) => key(n.id));
      // ★条つきの項目名(例「Proposition 4.4, (ii)」)はグラフの節点語彙ではないので、番号までで切る。
      const base = c.serves ? c.serves.item.split(',')[0].trim() : null;
      const sk = base ? `${c.serves.paper} / ${base}` : null;
      if (sk && adj.has(sk)) {
        for (const t of tops) if (!adj.get(sk).includes(t)) { adj.get(sk).push(t); edged++; }
        chains++;
      } else unreached++;
    }
    console.log(`  分解チェーンを展開: ${added} 節点 / 辺 ${edged}` +
      `(連結したチェーン ${chains} 本、奉じる項目が未到達 ${unreached} 件)`);
  }
}


// ── Tarjan SCC(反復版) ────────────────────────────────────
const index = new Map(), low = new Map(), onstk = new Set(), stk = [], comp = new Map();
let idx = 0, nc = 0;
for (const s of adj.keys()) {
  if (index.has(s)) continue;
  const work = [[s, 0]];
  while (work.length) {
    const fr = work[work.length - 1], v = fr[0];
    if (fr[1] === 0) { index.set(v, idx); low.set(v, idx); idx++; stk.push(v); onstk.add(v); }
    let rec = false;
    const es = adj.get(v) ?? [];
    for (let i = fr[1]; i < es.length; i++) {
      const w = es[i]; if (!adj.has(w)) continue;
      if (!index.has(w)) { fr[1] = i + 1; work.push([w, 0]); rec = true; break; }
      else if (onstk.has(w)) low.set(v, Math.min(low.get(v), index.get(w)));
    }
    if (rec) continue;
    if (low.get(v) === index.get(v)) { const c = nc++; while (true) { const w = stk.pop(); onstk.delete(w); comp.set(w, c); if (w === v) break; } }
    work.pop();
    if (work.length) { const p = work[work.length - 1][0]; low.set(p, Math.min(low.get(p), low.get(v))); }
  }
}
const members = new Map();
for (const [k, c] of comp) { if (!members.has(c)) members.set(c, []); members.get(c).push(k); }

// 凝縮 DAG と層(左 0 = 依存なし)
const cadj = new Map();
for (const [v, es] of adj) {
  const cv = comp.get(v); if (!cadj.has(cv)) cadj.set(cv, new Set());
  for (const w of es) { const cw = comp.get(w); if (cw !== undefined && cw !== cv) cadj.get(cv).add(cw); }
}
for (const c of members.keys()) if (!cadj.has(c)) cadj.set(c, new Set());
const memo = new Map();
const layerOf = (c, st = new Set()) => {
  if (memo.has(c)) return memo.get(c);
  if (st.has(c)) return 0;
  st.add(c); let d = 0;
  for (const w of cadj.get(c) ?? []) d = Math.max(d, layerOf(w, st) + 1);
  st.delete(c); memo.set(c, d); return d;
};


// ── 我々の位置(3層) ───────────────────────────────────────
function walk(d, a = []) { for (const f of readdirSync(d)) { const p = join(d, f); if (statSync(p).isDirectory()) walk(p, a); else if (p.endsWith('.lean')) a.push(p); } return a; }
const ours = new Map();
{
  const SRC_RE = /paper\s*:=\s*"([^"]*)"[\s\S]{0,400}?item\s*:=\s*"([^"]*)"/g;
  const EDGE_RE = /\.otherPaper\s+"\[?([A-Za-z]+)\]?"\s+"([^"]*)"/g;
  const ITEM_RE = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`);
  // ★`done` = `Found/` に `.src` がある(その原典項目を**完全に**実装したという主張)。
  //   `tools/frdi-progress.mjs` の分子と同じ規則である。部分実装には `.src` を付けない。
  //
  // ★★2026-08-16 修正: **項目名が丸ごと `Kind N.M` のものだけを `done` にする。**
  //   それまでは `ITEM_RE` が先頭一致だったので `item := "Lemma 3.1, (i)"` が
  //   キー `Lemma 3.1` に潰れ、★**条が 1 つでもあれば命題全体が完了として描かれていた。**
  //   これは上の 2 行が宣言している規則そのものと食い違っていた。
  //   ★`tools/frdi-progress.mjs` は**同じ欠陥を同日に直している**が、こちらは直っておらず、
  //   [GenEll] `Lemma 3.1`(条つき 3 個、(iv) 未実装)が `done` と表示されて発覚した。
  //   条つきは `partial` にする——**触れてはいるが完了ではない**。
  //   テンプレートは `o === 'done'` だけを ■(完了)として扱うので、`partial` は ●(触れた)になる。
  const EXACT_RE = new RegExp(`^\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\s*$`);
  const rank = { named: 0, skeleton: 1, landed: 2, partial: 3, done: 4 };
  const put = (k, kind, file) => { const c = ours.get(k); if (!c || rank[kind] > rank[c.kind]) ours.set(k, { kind, file }); };
  for (const f of walk(LEAN)) {
    const rel = f.slice(ROOT.length + 1).replace(/\\/g, '/'), t = readFileSync(f, 'utf8');
    const landed = /\.inProject|\.inMathlib/.test(t);
    const inFound = rel.includes('/Found/');   // ★実装が置かれる場所
    for (const m of t.matchAll(SRC_RE)) {
      const im = ITEM_RE.exec(m[2]);
      if (!im) continue;
      const kind = inFound ? (EXACT_RE.test(m[2]) ? 'done' : 'partial') : (landed ? 'landed' : 'skeleton');
      put(`${m[1]} / ${im[1]} ${im[2]}`, kind, rel);
    }
    for (const m of t.matchAll(EDGE_RE)) { const im = ITEM_RE.exec(m[2]); if (im) put(`${m[1]} / ${im[1]} ${im[2]}`, 'named', rel); }
  }
}

// ★義務節点の状態を `ours` に載せる(色付けのため)
for (const k of adj.keys()) {
  if (!k.startsWith('義務 / ')) continue;
  ours.set(k, { kind: k.includes('[埋まった]') ? 'done' : 'skeleton', file: 'lean/ABC3/Interface/' });
}

/** 節点キー `FrdI / Proposition 5.5` から節番号 `5` を取る。
 *  ★ボックスの名前を `[論文] §N` の一つの形に揃えるために使う(2026-09-03)。 */
const secOf = (k) => {
  const m = /\s(\d+)(?:\.\d+)*$/.exec(k.split(' / ')[1] ?? '');
  return m ? m[1] : '—';
};
/** ボックスの名前を `[論文] §N` の一つの形にする。
 *  ★番号を持たない節点(概念・`§0` の語)もここで同じ形に寄せる。 */
const tagOf = (k) => {
  const t = k.split(' / ')[0];
  if (t === 'CONCEPT') return '語彙(番号なし)';
  const z = /^§0:(.+)$/.exec(t);          // `§0:AbsTopI` は「その論文の §0」である
  if (z) return `${z[1]} §0`;
  if (t.startsWith('義務')) return '義務(Interface)';
  if (t === '基礎') return '基礎(mathlib の在庫)';
  if (t === '分解') return '分解(段への割り付け)';
  return `${t} §${secOf(k)}`;
};

// ── ボックス化: SCC(>1) はそのまま / 単独節点は (層, 論文, 節) でまとめる ──
const boxes = new Map();  // id -> {layer, label, tags, items[], scc}
const boxOfNode = new Map();
for (const [c, mem] of members) {
  const L = layerOf(c);
  if (mem.length > 1) {
    const id = `scc${c}`;
    const tg = new Map();
    for (const m of mem) { const t = tagOf(m); tg.set(t, (tg.get(t) ?? 0) + 1); }
    boxes.set(id, { layer: L, tags: [...tg].sort((a, b) => b[1] - a[1]), items: mem, scc: true });
    for (const m of mem) boxOfNode.set(m, id);
  } else {
    // ★ボックスの表記を揃えるため、(層, 論文) ではなく **(層, 論文, 節)** でまとめる
    //   (2026-09-03)。こうするとどのボックスも `[論文] §N` の一つの形になる。
    const k = mem[0], lab = tagOf(k);
    const id = `L${L}:${lab}`;
    if (!boxes.has(id)) boxes.set(id, { layer: L, tags: [[lab, 0]], items: [], scc: false });
    boxes.get(id).items.push(k);
    boxes.get(id).tags[0][1]++;
    boxOfNode.set(k, id);
  }
}
// ボックス間の辺
const bedge = new Map();
for (const [v, es] of adj) {
  const a = boxOfNode.get(v); if (!a) continue;
  for (const w of es) { const b = boxOfNode.get(w); if (!b || a === b) continue; const k = `${a} ${b}`; bedge.set(k, (bedge.get(k) ?? 0) + 1); }
}

// ── 配置 ─────────────────────────────────────────────────
const BW = 150, GAPX = 58, GAPY = 10, PADT = 62;
const maxL = Math.max(...[...boxes.values()].map((b) => b.layer));
const byLayer = new Map();
for (const [id, b] of boxes) { if (!byLayer.has(b.layer)) byLayer.set(b.layer, []); byLayer.get(b.layer).push(id); }
let maxH = 0;
for (const [L, ids] of byLayer) {
  ids.sort((a, b) => boxes.get(b).items.length - boxes.get(a).items.length);
  let y = PADT;
  for (const id of ids) {
    const b = boxes.get(id);
    b.h = Math.max(28, Math.min(170, 24 + b.items.length * 1.9));
    b.w = BW;
    // ★層 0(依存なし)が最も左、根(最大層)が最も右。
    //   以前 (maxL - b.layer) と書いており **左右が反転していた**(2026-08-15 訂正)。
    b.x = b.layer * (BW + GAPX);
    b.y = y; y += b.h + GAPY;
  }
  maxH = Math.max(maxH, y);
}
// x を反転(層 0 = 左)。上で (maxL - layer) にしているので根が x 最大 = 右。
const B = [...boxes.entries()].map(([id, b]) => ({
  id, layer: b.layer, x: Math.round(b.x), y: Math.round(b.y), w: b.w, h: Math.round(b.h),
  tags: b.tags, n: b.items.length, scc: b.scc ? 1 : 0,
  ours: b.items.filter((k) => ours.has(k)).length,
  done: b.items.filter((k) => ours.get(k)?.kind === 'done').length,   // ★Lean 形式化が完了した項目数
  landed: b.items.filter((k) => ours.get(k)?.kind === 'landed').length,
  root: b.items.includes(ROOTK) ? 1 : 0,
  wip: b.items.filter((k) => WIP.has(k)).length,   // ★現在作業中の節点数
  items: b.items.slice(0, 400).map((k) => ({ k, p: page.get(k) ?? 0, o: ours.get(k)?.kind ?? null })),
}));
const bidx = new Map(B.map((b, i) => [b.id, i]));
const BE = [...bedge.entries()].map(([k, w]) => { const [a, b] = k.split(' '); return [bidx.get(a), bidx.get(b), w]; })
  .filter(([a, b]) => a !== undefined && b !== undefined);

const layerStats = [...byLayer.keys()].sort((a, b) => a - b).map((L) => ({
  L, boxes: byLayer.get(L).length,
  nodes: byLayer.get(L).reduce((s, id) => s + boxes.get(id).items.length, 0),
  ours: byLayer.get(L).reduce((s, id) => s + boxes.get(id).items.filter((k) => ours.has(k)).length, 0),
}));

const data = {
  B, BE, layerStats, maxL,
  stats: { nodes: adj.size, sccs: nc, bigSccs: [...members.values()].filter((m) => m.length > 1).length,
    boxes: B.length, layers: maxL + 1,
    // ★被覆数から**義務節点を除く**——語彙が違う(原文の項目番号 対 我々の構造体)。
    //   混ぜると「触れた節点」が 28 だけ水増しされ、進捗の読みを狂わせる。
    ours: [...ours.keys()].filter((k) => adj.has(k) && !k.startsWith('義務 / ')).length,
    oblig: [...adj.keys()].filter((k) => k.startsWith('義務 / ')).length,
    obligDone: [...adj.keys()].filter((k) => k.startsWith('義務 / ') && k.includes('[埋まった]')).length },
  width: (maxL + 1) * (BW + GAPX), height: maxH + 40,
};
const html = readFileSync(join(ROOT, 'tools', 'graph-layers.template.html'), 'utf8').replace('/*__DATA__*/', JSON.stringify(data));
writeFileSync(OUT, html, 'utf8');
// ★統計を JSON でも出す——index.html がこれを読む。
// ★数字を手で書かないためである。手で書くと古びて、古びた数字は嘘になる。
writeFileSync(join(ROOT, 'ResearchPaper', 'graph-layers-stats.json'),
  JSON.stringify({ ...data.stats, maxLayer: maxL }, null, 2), 'utf8');
console.log(`書き出し: ${OUT}`);
console.log(`  節点 ${data.stats.nodes} / SCC ${data.stats.sccs}(サイズ>1 は ${data.stats.bigSccs})`);
console.log(`  ボックス ${data.stats.boxes} / 層 ${data.stats.layers}(左 0 = 依存なし、右 ${maxL} = 根)`);
console.log(`  我々が触れている節点 ${data.stats.ours}(★義務は含めない)`);
console.log(`  Interface の義務: ${data.stats.obligDone} / ${data.stats.oblig}`);
