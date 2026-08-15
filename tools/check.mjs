#!/usr/bin/env node
/**
 * check.mjs — 本プロジェクト唯一のゲート
 *
 * 検査するもの:
 *   [A] 1_Structured の構造化文書  … S1-S6(規約は 1_Structured/README.md §4)
 *   [B] Lean プロジェクト          … ビルド / sorry 台帳 / axiom 禁止(G4)
 *
 * 使い方:
 *   node tools/check.mjs                 全部
 *   node tools/check.mjs --structured    [A] のみ
 *   node tools/check.mjs --lean          [B] のみ
 *   node tools/check.mjs --selftest      器具自身の較正
 *
 * 終了コード: 0 = 全 PASS / 1 = NG あり
 *
 * ★原則: 器具は自分の欠陥を報告できなければならない。
 *   非ゲート指標(件数等)も必ず数値を印字し、「緑1語」で隠さない。
 *
 * ════════════════════════════════════════════════════════════════
 * ★★ この器具が検査して **いない** こと(必読)
 * ════════════════════════════════════════════════════════════════
 *
 * 「selftest が全部 PASS」は「欠陥が無い」ではない。**用意した項目を通る**
 * というだけで、用意していない失敗形については何も言っていない。
 * 2026-08-14 の監査で、20/20 PASS の状態で以下の2件が素通りしていた:
 *   - `.needs` が存在しない物理ページを指す(`.src` では検査していたのに)
 *   - `Found/` に `sorry` が残る(bucket の規則に明記してあるのに)
 * どちらも fixture を書いていなかったから見えなかった。
 *
 * 2026-08-14(同日、`ResidueCardinality` の discharge 時)にもう1件:
 *   - `Interface/` が `Found/` を import する(= Skeleton が実装に推移的に依存する)
 * これは 22/22 PASS の状態で素通りした。PLAN §3 の2トラック構成そのものを壊す向きなのに、
 * 規則がどこにも書かれておらず(暗黙の前提だった)、当然 fixture も無かった。
 * D23/D24 として機械化した。**「図に描いてあるから守られる」は成り立たない。**
 *
 * 2026-08-15、依存グラフを作る過程でさらに 3 件(いずれも 28/28 PASS の状態で素通り):
 *   - `.needs` の `otherPaper` は **別の論文**を指すのに、ページを **所有論文**の
 *     `pdfPages` と比べていた。つまり辺の先は事実上検査していなかった。
 *     現に `[IUTchI]` `[IUTchII]` が `papers.json` **未登記のまま**通っていた。
 *     → D26/D27/D28 として機械化(D28 が「所有論文では範囲内・辺の先では範囲外」の形)。
 *   - `.needs` の `page` の**意味が混在**していた。`cor_3_12.needs` の 6 本のうち
 *     `Theorem 3.11` だけが辺の先のページで、他 5 本は引用している側のページだった。
 *     三つ組が自己完結しないと辺の先は検査できない。`Meta/Claim.lean` に明記した。
 *   - `.needs` の本体抽出が `src.indexOf("foo.needs")` だったため、**docstring 中の
 *     散文の言及**を先に拾い、そこから最初の `[` までを本体と誤認していた。
 *     件数には数えられるのに中身は 1 件も集計されない、という壊れ方をする。
 *     → D29 として機械化。`.src` 側にも同じ脆さがあったので同時に直した。
 *
 * ★このうち 2 件目は「器具の穴」ではなく **データの穴**である。器具を足しても、
 *   書く側の規約が曖昧なままなら検査は空回りする。
 *
 * ── A. 原理的に検査できない(自己申告に依存する)
 *   A1. `data-notation-checked` の日付 — **実際に PDF を目視したかは検査不能**。
 *       形式(YYYY-MM-DD か "none")しか見ていない。
 *   A2. `LeanStatus` の `inMathlib` / `absent` — 実測したかは検査不能。
 *       ★★**「無い」と言う前の探索手順**(2026-08-14 制定。同じ失敗を2回してから決めた)。
 *       失敗の実例:
 *         1. `relatively compact` で grep して「Θ 側について 0 件」と報告した。
 *            原文は同じことを `compactness` という**別の語**で述べていた(IUTchIII p.175)。
 *         2. `subset` / `contains` で grep して「原文に無い」と報告した。
 *            原文は `contained in tensor packets of log-shells` と述べていた(同 p.131)。
 *       2回とも **Lean 側で付けた名前の語**で原文を探していた。原文は
 *       **原文側の語**(その対象を原文が何と呼んでいるか)を主語にして述べていた。
 *       ゆえに規則:
 *         (S1) **原文側の呼称を先に決める**。Lean の識別子ではなく、原文がその対象を
 *              何と呼んでいるか(例: “possible images”, “log-shell”, “holomorphic hull”)。
 *         (S2) その呼称の**全出現を列挙する**(`grep -n` の件数を数え、全部見る)。
 *              件数が多すぎて全部見られないなら、それは「無い」と言える状態ではない。
 *         (S3) そのうえで Lean 側の語(`subset`/`compact` 等)で絞る。
 *         (S4) `absent (searched)` には **(S1) の呼称と (S2) の件数**を書く。
 *              「`X` で grep して 0 件」だけでは不十分——`X` は我々の語かもしれない。
 *       ★この手順を踏んだかは**検査不能**(だから A 群にある)。書かせることしかできない。
 *   A3. `.needs` の内容が原文の証明文を反映しているか — 件数と形式しか見ない。
 *       原文が言っていない依存を書いても通る。
 *   A4. `.reading` に原文が混ざっていないか — `.verbatim` しか照合していない。
 *   A5. ★**statement が原典を忠実に表しているか** — Lean も check.mjs も判定しない。
 *       これが「事実1」の核心であり、この器具の存在理由でもある。器具が肩代わり
 *       できるのは**材料を揃えさせること**までで、判断そのものは残る。
 *   A7. ★★**docstring が「実在しない検査・宣言」を約束していないか**(2026-08-14 追加。
 *       **実際に1件すり抜けた**)。実例: `Found/IUTchIII/PadicLog.lean` の docstring に
 *       「(`logOneAdd = 0` という退化ではないことは、下の `Check/` で別に示す)」と
 *       書いたが、**その `Check/` は存在しなかった**。報告では「非退化性は未証明」と
 *       正直に書いていたのに、**コードの中には「別に示す」と書いてあった**。
 *       親セッションとの照合で発見。
 *       **機械化を検討し、実測したうえで却下した**:
 *         - 提案: docstring 内のバッククォート付き識別子が実在するか照合する。
 *         - 実測(2026-08-14): 識別子の形をしたトークン 602 件中 **444 件(74%)が
 *           プロジェクト内宣言に解決しない**。snake_case/dotted に絞っても
 *           196 件中 87 件(**44%**)。残りは正当なもの——mathlib の補題名
 *           (`WithTop.coe_ne_top` 等)・ファイル名(`check.mjs`)・数式の記号(`I_K`)。
 *         - ★**そもそもこの実例を捕まえられない**。約束は散文で書かれており、
 *           バッククォート内は `logOneAdd = 0`(式)と `Check/`(パス。しかも実在する)。
 *           **動機になった事例を捕まえない検査は、偽の安心そのもの**である。
 *         - ★★**決定打**: 本プロジェクトは「**訂正の記録を消さない**」を規律にしている。
 *           訂正記録は**削除した宣言の名前を必然的に含む**(実測: 未解決トークン上位に
 *           `thetaUnion_isCompact`・`logVol_hull_ne_top_of_isCompact`・
 *           `cor_3_12_refutable_under_current_interface` 等、意図的に削除したものが並ぶ)。
 *           「参照先はすべて実在せよ」という検査は、**この規律と正面から衝突する**。
 *       → 機械化しない。代わりに**書き方の規則**を置く:
 *         **docstring に「後で示す」「別に示す」類の約束を書かない。**
 *         書いてよいのは (a) 実在する宣言の名前、(b)「未証明である」の明示、のどちらかだけ。
 *         ★これを守ったかは検査不能(だから A 群にある)。
 *   A6. ★**`Interface` の仮説が原文より強くないか**(2026-08-14 追加)。
 *       G5 は結論の弱化を禁じるが、その鏡像である**仮説の強化**は誰も見ていない。
 *       実例: IUTchIII Cor 3.12 の `qLogVol_mem` は原文 p.184 の「The inclusion …」を
 *       写したものだが、原文が直前の文に付けている限定
 *       「subject to the condition」「perhaps only up to some sort of *approximation*」を
 *       落としている(意図的・`.needs` に `.implicitStep` で明示)。
 *       **仮説を強めれば結論はいくらでも証明できる**ので、これは「型が付くが
 *       何も言っていない」の第4の顔である。
 *       機械化できない理由: 「原文より強いか」の判定は逐語の**含意関係**の判定であり、
 *       逐語照合(S4)は部分文字列一致しか見ない。A5 と同じ壁。
 *       → 器具にできるのは、落とした限定を `.needs` に書かせることまで。
 *          書いたかどうかは**自己申告**である。
 *
 * ── B. 規則はあるが機械検査が無い
 *   B1. **G5(弱化禁止)** — 完全に規律のみ。証明できないときに原典の主張を
 *       導ける形へ弱めても、何も落ちない。
 *   B2. `waiting` に期限が無い — `Interface` は永久に waiting でいられる。
 *   B3. `*.legacy.*` は検査対象外のまま。再作成の期限が無い。
 *   B4. `loadBearing` の印は自己申告 — 付け忘れれば G3 は発火しない。
 *   B6. ★★**器具の逆インセンティブ**(2026-08-14 実測)。`Skeleton` の定理は
 *       `Interface` にフィールドを posit すれば証明できる。すると
 *       **見出しの `sorry` 件数は減り、実際に残っている仕事は増える**。
 *       実例: IUTchIII `cor_3_12` は原文 p.184 と p.175 の段を `Interface` へ
 *       輸入した結果 `sorry` が消えた(5 件 → 4 件)が、我々は1つも導出していない。
 *       「Interface 実装待ち」は `structure` 単位なので、フィールドが 5 増えても 2 のまま。
 *       **機械化の検討と却下**:
 *         - `Interface` の**フィールド数**を印字する案 → 却下。重みが違いすぎる
 *           (`qAbs : ℝ` と `thetaUnion_isCompact` が同じ 1)うえ、**束ねれば 1 に減る**
 *           ので、悪くする方向に指標が改善しうる。
 *         - 「証明が `Found/` の補題を使わない」で近似する案 → 却下。IUT 方向には
 *           `Found/` が存在しないので**恒真**になる。定数の出力は検査ではない(G7 の基準)。
 *         - 「証明が `Interface` の射影だけからなる」なら機械化できるが、
 *           mathlib の補題を1つ挟むだけで逃げられる**偽陰性**が容易。
 *       → 新しい指標は足さず、(a) `sorry` 件数の直後に「単独で読むな」と印字し、
 *          (b) 輸入した段を `.needs` の `.implicitStep` に書かせる(既存の
 *          「暗黙の段」件数が実際に 3 → 5 と動いた)。**ただし (b) は自己申告である。**
 *   B5. ★**G2 は退化 witness で満たせる、しかも名前で操作できる**(2026-08-14 実証)。
 *       `Nonempty (ResidueCardinality p)` は退化した実例で証明でき、Track B の
 *       実装を待たずに通る。さらに `X.nonvacuous` という**名前を付けるかどうか**で
 *       実装待ちキューに載るかが決まるので、退化 witness にその名を付ければ
 *       実装ゼロのままキューが空になる。現に `SubgroupCorrespondence` の
 *       `Nonempty` は証明済みだが、意図的にその名を付けていない(人手の判断)。
 *       → G2 の通過は**内容の保証ではない**。内容は G3(負の対照)の側で測る。
 *       この穴を機械で塞ぐ案は検討したが見送った: `Found/` に退化 witness を
 *       置かれる場合を機械で区別できず、**部分的な塞ぎが偽の安心を与える**ため。
 *
 * ── 運用規律
 *   **新しい失敗形を見つけたら、直すと同時に fixture を足す。**
 *   直すだけでは、次に同じ穴が開いたとき気づけない。
 */

import {
  readFileSync, readdirSync, statSync, existsSync,
  mkdtempSync, mkdirSync, writeFileSync, rmSync,
} from 'node:fs';
import { execFileSync } from 'node:child_process';
import { join, dirname, relative, basename } from 'node:path';
import { fileURLToPath } from 'node:url';
import { tmpdir } from 'node:os';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const SOURCE_DIR = join(ROOT, 'ResearchPaper', '0_Source');
const STRUCT_DIR = join(ROOT, 'ResearchPaper', '1_Structured');
const PAPERS_JSON = join(ROOT, 'ResearchPaper', 'papers.json');
/** 原文からの機械概算(規模の分母)。無くても動く。 */
const SCALE_JSON = join(ROOT, 'ResearchPaper', 'dependency-scale.json');
/** 葉の仕分け(自己申告)。無くても動く。 */
const BLOCKED_JSON = join(ROOT, 'ResearchPaper', 'blocked-leaves.json');
const LEAN_DIR = join(ROOT, 'lean');
const LEAN_SRC = join(LEAN_DIR, 'ABC3');

/** 較正デモだけは axiom を持つことを許す(それが実演の内容だから) */
const AXIOM_EXEMPT = ['Meta/Calibration.lean']; // lean/ABC3 からの相対

let NG = 0;
const ng = (where, msg) => { NG++; console.log(`  NG  ${where}\n      ${msg}`); };
const ok = (msg) => console.log(`  ok  ${msg}`);
const h1 = (s) => console.log(`\n=== ${s} ===`);

// ────────────────────────────────────────────────────────────────
// テキスト正規化 / HTML の最小パーサ
// ────────────────────────────────────────────────────────────────

const ENTITIES = {
  amp: '&', lt: '<', gt: '>', quot: '"', apos: "'", nbsp: ' ',
  ldquo: '“', rdquo: '”', lsquo: '‘', rsquo: '’',
  minus: '−', times: '×', ge: '≥', le: '≤',
  rarr: '→', larr: '←', harr: '↔',
  sube: '⊆', supe: '⊇', sub: '⊂', sup: '⊃',
  cong: '≅', or: '∨', and: '∧',
  Gamma: 'Γ', alpha: 'α', chi: 'χ', prime: '′',
  sect: '§', hellip: '…', mdash: '—', ndash: '–',
};

function decodeEntities(s) {
  return s
    .replace(/&#x([0-9a-fA-F]+);/g, (_, h) => String.fromCodePoint(parseInt(h, 16)))
    .replace(/&#(\d+);/g, (_, d) => String.fromCodePoint(parseInt(d, 10)))
    .replace(/&([a-zA-Z]+);/g, (m, name) => (name in ENTITIES ? ENTITIES[name] : m));
}

/**
 * 照合用の正規化。両側に同じものを適用する。
 *
 * ここで吸収するのは **レンダラ由来の差** だけ——合字・引用符・ダッシュ・
 * 記号の分解(pdftotext は ≅ を "∼=" や "=∼" として出力する)。
 * **論文由来の差**(装飾の有無)はここで吸収してはならない——それは
 * `data-txt` で明示するもので、自動で潰すと事実2の検出力が消える。
 */
function normalize(s) {
  return decodeEntities(s)
    .replace(/ﬀ/g, 'ff').replace(/ﬁ/g, 'fi').replace(/ﬂ/g, 'fl')
    .replace(/ﬃ/g, 'ffi').replace(/ﬄ/g, 'ffl')
    .replace(/­/g, '')
    .replace(/[‐-―]/g, '-')
    .replace(/[‘’]/g, "'").replace(/[“”]/g, '"')
    .replace(/∼=|=∼/g, '≅')
    .replace(/≝/g, 'd=ef')
    .replace(/\s+/g, ' ')
    .trim();
}

/**
 * Lean のコメント内で装飾を表す平文記法 `X[cls]` の照合射影。
 *
 * `1_Structured` の HTML クラスと **同じ語彙**を使う(1つの体系・2つの表現)。
 * 実測により、装飾は `pdftotext` で一律に脱落し **基底文字だけが残る**——
 * K̄→"K" / Ẑ→"Z" / 𝒪→"O" / **Z**→"Z" / Y̲̲→"Y" / K′→"K"。
 * ゆえに射影は「角括弧を落とす」だけでよい。
 *
 * `_`(下付き)・`^`(上付き)も落とす: 原文 "Γ_K" は pdftotext で "ΓK"。
 */
// `frak` は 2026-08-14 追加。IUTchIII は**書体だけで別対象を区別する**——
// 添字 MOD(ローマン)対 𝔪𝔬𝔡(フラクトゥール)、LGP 対 𝔩𝔤𝔭、log-volume の log 対
// 𝔩𝔬𝔤-link の 𝔩𝔬𝔤。原文 物理 p.156 は「MOD の方は両立するが 𝔪𝔬𝔡 の方は両立しない」と
// 明示的に対比しているので、この区別は装飾と同じ資格で対象の同定に効く。
// HTML 側のクラス表(1_Structured/README.md §3)と**同じ語彙**でなければならない。
const DECOR_CLASSES = 'ul1|ul2|bar|hat|tilde|dot1|dot2|bb|scr|frak|prime';
function leanQuoteProjection(text) {
  return squash(text
    .replace(new RegExp(`\\[(?:${DECOR_CLASSES})\\]`, 'g'), '')
    .replace(/[_^]/g, ''));
}

/**
 * 照合は空白を完全に無視して行う。
 *
 * 理由: pdftotext は上付き・下付きを別のテキスト run として出力するため、
 * 添字の直後に空白が入る(原文 "Γ_K." → 出力 "ΓK ." / "Γ_K-module" →
 * "ΓK -module")。この空白は組版の産物であって原文の一部ではないので、
 * 照合の対象にしない。空白を無視しても、数百文字の部分文字列一致という
 * 検査の識別力はほとんど落ちない。
 */
const squash = (s) => normalize(s).replace(/\s+/g, '');

/**
 * .verbatim の内側から「照合射影」を作る。
 * data-txt を持つ要素はその値に置換する(空文字なら削除)。それ以外はタグを剥がす。
 */
function matchProjection(html) {
  let s = html;
  const withTxt = /<(\w+)\b[^>]*\bdata-txt="([^"]*)"[^>]*>([\s\S]*?)<\/\1>/;
  for (let guard = 0; guard < 500 && withTxt.test(s); guard++) {
    s = s.replace(withTxt, (_m, _tag, txt) => txt);
  }
  s = s.replace(/<[^>]+>/g, '');
  return squash(s);
}

/** 極小の section 抽出。属性は key="value" 形式のみを想定する。 */
function extractSections(html) {
  const out = [];
  const re = /<section\b([^>]*)>([\s\S]*?)<\/section>/g;
  let m;
  while ((m = re.exec(html)) !== null) {
    const attrs = {};
    const ar = /([\w-]+)="([^"]*)"/g;
    let a;
    while ((a = ar.exec(m[1])) !== null) attrs[a[1]] = a[2];
    const vm = /<div class="verbatim">([\s\S]*?)<\/div>/.exec(m[2]);
    out.push({ attrs, verbatimHtml: vm ? vm[1] : null });
  }
  return out;
}

function walk(dir, acc = []) {
  if (!existsSync(dir)) return acc;
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) walk(p, acc);
    else acc.push(p);
  }
  return acc;
}

// ────────────────────────────────────────────────────────────────
// PDF ページのテキスト(キャッシュ付き)
// ────────────────────────────────────────────────────────────────

/**
 * pdftotext の抽出モード。
 *
 * 既定モードは、2段組でない本文でも **行の順序を入れ替える**ことがある
 * (実測: pGC 物理 p.3 の Corollary 1.3 が "Corollary 1.3: from ΓK ." の後に
 * 本体行が来る)。`-layout` は組版上の行順を保つのでこれを回避できる。
 * 逆に `-layout` が不利な版面もありうるため、3モードを順に試し、
 * **どのモードで一致したかを必ず印字する**(既定以外での一致は、その版面が
 * 素直でないことの情報になる)。
 */
const PDF_MODES = [
  ['layout', ['-layout']],
  ['default', []],
  ['raw', ['-raw']],
];

const pageCache = new Map();
/** @returns {Array<[string, string]>|null} [モード名, squash済みテキスト] の配列 */
function pdfPageTexts(pdfPath, page) {
  const key = `${pdfPath}#${page}`;
  if (pageCache.has(key)) return pageCache.get(key);
  const out = [];
  for (const [name, flags] of PDF_MODES) {
    try {
      out.push([name, squash(execFileSync('pdftotext',
        ['-enc', 'UTF-8', ...flags, '-f', String(page), '-l', String(page), pdfPath, '-'],
        { encoding: 'utf8', maxBuffer: 32 * 1024 * 1024 }))]);
    } catch { /* このモードは使えない */ }
  }
  const val = out.length ? out : null;
  pageCache.set(key, val);
  return val;
}

// ────────────────────────────────────────────────────────────────
// [A] 1_Structured の検査
// ────────────────────────────────────────────────────────────────

const DATE_RE = /^\d{4}-\d{2}-\d{2}$/;

/** `sorry` の検出。識別子の一部や `.sorry` は拾わない。 */
const SORRY_RE = /(^|[^\w.])sorry([^\w]|$)/;

/**
 * コメント(`/- ... -/` と `--` 行末)**と文字列リテラル**を、
 * 改行を保ったまま空白へ潰す。
 *
 * docstring 内で "sorry" に言及しただけのものを数えると台帳が信用できなくなる。
 * ★2026-08-15 追加: **文字列リテラルも潰す**。`.needs` の `LeanStatus.inProject` は
 * 「あちらの宣言は `sorry` 無しである」と**文字列で**書く場所なので、
 * 潰さないと台帳が誤って +1 される(実際 `Skeleton/AbsTopIII/LogShell.lean` で起きた)。
 * 文字列の中の `sorry` が本物の `sorry` であることは無いので、潰して安全。
 */
const stripLeanComments = (src) => src
  .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
  .replace(/--[^\n]*/g, (m) => ' '.repeat(m.length))
  .replace(/"(?:[^"\\\n]|\\.)*"/g, (m) => m.replace(/[^\n]/g, ' '));

/** `X.src : Source` の中身を読む正準形パターン(1行で書くこと)。 */
const SRC_RE =
  /\.src\b[\s\S]{0,400}?paper\s*:=\s*"([^"]*)"[\s\S]{0,300}?pdfPage\s*:=\s*(\d+)[\s\S]{0,300}?sectionId\s*:=\s*"([^"]*)"/;

function checkStructured({ files = null, papersPath = PAPERS_JSON, quiet = false } = {}) {
  const before = NG;
  const reg = JSON.parse(readFileSync(papersPath, 'utf8')).papers;

  const targets = files ?? walk(STRUCT_DIR)
    .filter((p) => p.endsWith('.html'))
    .filter((p) => !basename(p).includes('.legacy.'));

  let nSection = 0;
  let nVerbatim = 0;
  const nonPrimaryMode = [];
  const nLegacy = files ? 0 : walk(STRUCT_DIR).filter((p) => basename(p).includes('.legacy.')).length;

  for (const file of targets) {
    const rel = relative(ROOT, file);
    const html = readFileSync(file, 'utf8');
    const seenIds = new Set();

    for (const { attrs, verbatimHtml } of extractSections(html)) {
      if (!/\bstatement\b/.test(attrs['class'] ?? '')) continue;
      nSection++;
      const id = attrs['id'] ?? '(id無し)';
      const at = `${rel} #${id}`;

      const required = ['id', 'data-paper', 'data-pdf-page', 'data-item', 'data-notation-checked'];
      const missing = required.filter((k) => !(k in attrs));
      if (missing.length) { ng(at, `S1 必須属性が欠落: ${missing.join(', ')}`); continue; }

      if (seenIds.has(id)) ng(at, 'S6 id が重複している');
      seenIds.add(id);

      const tag = attrs['data-paper'];
      const paper = reg[tag];
      if (!paper) { ng(at, `S2 data-paper="${tag}" が papers.json に無い`); continue; }

      const page = Number(attrs['data-pdf-page']);
      if (!Number.isInteger(page) || page < 1 || page > paper.pdfPages) {
        ng(at, `S3 data-pdf-page=${attrs['data-pdf-page']} が範囲外(1..${paper.pdfPages})`);
        continue;
      }

      const nc = attrs['data-notation-checked'];
      if (!(DATE_RE.test(nc) || nc === 'none')) {
        ng(at, `S5 data-notation-checked="${nc}" は YYYY-MM-DD か "none" でなければならない`);
      }

      if (verbatimHtml === null) { ng(at, 'S4 .verbatim が無い'); continue; }
      nVerbatim++;

      const pdfPath = join(SOURCE_DIR, `${paper.file}.pdf`);
      if (!existsSync(pdfPath)) { ng(at, `S4 PDF が見つからない: ${pdfPath}`); continue; }

      const texts = pdfPageTexts(pdfPath, page);
      if (texts === null) { ng(at, `S4 pdftotext が失敗した(物理 p.${page})`); continue; }

      const proj = matchProjection(verbatimHtml);
      const hit = texts.find(([, t]) => t.includes(proj));
      if (hit) {
        if (hit[0] !== PDF_MODES[0][0]) nonPrimaryMode.push(`${at} (${hit[0]})`);
      } else {
        // 最も長く一致したモードで、どこまで一致したかを二分探索で示す
        let best = { mode: '', lo: -1 };
        for (const [mode, t] of texts) {
          let lo = 0;
          let hi = proj.length;
          while (lo < hi) {
            const mid = Math.ceil((lo + hi) / 2);
            if (t.includes(proj.slice(0, mid))) lo = mid; else hi = mid - 1;
          }
          if (lo > best.lo) best = { mode, lo };
        }
        ng(at,
          `S4 逐語が物理 p.${page} に見つからない(最良 ${best.mode} モードで先頭 ${best.lo}/${proj.length} 文字まで一致)\n` +
          `      一致した末尾: ...${JSON.stringify(proj.slice(Math.max(0, best.lo - 60), best.lo))}\n` +
          `      次に来るはず: ${JSON.stringify(proj.slice(best.lo, best.lo + 60))}`);
      }
    }
  }

  if (!quiet) {
    if (nonPrimaryMode.length) {
      console.log(`  -- layout 以外のモードで一致した逐語: ${nonPrimaryMode.length} 件`);
      for (const x of nonPrimaryMode) console.log(`     ${x}`);
    }
    console.log(`  -- 構造単位 ${nSection} 件 / 逐語照合 ${nVerbatim} 件 / 対象文書 ${targets.length} 件` +
      (files ? '' : ` / legacy 除外 ${nLegacy} 件`));
    if (NG === before) ok('1_Structured: S1-S6 すべて PASS');
  }
  return NG - before;
}

// ────────────────────────────────────────────────────────────────
// [B-2] Lean ソースの台帳検査(G1-Lean / G2 / G4 / Gap)
//
// ★分業: **型が付くかどうかは Lean が検査する**。ここが見るのは
//   「必要な宣言が存在するか」という帳簿だけ。
//   `Foo.nonvacuous` が本当に `Nonempty Foo` を証明しているかは `lake build`
//   が保証する——check.mjs はその存在しか見ていない。この限界を明示する。
// ────────────────────────────────────────────────────────────────

const DECL_RE = /^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*(theorem|lemma|def|abbrev|instance|structure|inductive|example)\s+([A-Za-z_][\w.'!?₀-₉]*)/;

/** Lean ソース木を走査して、宣言名の集合と structure 一覧を作る。 */
function scanLeanTree(dir) {
  const files = walk(dir).filter((p) => p.endsWith('.lean'));
  const names = new Set();
  const structures = [];   // { name, file, line, bucket }
  const decls = [];        // { kind, name, file, line, bucket }
  const axioms = [];       // { file, line }
  const texts = new Map();

  for (const f of files) {
    const src = readFileSync(f, 'utf8');
    texts.set(f, src);
    const relf = relative(dir, f).replace(/\\/g, '/');
    const bucket = relf.split('/')[0];               // Interface / Skeleton / Found / Gap / Meta
    // コメントを潰してから宣言を拾う(docstring 内の例を宣言と誤認しない)
    const stripped = src
      .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
      .replace(/--[^\n]*/g, (m) => ' '.repeat(m.length));
    stripped.split('\n').forEach((line, i) => {
      const m = DECL_RE.exec(line);
      if (m) {
        names.add(m[2]);
        decls.push({ kind: m[1], name: m[2], file: f, line: i + 1, bucket });
        if (m[1] === 'structure' || m[1] === 'inductive') {
          structures.push({ name: m[2], file: f, line: i + 1, bucket });
        }
      }
      if (/^\s*axiom\s+/.test(line)) axioms.push({ file: f, line: i + 1 });
    });
  }
  return { files, names, structures, decls, axioms, texts };
}

/** 1_Structured から「論文タグ → section id の集合」を集める。 */
function collectStructuredIds() {
  const map = new Map();
  const files = walk(STRUCT_DIR)
    .filter((p) => p.endsWith('.html'))
    .filter((p) => !basename(p).includes('.legacy.'));
  for (const f of files) {
    for (const { attrs } of extractSections(readFileSync(f, 'utf8'))) {
      if (!/\bstatement\b/.test(attrs['class'] ?? '')) continue;
      const paper = attrs['data-paper'];
      const id = attrs['id'];
      if (!paper || !id) continue;
      if (!map.has(paper)) map.set(paper, new Set());
      map.get(paper).add(id);
    }
  }
  return map;
}

/**
 * @param {object} o
 * @param {string} o.dir            走査するルート(本番は lean/ABC3、selftest は fixture)
 * @param {string[]} o.axiomExempt  bucket 相対パスで axiom を許すファイル
 */
function checkLeanLedger({ dir, axiomExempt = [], papersPath = PAPERS_JSON, quiet = false } = {}) {
  const before = NG;
  if (!existsSync(dir)) { ng(relative(ROOT, dir), 'ディレクトリが無い'); return NG - before; }

  const reg = JSON.parse(readFileSync(papersPath, 'utf8')).papers;
  const structuredIds = collectStructuredIds();
  const { names, structures, decls, axioms, texts } = scanLeanTree(dir);
  const at = (x) => `${relative(ROOT, x.file)}:${x.line}`;

  // ── G4: axiom 禁止
  for (const a of axioms) {
    const relf = relative(dir, a.file).replace(/\\/g, '/');
    if (!axiomExempt.includes(relf)) {
      ng(`${relative(ROOT, a.file)}:${a.line}`,
        'G4 axiom 宣言は禁止(未構築の基礎は Interface/ の structure で受ける)');
    }
  }

  // ── G2: Interface の structure は nonvacuous witness か waiting を持つ
  const waiting = [];
  for (const s of structures.filter((x) => x.bucket === 'Interface')) {
    const hasWitness = names.has(`${s.name}.nonvacuous`);
    const hasWaiting = names.has(`${s.name}.waiting`);
    if (!hasWitness && !hasWaiting) {
      ng(at(s), `G2 非空虚 witness が無い: \`${s.name}.nonvacuous : Nonempty ${s.name}\` を構成するか、` +
                `構成できないなら \`${s.name}.waiting : ABC3.Meta.WaitingFor\` で何を待っているか書く`);
    } else if (hasWaiting && !hasWitness) {
      const src = texts.get(s.file) ?? '';
      const m = new RegExp(`${s.name.replace(/[.*+?^${}()|[\\]\\\\]/g, '\\\\$&')}\\.waiting[\\s\\S]{0,400}?what\\s*:=\\s*"([^"]*)"[\\s\\S]{0,200}?trackB\\s*:=\\s*"([^"]*)"`).exec(src);
      waiting.push(m ? `${s.name} — ${m[1]} → ${m[2]}` : `${s.name} — (waiting の内容を読めなかった)`);
    }
  }

  // ── Gap: 各 structure は GapRecord を伴う(falsifier は型が必須にしている)
  const gaps = [];
  for (const s of structures.filter((x) => x.bucket === 'Gap')) {
    if (!names.has(`${s.name}.record`)) {
      ng(at(s), `Gap の記録が無い: \`${s.name}.record : ABC3.Meta.GapRecord\` を書く` +
                '(falsifier は GapRecord の必須フィールド——書けないうちは ③ ではない)');
      continue;
    }
    const src = texts.get(s.file) ?? '';
    const idx = src.indexOf(`${s.name}.record`);
    const m = /classification\s*:=\s*\.?(\w+)[\s\S]{0,400}?falsifier\s*:=\s*"([^"]*)"/
      .exec(idx >= 0 ? src.slice(idx) : src);
    gaps.push(m ? `${s.name} — ${m[1]}: ${m[2]}` : `${s.name} — (record の内容を読めなかった)`);
  }

  // ── G3: load-bearing と印を付けた宣言には負の対照が要る
  for (const d of decls.filter((x) => x.bucket === 'Skeleton')) {
    if (!names.has(`${d.name}.loadBearing`)) continue;
    if (!names.has(`${d.name}.negControl`)) {
      ng(at(d), `G3 負の対照が無い: \`${d.name}.negControl : ABC3.Meta.NegControl\` を書く` +
                '(性質を1つだけ落とした対照が破れることを確認する。破れないならその性質は効いていない)');
    }
  }

  // ── Found/ と Interface/ は sorry を残さない(各 bucket の docstring の規則)
  //    Skeleton/ の sorry は設計どおりなのでゲートしない。
  //    2026-08-14 の監査で、この規則が機械化されていないことが発覚したので追加。
  const SORRY_FREE_BUCKETS = ['Found', 'Interface'];
  for (const [f, src] of texts) {
    const bucket = relative(dir, f).replace(/\\/g, '/').split('/')[0];
    if (!SORRY_FREE_BUCKETS.includes(bucket)) continue;
    stripLeanComments(src).split('\n').forEach((line, i) => {
      if (SORRY_RE.test(line)) {
        ng(`${relative(ROOT, f)}:${i + 1}`,
          `${bucket}/ に sorry を残してはならない(規則は ABC3/${bucket}.lean の docstring)`);
      }
    });
  }

  // ── 主語の分離: Skeleton(原典の主張)は Check(我々のモデルの検査)に依存できない
  for (const [f, src] of texts) {
    const relf = relative(dir, f).replace(/\\/g, '/');
    if (relf.split('/')[0] !== 'Skeleton') continue;
    src.split('\n').forEach((line, i) => {
      if (/^\s*import\s+ABC3\.Check\b/.test(line)) {
        ng(`${relative(ROOT, f)}:${i + 1}`,
          '主語の分離: Skeleton(原典の主張)から Check(我々のモデルの検査)を import してはならない');
      }
    });
  }

  // ── 条件付き形式化の向き: Interface(まだ無い基礎の型)は Found(実装)に依存できない
  //
  //    Skeleton は Interface を import する。したがって Interface が Found を import すると
  //    **Skeleton が Found を推移的に import する**ことになり、
  //    「実装が無くても statement を書ける」という2トラック構成の要点(PLAN §3)が壊れる。
  //    実装が無い間 Skeleton がビルドできなくなり、Interface で受ける意味が消える。
  //
  //    2026-08-14 追加。G2 witness を structure と同じファイルに置こうとして実際に
  //    この向きを作ってしまったのが発端(witness は Found 側から Interface 名前空間へ足す)。
  for (const [f, src] of texts) {
    const relf = relative(dir, f).replace(/\\/g, '/');
    if (relf.split('/')[0] !== 'Interface') continue;
    src.split('\n').forEach((line, i) => {
      if (/^\s*import\s+ABC3\.Found\b/.test(line)) {
        ng(`${relative(ROOT, f)}:${i + 1}`,
          '条件付き形式化の向き: Interface(まだ無い基礎の型)から Found(実装)を import してはならない' +
          '(Skeleton が Found を推移的に引くことになる)。' +
          '非空虚 witness が実装を要するなら、Found 側のファイルで Interface 名前空間へ足す');
      }
    });
  }

  // ── Lean コメント内の引用も PDF に対して照合する
  //    形式:  原文 (<タグ> p.<物理ページ>):
  //           > ...
  //           > ...
  //    引用を2箇所に持つ以上、ズレたら落ちるようにしておく。
  const QUOTE_RE = /原文\s*\(([A-Za-z][\w-]*)\s+p\.(\d+)\)\s*[:：]\s*\n((?:[^\n]*>[^\n]*\n)+)/g;
  let nQuote = 0;
  for (const [f, src] of texts) {
    QUOTE_RE.lastIndex = 0;
    let q;
    while ((q = QUOTE_RE.exec(src)) !== null) {
      nQuote++;
      const line = src.slice(0, q.index).split('\n').length;
      const where = `${relative(ROOT, f)}:${line}`;
      const [, tag, pageStr, body] = q;
      const page = Number(pageStr);
      const paper = reg[tag];
      if (!paper) { ng(where, `引用照合: タグ "${tag}" が papers.json に無い`); continue; }
      if (page < 1 || page > paper.pdfPages) {
        ng(where, `引用照合: p.${page} が範囲外(1..${paper.pdfPages})`); continue;
      }
      const quoted = body.split('\n')
        .map((l) => l.replace(/^[^>]*>\s?/, ''))
        .join(' ');
      const proj = leanQuoteProjection(quoted);
      const pdfPath = join(SOURCE_DIR, `${paper.file}.pdf`);
      const pages = existsSync(pdfPath) ? pdfPageTexts(pdfPath, page) : null;
      if (!pages) { ng(where, `引用照合: PDF を読めない(${paper.file})`); continue; }
      if (!pages.some(([, t]) => t.includes(proj))) {
        let best = { mode: '', lo: -1 };
        for (const [mode, t] of pages) {
          let lo = 0; let hi = proj.length;
          while (lo < hi) {
            const mid = Math.ceil((lo + hi) / 2);
            if (t.includes(proj.slice(0, mid))) lo = mid; else hi = mid - 1;
          }
          if (lo > best.lo) best = { mode, lo };
        }
        ng(where,
          `引用照合: 逐語が ${tag} 物理 p.${page} に見つからない(${best.mode} で ${best.lo}/${proj.length} 文字まで一致)\n` +
          `      次に来るはず: ${JSON.stringify(proj.slice(best.lo, best.lo + 60))}`);
      }
    }
  }

  // ── G6: Skeleton の theorem/lemma は「証明が要求するもの」を持つ
  //    空リストは省略ではなく **主張**(原文の証明は外部依存を持たない)。
  const OBLIGATION_KINDS = ['citation', 'folklore', 'implicitStep', 'otherPaper', 'derivation'];
  const tally = Object.fromEntries(OBLIGATION_KINDS.map((k) => [k, 0]));
  const statusTally = { inMathlib: 0, inProject: 0, inProgress: 0, absent: 0, unmeasured: 0 };
  let nNeeds = 0;
  /** 依存グラフの辺(otherPaper)。`.needs` を辺として読む */
  const edges = [];
  for (const d of decls.filter((x) => x.bucket === 'Skeleton')) {
    // 台帳の付随宣言そのものは対象外
    if (/\.(src|needs|nonvacuous|waiting|record|loadBearing|negControl)$/.test(d.name)) continue;
    // `.needs` を **要求** するのは theorem/lemma だけ。ただし def/structure が
    // 書いている場合は読む——定義も別論文に依拠しうるので、辺はそこからも出る。
    const requiresNeeds = ['theorem', 'lemma'].includes(d.kind);
    if (!names.has(`${d.name}.needs`)) {
      if (requiresNeeds) {
        ng(at(d), `G6 「証明が要求するもの」が無い: \`${d.name}.needs : List ABC3.Meta.ProofObligation\` を書く` +
                  '(依存が無いと考えるなら `[]` と明記する——空欄は不可)');
      }
      continue;
    }
    nNeeds++;
    const src = texts.get(d.file) ?? '';
    // ★`def ` を付けて探す。付けないと docstring 中の言及(``foo.needs`` と書いた散文)を
    //   拾ってしまい、そこから最初の `[` までを本体と誤認する
    //   (2026-08-15 実測: `InitialThetaDataFieldPart` で実際に起きた——
    //    `.needs` は 7 件と数えられたのに中身は 1 件も集計されなかった)。
    let i0 = src.indexOf(`def ${d.name}.needs`);
    if (i0 < 0) i0 = src.indexOf(`${d.name}.needs`);
    const open = src.indexOf('[', i0);
    if (i0 < 0 || open < 0) continue;
    let depth = 0; let end = open;
    for (let i = open; i < src.length; i++) {
      if (src[i] === '[') depth++;
      else if (src[i] === ']') { depth--; if (depth === 0) { end = i; break; } }
    }
    const body = src.slice(open, end + 1);

    // 各 obligation の末尾の数値は物理ページ。範囲外なら NG
    // (`.src` では検査していたのに `.needs` では見ていなかった——2026-08-14 の監査で発覚)
    const srcIdx = src.indexOf(`${d.name}.src`);
    const srcM = srcIdx >= 0 ? SRC_RE.exec(src.slice(srcIdx)) : null;
    const paperTag = srcM ? srcM[1] : null;
    const maxPage = paperTag && reg[paperTag] ? reg[paperTag].pdfPages : null;

    // ── obligation を1件ずつ切り出す(先頭の `.kind` で区切る)。
    //    ★以前は body 全体から数値を拾って **所有論文** の範囲と比べていた。
    //    `otherPaper` は **別の論文** を指すのだから、それでは辺の先を検査していない
    //    (実際 [IUTchI] [IUTchII] が未登記のまま素通りしていた——2026-08-15 の監査)。
    const obs = [];
    {
      const marks = [];
      const kre = /\.(citation|folklore|implicitStep|otherPaper|derivation)\b/g;
      let km;
      while ((km = kre.exec(body))) marks.push({ kind: km[1], at: km.index });
      for (let i = 0; i < marks.length; i++) {
        obs.push({
          kind: marks[i].kind,
          text: body.slice(marks[i].at, i + 1 < marks.length ? marks[i + 1].at : body.length),
        });
      }
    }
    const stripStr = (t) => t.replace(/"[^"]*"/g, '""');
    const lastNum = (t) => {
      const ns = stripStr(t).match(/\d+/g);
      return ns ? Number(ns[ns.length - 1]) : null;
    };

    for (const ob of obs) {
      const pg = lastNum(ob.text);
      if (ob.kind === 'otherPaper') {
        // 辺は三つ組 (paper, item, page) なので papers.json と照合できる
        const strs = [...ob.text.matchAll(/"([^"]*)"/g)].map((x) => x[1]);
        const tag = (strs[0] ?? '').replace(/^\[|\]$/g, '');
        const item = strs[1] ?? '';
        const tp = reg[tag];
        if (!tp) {
          ng(at(d), `G6 \`${d.name}.needs\` の otherPaper が指す "[${tag}]" が papers.json に無い` +
                    '——辺の先を検査できないので、先に登記すること');
          continue;
        }
        if (pg === null || pg < 1 || pg > tp.pdfPages) {
          ng(at(d), `G6 \`${d.name}.needs\` の otherPaper "[${tag}] ${item}" が` +
                    `物理 p.${pg} を指すが [${tag}] の範囲外(1..${tp.pdfPages})`);
          continue;
        }
        // fromPaper を持たせるのは、同一論文内の辺と論文をまたぐ辺を分けて数えるため。
        // ★`otherPaper` という名前に反して、tag は **同じ論文** でよい(そう使われている)。
        edges.push({ from: d.name, fromPaper: paperTag, tag, item, page: pg });
      } else if (maxPage !== null && pg !== null && (pg < 1 || pg > maxPage)) {
        ng(at(d), `G6 \`${d.name}.needs\` が物理 p.${pg} を指しているが範囲外(1..${maxPage}, ${paperTag})`);
      }
    }

    for (const k of OBLIGATION_KINDS) {
      tally[k] += (body.match(new RegExp(`\\.${k}\\b`, 'g')) ?? []).length;
    }
    for (const k of Object.keys(statusTally)) {
      statusTally[k] += (body.match(new RegExp(`\\.${k}\\b`, 'g')) ?? []).length;
    }
  }

  // ── G1-Lean: Skeleton の宣言は Source を伴い、その locator が実在する
  let nSrcOk = 0;
  for (const d of decls.filter((x) => x.bucket === 'Skeleton')) {
    // 台帳の付随宣言そのものには出典を要求しない
    if (/\.(src|needs|nonvacuous|waiting|record|loadBearing|negControl)$/.test(d.name)) continue;
    if (!['theorem', 'lemma', 'def', 'abbrev', 'structure'].includes(d.kind)) continue;
    if (!names.has(`${d.name}.src`)) {
      ng(at(d), `G1 出典が無い: \`${d.name}.src : ABC3.Meta.Source\` を書く`);
      continue;
    }
    const src = texts.get(d.file) ?? '';
    let idx = src.indexOf(`def ${d.name}.src`);
    if (idx < 0) idx = src.indexOf(`${d.name}.src`);
    const m = SRC_RE.exec(idx >= 0 ? src.slice(idx) : src);
    if (!m) { ng(at(d), `G1 \`${d.name}.src\` の中身を読めなかった(1行の正準形で書くこと)`); continue; }
    const [, paper, pageStr, sectionId] = m;
    const page = Number(pageStr);
    const p = reg[paper];
    if (!p) { ng(at(d), `G1 paper="${paper}" が papers.json に無い`); continue; }
    if (page < 1 || page > p.pdfPages) { ng(at(d), `G1 pdfPage=${page} が範囲外(1..${p.pdfPages})`); continue; }
    const ids = structuredIds.get(paper);
    if (!ids || !ids.has(sectionId)) {
      ng(at(d), `G1 sectionId="${sectionId}" が 1_Structured(${paper})に無い` +
                `——先に構造化するか、id を直す`);
      continue;
    }
    nSrcOk++;
  }

  // ── 依存グラフ(指標。★ゲートではない)
  //    `.needs` の otherPaper を辺として推移閉包を取り、
  //    (a) 着地した葉 (b) 未展開の葉 (c) 循環 (d) 深さ を分けて印字する。
  const ITEM_RE = /^\s*(Theorem|Proposition|Corollary|Definition|Lemma|Remark|Example)\s+([0-9]+(?:\.[0-9]+)*)/;
  // ★番号付き項目でない「(L1)」「(Ob1)」型のラベルも節点キーにする。
  //   これが無いと、辺は説明つきの長い文字列、`.src` は `(L1)` だけ、で
  //   **同じものが別節点に見え**、張ったのに未展開の葉のまま残る(2026-08-15 実測)。
  const LABEL_RE = /^\s*\(([A-Za-z]+[0-9]*)\)/;
  const nodeKey = (tag, item) => {
    const m = ITEM_RE.exec(item);
    if (m) return `[${tag}] ${m[1]} ${m[2]}`;
    const l = LABEL_RE.exec(item);
    if (l) return `[${tag}] (${l[1]})`;
    return `[${tag}] ${item.trim().slice(0, 30)}`;
  };
  const SRC_ITEM_RE =
    /\.src[\s\S]{0,400}?paper\s*:=\s*"([^"]*)"[\s\S]{0,300}?item\s*:=\s*"([^"]*)"/;
  /** 節点キー → それを張っている宣言名(スケルトンがある = 展開済み) */
  const expanded = new Map();
  const declKey = new Map();
  for (const d of decls.filter((x) => x.bucket === 'Skeleton')) {
    if (/\.(src|needs|nonvacuous|waiting|record|loadBearing|negControl)$/.test(d.name)) continue;
    const s2 = texts.get(d.file) ?? '';
    const i2 = s2.indexOf(`${d.name}.src`);
    if (i2 < 0) continue;
    const m2 = SRC_ITEM_RE.exec(s2.slice(i2));
    if (!m2) continue;
    const k = nodeKey(m2[1], m2[2]);
    if (!expanded.has(k)) expanded.set(k, d.name);
    declKey.set(d.name, k);
  }
  const adj = new Map();
  const unexpanded = new Map();
  for (const e of edges) {
    const from = declKey.get(e.from);
    if (!from) continue;
    const to = nodeKey(e.tag, e.item);
    if (!adj.has(from)) adj.set(from, new Set());
    adj.get(from).add(to);
    if (!expanded.has(to) && !unexpanded.has(to)) unexpanded.set(to, e);
  }
  const incoming = new Set();
  for (const st of adj.values()) for (const tt of st) incoming.add(tt);
  const roots = [...expanded.keys()].filter((k) => !incoming.has(k));
  let maxDepth = 0;
  const cycles = [];
  for (const r of roots) {
    const stack = [[r, 0, new Set([r])]];
    while (stack.length) {
      const [node, dep, path] = stack.pop();
      if (dep > maxDepth) maxDepth = dep;
      for (const tt of adj.get(node) ?? []) {
        if (path.has(tt)) { cycles.push(`${node} → ${tt}`); continue; }
        stack.push([tt, dep + 1, new Set([...path, tt])]);
      }
    }
  }

  if (!quiet) {
    console.log(`  -- Lean 宣言 ${decls.length} / structure ${structures.length} / axiom ${axioms.length} 件` +
                `(免除 ${axiomExempt.length} ファイル)`);
    const nCheck = decls.filter((x) => x.bucket === 'Check').length;
    console.log(`  -- Skeleton の出典照合: ${nSrcOk} 件 OK / Check(我々のモデルの検査): ${nCheck} 宣言`);
    console.log(`  -- Interface 実装待ち: ${waiting.length} 件(= Track B の作業キュー)`);
    for (const w of waiting) console.log(`     ${w}`);
    console.log(`  -- Gap(飛躍): ${gaps.length} 件`);
    for (const g of gaps) console.log(`     ${g}`);
    console.log(`  -- 引用照合(Lean コメント内): ${nQuote} 件`);
    console.log('  -- 依存グラフ(`.needs` の otherPaper を辺として推移閉包。★指標であってゲートではない):');
    console.log(`     スケルトンのある節点  ${expanded.size} 件(うち根 ${roots.length})`);
    const sameEdges = edges.filter((e) => e.fromPaper === e.tag).length;
    console.log(`     辺(otherPaper)        ${edges.length} 本` +
                `(同一論文内 ${sameEdges} / 論文をまたぐ ${edges.length - sameEdges})` +
                ` / 最大深さ ${maxDepth} / 循環 ${cycles.length} 件`);
    for (const c of cycles) console.log(`       循環: ${c}`);
    // ★「循環 N 件」の意味(2026-08-15 に辺の意味を定めてから、この行は主張になった)。
    //   辺は「原文の証明・主張が実際に依拠しているもの」であり、前方参照・解説案内・
    //   導入部からの指しは辺ではない(`Meta/Claim.lean` の `otherPaper`)。
    //   したがって:
    //     循環 0 件 = **証明の依拠関係に循環が無い**(= 我々が写した範囲では、
    //                 原典の議論に順序が付く)。
    //     循環あり  = 次のどれか。**既定は①**:
    //       ① 我々の転写の誤り(前方参照や解説案内を辺として書いてしまった)。
    //          原文側のグラフで実測: 目印つきの辺を落とすと最大 SCC が 262 → 8 になり、
    //          残る循環はすべて1論文の中に閉じた(`ResearchPaper/cycle-analysis.md`)。
    //          すなわち**循環の大半は辺の定義の副作用**である。
    //       ② 同時再帰的な定義・主張(原典が意図して同時に立てているもの)。
    //          この場合は循環ではなく「束ねて1つの節点にする」のが正しい写し方。
    //       ③ 原典側の問題。★これを主張するには ① と ② を先に潰すこと。
    if (cycles.length === 0) {
      console.log('       ★循環 0 件 = 我々が写した範囲では、証明の依拠関係に順序が付く');
    } else {
      console.log('       ★循環あり——既定の疑いは**我々の転写の誤り**(前方参照や解説案内を');
      console.log('          辺にしていないか)。原典側の問題と判断する前に、それを潰すこと');
    }
    console.log(`     着地した葉            ${statusTally.inMathlib + statusTally.inProject} 件` +
                `(mathlib ${statusTally.inMathlib} / 公開プロジェクト ${statusTally.inProject})`);
    // ★葉を仕分ける。「次に張れる」と「中層が無いので今は張れない」は別物である。
    //   ここは**自己申告**の読み込みであって機械判定ではない(A 群と同じ性質)。
    let blockedMap = new Map(); let workMap = new Map();
    if (existsSync(BLOCKED_JSON)) {
      try {
        const bl = JSON.parse(readFileSync(BLOCKED_JSON, 'utf8'));
        for (const b of bl.blocked ?? []) blockedMap.set(b.key, b);
        for (const b of bl.expandableWithWork ?? []) workMap.set(b.key, b);
      } catch { /* 読めなければ仕分けなしで印字する */ }
    }
    const unsorted = [...unexpanded.keys()].filter((k) => !blockedMap.has(k) && !workMap.has(k));
    console.log(`     ★未展開の葉          ${unexpanded.size} 件` +
                `(中層待ち ${[...unexpanded.keys()].filter((k) => blockedMap.has(k)).length}` +
                ` / 手間だけ ${[...unexpanded.keys()].filter((k) => workMap.has(k)).length}` +
                ` / 未仕分け ${unsorted.length})`);
    for (const [k, e] of unexpanded) {
      const tagStr = blockedMap.has(k) ? '中層待ち' : workMap.has(k) ? '手間だけ' : '★未仕分け';
      console.log(`       [${tagStr}] ${k} — 物理 p.${e.page}(${e.from} から)`);
      const b = blockedMap.get(k);
      if (b) console.log(`                 解消条件: ${b.unblockedBy}`);
    }
    if (unsorted.length) {
      console.log(`     ★未仕分けの葉が ${unsorted.length} 件ある` +
                  `——${relative(ROOT, BLOCKED_JSON)} に「なぜ張れないか」を書くこと`);
    }
    console.log('     ★★この数は**張れば増える**——辺の先を張ると、その先の辺が新たに現れる。');
    console.log('        減ったことを進捗と読まないこと(辿るのをやめても減る)。');
    console.log('        進捗として読むなら「スケルトンのある節点」と「最大深さ」の側である。');
    // ── ★被覆率: 我々が **書けている** 量 対 原文の機械概算
    //    器具は自分の欠陥を報告できなければならない(冒頭の原則)。ここは
    //    「グラフがどれだけ足りていないか」を毎回目に入れるための行である。
    if (existsSync(SCALE_JSON)) {
      try {
        const sc = JSON.parse(readFileSync(SCALE_JSON, 'utf8'));
        const tr = sc.trustedTotals ?? {};
        const pct = (a, b) => (b ? `${((100 * a) / b).toFixed(1)}%` : '—');
        console.log(`     ★被覆率(対 原文の機械概算、${sc.measuredAt}、${sc.tool}):`);
        console.log(`        節点 ${expanded.size} / ${tr.reachable} = ${pct(expanded.size, tr.reachable)}` +
                    `   辺 ${edges.length} / ${tr.edges} = ${pct(edges.length, tr.edges)}`);
        console.log('        ★★分母は**下界でも上界でもない**——番号の無い依存を数え落とし、');
        console.log('           「cf.」の案内を依存として数え過ぎる。桁を見るためだけに使うこと。');
        console.log(`           限界の全文は ${relative(ROOT, SCALE_JSON)} の "★limits"。`);
      } catch {
        console.log('     ★被覆率: dependency-scale.json を読めなかった');
      }
    }
    const total = Object.values(tally).reduce((a, b) => a + b, 0);
    console.log(`  -- 規模(原文の証明文からの抽出、${nNeeds} 定理ぶん、★下界):`);
    console.log(`     引用 ${tally.citation} 件 — mathlib ${statusTally.inMathlib} / ` +
                `公開プロジェクト(完成) ${statusTally.inProject} / 同(作業中) ${statusTally.inProgress} / ` +
                `不在 ${statusTally.absent} / 未測定 ${statusTally.unmeasured}`);
    console.log(`     典拠なし(well-known 等) ${tally.folklore} 件 ← 大きさ未知`);
    console.log(`     暗黙の段(Gap 候補)      ${tally.implicitStep} 件`);
    console.log(`     別論文への枝            ${tally.otherPaper} 件`);
    console.log(`     原文内の導出            ${tally.derivation} 件`);
    console.log(`     合計 ${total} 件`);
    if (statusTally.unmeasured > 0) {
      console.log(`     ★ 未測定が ${statusTally.unmeasured} 件ある——集計は未確定`);
    }
    if (statusTally.inProgress > 0) {
      console.log(`     ★ 作業中の公開プロジェクトに依る項目が ${statusTally.inProgress} 件` +
                  `——独立に作ると重複投資(ResearchPaper/lean-ecosystem.json 参照)`);
    }
    console.log('  -- 注意: ここは宣言の**存在**しか見ていない。型の正しさは lake build が保証する');
    console.log('  -- 注意: 規模は原文が挙げた依存のみ。証明を書いて初めて要ると分かるものは写らない');
  }
  return NG - before;
}

// ────────────────────────────────────────────────────────────────
// [B] Lean の検査
// ────────────────────────────────────────────────────────────────

function checkLean() {
  const before = NG;
  if (!existsSync(join(LEAN_DIR, 'lakefile.toml'))) { ng('lean/', 'lakefile.toml が無い'); return 1; }

  try {
    execFileSync('lake', ['build'], { cwd: LEAN_DIR, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024 });
    ok('lake build 成功');
  } catch (e) {
    const out = `${e.stdout ?? ''}${e.stderr ?? ''}`;
    ng('lean/', `lake build 失敗\n${out.split('\n').slice(-15).join('\n')}`);
  }

  const leanFiles = walk(LEAN_SRC).filter((p) => p.endsWith('.lean'));

  // sorry 台帳(件数はゲートではない。必ず印字する)
  // bucket 別のゲート(Found/ Interface/ は sorry 禁止)は checkLeanLedger にある
  // ——そちらは fixture で較正できる場所だから。
  const sorries = [];
  for (const f of leanFiles) {
    stripLeanComments(readFileSync(f, 'utf8')).split('\n').forEach((line, i) => {
      if (SORRY_RE.test(line)) sorries.push(`${relative(ROOT, f)}:${i + 1}`);
    });
  }
  console.log(`  -- sorry: ${sorries.length} 件${sorries.length ? `\n     ${sorries.join('\n     ')}` : ''}`);
  // ★この数字を単独で読まないこと(2026-08-14。B6 参照)。
  //   `Interface` にフィールドを足せば `Skeleton` の sorry は消える——仕事は増えているのに。
  console.log('     ★この数字を単独で進捗と読まないこと: `Interface` に条件を posit すれば' +
              ' sorry は消える(実例: IUTchIII cor_3_12、2026-08-14)。下の「暗黙の段」と対で見ること');


  // G1-Lean / G2 / G4 / Gap の台帳検査
  checkLeanLedger({ dir: LEAN_SRC, axiomExempt: AXIOM_EXEMPT });

  if (NG === before) ok('Lean: ゲート PASS');
  return NG - before;
}

// ────────────────────────────────────────────────────────────────
// selftest — 器具自身の較正
// ────────────────────────────────────────────────────────────────

function selftest() {
  h1('selftest: 器具は壊れた入力を落とせるか');
  const tmp = mkdtempSync(join(tmpdir(), 'abc3-selftest-'));
  const papers = join(tmp, 'papers.json');
  writeFileSync(papers, readFileSync(PAPERS_JSON, 'utf8'));

  const mk = (name, body) => {
    const p = join(tmp, name);
    writeFileSync(p, `<!DOCTYPE html><html><body>${body}</body></html>`, 'utf8');
    return p;
  };
  const sec = (attrs, verbatim) =>
    `<section class="statement proposition" ${attrs}><div class="verbatim">${verbatim}</div></section>`;

  const cases = [
    ['D1 存在しないPDFページを指す locator',
      mk('d1.html', sec('id="d1" data-paper="pGC" data-pdf-page="999" data-item="X" data-notation-checked="2026-08-14"', 'whatever'))],
    ['D2 逐語が該当ページに無い',
      mk('d2.html', sec('id="d2" data-paper="pGC" data-pdf-page="3" data-item="X" data-notation-checked="2026-08-14"', 'This sentence does not appear anywhere in the paper.'))],
    ['D3 必須属性の欠落(data-notation-checked 無し)',
      mk('d3.html', sec('id="d3" data-paper="pGC" data-pdf-page="3" data-item="X"', 'x'))],
    ['D4 未登記の論文タグ',
      mk('d4.html', sec('id="d4" data-paper="NOPE" data-pdf-page="1" data-item="X" data-notation-checked="none"', 'x'))],
    ['D5 目視確認日が日付でない',
      mk('d5.html', sec('id="d5" data-paper="pGC" data-pdf-page="3" data-item="Proposition 1.1" data-notation-checked="yes"', 'The cyclotomic character'))],
  ];

  let passed = 0;
  for (const [label, file] of cases) {
    const before = NG;
    checkStructured({ files: [file], papersPath: papers, quiet: true });
    const caught = NG > before;
    NG = before; // selftest の結果は本体の NG に混ぜない
    console.log(`  ${caught ? 'ok ' : 'NG '} ${label} → ${caught ? '落とせた' : '★素通りした'}`);
    if (caught) passed++;
  }

  // D6: 正しい入力が通ること(偽陽性の検査)
  const real = join(STRUCT_DIR,
    'A Version of the Grothendieck Conjecture for p-adic Local Fields', 'section-1.html');
  let d6 = false;
  if (existsSync(real)) {
    const before = NG;
    checkStructured({ files: [real], quiet: true });
    d6 = NG === before;
    NG = before;
  }
  console.log(`  ${d6 ? 'ok ' : 'NG '} D6 正しい入力(pGC §1)は通る → ${d6 ? '通った' : '★落ちた(偽陽性)'}`);
  if (d6) passed++;

  // ── Lean 側の較正。fixture を Interface/ Gap/ に配置した木を組み立てて叩く。
  const leanCases = [
    ['D7 Interface に非空虚 witness も waiting も無い', 'Interface', 'd7-vacuous-interface.lean', true],
    ['D8 プロジェクト内の axiom 宣言', 'Interface', 'd8-axiom.lean', true],
    ['D9 Gap に GapRecord(falsifier)が無い', 'Gap', 'd9-gap-no-record.lean', true],
    ['D10 非空虚 witness を持つ Interface は通る', 'Interface', 'd10-good-witness.lean', false],
    ['D11 waiting を書いた Interface は通る', 'Interface', 'd11-waiting.lean', false],
    ['D12 load-bearing なのに負の対照が無い', 'Skeleton', 'd12-loadbearing-no-negcontrol.lean', true],
    ['D13 load-bearing + 負の対照ありは通る', 'Skeleton', 'd13-loadbearing-ok.lean', false],
    ['D14 Skeleton が Check を import している', 'Skeleton', 'd14-skeleton-imports-check.lean', true],
    ['D15 Skeleton の theorem に .needs が無い', 'Skeleton', 'd15-no-needs.lean', true],
    ['D16 .needs があれば通る(空リストも主張)', 'Skeleton', 'd16-needs-ok.lean', false],
    ['D17 コメント内の引用が該当ページに無い', 'Skeleton', 'd17-bad-quote.lean', true],
    ['D18 装飾記法つきの正しい引用は通る', 'Skeleton', 'd18-good-quote.lean', false],
    ['D19 .needs が範囲外のページを指す', 'Skeleton', 'd19-needs-bad-page.lean', true],
    ['D20 .needs のページが範囲内なら通る', 'Skeleton', 'd20-needs-good-page.lean', false],
    ['D21 Found/ に sorry が残っている', 'Found', 'd21-found-sorry.lean', true],
    ['D22 Found/ が sorry 無し(docstring の言及は誤検出しない)', 'Found', 'd22-found-clean.lean', false],
    ['D23 Interface が Found を import している', 'Interface', 'd23-interface-imports-found.lean', true],
    ['D24 Interface が Found を import していない', 'Interface', 'd24-interface-no-found-import.lean', false],
    ['D25 フラクトゥール記法つきの正しい引用は通る', 'Skeleton', 'd25-frak-quote.lean', false],
    ['D26 otherPaper の辺が未登記の論文を指す', 'Skeleton', 'd26-edge-unregistered.lean', true],
    ['D27 otherPaper の辺が登記済み・範囲内なら通る', 'Skeleton', 'd27-edge-ok.lean', false],
    ['D28 辺のページが引用元では範囲内・辺の先では範囲外', 'Skeleton',
      'd28-edge-page-of-wrong-paper.lean', true],
    ['D29 docstring の言及を .needs 本体と取り違えない', 'Skeleton',
      'd29-needs-mentioned-in-docstring.lean', true],
    ['D30 同一論文内の辺も範囲検査される', 'Skeleton', 'd30-edge-same-paper.lean', true],
    ['D31 文字列内の sorry を台帳に数えない', 'Skeleton', 'd31-sorry-in-string.lean', false],
    ['D32 文字列を潰しても本物の sorry は落とせる', 'Found', 'd32-real-sorry-in-found.lean', true],
  ];
  const FIXTURES = join(ROOT, 'tools', 'selftest-fixtures');
  for (const [label, bucket, fixture, shouldFail] of leanCases) {
    const root = join(tmp, `lean-${fixture}`);
    mkdirSync(join(root, bucket), { recursive: true });
    writeFileSync(join(root, bucket, fixture), readFileSync(join(FIXTURES, fixture), 'utf8'), 'utf8');
    const before = NG;
    checkLeanLedger({ dir: root, quiet: true });
    const failed = NG > before;
    NG = before;
    const good = failed === shouldFail;
    console.log(`  ${good ? 'ok ' : 'NG '} ${label} → ` +
      (shouldFail
        ? (failed ? '落とせた' : '★素通りした')
        : (failed ? '★落ちた(偽陽性)' : '通った')));
    if (good) passed++;
  }

  rmSync(tmp, { recursive: true, force: true });
  const total = cases.length + 1 + leanCases.length;
  console.log(`\n  selftest: ${passed}/${total} PASS`);
  if (passed !== total) NG++;
  return passed === total;
}

// ────────────────────────────────────────────────────────────────

const args = process.argv.slice(2);
const only = (f) => args.includes(f);
const all = args.length === 0;

console.log(`Math_ABC3 check — ${new Date().toISOString().slice(0, 19).replace('T', ' ')}`);

if (all || only('--selftest')) selftest();
if (all || only('--structured')) { h1('1_Structured'); checkStructured(); }
if (all || only('--lean')) { h1('Lean'); checkLean(); }

console.log(`\n${NG === 0 ? 'PASS' : `NG ${NG} 件`}`);
process.exit(NG === 0 ? 0 : 1);
