import ABC3.Found.FrdI.Prop21

/-!
# Gap — [FrdI] §2 で原典の量化子が合わなかった 1 段

★**`Found/FrdI/Prop21.lean` で `Proposition 2.1` の (i)(ii)(iii) をすべて実装した。**
★**ただし (iii) の `⟸` は、原文どおりの量化子では出ない。**

## ★原文の読み

原文 (FrdI p.44):
> (The Naive Frobenius Functor) Let d ∈N≥1. Then:

原文 (FrdI p.44):
> (iii) C is of perfect type if and only if Ψ is an equivalence of categories.

★**`d` は命題の冒頭で固定されている。**したがって (iii) を字義どおり読むと
「**ある固定した `d`** について `Ψ_d` が圏同値 ⟺ `𝒞` が perfect 型」である。

## ★★不整合

`perfect object`(`Definition 1.2, (iv)`、`MorphismTypes.lean:422`)は
★**すべての `n ∈ ℕ≥1`** を量化する:

- (a) 各次数 `n` について、`A` に base-isomorphic な対象は
  **次数 `n` の Frobenius 型射の終域として現れる**
- (b) その上の pre-step が一意に降りる

★★**`Ψ_d` が圏同値であることから出るのは `n = d` の場合だけ**である ——
本質的全射性は「次数 `d` の Frobenius 型射の終域」しか与えず、
充満性・忠実性も次数 `d` の四角形しか与えない。
★**次数は乗法的なので、1 つの `d` から他の `n` は作れない**
(`d` の冪しか出ず、しかも `perfect` は素数次数も要求する)。

## ★我々の実装

- `prop_2_1_iii_mp` —— ★**`⟹` は `d` を固定したまま成り立つ**(原文どおり)
- `prop_2_1_iii_mpr` —— ★**`⟸` は「すべての `n`」を仮定して証明した**
- `prop_2_1_iii` —— 両者を合わせた `IsOfPerfectType P ↔ ∀ n, (naiveFrob P F n).IsEquivalence`

★★**2026-08-17、利用者の判断により決着を改めた。**

以前ここには「`.src` は付けない —— 原文の (iii) の `⟸` を字義どおりの量化子で
実装したわけではないから」と書いていた。★**方針を次のように改めた**:

- ★**指標が測るのは「原典項目の**数学的内容**を完全に実装したか」**である。
  逐語の文が偽である場合、その文を実装することは**原理的に不可能**であり、
  「付けない」を貫くと**その項目は永久に数えられない**。
- ★したがって `Found/FrdI/Prop21.lean` に `prop_2_1.src`(条なし)を付けた。
- ★★**代わりに逸脱を 3 箇所で開示する**: (1) `prop_2_1.src` の docstring、
  (2) この `GapRecord`(分類 ③・反例つき、存置)、
  (3) ★**`index.html` の「原文の記載不備と、我々が足した前提」**
  (`tools/index-html.mjs` が生成)——**原文の訂正が必要**であることが
  後から誰にでも分かるようにするため。
-/

namespace ABC3.Gap.FrdI

open CategoryTheory ABC3.Found.FrdI

universe v u w u2 v2

/-- ★**`Proposition 2.1, (iii)` の `⟸` に不足しているもの**。

★原典の語彙で書けば「**1 つの次数 `d` で `Ψ_d` が圏同値なら、
すべての次数で Frobenius 型射の終域が取れる**」。

★★**「`d` 固定」の読みは機械検証で偽になった**(`Check/FrdI/Prop21QuantifierGap.lean`)。
★だから分類は ③(`sourceGap`)である。
★**ただし数学の誤りではなく書き方の不備と見るのが妥当**で、
意図された主張(すべての `d`)は `prop_2_1_iii` で証明済みである。 -/
structure Gap_2_1_iii {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
    {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P) (d : ℕ+) :
    Prop where
  /-- 不足: **1 つの次数から全次数へ渡ること**。 -/
  allDegrees : (naiveFrob P F d).IsEquivalence → ∀ n : ℕ+, (naiveFrob P F n).IsEquivalence

def Gap_2_1_iii.record : ABC3.Meta.GapRecord :=
  { source :=
      { paper := "FrdI", pdfPage := 44, item := "Proposition 2.1, (iii)",
        sectionId := "frdi-prop-2-1" },
    classification := ABC3.Meta.GapClass.sourceGap,
    falsifier :=
      "★★**機械検証済み**(`Check/FrdI/Prop21QuantifierGap.lean`)。" ++
      "★原文は `Let d ∈N≥1.` と冒頭で `d` を固定したうえで (iii) を述べるが、" ++
      "**その読みでは偽**である: " ++
      "(1) `prop_2_1_iii_fixed_degree_false_general` —— d = 1 では " ++
      "「素朴 Frobenius 函手が恒等函手と同型」がつねに成り立つので、" ++
      "perfect 型でない Frobenioid **すべて**で反例になる。" ++
      "(2) `prop_2_1_iii_fixed_degree_false_deg_three` —— d ≥ 2 でも同じ: " ++
      "Φ = Z/2(定数)では 3 倍が全単射だが 2 倍は非全射なので、" ++
      "次数 3 の素朴 Frobenius 函手は圏同値だが n = 2 で perfect の条件が破れる。" ++
      "★原因は `perfect object`(Definition 1.2, (iv))が**すべての n ∈ N≥1** を量化すること。" ++
      "★★**ただしこれは書き方の不備であって、数学の誤りではない可能性が高い** —— " ++
      "意図された主張「すべての d で Ψ_d が圏同値 ⟺ perfect 型」は " ++
      "`Found/FrdI/Prop21.lean` の `prop_2_1_iii` で**証明済み**である。" ++
      "★`⟹` の側は原文どおり d を固定したまま成り立つ(`prop_2_1_iii_mp`)。" ++
      "★これが覆るのは、原典の `d` が実は全称だと確定した場合のみである。" ++
      "★2026-08-17、利用者の判断により `.src`(条なし)を付けた。" ++
      "逸脱は index.html の「原文の記載不備と、我々が足した前提」で開示している。" }

/-! ## ★★`Definition 2.8, (ii)` —— **§2 の最後の 1 件が止まっている所**(2026-08-17 に測定)

原文 (FrdI p.52):
> (ii) Suppose that M is a topologically finitely generated profinite abelian group.

原文 (FrdI p.52):
> Thus, M decomposes as a direct product of pro-l groups M [l], where l varies over

★★**`Definition 2.8` は §2 に残る唯一の未実装項目**である(他の 7 件は実装済み)。
★★★**律速は「位相的有限生成な副有限アーベル群の pro-`l` 分解」であり、
これは mathlib に無い。**

## ★測定(2026-08-17、探索範囲つき)

| 要るもの | mathlib | 判定 |
|---|---|---|
| `ProfiniteGrp`(圏・極限表示) | `Topology/Algebra/Category/ProfiniteGrp/`(Basic/Completion/Limits) | ★**ある**。`continuousMulEquivLimittoFiniteQuotientFunctor` で開正規部分群による極限表示が取れる |
| 有限アーベル群の準素分解 | `Algebra/Module/Torsion/PrimaryComponent.lean`(`iSup_primaryComponent_eq_top` / `iSupIndep_primaryComponent`) | ★**ある**(捻れ加群について) |
| **pro-`l` 群** | ★★**実質不在**(`pro-p` で 2 件、いずれも別物。2026-08-16 の `Found/GenEll/Lemma31.lean` の測定と一致) | ★★**無い** |
| 位相的有限生成 | 該当なし | ★**無い**(定義するのは容易) |

★**道はある**: `M ≅ lim_U M/U`(mathlib)＋ 各 `M/U` の準素分解(mathlib)から
`M[l] := lim_U (M/U)[l]` を作り、極限と積の交換で `M ≅ ∏_l M[l]` を出す。

## ★★2026-08-18 —— 「壁」と呼ぶのをやめ、statement を型で固定して道を割った

★CLAUDE.md の**姿勢**(「工数の山を『壁』と呼ばない。既知数学の person-years は
壁でなく道」)に従い、次を行った:

- ★statement を **`Skeleton/FrdI/Def28ProL.lean`** に型で固定した
  (`IsProL` の定義と `def_2_8_ii`)。`sorry` は `Skeleton/` では正しい状態である。
- ★道を **`ResearchPaper/frdi-decomposition.json`** の `prol` チェーン
  (7 節点・4 層・葉 4 個)に割った。`node tools/frdi-newleaves.mjs` が層と葉を印字する。

★★測り直すと、**部品の 4/5 は既に在庫にある**:
`ProfiniteGrp.continuousMulEquivLimittoFiniteQuotientFunctor`(極限表示)・
`Ideal.iSup_primaryComponent_eq_top` / `Ideal.iSupIndep_primaryComponent`(準素分解)・
`Ideal.primaryComponent.map`(その関手性)・
`CategoryTheory.Limits.limitFlipCompLimIsoLimitCompLim`(極限と積の交換)。
★足りないのは **pro-`l` の語彙**(我々が定義した)と、それらを繋ぐ作業である。

## ★なぜ迂回できないか

`Definition 2.8, (ii)` の分解は**下流で本質的に使われる**:
原文 p.106–107 は「`∏_p 𝒪^×(A)[p] ≅ 𝒪^×(A)`(cf. Definition 2.8, (ii))」を
**そのまま等式として**使い、各 `p` 成分ごとに `u_p` を作って貼り合わせる。
★したがって「分解を仮定に置く」形では下流が posit に依ることになり、
`check.mjs` の B5(条件を posit して `sorry` を消しても進捗ではない)に触れる。

★**(i) と (iii) の `co-prime type` は分解を要さない**(定義だけ)。
★しかし (iii) の「λ 乗写像」は**分解が無いと定義できない**——
`M[l]` ごとに `λ(l)` 乗する、という定義だからである。
-/

/-- ★★**`Definition 2.8, (ii)` に不足しているもの**。

★分類は **② `missingMath`** —— 原文の主張は標準的な数学であり真である。
**足りないのは mathlib の在庫**であって、原典の飛躍ではない。 -/
structure Gap_2_8_ii : Prop where
  /-- 不足: **位相的有限生成な副有限アーベル群の pro-`l` 分解**。 -/
  proLDecomposition : True

def Gap_2_8_ii.record : ABC3.Meta.GapRecord :=
  { source :=
      { paper := "FrdI", pdfPage := 52, item := "Definition 2.8, (ii)",
        sectionId := "frdi-def-2-8" },
    classification := ABC3.Meta.GapClass.missingMath,
    falsifier :=
      "★**これが ① に落ちる条件**: mathlib(または我々)が" ++
      "「位相的有限生成な副有限アーベル群 M は pro-l 群の直積に分解する」を" ++
      "実装すること。★道は測定済みで、`ProfiniteGrp` の極限表示" ++
      "(`continuousMulEquivLimittoFiniteQuotientFunctor`)と" ++
      "有限アーベル群の準素分解(`Algebra/Module/Torsion/PrimaryComponent.lean`)から" ++
      "`M[l] := lim_U (M/U)[l]` を作り、極限と積の交換で出す。" ++
      "★**③ に上がることはない** —— 原文の主張は標準的な数学であり、" ++
      "反例は存在しない。★規模は mathlib の PR 数本ぶんと見積もった(2026-08-17)。" ++
      "★これが §2 を 7/8 で止めている唯一の原因である。" }

end ABC3.Gap.FrdI
