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

★★**`.src`(＝原典項目を完全に実装したという主張)は付けない。**
★原文の (iii) の `⟸` を字義どおりの量化子で実装したわけではないからである。
-/

namespace ABC3.Gap.FrdI

open CategoryTheory ABC3.Found.FrdI

universe v u w u2 v2

/-- ★**`Proposition 2.1, (iii)` の `⟸` に不足しているもの**。

★原典の語彙で書けば「**1 つの次数 `d` で `Ψ_d` が圏同値なら、
すべての次数で Frobenius 型射の終域が取れる**」。

★★**これは量化子の問題であって、数学の穴とは限らない** ——
原文が `d` を固定して書いたのは**書き方の都合**で、
著者は「すべての `d` について」を意図していた可能性が高い。
★だから分類は ①(`modelError`)でも ③(`sourceGap`)でもなく ②(`missingMath`)とする。 -/
structure Gap_2_1_iii {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
    {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P) (d : ℕ+) :
    Prop where
  /-- 不足: **1 つの次数から全次数へ渡ること**。 -/
  allDegrees : (naiveFrob P F d).IsEquivalence → ∀ n : ℕ+, (naiveFrob P F n).IsEquivalence

def Gap_2_1_iii.record : ABC3.Meta.GapRecord :=
  { source :=
      { paper := "FrdI", pdfPage := 44, item := "Proposition 2.1, (iii)",
        sectionId := "frdi-prop-2-1" },
    classification := ABC3.Meta.GapClass.missingMath,
    falsifier :=
      "★原文は `Let d ∈N≥1.` と冒頭で `d` を固定したうえで (iii) を述べるが、" ++
      "`perfect object`(Definition 1.2, (iv))は**すべての n ∈ N≥1** を量化する。" ++
      "`Ψ_d` が圏同値であることから出るのは n = d の場合だけで、" ++
      "次数は乗法的なので 1 つの d から他の n は作れない(d の冪しか出ない)。" ++
      "★これが覆るのは (a) 1 つの d から全次数の Frobenius 型射の終域が作れることが" ++
      "示された場合(→ ①、我々の読み落とし)、または (b) 原文の意図が" ++
      "「すべての d について」だと確定した場合(→ 実質 ①、書き方の問題)である。" ++
      "★★我々の実装 `prop_2_1_iii` は `⟸` に「すべての n」を仮定しており、" ++
      "`⟹` の側は原文どおり d を固定したまま成り立つ(`prop_2_1_iii_mp`)。" ++
      "★`.src` は付けていない。" }

end ABC3.Gap.FrdI
