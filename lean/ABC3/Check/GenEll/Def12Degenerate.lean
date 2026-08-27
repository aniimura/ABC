import ABC3.Found.GenEll.Def12Height

/-!
# `ht` は算術直線束の `Pic` 類を見ていない —— **`Definition 1.2` の主張を下げる検査**

`Found/GenEll/Def12Height.lean` は `Definition 1.2` に**項目全体の `.src`** を置いていた。
★**それは過剰な主張だった**——本ファイルはその理由を機械にかける。

**これは原典の主張ではない**（我々のモデルについての事実）ので `.src` を持たない。

## ★★★何が起きているか

`Found/Arakelov/ADegEmb.lean` の正規化次数は

    degFOf F x = -(∑_{σ : F →+* ℂ} x.2.fn (embPoint F σ)) / [F:ℚ]

と定義されている。★**`x.1`（`Pic` 類）が式に現れない。**
したがって `htOf`・`ht` も `Pic` 類を見ない（`ht_indep_pic`）。

## ★★★★なぜ Interface が止められなかったか

`Interface/Arakelov/APic.lean` の `APicSpecData` は `deg_F` について
`degF_mul` / `degF_baseChange` / `degF_scale` の 3 欄しか持たない。
★**`deg_F` を有限素点（`ADiv(F)` の非アルキメデス側）と結ぶ欄が無い。**
だからアルキメデス側だけを見る `deg_F` が全欄を満たしてしまう。

★★これは 2026-08-27 に見つけた `APicData` の 2 つの穴
（`ι_X` 両立の欄が無い・`pullback_comp` の欄が無い）と**同じ種類**である。

## ★★★★★原文が言っていること

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if

★算術因子は**有限素点の整数係数**と**アルキメデス素点の実係数**の対である。
`deg_F` はその両方を測る。★★アルキメデス側だけでは `deg_F` ではない。

## ★★★★★★残っているもの・残っていないもの

| `Definition 1.2` の中身 | 状態 |
|---|---|
| (i) `ht` が `X(ℚ̄)` の関数として well-defined（底変換不変） | ★**真**（`htOf_baseChange`）——本検査は揺るがさない |
| (i) `ht_M̄` が算術直線束 `M̄` の高さであること | ★★**未**——`deg_F` が `Pic` 類を見ていない |
| (ii) BD-類 | ★**真**（`BDClass.lean`、`Pic` に依らない） |

★★★したがって `Definition 1.2` は**条つき**に下げる。
well-defined 性（原文の "any morphism"）は取れている——そこは失われていない。
-/

namespace ABC3.Check.GenEll

open AlgebraicGeometry CategoryTheory NumberField
open ABC3.Found.Arakelov ABC3.Interface.Arakelov ABC3.Found.GenEll

/-- ★★★★★★**退化** —— `deg_F` は `Pic` 類を見ていない。

★`rfl` で落ちる——式に `x.1` が現れないからである。 -/
theorem degFOf_indep_pic (F : Type) [Field F] [NumberField F]
    (a b : picardDataWitness.Pic (Spec (CommRingCat.of (𝓞 F))))
    (m : Multiplicative (arcCM (Spec (CommRingCat.of (𝓞 F))))) :
    degFOf F (a, m) = degFOf F (b, m) := rfl

/-- ★★★★★★**その帰結** —— 代表元の水準の高さも `Pic` 類を見ていない。 -/
theorem htOf_indep_pic (X : Scheme.{0})
    (A B : picardDataWitness.Pic X) (m : Multiplicative (arcCM X)) (x : NFPoint X) :
    htOf X (A, m) x = htOf X (B, m) x := rfl

/-- ★★★★★★★★**`ht` は算術直線束の `Pic` 類を見ていない**。

★★★これが `Definition 1.2` の項目全体の `.src` を下げる理由である。
原文の `ht_M̄` は `M̄` の高さでなければならないが、
本実装の `ht` は**計量だけ**で決まってしまう。

★★well-defined 性（`htOf_baseChange`）は**別の話であり、そちらは真である**
——`Definition 1.2, (i)` の原文「any morphism」の要求はそこにある。 -/
theorem ht_indep_pic (X : Scheme.{0})
    (A B : picardDataWitness.Pic X) (m : Multiplicative (arcCM X)) (x : NFAlgPoint X) :
    ht X (A, m) x = ht X (B, m) x := by
  induction x using Quot.inductionOn with
  | _ y => exact htOf_indep_pic X A B m y

end ABC3.Check.GenEll
