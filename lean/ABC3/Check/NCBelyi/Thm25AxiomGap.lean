import ABC3.Skeleton.NCBelyi.Theorem25

/-!
# ★★★[NCBelyi] `Theorem 2.5` も現在の `Interface` の下では**偽**である

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

## ★★★★存在主張は「射の型が空」で破れる

`Prop17AxiomGap.lean`(`hyp := True` で仮定を無力化)と
`Thm21AxiomGap.lean`(`Reduced := False` で片側を空虚化)に続く**第 3 の形**である。

★`BelyiSetup` の `Mor : Curve → Type` は**公理を持たない型のフィールド**である。
`Mor X := Empty` と置けば `∃ φ : B.Mor X, …` は**必ず偽**になる。
一方、仮定(`S.Finite`・`T.Finite`・`S ∩ T = ∅`)は `S = T = ∅` で満たされる。

★★★**存在主張を posit した型の上で書くと、その型が空でないことすら保証されない。**
`HeightTheoryData` のときは値が縛られないことが問題だったが、
ここでは**型そのものが空になれる**——同じ病の別の面である。

## ★★★★★では何をすべきか

`Remark 1.5.1` と同じ治療である——**構成に載せ替える**。
`Found/NCBelyi/` 側には既に

| 部品 | 場所 |
|---|---|
| `Lemma 2.1`(計算核) | `Lemma21.lean` |
| `Lemma 2.2`(有理点版、★`sorry` 無しで `Skeleton` にある) | `Lemma22.lean` |
| `Lemma 2.3`(分離、★`sorry` 無しで `Skeleton` にある) | `Separation.lean` |
| `Lemma 2.4` の (a)(c) | `DescendData.lean`(第 417-418) |

がある。★残るのは `Lemma 2.4` の (b) と、`X = ℙ¹` への帰着(曲線の Riemann–Roch)である。
-/

namespace ABC3.Check.NCBelyi

open ABC3.Interface.NCBelyi ABC3.Skeleton.NCBelyi

/-! ## ★射の型が空な setup -/

/-- ★**`Mor X := Empty` の `BelyiSetup`** —— 存在主張を破る。 -/
def badBelyi : BelyiSetup where
  Curve := Unit
  Point := fun _ => Unit
  ProjPoint := Unit
  three := ∅
  Mor := fun _ => Empty
  app := fun {_} e _ => e.elim
  UnramifiedOutsideThree := fun {_} _ => True
  NumField := Unit
  CurveOver := fun _ _ => True
  PointsOver := fun {_} _ _ => True
  MorOver := fun {_} _ _ => True

/-! ## ★★★反例 -/

/-- ★★★★**`Theorem 2.5` の statement は任意の `BelyiSetup` については偽**。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

`badBelyi` では `Mor X = Empty` なので `∃ φ : B.Mor X, …` は偽であり、
仮定は `S = T = ∅` で満たされる。

★★**この `sorry` も「まだ証明していない」ではなく「証明できない」である**。 -/
theorem theorem_2_5_false :
    ¬ (∀ (B : BelyiSetup) (X : B.Curve) (S T : Set (B.Point X)),
        S.Finite → T.Finite → S ∩ T = ∅ →
        (∃ φ : B.Mor X, B.UnramifiedOutsideThree φ
            ∧ (∀ s ∈ S, B.app φ s ∈ B.three)
            ∧ (∀ t ∈ T, B.app φ t ∉ B.three))
      ∧ (∀ F : B.NumField, B.CurveOver X F → B.PointsOver S F → B.PointsOver T F →
            ∃ φ : B.Mor X, B.MorOver φ F ∧ B.UnramifiedOutsideThree φ
              ∧ (∀ s ∈ S, B.app φ s ∈ B.three)
              ∧ (∀ t ∈ T, B.app φ t ∉ B.three))) := by
  intro h
  obtain ⟨⟨φ, -⟩, -⟩ :=
    h badBelyi () ∅ ∅ Set.finite_empty Set.finite_empty (by simp)
  exact φ.elim

end ABC3.Check.NCBelyi
