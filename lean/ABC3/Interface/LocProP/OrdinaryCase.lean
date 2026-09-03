import Mathlib.GroupTheory.QuotientGroup.Defs

/-!
# [LocProP] §1 —— The Ordinary Case の posit(`Interface`)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、物理 p.17-19。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

`Definition 1.1`(ordinary な `X_K`)・`Lemma 1.2`(節 `θ_X` の一意性)・
`Lemma 1.3`(`Δ_X → Δ_X^et` の核の群論的特徴づけ)・`Lemma 1.4`(`α_x^et = θ_X`)は
すべて、`Δ_X ⊴ Π_X_K → Γ_K` から作られる商 `Δ_X^et`・`Π_X_K^et` とその上の節 `θ_X` を
土台にする。★★**その土台(Jacobian・ordinary abelian variety・étale cohomology 経由の
構成)は mathlib に丸ごと無い**ので、`Interface` として posit する
(`Definition 0.3`・`Lemma 0.4` と同じ扱い)。

★★★**逸脱**: 位相(`Δ_X`・`Π_X` が profinite であること)はこの 4 項目の主張には
使わないので posit から省いた(`Remark Concerning the Term "Group-Theoretic"`(p.19)が
「この節の結果は本論文の主定理の証明では使わない」と明言している箇所と対応する)。
`Lemma 1.3` の条件 `(∗)^et` も抽象化(`satisfiesConditionEt`)して posit する。
-/

namespace ABC3.Interface.LocProP

universe u v w x

/-- 位相を持たない自明群(non-vacuous witness 用)。 -/
instance trivialGroup : Group Unit where
  mul _ _ := ()
  one := ()
  inv _ := ()
  mul_assoc := by intros; rfl
  one_mul := by intros; rfl
  mul_one := by intros; rfl
  inv_mul_cancel := by intros; rfl

/-- ★posit —— §1(The Ordinary Case)の骨組み。 -/
structure OrdinaryCaseSetup where
  Delta : Type u
  [deltaGroup : Group Delta]
  Gamma : Type w
  [gammaGroup : Group Gamma]
  /-- **[LocProP] Definition 1.1**(posit)。 -/
  isOrdinary : Prop
  /-- `Π_X_K^et`。 -/
  PiEt : Type v
  [piEtGroup : Group PiEt]
  /-- `Δ_X^et`。 -/
  DeltaEt : Type u
  [deltaEtGroup : Group DeltaEt]
  piEtToGamma : PiEt →* Gamma
  deltaEtIncl : DeltaEt →* PiEt
  /-- `θ_X : Γ_K → Π_X_K^et`(節)。 -/
  theta : Gamma →* PiEt
  theta_section : ∀ g, piEtToGamma (theta g) = g
  theta_commutes : ∀ g d, theta g * deltaEtIncl d = deltaEtIncl d * theta g
  /-- **[LocProP] Lemma 1.2**。 -/
  theta_unique : ∀ (θ' : Gamma →* PiEt), (∀ g, piEtToGamma (θ' g) = g) →
      (∀ g d, θ' g * deltaEtIncl d = deltaEtIncl d * θ' g) → θ' = theta
  Delta_toDeltaEt : Delta →* DeltaEt
  /-- 条件 `(∗)^et`(posit、抽象化)。 -/
  satisfiesConditionEt : Subgroup Delta → Prop
  /-- **[LocProP] Lemma 1.3**。 -/
  ker_eq_iInter : Delta_toDeltaEt.ker
      = ⨅ H : {H : Subgroup Delta // satisfiesConditionEt H}, (H : Subgroup Delta)
  /-- `X_K(K)` の点の集合(添字だけ posit する)。 -/
  PtX : Type x
  /-- `α_x^et : Γ_K → Π_X_K^et`(点 `x` ごと)。 -/
  alphaEt : PtX → Gamma →* PiEt
  /-- **[LocProP] Lemma 1.4**。 -/
  alphaEt_eq_theta : ∀ x, alphaEt x = theta

/-- `Subgroup Unit` は(自明群の部分群は 1 つしか無いので)`Subsingleton`。
`ker_eq_iInter` を機械的な等式展開ではなく一発で潰すために使う。 -/
instance subgroupUnitSubsingleton : Subsingleton (Subgroup Unit) :=
  ⟨fun H H' => by
    apply SetLike.ext
    intro x
    simp [Subsingleton.elim x (1 : Unit)]⟩

/-- ★★**非退化性の witness**(具体項、`Classical.choice` を経由しない)——
すべて自明群に潰す(θ_unique・ker_eq_iInter とも `Unit` 上で機械的に成り立つ)。 -/
@[reducible] def OrdinaryCaseSetup.example : OrdinaryCaseSetup.{0,0,0,0} where
  Delta := Unit; deltaGroup := inferInstance
  Gamma := Unit; gammaGroup := inferInstance
  isOrdinary := True
  PiEt := Unit; piEtGroup := inferInstance
  DeltaEt := Unit; deltaEtGroup := inferInstance
  piEtToGamma := 1; deltaEtIncl := 1; theta := 1
  theta_section := fun _ => rfl
  theta_commutes := fun _ _ => rfl
  theta_unique := fun _ _ _ => Subsingleton.elim _ _
  Delta_toDeltaEt := 1
  satisfiesConditionEt := fun _ => True
  ker_eq_iInter := Subsingleton.elim _ _
  PtX := Unit
  alphaEt := fun _ => 1
  alphaEt_eq_theta := fun _ => rfl

@[reducible] def OrdinaryCaseSetup.nonvacuous : Nonempty (OrdinaryCaseSetup.{0,0,0,0}) :=
  ⟨OrdinaryCaseSetup.example⟩

end ABC3.Interface.LocProP
