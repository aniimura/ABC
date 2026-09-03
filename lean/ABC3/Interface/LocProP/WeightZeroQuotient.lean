import ABC3.Meta.Claim
import ABC3.Interface.LocProP.GaloisCohomologyReview
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs

/-!
# [LocProP] §3 —— The Weight Zero Quotient の posit(`Interface`)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、物理 p.22-27。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

`Lemma 3.1`(`c_X` の核)・`Definition 3.2`/`3.4`(weight zero quotient `Z_X` の命名)・
`Proposition 3.3`/`3.5`(拡大類の消滅・完全列の分裂)は、Malčev completion・
Tannakian category・unipotent algebraic group の理論を土台にする。
★★**mathlib にこの理論は丸ごと無い**ので posit する(§0-§2 と同じ扱い)。
-/

namespace ABC3.Interface.LocProP

universe u

/-- ★posit —— §3(The Weight Zero Quotient)の骨組み。 -/
structure WeightZeroQuotientSetup where
  /-- `H_{Qp}`(§2 末尾の `H` を `Q_p` 化したもの)。 -/
  HQp : Type u
  HQpGrp : AddCommGroup HQp
  /-- `𝔐_X[1]`。 -/
  M1 : Type u
  M1Grp : AddCommGroup M1
  Wedge2HQp : Type u
  Wedge2HQpGrp : AddCommGroup Wedge2HQp
  /-- `c_X : ∧² H_{Qp} → 𝔐_X[1]`(交換子写像)。 -/
  cX : Wedge2HQp →+ M1
  r : ℕ
  IntersectionDual : Type u
  IntersectionDualGrp : AddCommGroup IntersectionDual
  intersectionMap : IntersectionDual →+ Wedge2HQp
  /-- **[LocProP] Lemma 3.1**。 -/
  lemma_3_1 : (r > 0 → ∀ x, cX x = 0 → x = 0) ∧
      (r = 0 → ∀ x, cX x = 0 ↔ ∃ y, intersectionMap y = x)
  /-- **[LocProP] Definition 3.2 / Definition 3.4** —— `Z_X`(weight zero quotient)。
  両項目とも同じ対象を(パラメータ付き版・単一曲線版として)命名するだけなので
  同じ posited field で受ける。 -/
  ZX : Type u
  ZXGrp : AddCommGroup ZX
  /-- **[LocProP] Proposition 3.3** —— `η_α = 0`(拡大類の消滅)。 -/
  EtaGroup : Type u
  EtaGroupGrp : AddCommGroup EtaGroup
  etaAlpha : EtaGroup
  prop_3_3 : etaAlpha = 0
  /-- `0 → 𝔷_X[1] → 𝔷_X → 𝔷_X[0] → 0` の completion。 -/
  Z1 : Type u
  Z1Grp : AddCommGroup Z1
  Z0 : Type u
  Z0Grp : AddCommGroup Z0
  zIncl : Z1 →+ ZX
  zProj : ZX →+ Z0
  zSplit : Z0 →+ ZX
  /-- **[LocProP] Proposition 3.5** —— 完全列が分裂する。 -/
  prop_3_5 : ∀ z, zProj (zSplit z) = z

/-- ★★**非退化性の witness**(具体項)。 -/
@[reducible] def WeightZeroQuotientSetup.example : WeightZeroQuotientSetup.{0} where
  HQp := ℤ; HQpGrp := inferInstance
  M1 := ℤ; M1Grp := inferInstance
  Wedge2HQp := ℤ; Wedge2HQpGrp := inferInstance
  cX := AddMonoidHom.id ℤ
  r := 1
  IntersectionDual := Unit
  IntersectionDualGrp := inferInstance
  intersectionMap := 0
  lemma_3_1 := ⟨fun _ _ h => h, fun hr => absurd hr (by omega)⟩
  ZX := ℤ; ZXGrp := inferInstance
  EtaGroup := ℤ; EtaGroupGrp := inferInstance
  etaAlpha := 0
  prop_3_3 := rfl
  Z1 := ℤ; Z1Grp := inferInstance
  Z0 := ℤ; Z0Grp := inferInstance
  zIncl := AddMonoidHom.id ℤ
  zProj := AddMonoidHom.id ℤ
  zSplit := AddMonoidHom.id ℤ
  prop_3_5 := fun _ => rfl

@[reducible] def WeightZeroQuotientSetup.nonvacuous : Nonempty (WeightZeroQuotientSetup.{0}) :=
  ⟨WeightZeroQuotientSetup.example⟩

end ABC3.Interface.LocProP
