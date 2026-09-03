import ABC3.Meta.Claim
import ABC3.Interface.Falt1.AlmostMath
import Mathlib.Algebra.Exact.Basic
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs

/-!
# [Falt1] Chapter I §4, Theorem 4.1 の posit(`Interface`)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)。物理 p.13-14(印字 p.266-267)。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

★注: OCR 層は数式記号を激しく壊すため地の文で写す
(Interface/Falt1/Ramification.lean 冒頭の注記と同じ運用)。

微分加群 `Ω_{R/V}`・`Ω̄_{R/V}` 等そのものは posit する(§1 の微分加群論の
延長で、almost mathematics の枠組みに乗せた具体的な計算——mathlib への
対応物の有無は Found 着手時に調査)。
-/

namespace ABC3.Interface.Falt1

universe u

structure DifferentialsSetup where
  -- Theorem 4.1(i) の1つ目の almost isomorphism
  -- `Ω_{R/V} ⊗_R R̄[1/p] ≅ Ω̄_{R/V}`(核・余核が m で零化される)。
  Ker41i1 : Type u
  Ker41i1Grp : AddCommGroup Ker41i1
  Coker41i1 : Type u
  Coker41i1Grp : AddCommGroup Coker41i1
  thm41i1_ker : MTorsionWitness Ker41i1
  thm41i1_coker : MTorsionWitness Coker41i1
  -- Theorem 4.1(i) の2つ目の almost isomorphism
  -- `Ω_{R/V} ⊗_R (R̄[1/p]/R̄) ≅ Ω̄_{R/R⊗_V V̄}`。
  Ker41i2 : Type u
  Ker41i2Grp : AddCommGroup Ker41i2
  Coker41i2 : Type u
  Coker41i2Grp : AddCommGroup Coker41i2
  thm41i2_ker : MTorsionWitness Ker41i2
  thm41i2_coker : MTorsionWitness Coker41i2
  -- Theorem 4.1(ii): 完全列
  -- `0 → Ω_{V̄/V} ⊗_V̄ R̄ → Ω̄_{R̄/R} → Ω̄_{R̄/R⊗_V V̄} → 0`。
  A41ii : Type u
  A41iiGrp : AddCommGroup A41ii
  B41ii : Type u
  B41iiGrp : AddCommGroup B41ii
  C41ii : Type u
  C41iiGrp : AddCommGroup C41ii
  thm41ii_f : A41ii →+ B41ii
  thm41ii_g : B41ii →+ C41ii
  thm41ii_injective : Function.Injective thm41ii_f
  thm41ii_exact : Function.Exact thm41ii_f thm41ii_g
  thm41ii_surjective : Function.Surjective thm41ii_g

@[reducible] def DifferentialsSetup.example : DifferentialsSetup.{0} where
  Ker41i1 := ℤ; Ker41i1Grp := inferInstance
  Coker41i1 := ℤ; Coker41i1Grp := inferInstance
  thm41i1_ker := MTorsionWitness.trivial
  thm41i1_coker := MTorsionWitness.trivial
  Ker41i2 := ℤ; Ker41i2Grp := inferInstance
  Coker41i2 := ℤ; Coker41i2Grp := inferInstance
  thm41i2_ker := MTorsionWitness.trivial
  thm41i2_coker := MTorsionWitness.trivial
  A41ii := ℤ; A41iiGrp := inferInstance
  B41ii := ℤ; B41iiGrp := inferInstance
  C41ii := PUnit; C41iiGrp := inferInstance
  thm41ii_f := AddMonoidHom.id ℤ
  thm41ii_g := 0
  thm41ii_injective := fun _ _ h => h
  thm41ii_exact := by
    intro b
    constructor
    · intro _; exact ⟨b, rfl⟩
    · intro _; trivial
  thm41ii_surjective := fun _ => ⟨0, Subsingleton.elim _ _⟩

@[reducible] def DifferentialsSetup.nonvacuous : Nonempty (DifferentialsSetup.{0}) :=
  ⟨DifferentialsSetup.example⟩

end ABC3.Interface.Falt1
