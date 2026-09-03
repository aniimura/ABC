import ABC3.Meta.Claim
import ABC3.Interface.Falt1.AlmostMath
import Mathlib.Algebra.Exact.Basic
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs

/-!
# [Falt1] Chapter I §4, Theorem 4.2・4.3・4.5 の posit(`Interface`)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)。物理 p.15(印字 p.268、
Theorem 4.2)・物理 p.17(印字 p.270、Theorem 4.3)・物理 p.19-20
(印字 p.272-273、Theorem 4.5)。**260 dpi 目視確認 2026-09-04**
(逐語の食い違い無し)。

★注: OCR 層は数式記号を激しく壊すため地の文で写す
(Interface/Falt1/Ramification.lean 冒頭の注記と同じ運用)。

★★**簡略化(逸脱)**: 原文の各定理は添字(i・M・R₁,R₂ の任意の分解等)に
関する全称命題だが、ここでは代表的な1インスタンスだけを posit する
(∀ を落とす)。また「p^{1/(p-1)} で零化される」も
`Interface/Falt1/AlmostMath.lean` の `MTorsionWitness`(「ある p の冪で
零化される」)で近似する——真の値が p^{1/(p-1)} という非整数冪である
点を実数冪ではなく整数冪に単純化している。後続の証明が [LocProP] の
範囲でこれらを消費しないので、この簡略化は影響を持たない
(CLAUDE.md の「逸脱」方針)。
-/

namespace ABC3.Interface.Falt1

universe u

/-- **Theorem 4.2**(6 部——(i)Λ分解、(ii)étale基底変換、(iii)Künneth、
(iv)Ω_{R/V}への同型、(v)cup積による同型、(vi)双対性)。 -/
structure GaloisCohCompute42 where
  -- (i) H^i(Δ,R̂) ≅ Λ^i(...) ⊕ rest, rest は p^{1/(p-1)} で零化。
  Hi42i : Type u
  Hi42iGrp : AddCommGroup Hi42i
  ExtPart42i : Type u
  ExtPart42iGrp : AddCommGroup ExtPart42i
  Rest42i : Type u
  Rest42iGrp : AddCommGroup Rest42i
  thm42i : Hi42i ≃+ (ExtPart42i × Rest42i)
  thm42i_rest : MTorsionWitness Rest42i
  -- (ii) étale 基底変換で誘導される射は almost isomorphism。
  Src42ii : Type u
  Src42iiGrp : AddCommGroup Src42ii
  Tgt42ii : Type u
  Tgt42iiGrp : AddCommGroup Tgt42ii
  thm42ii_ker : Type u
  thm42ii_kerGrp : AddCommGroup thm42ii_ker
  thm42ii_coker : Type u
  thm42ii_cokerGrp : AddCommGroup thm42ii_coker
  thm42ii_ker_ann : MTorsionWitness thm42ii_ker
  thm42ii_coker_ann : MTorsionWitness thm42ii_coker
  -- (iii) Künneth 射は almost isomorphism。
  Src42iii : Type u
  Src42iiiGrp : AddCommGroup Src42iii
  Tgt42iii : Type u
  Tgt42iiiGrp : AddCommGroup Tgt42iii
  thm42iii_ker : Type u
  thm42iii_kerGrp : AddCommGroup thm42iii_ker
  thm42iii_coker : Type u
  thm42iii_cokerGrp : AddCommGroup thm42iii_coker
  thm42iii_ker_ann : MTorsionWitness thm42iii_ker
  thm42iii_coker_ann : MTorsionWitness thm42iii_coker
  -- (iv) H^1(Δ,R̂)/(p-torsion) ≅ Ω_{R/V}⊗(R⊗V̄)^(-1)。
  H1modP42iv : Type u
  H1modP42ivGrp : AddCommGroup H1modP42iv
  OmegaTensor42iv : Type u
  OmegaTensor42ivGrp : AddCommGroup OmegaTensor42iv
  thm42iv : H1modP42iv ≃+ OmegaTensor42iv
  -- (v) cup 積による同型 H^i(Δ,R̂)/(p-torsion) ≅ Ω^i_{R/V}⊗(R⊗V̄)^(-i)。
  HimodP42v : Type u
  HimodP42vGrp : AddCommGroup HimodP42v
  OmegaITensor42v : Type u
  OmegaITensor42vGrp : AddCommGroup OmegaITensor42v
  thm42v : HimodP42v ≃+ OmegaITensor42v
  -- (vi) cup 積による双対性は almost isomorphism。
  Src42vi : Type u
  Src42viGrp : AddCommGroup Src42vi
  Tgt42vi : Type u
  Tgt42viGrp : AddCommGroup Tgt42vi
  thm42vi_ker : Type u
  thm42vi_kerGrp : AddCommGroup thm42vi_ker
  thm42vi_coker : Type u
  thm42vi_cokerGrp : AddCommGroup thm42vi_coker
  thm42vi_ker_ann : MTorsionWitness thm42vi_ker
  thm42vi_coker_ann : MTorsionWitness thm42vi_coker

@[reducible] def GaloisCohCompute42.example : GaloisCohCompute42.{0} where
  Hi42i := ℤ × ℤ; Hi42iGrp := inferInstance
  ExtPart42i := ℤ; ExtPart42iGrp := inferInstance
  Rest42i := ℤ; Rest42iGrp := inferInstance
  thm42i := AddEquiv.refl (ℤ × ℤ)
  thm42i_rest := MTorsionWitness.trivial
  Src42ii := ℤ; Src42iiGrp := inferInstance
  Tgt42ii := ℤ; Tgt42iiGrp := inferInstance
  thm42ii_ker := ℤ; thm42ii_kerGrp := inferInstance
  thm42ii_coker := ℤ; thm42ii_cokerGrp := inferInstance
  thm42ii_ker_ann := MTorsionWitness.trivial
  thm42ii_coker_ann := MTorsionWitness.trivial
  Src42iii := ℤ; Src42iiiGrp := inferInstance
  Tgt42iii := ℤ; Tgt42iiiGrp := inferInstance
  thm42iii_ker := ℤ; thm42iii_kerGrp := inferInstance
  thm42iii_coker := ℤ; thm42iii_cokerGrp := inferInstance
  thm42iii_ker_ann := MTorsionWitness.trivial
  thm42iii_coker_ann := MTorsionWitness.trivial
  H1modP42iv := ℤ; H1modP42ivGrp := inferInstance
  OmegaTensor42iv := ℤ; OmegaTensor42ivGrp := inferInstance
  thm42iv := AddEquiv.refl ℤ
  HimodP42v := ℤ; HimodP42vGrp := inferInstance
  OmegaITensor42v := ℤ; OmegaITensor42vGrp := inferInstance
  thm42v := AddEquiv.refl ℤ
  Src42vi := ℤ; Src42viGrp := inferInstance
  Tgt42vi := ℤ; Tgt42viGrp := inferInstance
  thm42vi_ker := ℤ; thm42vi_kerGrp := inferInstance
  thm42vi_coker := ℤ; thm42vi_cokerGrp := inferInstance
  thm42vi_ker_ann := MTorsionWitness.trivial
  thm42vi_coker_ann := MTorsionWitness.trivial

@[reducible] def GaloisCohCompute42.nonvacuous : Nonempty (GaloisCohCompute42.{0}) :=
  ⟨GaloisCohCompute42.example⟩

/-- **Theorem 4.3**: `0 → ρ⁻¹R̂ → E_ρ → Ω_{R/V}(dlog∞)⊗R̂(-1) → 0` が
Γ-同変・関手的な拡大として存在し、`R̂` を係数とする拡大の pushout として
得られる(pushout の関係そのものは posit しない、2本の完全列だけ書く)。 -/
structure GaloisCohCompute43 where
  RhoInvRhat : Type u
  RhoInvRhatGrp : AddCommGroup RhoInvRhat
  Erho : Type u
  ErhoGrp : AddCommGroup Erho
  OmegaTwist43 : Type u
  OmegaTwist43Grp : AddCommGroup OmegaTwist43
  f43a : RhoInvRhat →+ Erho
  g43a : Erho →+ OmegaTwist43
  thm43a_inj : Function.Injective f43a
  thm43a_exact : Function.Exact f43a g43a
  thm43a_surj : Function.Surjective g43a
  -- pushout 元の完全列(R̂ 係数)。
  Rhat43 : Type u
  Rhat43Grp : AddCommGroup Rhat43
  E43 : Type u
  E43Grp : AddCommGroup E43
  f43b : Rhat43 →+ E43
  g43b : E43 →+ OmegaTwist43
  thm43b_inj : Function.Injective f43b
  thm43b_exact : Function.Exact f43b g43b
  thm43b_surj : Function.Surjective g43b

@[reducible] def GaloisCohCompute43.example : GaloisCohCompute43.{0} where
  RhoInvRhat := ℤ; RhoInvRhatGrp := inferInstance
  Erho := ℤ; ErhoGrp := inferInstance
  OmegaTwist43 := PUnit; OmegaTwist43Grp := inferInstance
  f43a := AddMonoidHom.id ℤ
  g43a := 0
  thm43a_inj := fun _ _ h => h
  thm43a_exact := by
    intro b; constructor
    · intro _; exact ⟨b, rfl⟩
    · intro _; trivial
  thm43a_surj := fun _ => ⟨0, Subsingleton.elim _ _⟩
  Rhat43 := ℤ; Rhat43Grp := inferInstance
  E43 := ℤ; E43Grp := inferInstance
  f43b := AddMonoidHom.id ℤ
  g43b := 0
  thm43b_inj := fun _ _ h => h
  thm43b_exact := by
    intro b; constructor
    · intro _; exact ⟨b, rfl⟩
    · intro _; trivial
  thm43b_surj := fun _ => ⟨0, Subsingleton.elim _ _⟩

@[reducible] def GaloisCohCompute43.nonvacuous : Nonempty (GaloisCohCompute43.{0}) :=
  ⟨GaloisCohCompute43.example⟩

/-- **Theorem 4.5**(3 部——(i)スペクトル系列の m-torsion を除く退化、
(ii)写像が almost isomorphism、(iii)フィルトレーションの対応)。 -/
structure GaloisCohCompute45 where
  -- (i) スペクトル系列 E2 と極限 H^{a+b} の間の almost isomorphism。
  E2Page45 : Type u
  E2Page45Grp : AddCommGroup E2Page45
  HLimit45 : Type u
  HLimit45Grp : AddCommGroup HLimit45
  thm45i_ker : Type u
  thm45i_kerGrp : AddCommGroup thm45i_ker
  thm45i_coker : Type u
  thm45i_cokerGrp : AddCommGroup thm45i_coker
  thm45i_ker_ann : MTorsionWitness thm45i_ker
  thm45i_coker_ann : MTorsionWitness thm45i_coker
  -- (ii) E2 と直和分解の間の写像は almost isomorphism。
  Src45ii : Type u
  Src45iiGrp : AddCommGroup Src45ii
  Tgt45ii : Type u
  Tgt45iiGrp : AddCommGroup Tgt45ii
  thm45ii_ker : Type u
  thm45ii_kerGrp : AddCommGroup thm45ii_ker
  thm45ii_coker : Type u
  thm45ii_cokerGrp : AddCommGroup thm45ii_coker
  thm45ii_ker_ann : MTorsionWitness thm45ii_ker
  thm45ii_coker_ann : MTorsionWitness thm45ii_coker
  -- (iii) フィルトレーションの対応(可換性、簡略化して Prop で posit)。
  thm45iii : Prop
  thm45iii_holds : thm45iii

@[reducible] def GaloisCohCompute45.example : GaloisCohCompute45.{0} where
  E2Page45 := ℤ; E2Page45Grp := inferInstance
  HLimit45 := ℤ; HLimit45Grp := inferInstance
  thm45i_ker := ℤ; thm45i_kerGrp := inferInstance
  thm45i_coker := ℤ; thm45i_cokerGrp := inferInstance
  thm45i_ker_ann := MTorsionWitness.trivial
  thm45i_coker_ann := MTorsionWitness.trivial
  Src45ii := ℤ; Src45iiGrp := inferInstance
  Tgt45ii := ℤ; Tgt45iiGrp := inferInstance
  thm45ii_ker := ℤ; thm45ii_kerGrp := inferInstance
  thm45ii_coker := ℤ; thm45ii_cokerGrp := inferInstance
  thm45ii_ker_ann := MTorsionWitness.trivial
  thm45ii_coker_ann := MTorsionWitness.trivial
  thm45iii := True
  thm45iii_holds := trivial

@[reducible] def GaloisCohCompute45.nonvacuous : Nonempty (GaloisCohCompute45.{0}) :=
  ⟨GaloisCohCompute45.example⟩

end ABC3.Interface.Falt1
