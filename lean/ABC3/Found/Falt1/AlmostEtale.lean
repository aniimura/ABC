import ABC3.Meta.Claim
import Mathlib.RingTheory.Etale.Basic
import Mathlib.RingTheory.Unramified.Finite
import Mathlib.RingTheory.Unramified.Pi
import Mathlib.RingTheory.Smooth.Pi
import Mathlib.RingTheory.RingHom.Etale
import Mathlib.RingTheory.Flat.Basic
import Mathlib.RingTheory.Finiteness.FinitePresentationLocal
import Mathlib.RingTheory.Trace.Defs
import Mathlib.RingTheory.Localization.Away.Basic
import Mathlib.RingTheory.Localization.NormTrace
import Mathlib.RingTheory.Localization.Finiteness
import Mathlib.RingTheory.LocalProperties.Projective
import Mathlib.RingTheory.TensorProduct.Basic
import Mathlib.RingTheory.Ideal.Norm.RelNorm
import Mathlib.Algebra.Homology.DerivedCategory.Ext.EnoughProjectives
import Mathlib.Algebra.Category.ModuleCat.Projective

/-!
# [Falt1] Definition 2.1(almost étale covering の定義)、Found(2026-09-04)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §2、物理 p.6
(印字 p.259)。`ABC3.Interface.Falt1.AlmostEtaleSetup.isAlmostEtale` が
posit していた `Prop` を、ここで mathlib の実在する道具
(`Module.Free`・`Module.Finite`・`Algebra.Etale`・`Algebra.trace`・
`Algebra.FormallyUnramified.elem`)だけを使って**実際に定義する**。

内容 (Falt1 p.6、260dpi 目視。OCR 層は数式記号を激しく壊すため地の文で
写す): Suppose A is a ring, B an A-algebra. B is called an almost étale
covering of A if (i) B[1/p] is a projective A[1/p]-module of finite rank
and an étale A[1/p]-algebra; (ii) the trace map tr_{B/A}: B[1/p]→A[1/p]
maps B into A; (iii) for the idempotent e_{B/A} corresponding to the
diagonal, p^ε e_{B/A} lies in the image of B⊗_A B for any ε>0.

## 逸脱(記録)

**(i) の "projective" を "free" に強めている**——mathlib の
`Algebra.trace`(`RingTheory/Trace/Defs.lean`)は `[Module.Free R S]` を
要求し(`if ∃ 基底 then ... else 0` という `LinearMap.trace` の定義に
由来)、有限射影加群(free とは限らない)の trace 構成は mathlib に無い
(2026-09-04 実測、`ResearchPaper/mathlib-gap.json` の
`falt1-almost-mathematics` に詳細記録)。`Theorem 2.2`-`2.4`(まだ
`Interface.AlmostEtaleSetup` の posit のまま、この定義を消費していない)
に影響しないため、CLAUDE.md の「逸脱」条項に従い許容する。

**(iii) の "B⊗_A B の像" は `diagonalCompare` という明示的な写像
(`TensorProduct.lift` 経由、局所化の可換性の補題は不要——ここが鍵の
発見)で捉える**。idempotent `e_{B/A}` そのものは mathlib の
`Algebra.FormallyUnramified.elem`(`RingTheory/Unramified/Finite.lean`)
がまさに対応する元を与える。
-/

namespace ABC3.Found.Falt1

/-- `A[1/p]` から `B[1/p]` への自然な `Algebra` 構造(`Localization.awayMap`
を `RingHom.toAlgebra` で配線)。`abbrev`/`instance` にはしない
(具体的な `TensorProduct` に対する既存のグローバルインスタンスと
衝突しうるため——`pushoutKaehlerSplit3` で確認した `Algebra.TensorProduct.
rightAlgebra` と同じ理由、tools/lean-idioms.md #23)。 -/
@[reducible] noncomputable def awayAlgebra {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] (p : A) :
    Algebra (Localization.Away p) (Localization.Away (algebraMap A B p)) :=
  (Localization.awayMap (algebraMap A B) p).toAlgebra

/-- `A → A[1/p] → B[1/p]` が scalar tower であること(`Localization.awayMap`
の自然性、`IsLocalization.map_eq` から)。 -/
theorem awayScalarTower {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] (p : A) :
    letI := awayAlgebra p (A := A) (B := B)
    IsScalarTower A (Localization.Away p) (Localization.Away (algebraMap A B p)) := by
  letI := awayAlgebra p (A := A) (B := B)
  apply IsScalarTower.of_algebraMap_eq
  intro x
  show algebraMap A (Localization.Away (algebraMap A B p)) x
    = Localization.awayMap (algebraMap A B) p (algebraMap A (Localization.Away p) x)
  rw [show Localization.awayMap (algebraMap A B) p (algebraMap A (Localization.Away p) x)
    = algebraMap A (Localization.Away (algebraMap A B p)) x from by
      unfold Localization.awayMap IsLocalization.Away.map
      rw [IsLocalization.map_eq]
      rw [← IsScalarTower.algebraMap_apply]]

/-- **`B⊗_A B → (A[1/p]の下で) B[1/p]⊗_{A[1/p]} B[1/p]` への比較写像**。
条件(iii)の「`B⊗_A B` の像」をこれで捉える。`TensorProduct.isPushout`
型の局所化と可換性の同一視は一切不要——`Algebra.TensorProduct.lift` で
2本の `AlgHom`(`includeLeft`・`includeRight` に `B → B[1/p]` を合成した
もの)を直接貼り合わせるだけで構成できることを発見した(局所化を経由
した込み入った同一視より遥かに単純)。 -/
noncomputable def diagonalCompare {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] (p : A) :
    letI := awayAlgebra p (A := A) (B := B)
    haveI := awayScalarTower p (A := A) (B := B)
    TensorProduct A B B →ₐ[A] TensorProduct (Localization.Away p) (Localization.Away (algebraMap A B p))
      (Localization.Away (algebraMap A B p)) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  set Bp := Localization.Away (algebraMap A B p)
  set Ap := Localization.Away p
  set f0 : B →ₐ[A] Bp := IsScalarTower.toAlgHom A B Bp
  set f : B →ₐ[A] TensorProduct Ap Bp Bp :=
    (Algebra.TensorProduct.includeLeft (R := Ap) (S := A) (A := Bp) (B := Bp)).comp f0
  set g : B →ₐ[A] TensorProduct Ap Bp Bp :=
    ((Algebra.TensorProduct.includeRight (R := Ap) (A := Bp) (B := Bp)).restrictScalars A).comp f0
  exact Algebra.TensorProduct.lift f g (fun x y => mul_comm (f x) (g y))

/-- **`[Falt1] Definition 2.1`**(`Module.Free` への逸脱つき、上記参照)。 -/
noncomputable def IsAlmostEtaleCovering {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] (p : A) : Prop :=
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  ∃ (_ : Module.Free (Localization.Away p) (Localization.Away (algebraMap A B p)))
    (_ : Module.Finite (Localization.Away p) (Localization.Away (algebraMap A B p)))
    (_ : Algebra.Etale (Localization.Away p) (Localization.Away (algebraMap A B p))),
    -- (ii) trace map maps B into A
    (∀ b : B, ∃ a : A, Algebra.trace (Localization.Away p) (Localization.Away (algebraMap A B p))
        (algebraMap B (Localization.Away (algebraMap A B p)) b)
      = algebraMap A (Localization.Away p) a) ∧
    -- (iii) p^ε e_{B/A} lies in the image of B ⊗_A B for any ε > 0
    (∀ n : ℕ, 0 < n → ∃ e : TensorProduct A B B, diagonalCompare p e
      = p^n • Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p)))

def isAlmostEtaleCovering.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 6, item := "Definition 2.1",
    sectionId := "falt1-def-2-1" }

/-- 恒等拡大(`B=R`、`A=R`ではなく一般の `R`)の場合、自明な `p:=1`
(`Localization.Away 1 ≃ R` なので実質 `Ap=Bp=R` に帰着)で
`Algebra.FormallyUnramified.elem R R = 1⊗ₜ1` であることの補助補題。
`isAlmostEtaleCovering` の非空虚性の witness を作る際に使う想定
(★2026-09-04 時点では `A=B`(恒等拡大)を witness に選ぶと `awayAlgebra`
(局所化の自己写像)と mathlib の標準的な自己代数(`Algebra.id`)が
衝突し行き詰まっていた)。

★2026-09-05: witness を `A:=R`・`B:=Fin 2 → R`・`p:=1` に取り換える
(`A≠B` なので上記の衝突する標準インスタンスがそもそも存在しない)
ことで、条件 (i) の 3 点(`Module.Free`・`Module.Finite`・
`Algebra.Etale`)を以下 `awayOne_fin2_etale`・`awayOne_fin2_freeFinite`
で sorry 無しに完成させた。鍵となったのは `AlgEquiv`/インスタンス
レベルの `▸`・`convert`・`Module.Finite.equiv` を諦め、
`RingHom.Etale`(bare ring hom の性質、instance 衝突が起きようが
ない)のレベルで `IsLocalization.ringHom_ext` により環準同型の等式を
直接示し、`RingHom.Etale.stableUnderComposition`/`of_bijective`/
`toAlgebra` で組み立てる戦略、および `Module.Finite.of_equiv_equiv`
(`≃ₗ` を経由しない、環準同型の可換四角形からの推移)。条件 (ii)
(trace)・(iii)(idempotent)は引き続き未着手。 -/
theorem elem_self {R : Type*} [CommRing R] : Algebra.FormallyUnramified.elem R R = (1:R) ⊗ₜ[R] (1:R) := by
  have hinj : Function.Injective (Algebra.TensorProduct.lmul' R : TensorProduct R R R →ₐ[R] R) := by
    have h : (Algebra.TensorProduct.lmul' R : TensorProduct R R R →ₐ[R] R).toLinearMap
        = (TensorProduct.lid R R : TensorProduct R R R ≃ₗ[R] R).toLinearMap := by
      apply TensorProduct.ext'
      intro x y
      simp [Algebra.TensorProduct.lmul'_apply_tmul]
    intro a b hab
    have hab' : (Algebra.TensorProduct.lmul' R : TensorProduct R R R →ₐ[R] R).toLinearMap a
        = (Algebra.TensorProduct.lmul' R : TensorProduct R R R →ₐ[R] R).toLinearMap b := hab
    rw [h] at hab'
    exact (TensorProduct.lid R R).injective hab'
  apply hinj
  rw [Algebra.FormallyUnramified.lmul_elem, Algebra.TensorProduct.lmul'_apply_tmul, mul_one]

/-- `Fin 2 → R` は `R` 上 étale(`Algebra.FormallyUnramified.pi_iff`・
`Algebra.FormallySmooth.pi_iff` で各成分 `R`(自明に unramified・smooth)
に帰着)。`p:=1` の witness で `B := Fin 2 → R` として使う補題。 -/
theorem etale_fin2 {R : Type*} [CommRing R] : Algebra.Etale R (Fin 2 → R) := by
  rw [Algebra.Etale.iff_formallyUnramified_and_smooth]
  refine ⟨?_, ?_, ?_⟩
  · rw [Algebra.FormallyUnramified.pi_iff]; intro i; infer_instance
  · rw [Algebra.FormallySmooth.pi_iff]; intro i; infer_instance
  · infer_instance

/-- **`A:=R`・`B:=Fin 2 → R`・`p:=1` での witness、条件 (i) の
`Algebra.Etale` 部分**。`p=1` は単元なので `algebraMap R (Localization.
Away 1)` と `algebraMap (Fin 2→R) (Localization.Away (algebraMap R
(Fin2→R) 1))` はどちらも全単射(`IsLocalization.atUnit`)。この2つの
`RingEquiv`(`ιR`・`ιB`)を用いて `awayAlgebra 1` の台になる環準同型
`Localization.awayMap (algebraMap R (Fin2→R)) 1` を `ιB ∘ algebraMap R
(Fin2→R) ∘ ιR.symm` に分解する式 `heqmap` を`IsLocalization.
ringHom_ext`(局所化からの環準同型の一意性)で示し、`RingHom.Etale`
(bare ring hom の étale 性、`AlgEquiv`/インスタンスの衝突が起こり
ようがない)の合成安定性(`stableUnderComposition`)と全単射からの
étale 性(`of_bijective`)を貼り合わせて結論する。 -/
theorem awayOne_fin2_etale {R : Type*} [CommRing R] :
    letI := awayAlgebra (1:R) (A := R) (B := Fin 2 → R)
    Algebra.Etale (Localization.Away (1:R)) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) := by
  have hbijR : Function.Bijective (algebraMap R (Localization.Away (1:R))) := by
    have e : R ≃ₐ[R] Localization.Away (1:R) := IsLocalization.atUnit R (Localization.Away (1:R)) (1:R) isUnit_one
    have heq : (e : R →+* Localization.Away (1:R)) = algebraMap R (Localization.Away (1:R)) := by
      ext x; exact e.commutes x
    rw [← heq]; exact e.bijective
  have hunitB : IsUnit (algebraMap R (Fin 2 → R) (1:R)) := by rw [map_one]; exact isUnit_one
  have hbijB : Function.Bijective (algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R)))) := by
    have e : (Fin 2 → R) ≃ₐ[Fin 2 → R] Localization.Away (algebraMap R (Fin 2 → R) (1:R)) :=
      IsLocalization.atUnit (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) (algebraMap R (Fin 2 → R) (1:R)) hunitB
    have heq : (e : (Fin 2 → R) →+* Localization.Away (algebraMap R (Fin 2 → R) (1:R)))
        = algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) :=
      RingHom.ext (fun x => e.commutes x)
    rw [← heq]; exact e.bijective
  set ιR := RingEquiv.ofBijective (algebraMap R (Localization.Away (1:R))) hbijR with hιRdef
  set ιB := RingEquiv.ofBijective (algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R)))) hbijB with hιBdef
  have heqmap : Localization.awayMap (algebraMap R (Fin 2 → R)) (1:R)
      = ιB.toRingHom.comp ((algebraMap R (Fin 2 → R)).comp ιR.symm.toRingHom) := by
    apply IsLocalization.ringHom_ext (Submonoid.powers (1:R))
    ext x
    have hιRx : ιR.symm (ιR x) = x := ιR.symm_apply_apply x
    show Localization.awayMap (algebraMap R (Fin 2 → R)) (1:R) (algebraMap R (Localization.Away (1:R)) x)
      = ιB (algebraMap R (Fin 2 → R) (ιR.symm (algebraMap R (Localization.Away (1:R)) x)))
    have hιRxeq : (algebraMap R (Localization.Away (1:R)) x) = ιR x := by rw [hιRdef]; rfl
    rw [hιRxeq, hιRx]
    show Localization.awayMap (algebraMap R (Fin 2 → R)) (1:R) (algebraMap R (Localization.Away (1:R)) x)
      = algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) (algebraMap R (Fin 2 → R) x)
    unfold Localization.awayMap IsLocalization.Away.map
    rw [IsLocalization.map_eq]
  have hgEtale : RingHom.Etale (algebraMap R (Fin 2 → R)) := RingHom.etale_algebraMap.mpr etale_fin2
  have hιRinvEtale : RingHom.Etale (ιR.symm.toRingHom) :=
    RingHom.Etale.of_bijective ιR.symm.bijective
  have hιBEtale : RingHom.Etale (ιB.toRingHom) :=
    RingHom.Etale.of_bijective ιB.bijective
  have hcompEtale : RingHom.Etale ((algebraMap R (Fin 2 → R)).comp ιR.symm.toRingHom) :=
    RingHom.Etale.stableUnderComposition ιR.symm.toRingHom (algebraMap R (Fin 2 → R)) hιRinvEtale hgEtale
  have hfullEtale : RingHom.Etale (ιB.toRingHom.comp ((algebraMap R (Fin 2 → R)).comp ιR.symm.toRingHom)) :=
    RingHom.Etale.stableUnderComposition ((algebraMap R (Fin 2 → R)).comp ιR.symm.toRingHom) ιB.toRingHom hcompEtale hιBEtale
  rw [← heqmap] at hfullEtale
  exact RingHom.Etale.toAlgebra hfullEtale

/-- **`A:=R`・`B:=Fin 2 → R`・`p:=1` での witness、条件 (i) の
`Module.Free`・`Module.Finite` 部分**。`awayOne_fin2_etale` と同じ
`ιR`・`ιB`・`heqmap` を再構成し、`ιB` を `R`-加群としての半線形同値
`(Fin 2→R) ≃ₛₗ[ιR.toRingHom] Localization.Away (algebraMap R (Fin2→R)
1)` に仕立てて `Module.Free.of_equiv`(半線形移送)、`ιR`・`ιB` からの
環準同型の可換四角形 `he` を `heqmap` から直接示して
`Module.Finite.of_equiv_equiv`(`≃ₗ` を経由しない代数レベルの移送)に
渡す。 -/
theorem awayOne_fin2_freeFinite {R : Type*} [CommRing R] :
    letI := awayAlgebra (1:R) (A := R) (B := Fin 2 → R)
    Module.Free (Localization.Away (1:R)) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) ∧
    Module.Finite (Localization.Away (1:R)) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) := by
  letI := awayAlgebra (1:R) (A := R) (B := Fin 2 → R)
  have hbijR : Function.Bijective (algebraMap R (Localization.Away (1:R))) := by
    have e : R ≃ₐ[R] Localization.Away (1:R) := IsLocalization.atUnit R (Localization.Away (1:R)) (1:R) isUnit_one
    have heq : (e : R →+* Localization.Away (1:R)) = algebraMap R (Localization.Away (1:R)) := by
      ext x; exact e.commutes x
    rw [← heq]; exact e.bijective
  have hunitB : IsUnit (algebraMap R (Fin 2 → R) (1:R)) := by rw [map_one]; exact isUnit_one
  have hbijB : Function.Bijective (algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R)))) := by
    have e : (Fin 2 → R) ≃ₐ[Fin 2 → R] Localization.Away (algebraMap R (Fin 2 → R) (1:R)) :=
      IsLocalization.atUnit (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) (algebraMap R (Fin 2 → R) (1:R)) hunitB
    have heq : (e : (Fin 2 → R) →+* Localization.Away (algebraMap R (Fin 2 → R) (1:R)))
        = algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) :=
      RingHom.ext (fun x => e.commutes x)
    rw [← heq]; exact e.bijective
  set ιR := RingEquiv.ofBijective (algebraMap R (Localization.Away (1:R))) hbijR with hιRdef
  set ιB := RingEquiv.ofBijective (algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R)))) hbijB with hιBdef
  have heqmap : Localization.awayMap (algebraMap R (Fin 2 → R)) (1:R)
      = ιB.toRingHom.comp ((algebraMap R (Fin 2 → R)).comp ιR.symm.toRingHom) := by
    apply IsLocalization.ringHom_ext (Submonoid.powers (1:R))
    ext x
    have hιRx : ιR.symm (ιR x) = x := ιR.symm_apply_apply x
    show Localization.awayMap (algebraMap R (Fin 2 → R)) (1:R) (algebraMap R (Localization.Away (1:R)) x)
      = ιB (algebraMap R (Fin 2 → R) (ιR.symm (algebraMap R (Localization.Away (1:R)) x)))
    have hιRxeq : (algebraMap R (Localization.Away (1:R)) x) = ιR x := by rw [hιRdef]; rfl
    rw [hιRxeq, hιRx]
    show Localization.awayMap (algebraMap R (Fin 2 → R)) (1:R) (algebraMap R (Localization.Away (1:R)) x)
      = algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) (algebraMap R (Fin 2 → R) x)
    unfold Localization.awayMap IsLocalization.Away.map
    rw [IsLocalization.map_eq]
  haveI : RingHomInvPair ιR.toRingHom ιR.symm.toRingHom := ⟨ιR.symm_toRingHom_comp_toRingHom, ιR.toRingHom_comp_symm_toRingHom⟩
  haveI : RingHomInvPair ιR.symm.toRingHom ιR.toRingHom := ⟨ιR.toRingHom_comp_symm_toRingHom, ιR.symm_toRingHom_comp_toRingHom⟩
  have hsmul : ∀ (r : R) (m : Fin 2 → R), ιB (r • m) = ιR r • ιB m := by
    intro r m
    show ιB ((algebraMap R (Fin 2 → R)) r * m) = (Localization.awayMap (algebraMap R (Fin 2 → R)) (1:R)) (ιR r) * ιB m
    rw [heqmap]
    show ιB ((algebraMap R (Fin 2 → R)) r * m) = ιB ((algebraMap R (Fin 2 → R)) (ιR.symm (ιR r))) * ιB m
    rw [ιR.symm_apply_apply]
    rw [← map_mul]
  set e : (Fin 2 → R) ≃ₛₗ[ιR.toRingHom] Localization.Away (algebraMap R (Fin 2 → R) (1:R)) :=
    { ιB with map_smul' := hsmul }
  have he : (algebraMap (Localization.Away (1:R)) (Localization.Away (algebraMap R (Fin 2 → R) (1:R)))).comp ιR.toRingHom
      = ιB.toRingHom.comp (algebraMap R (Fin 2 → R)) := by
    show (Localization.awayMap (algebraMap R (Fin 2 → R)) (1:R)).comp ιR.toRingHom
      = ιB.toRingHom.comp (algebraMap R (Fin 2 → R))
    rw [heqmap]
    ext x
    simp
  refine ⟨?_, Module.Finite.of_equiv_equiv ιR ιB he⟩
  exact Module.Free.of_equiv e

/-- **`A:=R`・`B:=Fin 2 → R`・`p:=1` での witness、条件 (ii)(trace が
`B` を `A` へ写す)**。`p=1` の場合は自明——`algebraMap R (Localization.
Away 1)` 自体が全単射(`IsLocalization.atUnit`)なので、**どんな元でも**
その原像が `R` に存在する。trace の値そのものを計算する必要すら無い。 -/
theorem awayOne_fin2_trace {R : Type*} [CommRing R] :
    letI := awayAlgebra (1:R) (A := R) (B := Fin 2 → R)
    ∀ b : Fin 2 → R, ∃ a : R,
      Algebra.trace (Localization.Away (1:R)) (Localization.Away (algebraMap R (Fin 2 → R) (1:R)))
          (algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) b)
        = algebraMap R (Localization.Away (1:R)) a := by
  letI := awayAlgebra (1:R) (A := R) (B := Fin 2 → R)
  intro b
  have hbijR : Function.Bijective (algebraMap R (Localization.Away (1:R))) := by
    have e : R ≃ₐ[R] Localization.Away (1:R) := IsLocalization.atUnit R (Localization.Away (1:R)) (1:R) isUnit_one
    have heq : (e : R →+* Localization.Away (1:R)) = algebraMap R (Localization.Away (1:R)) := by
      ext x; exact e.commutes x
    rw [← heq]; exact e.bijective
  obtain ⟨a, ha⟩ := hbijR.2
    (Algebra.trace (Localization.Away (1:R)) (Localization.Away (algebraMap R (Fin 2 → R) (1:R)))
      (algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) b))
  exact ⟨a, ha.symm⟩

/-- `diagonalCompare p` の純テンソルでの値。`f0 := algebraMap B Bp` を
両成分に施すだけ(`f0` は `B` の局所化への自然な写像そのもの——
`IsScalarTower.toAlgHom A B Bp` を `IsScalarTower.toAlgHom_apply` で
`algebraMap B Bp` に戻す)。条件 (iii) の witness 構成の土台。 -/
theorem diagonalCompare_tmul {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] (p : A) (b1 b2 : B) :
    letI := awayAlgebra p (A := A) (B := B)
    haveI := awayScalarTower p (A := A) (B := B)
    diagonalCompare p (b1 ⊗ₜ[A] b2)
      = (algebraMap B (Localization.Away (algebraMap A B p)) b1) ⊗ₜ[Localization.Away p]
        (algebraMap B (Localization.Away (algebraMap A B p)) b2) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  unfold diagonalCompare
  simp only [Algebra.TensorProduct.lift_tmul, AlgHom.comp_apply, Algebra.TensorProduct.includeLeft_apply,
    IsScalarTower.toAlgHom_apply, AlgHom.coe_restrictScalars', Algebra.TensorProduct.includeRight_apply]
  rw [Algebra.TensorProduct.tmul_mul_tmul]
  simp

/-- **`algebraMap B Bp` が全射ならば `diagonalCompare p` も全射**。
`Bp ⊗_{Ap} Bp` の任意の元は純テンソルの有限和(`TensorProduct.
induction_on`)なので、各純テンソル `x⊗y` を `algebraMap B Bp` の全射性
で `f0(b1)⊗f0(b2)`(`diagonalCompare_tmul` より `diagonalCompare p
(b1⊗ₜb2)` に等しい)の形に引き戻し、和で閉じる。条件 (iii) は
`p^n • elem Ap Bp` という**特定の1点**への到達可能性を要求するだけ
なので、全射性さえ言えれば `elem` の値そのものを計算する必要が無い
(`elem` が `Exists.choose` で非構成的に定義されているため、これは
本質的な簡略化——値を計算せず全射性だけで押し切る)。 -/
theorem diagonalCompare_surjective_of_algebraMap_surjective {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    (p : A) (hsurj : Function.Surjective (algebraMap B (Localization.Away (algebraMap A B p)))) :
    letI := awayAlgebra p (A := A) (B := B)
    Function.Surjective (diagonalCompare (A := A) (B := B) p) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  intro z
  refine TensorProduct.induction_on z ⟨0, map_zero _⟩ ?_ ?_
  · intro x y
    obtain ⟨b1, hb1⟩ := hsurj x
    obtain ⟨b2, hb2⟩ := hsurj y
    exact ⟨b1 ⊗ₜ[A] b2, by rw [diagonalCompare_tmul, hb1, hb2]⟩
  · rintro x y ⟨ex, hex⟩ ⟨ey, hey⟩
    exact ⟨ex + ey, by rw [map_add, hex, hey]⟩

/-- **`A:=R`・`B:=Fin 2 → R`・`p:=1` での witness、条件 (iii)
(idempotent `p^n e_{B/A}` が `B⊗_AB` の像に入る)**。`p=1` なので
`p^n=1`・`1•x=x`——`elem Ap Bp` の値そのものを問わず、
`algebraMap B Bp` の全射性(`p=1`単元なので `hbijB`)から
`diagonalCompare_surjective_of_algebraMap_surjective` で`elem Ap Bp`
自体の原像を直接引けば良い。 -/
theorem awayOne_fin2_idempotent {R : Type*} [CommRing R] :
    letI := awayAlgebra (1:R) (A := R) (B := Fin 2 → R)
    haveI := awayOne_fin2_etale (R := R)
    haveI := (awayOne_fin2_freeFinite (R := R)).2
    ∀ n : ℕ, 0 < n → ∃ e : TensorProduct R (Fin 2 → R) (Fin 2 → R), diagonalCompare (1:R) e
      = (1:R)^n • Algebra.FormallyUnramified.elem (Localization.Away (1:R))
          (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) := by
  letI := awayAlgebra (1:R) (A := R) (B := Fin 2 → R)
  haveI := awayOne_fin2_etale (R := R)
  haveI := (awayOne_fin2_freeFinite (R := R)).2
  intro n _
  have hunitB : IsUnit (algebraMap R (Fin 2 → R) (1:R)) := by rw [map_one]; exact isUnit_one
  have hbijB : Function.Bijective (algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R)))) := by
    have e : (Fin 2 → R) ≃ₐ[Fin 2 → R] Localization.Away (algebraMap R (Fin 2 → R) (1:R)) :=
      IsLocalization.atUnit (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) (algebraMap R (Fin 2 → R) (1:R)) hunitB
    have heq : (e : (Fin 2 → R) →+* Localization.Away (algebraMap R (Fin 2 → R) (1:R)))
        = algebraMap (Fin 2 → R) (Localization.Away (algebraMap R (Fin 2 → R) (1:R))) :=
      RingHom.ext (fun x => e.commutes x)
    rw [← heq]; exact e.bijective
  have hsurj := diagonalCompare_surjective_of_algebraMap_surjective (A := R) (B := Fin 2 → R) (1:R) hbijB.2
  obtain ⟨e, he⟩ := hsurj (Algebra.FormallyUnramified.elem (Localization.Away (1:R))
      (Localization.Away (algebraMap R (Fin 2 → R) (1:R))))
  exact ⟨e, by rw [he, one_pow, one_smul]⟩

/-- **`Definition 2.1`(`isAlmostEtaleCovering`)の non-vacuous witness、
完成**。`A:=R`・`B:=Fin 2 → R`・`p:=1` を選ぶと条件 (i)(ii)(iii) が
全て成立する——`awayOne_fin2_freeFinite`・`awayOne_fin2_etale`
(条件 (i))・`awayOne_fin2_trace`(条件 (ii))・`awayOne_fin2_idempotent`
(条件 (iii))を貼り合わせるだけ。 -/
theorem awayOne_fin2_isAlmostEtaleCovering {R : Type*} [CommRing R] :
    IsAlmostEtaleCovering (A := R) (B := Fin 2 → R) (1:R) := by
  unfold IsAlmostEtaleCovering
  refine ⟨(awayOne_fin2_freeFinite (R := R)).1, (awayOne_fin2_freeFinite (R := R)).2,
    awayOne_fin2_etale (R := R), ?_, ?_⟩
  · exact awayOne_fin2_trace (R := R)
  · exact awayOne_fin2_idempotent (R := R)

/-- 非空虚性の具体例(`R := ℤ`)。`Fin 2 → ℤ` が `ℤ` 上 `p:=1` で
almost étale covering になる、という具体的なインスタンスが実在する。 -/
example : IsAlmostEtaleCovering (A := ℤ) (B := Fin 2 → ℤ) (1:ℤ) :=
  awayOne_fin2_isAlmostEtaleCovering (R := ℤ)

/-! ## `p:=1` witness の一般化——`B=Fin 2 → R` に限らず、`Etale`・
`Finite`・`Free` を満たす**任意の**`B`について同じ議論が成り立つ
(2026-09-05)。特に `B := A`(恒等拡大)の場合が今まで解けなかった——
`awayAlgebra`(局所化の自己写像として人工的に配線した instance)と
mathlib の標準的な自己代数(`Algebra.id`)が衝突していたため
(`elem_self` の docstring 参照)。上の `awayOne_fin2_*` 系列で確立した
「`AlgEquiv` を経由せず `RingHom.Etale` レベルで押し切る」戦略は
`B=Fin2→R` という具体的な型に依存しない——**`A`・`B` を同じ universe
に揃えさえすれば**(`RingHom.Etale.stableUnderComposition` が
`{R S T : Type u}` と同一 universe を要求するため、`Type*` を素朴に
2回書くと別々の universe metavariable になり失敗する——`tools/
lean-idioms.md` #46 参照)、一般の `B` でそのまま動くことが分かった。 -/

universe u

/-- `awayOne_fin2_etale` の一般化。`B` が `Fin 2 → R` である必要は
無く、`Algebra.Etale A B` でありさえすれば良い。 -/
theorem awayOne_etale_of_etale {A B : Type u} [CommRing A] [CommRing B] [Algebra A B] [Algebra.Etale A B] :
    letI := awayAlgebra (1:A) (A := A) (B := B)
    Algebra.Etale (Localization.Away (1:A)) (Localization.Away (algebraMap A B (1:A))) := by
  have hbijA : Function.Bijective (algebraMap A (Localization.Away (1:A))) := by
    have e : A ≃ₐ[A] Localization.Away (1:A) := IsLocalization.atUnit A (Localization.Away (1:A)) (1:A) isUnit_one
    have heq : (e : A →+* Localization.Away (1:A)) = algebraMap A (Localization.Away (1:A)) := by
      ext x; exact e.commutes x
    rw [← heq]; exact e.bijective
  have hunitB : IsUnit (algebraMap A B (1:A)) := by rw [map_one]; exact isUnit_one
  have hbijB : Function.Bijective (algebraMap B (Localization.Away (algebraMap A B (1:A)))) := by
    have e : B ≃ₐ[B] Localization.Away (algebraMap A B (1:A)) :=
      IsLocalization.atUnit B (Localization.Away (algebraMap A B (1:A))) (algebraMap A B (1:A)) hunitB
    have heq : (e : B →+* Localization.Away (algebraMap A B (1:A)))
        = algebraMap B (Localization.Away (algebraMap A B (1:A))) :=
      RingHom.ext (fun x => e.commutes x)
    rw [← heq]; exact e.bijective
  set ιA := RingEquiv.ofBijective (algebraMap A (Localization.Away (1:A))) hbijA with hιAdef
  set ιB := RingEquiv.ofBijective (algebraMap B (Localization.Away (algebraMap A B (1:A)))) hbijB with hιBdef
  have heqmap : Localization.awayMap (algebraMap A B) (1:A)
      = ιB.toRingHom.comp ((algebraMap A B).comp ιA.symm.toRingHom) := by
    apply IsLocalization.ringHom_ext (Submonoid.powers (1:A))
    ext x
    have hιAx : ιA.symm (ιA x) = x := ιA.symm_apply_apply x
    show Localization.awayMap (algebraMap A B) (1:A) (algebraMap A (Localization.Away (1:A)) x)
      = ιB (algebraMap A B (ιA.symm (algebraMap A (Localization.Away (1:A)) x)))
    have hιAxeq : (algebraMap A (Localization.Away (1:A)) x) = ιA x := by rw [hιAdef]; rfl
    rw [hιAxeq, hιAx]
    show Localization.awayMap (algebraMap A B) (1:A) (algebraMap A (Localization.Away (1:A)) x)
      = algebraMap B (Localization.Away (algebraMap A B (1:A))) (algebraMap A B x)
    unfold Localization.awayMap IsLocalization.Away.map
    rw [IsLocalization.map_eq]
  have hgEtale : RingHom.Etale (algebraMap A B) := RingHom.etale_algebraMap.mpr ‹Algebra.Etale A B›
  have hιAinvEtale : RingHom.Etale (ιA.symm.toRingHom) :=
    RingHom.Etale.of_bijective ιA.symm.bijective
  have hιBEtale : RingHom.Etale (ιB.toRingHom) :=
    RingHom.Etale.of_bijective ιB.bijective
  have hcompEtale : RingHom.Etale ((algebraMap A B).comp ιA.symm.toRingHom) :=
    RingHom.Etale.stableUnderComposition ιA.symm.toRingHom (algebraMap A B) hιAinvEtale hgEtale
  have hfullEtale : RingHom.Etale (ιB.toRingHom.comp ((algebraMap A B).comp ιA.symm.toRingHom)) :=
    RingHom.Etale.stableUnderComposition ((algebraMap A B).comp ιA.symm.toRingHom) ιB.toRingHom hcompEtale hιBEtale
  rw [← heqmap] at hfullEtale
  exact RingHom.Etale.toAlgebra hfullEtale

/-- `awayOne_fin2_freeFinite` の一般化。`Algebra.Etale` は module の
有限性・自由性を含意しない(例: 開はめ込みは étale だが加群として
有限とは限らない)ため、`[Module.Finite A B] [Module.Free A B]` を
別途仮定として要求する(Definition 2.1 原文の「projective A[1/p]-module
of finite rank」に対応、`Module.Free` への逸脱は冒頭の記録通り)。 -/
theorem awayOne_freeFinite_of_etale {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Finite A B] [Module.Free A B] :
    letI := awayAlgebra (1:A) (A := A) (B := B)
    Module.Free (Localization.Away (1:A)) (Localization.Away (algebraMap A B (1:A))) ∧
    Module.Finite (Localization.Away (1:A)) (Localization.Away (algebraMap A B (1:A))) := by
  letI := awayAlgebra (1:A) (A := A) (B := B)
  have hbijA : Function.Bijective (algebraMap A (Localization.Away (1:A))) := by
    have e : A ≃ₐ[A] Localization.Away (1:A) := IsLocalization.atUnit A (Localization.Away (1:A)) (1:A) isUnit_one
    have heq : (e : A →+* Localization.Away (1:A)) = algebraMap A (Localization.Away (1:A)) := by
      ext x; exact e.commutes x
    rw [← heq]; exact e.bijective
  have hunitB : IsUnit (algebraMap A B (1:A)) := by rw [map_one]; exact isUnit_one
  have hbijB : Function.Bijective (algebraMap B (Localization.Away (algebraMap A B (1:A)))) := by
    have e : B ≃ₐ[B] Localization.Away (algebraMap A B (1:A)) :=
      IsLocalization.atUnit B (Localization.Away (algebraMap A B (1:A))) (algebraMap A B (1:A)) hunitB
    have heq : (e : B →+* Localization.Away (algebraMap A B (1:A)))
        = algebraMap B (Localization.Away (algebraMap A B (1:A))) :=
      RingHom.ext (fun x => e.commutes x)
    rw [← heq]; exact e.bijective
  set ιA := RingEquiv.ofBijective (algebraMap A (Localization.Away (1:A))) hbijA with hιAdef
  set ιB := RingEquiv.ofBijective (algebraMap B (Localization.Away (algebraMap A B (1:A)))) hbijB with hιBdef
  have heqmap : Localization.awayMap (algebraMap A B) (1:A)
      = ιB.toRingHom.comp ((algebraMap A B).comp ιA.symm.toRingHom) := by
    apply IsLocalization.ringHom_ext (Submonoid.powers (1:A))
    ext x
    have hιAx : ιA.symm (ιA x) = x := ιA.symm_apply_apply x
    show Localization.awayMap (algebraMap A B) (1:A) (algebraMap A (Localization.Away (1:A)) x)
      = ιB (algebraMap A B (ιA.symm (algebraMap A (Localization.Away (1:A)) x)))
    have hιAxeq : (algebraMap A (Localization.Away (1:A)) x) = ιA x := by rw [hιAdef]; rfl
    rw [hιAxeq, hιAx]
    show Localization.awayMap (algebraMap A B) (1:A) (algebraMap A (Localization.Away (1:A)) x)
      = algebraMap B (Localization.Away (algebraMap A B (1:A))) (algebraMap A B x)
    unfold Localization.awayMap IsLocalization.Away.map
    rw [IsLocalization.map_eq]
  haveI : RingHomInvPair ιA.toRingHom ιA.symm.toRingHom := ⟨ιA.symm_toRingHom_comp_toRingHom, ιA.toRingHom_comp_symm_toRingHom⟩
  haveI : RingHomInvPair ιA.symm.toRingHom ιA.toRingHom := ⟨ιA.toRingHom_comp_symm_toRingHom, ιA.symm_toRingHom_comp_toRingHom⟩
  have hsmul : ∀ (r : A) (m : B), ιB (r • m) = ιA r • ιB m := by
    intro r m
    rw [Algebra.smul_def, Algebra.smul_def]
    show ιB ((algebraMap A B) r * m) = (Localization.awayMap (algebraMap A B) (1:A)) (ιA r) * ιB m
    rw [heqmap]
    show ιB ((algebraMap A B) r * m) = ιB ((algebraMap A B) (ιA.symm (ιA r))) * ιB m
    rw [ιA.symm_apply_apply]
    rw [← map_mul]
  set e : B ≃ₛₗ[ιA.toRingHom] Localization.Away (algebraMap A B (1:A)) :=
    { ιB with map_smul' := hsmul }
  have he : (algebraMap (Localization.Away (1:A)) (Localization.Away (algebraMap A B (1:A)))).comp ιA.toRingHom
      = ιB.toRingHom.comp (algebraMap A B) := by
    show (Localization.awayMap (algebraMap A B) (1:A)).comp ιA.toRingHom
      = ιB.toRingHom.comp (algebraMap A B)
    rw [heqmap]
    ext x
    simp
  refine ⟨?_, Module.Finite.of_equiv_equiv ιA ιB he⟩
  exact Module.Free.of_equiv e

/-- `awayOne_fin2_trace` の一般化。`B` に関する仮定は不要——`p=1` の
全単射性だけから従う。 -/
theorem awayOne_trace_of_unit {A B : Type u} [CommRing A] [CommRing B] [Algebra A B] :
    letI := awayAlgebra (1:A) (A := A) (B := B)
    ∀ b : B, ∃ a : A,
      Algebra.trace (Localization.Away (1:A)) (Localization.Away (algebraMap A B (1:A)))
          (algebraMap B (Localization.Away (algebraMap A B (1:A))) b)
        = algebraMap A (Localization.Away (1:A)) a := by
  letI := awayAlgebra (1:A) (A := A) (B := B)
  intro b
  have hbijA : Function.Bijective (algebraMap A (Localization.Away (1:A))) := by
    have e : A ≃ₐ[A] Localization.Away (1:A) := IsLocalization.atUnit A (Localization.Away (1:A)) (1:A) isUnit_one
    have heq : (e : A →+* Localization.Away (1:A)) = algebraMap A (Localization.Away (1:A)) := by
      ext x; exact e.commutes x
    rw [← heq]; exact e.bijective
  obtain ⟨a, ha⟩ := hbijA.2
    (Algebra.trace (Localization.Away (1:A)) (Localization.Away (algebraMap A B (1:A)))
      (algebraMap B (Localization.Away (algebraMap A B (1:A))) b))
  exact ⟨a, ha.symm⟩

/-- `awayOne_fin2_idempotent` の一般化。 -/
theorem awayOne_idempotent_of_etale {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Algebra.Etale A B] [Module.Finite A B] [Module.Free A B] :
    letI := awayAlgebra (1:A) (A := A) (B := B)
    haveI := awayOne_etale_of_etale (A := A) (B := B)
    haveI := (awayOne_freeFinite_of_etale (A := A) (B := B)).2
    ∀ n : ℕ, 0 < n → ∃ e : TensorProduct A B B, diagonalCompare (1:A) e
      = (1:A)^n • Algebra.FormallyUnramified.elem (Localization.Away (1:A))
          (Localization.Away (algebraMap A B (1:A))) := by
  letI := awayAlgebra (1:A) (A := A) (B := B)
  haveI := awayOne_etale_of_etale (A := A) (B := B)
  haveI := (awayOne_freeFinite_of_etale (A := A) (B := B)).2
  intro n _
  have hunitB : IsUnit (algebraMap A B (1:A)) := by rw [map_one]; exact isUnit_one
  have hbijB : Function.Bijective (algebraMap B (Localization.Away (algebraMap A B (1:A)))) := by
    have e : B ≃ₐ[B] Localization.Away (algebraMap A B (1:A)) :=
      IsLocalization.atUnit B (Localization.Away (algebraMap A B (1:A))) (algebraMap A B (1:A)) hunitB
    have heq : (e : B →+* Localization.Away (algebraMap A B (1:A)))
        = algebraMap B (Localization.Away (algebraMap A B (1:A))) :=
      RingHom.ext (fun x => e.commutes x)
    rw [← heq]; exact e.bijective
  have hsurj := diagonalCompare_surjective_of_algebraMap_surjective (A := A) (B := B) (1:A) hbijB.2
  obtain ⟨e, he⟩ := hsurj (Algebra.FormallyUnramified.elem (Localization.Away (1:A))
      (Localization.Away (algebraMap A B (1:A))))
  exact ⟨e, by rw [he, one_pow, one_smul]⟩

/-- **`Definition 2.1` の non-vacuous witness、`B := Fin 2 → R` から
「`Etale`・`Finite`・`Free` な任意の `B`」への一般化、完成**。 -/
theorem awayOne_isAlmostEtaleCovering_of_etale {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Algebra.Etale A B] [Module.Finite A B] [Module.Free A B] :
    IsAlmostEtaleCovering (A := A) (B := B) (1:A) := by
  unfold IsAlmostEtaleCovering
  refine ⟨(awayOne_freeFinite_of_etale (A := A) (B := B)).1, (awayOne_freeFinite_of_etale (A := A) (B := B)).2,
    awayOne_etale_of_etale (A := A) (B := B), ?_, ?_⟩
  · exact awayOne_trace_of_unit (A := A) (B := B)
  · exact awayOne_idempotent_of_etale (A := A) (B := B)

/-- **`elem_self` の docstring・冒頭の記録が「次のセッションへ持ち越す」
としていた `A=B`(恒等拡大)の具体例、ついに解消**。`Algebra.Etale A A`・
`Module.Finite A A`・`Module.Free A A` は全て自明なインスタンス
(`Algebra.Etale.id` 等)なので、上の一般化定理をそのまま適用するだけで
良い——`awayAlgebra`(人工的な局所化の自己写像 instance)と mathlib の
標準的な自己代数(`Algebra.id`)の衝突は、`RingHom.Etale` レベルで
押し切る戦略には最初から関係が無かった、という事後的な確認。 -/
example {A : Type u} [CommRing A] : IsAlmostEtaleCovering (A := A) (B := A) (1:A) :=
  awayOne_isAlmostEtaleCovering_of_etale (A := A) (B := A)

/-! ## `p` を単元に限らない一般化(2026-09-05、大幅な前進)。

上の `awayOne_*` 系列は全て `p:=1`(単元)に依存していた——`Bp≅B`・
`Ap≅A`(局所化が実質恒等)という**退化した**特殊ケースで、Faltings の
理論が本来扱いたい「`p` が真の(単元でない)素元」の場合には触れて
いなかった。ここでは、**`p` が単元かどうかに一切依存しない**形で、
`B` が(古典的な意味で)`A` 上 étale・finite・free でありさえすれば
`IsAlmostEtaleCovering A B p` が**任意の `p`** について成り立つことを
示す——「非分岐拡大は(古典的な意味で既に)almost étale である」という
Faltings の理論の健全性チェックであり、かつ非退化な non-vacuous witness
そのもの。

鍵となった発見は2つ:
1. **`Algebra.FormallyUnramified.elem` の一意性**(`elem_unique_of_props`)
   ——`Exists.choose` で非構成的に定義された idempotent だが、
   その定義性質(`one_tmul_sub_tmul_one_mul_elem`・`lmul_elem`)を
   満たす元は実は**一意**であることを示せる(`S⊗_RS`の元`t`が
   `∀s,(1⊗s-s⊗1)*t=0`を満たすなら、任意の`x`について`x*t=(μx⊗1)*t`
   という「吸収」性質を持つ——これだけから`t*t'=t=t'`が出る)。
2. **`elem` の局所化に関する自然性**(`diagonalCompare_elem_eq`)
   ——一意性を武器に、`diagonalCompare p (elem A B) = elem Ap Bp`
   を**任意の `p`**(単元である必要が無い!)について示せる。
   条件(iii)の「`s'∈Bp` 全体」への拡張は、`f0(B)` が(1)を満たす
   ことから出発し、`Z:={s'|(1⊗s'-s'⊗1)*t=0}` が**部分環**になる
   こと(`Zclosed_add`・`Zclosed_mul`)、`π:=f0(p)`(単元、`Away`
   局所化の定義から)の逆元も`Z`に入ること(`Zclosed_inv`、
   `(π⊗π)`で割ってから可逆性でキャンセルする論法)、`Bp`の任意の元が
   `f0(a)*π⁻ⁿ`の形(`IsLocalization.Away.surj`)であることを組み合わせて
   全射的に拡張する。 -/

/-- `t` が `elem` の定義性質(1)を満たすなら、任意の `x` について
`x*t` は `μ(x)⊗1` の形の元と `t` の積に等しい——`elem` 一意性の鍵。 -/
theorem elem_absorb_of_prop {R S : Type u} [CommRing R] [CommRing S] [Algebra R S]
    (t : TensorProduct R S S) (ht1 : ∀ s : S, (1 ⊗ₜ[R] s - s ⊗ₜ[R] 1) * t = 0) :
    ∀ x : TensorProduct R S S, x * t = ((Algebra.TensorProduct.lmul' R) x ⊗ₜ[R] (1:S)) * t := by
  have hswap : ∀ s : S, (1 ⊗ₜ[R] s : TensorProduct R S S) * t = (s ⊗ₜ[R] 1) * t := by
    intro s
    have h := ht1 s
    rw [sub_mul] at h
    exact sub_eq_zero.mp h
  intro x
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul a b =>
    have step1 : (a ⊗ₜ[R] b : TensorProduct R S S) = (a ⊗ₜ[R] (1:S)) * ((1:S) ⊗ₜ[R] b) := by
      rw [Algebra.TensorProduct.tmul_mul_tmul, mul_one, one_mul]
    rw [step1, mul_assoc, hswap b, ← mul_assoc, Algebra.TensorProduct.tmul_mul_tmul, mul_one]
    simp [Algebra.TensorProduct.lmul'_apply_tmul]
  | add x y hx hy =>
    rw [add_mul, hx, hy, map_add, TensorProduct.add_tmul, add_mul]

/-- **`Algebra.FormallyUnramified.elem` の一意性**(mathlib に無かった
事実)。定義性質(annihilate `1⊗s-s⊗1`・augment to `1`)を満たす元は
一意——`elem_absorb_of_prop` を `1-t'`(定義性質(2)から `μ(1-t')=0`)
に適用すると `(1-t')*t=0` すなわち `t=t'*t`、対称に `t'=t*t'`、
可換性 `t*t'=t'*t` と合わせて `t=t'`。 -/
theorem elem_unique_of_props {R S : Type u} [CommRing R] [CommRing S] [Algebra R S]
    (t t' : TensorProduct R S S)
    (ht1 : ∀ s : S, (1 ⊗ₜ[R] s - s ⊗ₜ[R] 1) * t = 0) (ht2 : (Algebra.TensorProduct.lmul' R) t = 1)
    (ht1' : ∀ s : S, (1 ⊗ₜ[R] s - s ⊗ₜ[R] 1) * t' = 0) (ht2' : (Algebra.TensorProduct.lmul' R) t' = 1) :
    t = t' := by
  have h1 : (1 - t') * t = 0 := by
    rw [elem_absorb_of_prop t ht1 (1 - t')]
    rw [map_sub, map_one, ht2']
    simp
  have h2 : (1 - t) * t' = 0 := by
    rw [elem_absorb_of_prop t' ht1' (1 - t)]
    rw [map_sub, map_one, ht2]
    simp
  have e1 : t = t' * t := by
    have := h1
    rw [sub_mul, one_mul, sub_eq_zero] at this
    exact this
  have e2 : t' = t * t' := by
    have := h2
    rw [sub_mul, one_mul, sub_eq_zero] at this
    exact this
  rw [e1, mul_comm]
  exact e2.symm

/-- `Z := {s | (1⊗s-s⊗1)*t=0}` は加法について閉じる。 -/
theorem Zclosed_add {Ap Bp : Type u} [CommRing Ap] [CommRing Bp] [Algebra Ap Bp] (t'' : TensorProduct Ap Bp Bp)
    (s1 s2 : Bp)
    (h1 : ((1:Bp) ⊗ₜ[Ap] s1 - s1 ⊗ₜ[Ap] (1:Bp)) * t'' = 0)
    (h2 : ((1:Bp) ⊗ₜ[Ap] s2 - s2 ⊗ₜ[Ap] (1:Bp)) * t'' = 0) :
    ((1:Bp) ⊗ₜ[Ap] (s1+s2) - (s1+s2) ⊗ₜ[Ap] (1:Bp)) * t'' = 0 := by
  have expand : (1:Bp) ⊗ₜ[Ap] (s1+s2) - (s1+s2) ⊗ₜ[Ap] (1:Bp)
      = ((1:Bp) ⊗ₜ[Ap] s1 - s1 ⊗ₜ[Ap] (1:Bp)) + ((1:Bp) ⊗ₜ[Ap] s2 - s2 ⊗ₜ[Ap] (1:Bp)) := by
    rw [TensorProduct.tmul_add, TensorProduct.add_tmul]; ring
  rw [expand, add_mul, h1, h2, add_zero]

/-- `Z` は乗法についても閉じる(`1⊗(s1s2)-(s1s2)⊗1` を
`(1⊗s1)*[(1⊗s2)-(s2⊗1)] + [(1⊗s1)-(s1⊗1)]*(s2⊗1)` に分解し、
可換性で `t` を中に押し込んで両方消す)。 -/
theorem Zclosed_mul {Ap Bp : Type u} [CommRing Ap] [CommRing Bp] [Algebra Ap Bp] (t'' : TensorProduct Ap Bp Bp)
    (s1 s2 : Bp)
    (h1 : ((1:Bp) ⊗ₜ[Ap] s1 - s1 ⊗ₜ[Ap] (1:Bp)) * t'' = 0)
    (h2 : ((1:Bp) ⊗ₜ[Ap] s2 - s2 ⊗ₜ[Ap] (1:Bp)) * t'' = 0) :
    ((1:Bp) ⊗ₜ[Ap] (s1*s2) - (s1*s2) ⊗ₜ[Ap] (1:Bp)) * t'' = 0 := by
  have expand : (1:Bp) ⊗ₜ[Ap] (s1*s2) - (s1*s2) ⊗ₜ[Ap] (1:Bp)
      = ((1:Bp) ⊗ₜ[Ap] s1) * ((1:Bp) ⊗ₜ[Ap] s2 - s2 ⊗ₜ[Ap] (1:Bp))
        + ((1:Bp) ⊗ₜ[Ap] s1 - s1 ⊗ₜ[Ap] (1:Bp)) * (s2 ⊗ₜ[Ap] (1:Bp)) := by
    simp only [mul_sub, sub_mul, Algebra.TensorProduct.tmul_mul_tmul, mul_one, one_mul]
    abel
  rw [expand, add_mul]
  have e1 : (1:Bp) ⊗ₜ[Ap] s1 * (((1:Bp) ⊗ₜ[Ap] s2 - s2 ⊗ₜ[Ap] (1:Bp)) * t'') = 0 := by rw [h2, mul_zero]
  have e2 : (((1:Bp) ⊗ₜ[Ap] s1 - s1 ⊗ₜ[Ap] (1:Bp)) * t'') * (s2 ⊗ₜ[Ap] (1:Bp)) = 0 := by rw [h1, zero_mul]
  rw [← mul_assoc] at e1
  rw [mul_assoc] at e2
  rw [e1, zero_add]
  rw [show ((1:Bp) ⊗ₜ[Ap] s1 - s1 ⊗ₜ[Ap] (1:Bp)) * (s2 ⊗ₜ[Ap] (1:Bp)) * t''
        = ((1:Bp) ⊗ₜ[Ap] s1 - s1 ⊗ₜ[Ap] (1:Bp)) * (t'' * (s2 ⊗ₜ[Ap] (1:Bp))) from by ring]
  rw [e2]

/-- `π` が単元で `Z` に入っているなら、`π⁻¹` も `Z` に入る
(`(π⊗π)` で割ってから可逆性でキャンセル)。`Away` 局所化の
`π:=f0(p)` は定義から単元なので、これで `Z` の生成元が揃う。 -/
theorem Zclosed_inv {Ap Bp : Type u} [CommRing Ap] [CommRing Bp] [Algebra Ap Bp] (t'' : TensorProduct Ap Bp Bp)
    (π : Bp) (hπunit : IsUnit π)
    (hπZ : ((1:Bp) ⊗ₜ[Ap] π - π ⊗ₜ[Ap] (1:Bp)) * t'' = 0) :
    ((1:Bp) ⊗ₜ[Ap] (↑hπunit.unit⁻¹ : Bp) - (↑hπunit.unit⁻¹ : Bp) ⊗ₜ[Ap] (1:Bp)) * t'' = 0 := by
  set πinv : Bp := (↑hπunit.unit⁻¹ : Bp) with hπinvdef
  have hmulinv : π * πinv = 1 := by
    rw [hπinvdef]
    exact_mod_cast hπunit.unit.mul_inv
  have step : ((π:Bp) ⊗ₜ[Ap] π) * ((1:Bp) ⊗ₜ[Ap] πinv - πinv ⊗ₜ[Ap] (1:Bp))
      = -(((1:Bp) ⊗ₜ[Ap] π - π ⊗ₜ[Ap] (1:Bp))) := by
    simp only [mul_sub, Algebra.TensorProduct.tmul_mul_tmul, mul_one, hmulinv]
    ring
  have key : ((π:Bp) ⊗ₜ[Ap] π) * (((1:Bp) ⊗ₜ[Ap] πinv - πinv ⊗ₜ[Ap] (1:Bp)) * t'') = 0 := by
    rw [← mul_assoc, step, neg_mul, hπZ, neg_zero]
  have hunit2 : IsUnit ((π:Bp) ⊗ₜ[Ap] π : TensorProduct Ap Bp Bp) := by
    refine IsUnit.of_mul_eq_one ((πinv:Bp) ⊗ₜ[Ap] πinv) ?_
    rw [Algebra.TensorProduct.tmul_mul_tmul, hmulinv]; rfl
  rw [← mul_zero ((π:Bp) ⊗ₜ[Ap] π : TensorProduct Ap Bp Bp)] at key
  exact hunit2.mul_left_cancel key

/-- `Algebra.TensorProduct.lmul'` と `diagonalCompare` の可換性
(条件(iii)の性質(2)を移送するのに使う)。 -/
theorem lmul'_diagonalCompare {A B : Type u} [CommRing A] [CommRing B] [Algebra A B] (p : A)
    (x : TensorProduct A B B) :
    letI := awayAlgebra p (A := A) (B := B)
    haveI := awayScalarTower p (A := A) (B := B)
    Algebra.TensorProduct.lmul' (Localization.Away p) (diagonalCompare p x)
      = algebraMap B (Localization.Away (algebraMap A B p)) (Algebra.TensorProduct.lmul' A x) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul b1 b2 =>
    rw [diagonalCompare_tmul, Algebra.TensorProduct.lmul'_apply_tmul, Algebra.TensorProduct.lmul'_apply_tmul, map_mul]
  | add x y hx hy => rw [map_add, map_add, hx, hy, map_add, map_add]

/-- **`elem` の局所化に関する自然性、任意の `p`(単元である必要は無い)
について**。`diagonalCompare p` は idempotent `elem A B` を idempotent
`elem Ap Bp` へ**厳密に**(`p^n` を掛けることなく)送る——`p=1` の
witness で使った「全射性だけで押し切る」トリックとは異なり、こちらは
`elem` の**一意性**(`elem_unique_of_props`)を武器に、`diagonalCompare
p (elem A B)` が `Ap`・`Bp` に対する `elem` の定義性質を実際に満たす
ことを示して結論する。性質(1)(`∀s'∈Bp`)は `f0(B)` から出発し `Z` が
部分環であること・`π⁻¹∈Z`(`Zclosed_inv`)・`IsLocalization.Away.surj`
(`Bp`の任意の元は`f0(a)*π⁻ⁿ`の形)を組み合わせて拡張する。 -/
theorem diagonalCompare_elem_eq {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Algebra.FormallyUnramified A B] [Algebra.EssFiniteType A B] (p : A) :
    letI := awayAlgebra p (A := A) (B := B)
    haveI := awayScalarTower p (A := A) (B := B)
    [Algebra.FormallyUnramified (Localization.Away p) (Localization.Away (algebraMap A B p))] →
    [Algebra.EssFiniteType (Localization.Away p) (Localization.Away (algebraMap A B p))] →
    diagonalCompare p (Algebra.FormallyUnramified.elem A B)
      = Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p)) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  intro _ _
  set Bp := Localization.Away (algebraMap A B p)
  set Ap := Localization.Away p
  set t'' := diagonalCompare p (Algebra.FormallyUnramified.elem A B) with ht''def
  set π : Bp := algebraMap B Bp (algebraMap A B p) with hπdef
  have hZB : ∀ b : B, ((1:Bp) ⊗ₜ[Ap] (algebraMap B Bp b) - (algebraMap B Bp b) ⊗ₜ[Ap] (1:Bp)) * t'' = 0 := by
    intro b
    have h1 : (1 ⊗ₜ[A] b - b ⊗ₜ[A] 1) * Algebra.FormallyUnramified.elem A B = 0 :=
      Algebra.FormallyUnramified.one_tmul_sub_tmul_one_mul_elem b
    have h2 := congrArg (diagonalCompare p) h1
    rw [map_zero, map_mul, map_sub, diagonalCompare_tmul, diagonalCompare_tmul, map_one] at h2
    rw [ht''def]; exact h2
  have hπunit : IsUnit π := IsLocalization.Away.algebraMap_isUnit (algebraMap A B p)
  have hπZ : ((1:Bp) ⊗ₜ[Ap] π - π ⊗ₜ[Ap] (1:Bp)) * t'' = 0 := hZB (algebraMap A B p)
  have hZpow : ∀ n : ℕ, ((1:Bp) ⊗ₜ[Ap] (π^n) - (π^n) ⊗ₜ[Ap] (1:Bp)) * t'' = 0 := by
    intro n
    induction n with
    | zero => simp
    | succ k ih => rw [pow_succ]; exact Zclosed_mul t'' (π^k) π ih hπZ
  have hprop1 : ∀ s' : Bp, ((1:Bp) ⊗ₜ[Ap] s' - s' ⊗ₜ[Ap] (1:Bp)) * t'' = 0 := by
    intro s'
    obtain ⟨n, a, ha⟩ := IsLocalization.Away.surj (algebraMap A B p) s'
    have hun : IsUnit (π ^ n) := hπunit.pow n
    have hZinv : ((1:Bp) ⊗ₜ[Ap] (↑hun.unit⁻¹ : Bp) - (↑hun.unit⁻¹ : Bp) ⊗ₜ[Ap] (1:Bp)) * t'' = 0 :=
      Zclosed_inv t'' (π^n) hun (hZpow n)
    have hZa : ((1:Bp) ⊗ₜ[Ap] (algebraMap B Bp a) - (algebraMap B Bp a) ⊗ₜ[Ap] (1:Bp)) * t'' = 0 :=
      hZB a
    have hs'eq : s' = (algebraMap B Bp a) * (↑hun.unit⁻¹ : Bp) := by
      have hmulinv : (π^n) * (↑hun.unit⁻¹ : Bp) = 1 := by exact_mod_cast hun.unit.mul_inv
      calc s' = s' * (π^n) * (↑hun.unit⁻¹ : Bp) := by rw [mul_assoc, hmulinv, mul_one]
        _ = (algebraMap B Bp a) * (↑hun.unit⁻¹ : Bp) := by rw [ha]
    rw [hs'eq]
    exact Zclosed_mul t'' (algebraMap B Bp a) (↑hun.unit⁻¹ : Bp) hZa hZinv
  have hprop2 : Algebra.TensorProduct.lmul' Ap t'' = 1 := by
    rw [ht''def, lmul'_diagonalCompare, Algebra.FormallyUnramified.lmul_elem, map_one]
  exact elem_unique_of_props t'' (Algebra.FormallyUnramified.elem Ap Bp)
    hprop1 hprop2
    (Algebra.FormallyUnramified.one_tmul_sub_tmul_one_mul_elem) (Algebra.FormallyUnramified.lmul_elem)

/-- **`Definition 2.1`(`IsAlmostEtaleCovering`)の non-vacuous witness、
`p` を単元に限らず一般化して完成**。`B` が `A` 上(古典的な意味で)
étale・finite・free でありさえすれば、**任意の `p : A`**(単元である
必要は無い、真の非単元素元でも良い)について成り立つ。条件(i)は
`RingHom.Etale.propertyIsLocal.localizationAwayPreserves`(étale性の
`Away`局所化保存性、mathlibの「局所的な性質」framework)・`Module.
Finite.of_isLocalization`・`Module.free_of_isLocalizedModule`という、
p=1 witness の `IsLocalization.atUnit`/`RingEquiv` 経由の議論より
遥かに直接的な道具で閉じる。条件(ii)は `Algebra.trace_localization`
(mathlib既存、trace の局所化との可換性)一発。条件(iii)は`p^n •`を
`diagonalCompare`の`A`-線形性で外に出し、`diagonalCompare_elem_eq`
(`n`に依らない、`elem`の完全な自然性)を適用するだけ。 -/
theorem isAlmostEtaleCovering_of_etale_general {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Algebra.Etale A B] [Module.Finite A B] [Module.Free A B] (p : A) :
    IsAlmostEtaleCovering (A := A) (B := B) p := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  unfold IsAlmostEtaleCovering
  have hFree : Module.Free (Localization.Away p) (Localization.Away (algebraMap A B p)) := by
    have hM : Algebra.algebraMapSubmonoid B (Submonoid.powers p) = Submonoid.powers (algebraMap A B p) := by
      rw [Algebra.algebraMapSubmonoid, Submonoid.map_powers]
    have hIsLoc : IsLocalization (Algebra.algebraMapSubmonoid B (Submonoid.powers p)) (Localization.Away (algebraMap A B p)) := by
      rw [hM]; infer_instance
    have hIsLocMod : IsLocalizedModule (Submonoid.powers p)
        (IsScalarTower.toAlgHom A B (Localization.Away (algebraMap A B p))).toLinearMap :=
      instIsLocalizedModuleToLinearMapToAlgHomOfIsLocalizationAlgebraMapSubmonoid (Submonoid.powers p)
    exact Module.free_of_isLocalizedModule (Submonoid.powers p)
      (IsScalarTower.toAlgHom A B (Localization.Away (algebraMap A B p))).toLinearMap
  have hFinite : Module.Finite (Localization.Away p) (Localization.Away (algebraMap A B p)) := by
    have hM : Algebra.algebraMapSubmonoid B (Submonoid.powers p) = Submonoid.powers (algebraMap A B p) := by
      rw [Algebra.algebraMapSubmonoid, Submonoid.map_powers]
    have : IsLocalization (Algebra.algebraMapSubmonoid B (Submonoid.powers p)) (Localization.Away (algebraMap A B p)) := by
      rw [hM]; infer_instance
    exact Module.Finite.of_isLocalization A B (Submonoid.powers p)
  have hEtale : Algebra.Etale (Localization.Away p) (Localization.Away (algebraMap A B p)) := by
    have hf : RingHom.Etale (algebraMap A B) := RingHom.etale_algebraMap.mpr ‹Algebra.Etale A B›
    have h2 : RingHom.Etale (IsLocalization.Away.map (Localization.Away p) (Localization.Away (algebraMap A B p)) (algebraMap A B) p) :=
      RingHom.Etale.propertyIsLocal.localizationAwayPreserves (algebraMap A B) p (Localization.Away p) (Localization.Away (algebraMap A B p)) hf
    have heq : IsLocalization.Away.map (Localization.Away p) (Localization.Away (algebraMap A B p)) (algebraMap A B) p
        = Localization.awayMap (algebraMap A B) p := rfl
    rw [heq] at h2
    exact RingHom.Etale.toAlgebra h2
  refine ⟨hFree, hFinite, hEtale, ?_, ?_⟩
  · intro b
    exact ⟨Algebra.trace A B b, Algebra.trace_localization A (Submonoid.powers p) b⟩
  · intro n _
    haveI := hEtale
    haveI := hFinite
    refine ⟨p^n • Algebra.FormallyUnramified.elem A B, ?_⟩
    rw [map_smul, diagonalCompare_elem_eq]

/-- 非空虚性の具体例——**真の非単元素元 `p:=5` で**、`B := Fin 2 → ℤ` が
`ℤ` 上 almost étale covering になる。`5` はもちろん `ℤ` の単元ではない
(`awayOne_*` 系列の `p:=1` 退化ケースとは質的に異なる)。 -/
example : IsAlmostEtaleCovering (A := ℤ) (B := Fin 2 → ℤ) (5:ℤ) :=
  isAlmostEtaleCovering_of_etale_general (A := ℤ) (B := Fin 2 → ℤ) 5

/-! ## `Theorem 2.4` への足場(2026-09-05、部分的な前進)。

`Theorem 2.4`(物理p.7-8=印字p.260-261、260dpi精読)は (i)(ii) で難度が
質的に違う——(i)(`Ω_A⊗B→Ω_B` が almost isomorphism)は`Theorem 2.2`
と同じ「Hochschild cohomology を使う」の一文で済まされ未着手のまま
だが、(ii)(有限群`G`の semilinear 作用に関する`H^i(G,M)`の almost
消滅)の証明は「`m` が `A/tr_{B/A}(B)` を零化する」ことに帰着し、これは
Definition 2.1 直後の remark(iii)(物理p.6末尾)そのもの:「`p^ε・
e_{B/A}=Σxᵢ⊗yᵢ`(条件(iii)の witness)なら**任意の`b∈B`について
`p^ε・b=Σtr_{B/A}(b・xᵢ)・yᵢ`**」という trace 恒等式一発で閉じる。

この恒等式を、`elem`の一意性の議論(`elem_absorb_of_prop`)と同じ
「annihilation性質だけから出発する」流儀で**部分的に**証明した——
`Ψ_s(t):=Σtr(s・xᵢ)yᵢ`(`t=Σxᵢ⊗yᵢ`の`Tr1map`版)が`s・Tr1map(t)`に
等しいことまでは annihilation 性質だけから閉じる(`Ψ_eq_smul_Tr1`)。
残るのは**`Tr1map(elem R S) = 1`という1個のスカラー等式**のみ——
これは「trace form が非退化」という分離拡大論の古典的事実の系だが、
mathlib の `traceForm_nondegenerate` 系は体の拡大専用(`Algebra.
IsSeparable K L`)で、一般の有限自由な非分岐環拡大向けの版が無い
ため、次回への未解決課題として残す。 -/

/-- `Tr1map R S t := Σᵢtr(xᵢ)•yᵢ` for `t=Σxᵢ⊗yᵢ∈S⊗_RS`(`R`-線形に
延長)。`Theorem 2.4(ii)`の remark(iii)恒等式の土台。 -/
noncomputable def Tr1map (R S : Type u) [CommRing R] [CommRing S] [Algebra R S] :
    TensorProduct R S S →ₗ[R] S :=
  TensorProduct.lift
    { toFun := fun x => (Algebra.trace R S x) • LinearMap.id (M := S)
      map_add' := by intro x y; simp [add_smul]
      map_smul' := by intro r x; simp [mul_smul] }

theorem Tr1map_tmul {R S : Type u} [CommRing R] [CommRing S] [Algebra R S] (x y : S) :
    Tr1map R S (x ⊗ₜ[R] y) = (Algebra.trace R S x) • y := by
  simp [Tr1map]

/-- `t`の annihilation 性質(`elem`の定義性質(1))から、`(1⊗s)*t`への
`Tr1map`は`s`をそのまま外に出せる(`Tr1map`は第2因子への`s`倍と
可換)。 -/
theorem Tr1_swap_smul {R S : Type u} [CommRing R] [CommRing S] [Algebra R S]
    (t : TensorProduct R S S) (s : S) :
    Tr1map R S ((1:S) ⊗ₜ[R] s * t) = s • Tr1map R S t := by
  induction t using TensorProduct.induction_on with
  | zero => simp
  | tmul x y =>
    rw [Algebra.TensorProduct.tmul_mul_tmul, one_mul, Tr1map_tmul, Tr1map_tmul]
    simp [Algebra.smul_def]; ring
  | add x y hx hy => rw [mul_add, map_add, map_add, hx, hy, smul_add]

/-- **`Theorem 2.4(ii)`の remark(iii)恒等式、核心部分**。`t`が
annihilation 性質(`elem`の定義性質(1))を満たすなら、`Ψ_s(t):=
Tr1map((s⊗1)*t)`(`t=Σxᵢ⊗yᵢ`なら`Σtr(s・xᵢ)yᵢ`に等しい)は
`s・Tr1map(t)`に等しい——annihilation の swap 性質(`(1⊗s)*t=(s⊗1)*t`)
と`Tr1_swap_smul`を貼り合わせるだけ。`t:=p^ε・elem A B`・`Tr1map
(elem A B)=1`(未証明、上記コメント参照)が分かれば、`s・p^ε =
Σtr(s・xᵢ)yᵢ`という remark(iii)の恒等式そのものが従う。 -/
theorem Ψ_eq_smul_Tr1 {R S : Type u} [CommRing R] [CommRing S] [Algebra R S]
    (t : TensorProduct R S S) (ht1 : ∀ s : S, (1 ⊗ₜ[R] s - s ⊗ₜ[R] 1) * t = 0) (s : S) :
    Tr1map R S ((s ⊗ₜ[R] (1:S)) * t) = s • Tr1map R S t := by
  have hs := ht1 s
  rw [sub_mul, sub_eq_zero] at hs
  rw [← hs, Tr1_swap_smul]

/-! ## `Tr1map(elem R S) = 1`、ついに証明(2026-09-05、続)。

上のコメントで「未解決課題」としていたスカラー等式を、基底を用いた
添字計算で完全に証明した。`S`(`R`上有限自由)の基底`{bᵢ}`を選び、
`elem R S`をその基底で`Σᵢbᵢ⊗aᵢ`と分解する(`tensorDecompEquiv`)と、
annihilation性質から構造定数を経由した関係式
`bⱼ・aₖ=Σᵢ(b.repr(bⱼ・bᵢ)k)・aᵢ`が出る(`hextract`)。これと
augmentation(`Σᵢbᵢaᵢ=1`)・trace の行列表示(`tr(bⱼ)=Σᵢb.repr(bⱼbᵢ)i`、
`Algebra.trace_eq_matrix_trace`)を、`S`の可換性(構造定数の対称性
`b.repr(bⱼbᵢ)k=b.repr(bᵢbⱼ)k`)で束ねる(`Finset.sum_comm`)と、
`Tr1map(elem R S)`と`1(=augmentation)`が**同じ二重和の並べ替え**である
ことが分かる。これは分離拡大論の古典的事実(trace form の非退化性の
系)だが、mathlibには一般の環拡大向けの版が無かった。 -/

/-- 基底`b`を経由した`S⊗_RS ≃ₗ[R] (ι→₀S)`(第1因子をbの座標で展開)。 -/
noncomputable def tensorDecompEquiv {R S : Type u} [CommRing R] [CommRing S] [Algebra R S] {ι : Type*}
    [Fintype ι] [DecidableEq ι] (b : Module.Basis ι R S) :
    TensorProduct R S S ≃ₗ[R] (ι →₀ S) :=
  (TensorProduct.congr b.repr (LinearEquiv.refl R S)).trans (TensorProduct.finsuppScalarLeft R S ι)

theorem tensorDecompEquiv_tmul {R S : Type u} [CommRing R] [CommRing S] [Algebra R S] {ι : Type*}
    [Fintype ι] [DecidableEq ι] (b : Module.Basis ι R S) (i : ι) (z : S) :
    tensorDecompEquiv b (b i ⊗ₜ[R] z) = Finsupp.single i z := by
  show TensorProduct.finsuppScalarLeft R S ι (TensorProduct.congr b.repr (LinearEquiv.refl R S) (b i ⊗ₜ[R] z))
    = Finsupp.single i z
  rw [show TensorProduct.congr b.repr (LinearEquiv.refl R S) (b i ⊗ₜ[R] z) = Finsupp.single i 1 ⊗ₜ[R] z from by
    simp [TensorProduct.congr, Module.Basis.repr_self]]
  ext j
  rw [TensorProduct.finsuppScalarLeft_apply]
  simp [Finsupp.lapply_apply, Finsupp.single_apply]

theorem tensorDecompEquiv_eq {R S : Type u} [CommRing R] [CommRing S] [Algebra R S] {ι : Type*}
    [Fintype ι] [DecidableEq ι] (b : Module.Basis ι R S) (t : TensorProduct R S S) :
    t = Finset.univ.sum (fun i => b i ⊗ₜ[R] (tensorDecompEquiv b t i)) := by
  apply (tensorDecompEquiv b).injective
  rw [map_sum]
  simp only [tensorDecompEquiv_tmul]
  ext j
  simp

/-- annihilation性質から出る、基底の構造定数を経由した関係式
`bⱼ・aₖ=Σᵢ(b.repr(bⱼ・bᵢ)k)・aᵢ`(`t=Σᵢbᵢ⊗aᵢ`の座標`aᵢ:=tensorDecompEquiv
b t i`)。`Tr1map(elem R S)=1`の証明の核心部品。 -/
theorem tensorDecompEquiv_extract {R S : Type u} [CommRing R] [CommRing S] [Algebra R S] {ι : Type*}
    [Fintype ι] [DecidableEq ι] (b : Module.Basis ι R S) (t : TensorProduct R S S)
    (ht1 : ∀ s : S, (1 ⊗ₜ[R] s - s ⊗ₜ[R] 1) * t = 0) (j k : ι) :
    b j * (tensorDecompEquiv b t k) =
      Finset.univ.sum (fun i => (b.repr (b j * b i) k) • (tensorDecompEquiv b t i)) := by
  set a := tensorDecompEquiv b t with hadef
  have hdecomp : t = Finset.univ.sum (fun i => b i ⊗ₜ[R] a i) := tensorDecompEquiv_eq b t
  have hswap : (1 ⊗ₜ[R] (b j) : TensorProduct R S S) * t = (b j ⊗ₜ[R] 1) * t := by
    have h := ht1 (b j)
    rw [sub_mul, sub_eq_zero] at h
    exact h
  have hLHS : tensorDecompEquiv b ((1 ⊗ₜ[R] (b j) : TensorProduct R S S) * t) k = b j * a k := by
    rw [hdecomp, Finset.mul_sum]
    rw [show (Finset.univ.sum fun i => (1:S) ⊗ₜ[R] (b j) * (b i ⊗ₜ[R] a i))
        = Finset.univ.sum (fun i => b i ⊗ₜ[R] (b j * a i)) from by
      apply Finset.sum_congr rfl; intro i _
      rw [Algebra.TensorProduct.tmul_mul_tmul, one_mul]]
    rw [map_sum]
    simp only [tensorDecompEquiv_tmul]
    simp [Finsupp.finsetSum_apply, Finsupp.single_apply]
  have hRHS : tensorDecompEquiv b ((b j ⊗ₜ[R] (1:S)) * t) k =
      Finset.univ.sum (fun i => (b.repr (b j * b i) k) • a i) := by
    rw [hdecomp, Finset.mul_sum]
    rw [show (Finset.univ.sum fun i => (b j : S) ⊗ₜ[R] (1:S) * (b i ⊗ₜ[R] a i))
        = Finset.univ.sum (fun i => (b j * b i) ⊗ₜ[R] a i) from by
      apply Finset.sum_congr rfl; intro i _
      rw [Algebra.TensorProduct.tmul_mul_tmul, one_mul]]
    rw [show (Finset.univ.sum fun i => (b j * b i) ⊗ₜ[R] a i)
        = Finset.univ.sum (fun i => Finset.univ.sum (fun m => (b.repr (b j * b i) m) • (b m ⊗ₜ[R] a i))) from by
      apply Finset.sum_congr rfl; intro i _
      conv_lhs => rw [← b.sum_repr (b j * b i)]
      rw [TensorProduct.sum_tmul]
      simp [TensorProduct.smul_tmul]]
    rw [map_sum, Finsupp.finsetSum_apply]
    apply Finset.sum_congr rfl
    intro i _
    rw [map_sum, Finsupp.finsetSum_apply]
    simp only [map_smul, tensorDecompEquiv_tmul, Finsupp.smul_apply, Finsupp.single_apply]
    simp
  rw [← hLHS, hswap, hRHS]

/-- **`Tr1map(elem R S)=1`、Faltings remark(iii)の核心のスカラー等式**
(mathlibに一般形が無かった、分離拡大論の古典的事実)。基底
`{bᵢ}`での`elem R S=Σᵢbᵢ⊗aᵢ`分解を使い、augmentation(`Σbᵢaᵢ=1`)と
`tensorDecompEquiv_extract`(annihilationから来る構造定数関係式)を
`Finset.sum_comm`+可換性(構造定数の対称性)で束ねると、`Tr1map(elem)`
と`1`が**同じ二重和**であることが分かる。 -/
theorem Tr1map_elem_eq_one {R S : Type u} [CommRing R] [CommRing S] [Algebra R S] {ι : Type*}
    [Fintype ι] [DecidableEq ι] (b : Module.Basis ι R S) (t : TensorProduct R S S)
    (ht1 : ∀ s : S, (1 ⊗ₜ[R] s - s ⊗ₜ[R] 1) * t = 0) (ht2 : Algebra.TensorProduct.lmul' R t = 1) :
    Tr1map R S t = 1 := by
  set a := tensorDecompEquiv b t with hadef
  have hdecomp : t = Finset.univ.sum (fun i => b i ⊗ₜ[R] a i) := tensorDecompEquiv_eq b t
  have hAug : Finset.univ.sum (fun i => b i * a i) = 1 := by
    rw [← ht2, hdecomp, map_sum]
    apply Finset.sum_congr rfl
    intro i _
    rw [Algebra.TensorProduct.lmul'_apply_tmul]
  have hTr1 : Tr1map R S t = Finset.univ.sum (fun j => (Algebra.trace R S (b j)) • a j) := by
    rw [hdecomp, map_sum]
    apply Finset.sum_congr rfl
    intro j _
    exact Tr1map_tmul (b j) (a j)
  have htrace : ∀ j, Algebra.trace R S (b j) =
      Finset.univ.sum (fun i => b.repr (b j * b i) i) := by
    intro j
    rw [Algebra.trace_eq_matrix_trace b]
    unfold Matrix.trace Matrix.diag
    apply Finset.sum_congr rfl
    intro i _
    rw [Algebra.leftMulMatrix_eq_repr_mul]
  have hSUM1 : Finset.univ.sum (fun j : ι => Finset.univ.sum
      (fun i => (b.repr (b j * b i) j) • a i)) = 1 := by
    rw [← hAug]
    apply Finset.sum_congr rfl
    intro j _
    rw [← tensorDecompEquiv_extract b t ht1 j j]
  have hTr1' : Tr1map R S t = Finset.univ.sum (fun j : ι => Finset.univ.sum
      (fun i => (b.repr (b j * b i) i) • a j)) := by
    rw [hTr1]
    apply Finset.sum_congr rfl
    intro j _
    rw [htrace j, Finset.sum_smul]
  rw [hTr1', ← hSUM1]
  rw [Finset.sum_comm (s := (Finset.univ : Finset ι)) (t := (Finset.univ : Finset ι))]
  apply Finset.sum_congr rfl
  intro i _
  apply Finset.sum_congr rfl
  intro j _
  rw [show b i * b j = b j * b i from mul_comm _ _]

/-- `Tr1map`と`diagonalCompare`の可換性(`lmul'_diagonalCompare`と同型の
議論、`Algebra.trace_localization`経由)。remark(iii)を局所化した先
(`Ap`・`Bp`、honestに étale)から`A`・`B`側へ引き戻すのに使う。 -/
theorem Tr1map_diagonalCompare {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Free A B] [Module.Finite A B] (p : A)
    (x : TensorProduct A B B) :
    letI := awayAlgebra p (A := A) (B := B)
    haveI := awayScalarTower p (A := A) (B := B)
    Tr1map (Localization.Away p) (Localization.Away (algebraMap A B p)) (diagonalCompare p x)
      = algebraMap B (Localization.Away (algebraMap A B p)) (Tr1map A B x) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul b1 b2 =>
    rw [diagonalCompare_tmul, Tr1map_tmul, Tr1map_tmul,
      Algebra.trace_localization A (Submonoid.powers p) b1, Algebra.smul_def, Algebra.smul_def,
      ← IsScalarTower.algebraMap_apply A (Localization.Away p) (Localization.Away (algebraMap A B p)),
      IsScalarTower.algebraMap_apply A B (Localization.Away (algebraMap A B p)), map_mul]
  | add x y hx hy => rw [map_add, map_add, hx, hy, map_add, map_add]

/-- **`Theorem 2.4(ii)`の remark(iii)恒等式、完全に証明**(2026-09-05)。
`IsAlmostEtaleCovering`の条件(iii)の witness `e`(`diagonalCompare p e
= p^n・elem Ap Bp`)から、Faltings が主張する trace 恒等式
`p^n・b = Σtr_{B/A}(b・xᵢ)yᵢ`(`Tr1map A B ((b⊗1)*e)`と同じ)を
導く——`Tr1map_elem_eq_one`(局所化先`Ap,Bp`で)・`Ψ_eq_smul_Tr1`・
`Tr1map_diagonalCompare`(局所化の可換性で`A,B`側へ引き戻す)を
組み合わせる。`algebraMap B Bp`の単射性(Faltings 自身も「全ての環で
`p`による乗法は単射」と明記している標準仮定)を経由する。 -/
theorem remark_iii_trace_identity {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Algebra.Etale A B] [Module.Finite A B] [Module.Free A B] (p : A) :
    letI := awayAlgebra p (A := A) (B := B)
    haveI := isAlmostEtaleCovering_of_etale_general (A := A) (B := B) p |>.2.2.1
    haveI := (isAlmostEtaleCovering_of_etale_general (A := A) (B := B) p |>.2.1 : Module.Finite _ _)
    Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))) →
    ∀ (n : ℕ) (e : TensorProduct A B B) (b : B),
    (diagonalCompare p e
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p))) →
    algebraMap A B (p ^ n) * b = Tr1map A B ((b ⊗ₜ[A] (1:B)) * e) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  haveI hEtale := isAlmostEtaleCovering_of_etale_general (A := A) (B := B) p |>.2.2.1
  haveI hFinite := (isAlmostEtaleCovering_of_etale_general (A := A) (B := B) p |>.2.1 : Module.Finite _ _)
  intro hf0inj n e b he
  set Ap := Localization.Away p
  set Bp := Localization.Away (algebraMap A B p)
  set f0 := algebraMap B Bp with hf0def
  haveI hFreeAp : Module.Free Ap Bp := isAlmostEtaleCovering_of_etale_general (A := A) (B := B) p |>.1
  set bB := Module.Free.chooseBasis Ap Bp
  have key1 : Tr1map Ap Bp ((f0 b ⊗ₜ[Ap] (1:Bp)) * (diagonalCompare p e)) = f0 (algebraMap A B (p^n) * b) := by
    rw [he]
    have hscale : (f0 b ⊗ₜ[Ap] (1:Bp)) * (p^n • Algebra.FormallyUnramified.elem Ap Bp)
        = p^n • ((f0 b ⊗ₜ[Ap] (1:Bp)) * Algebra.FormallyUnramified.elem Ap Bp) := by
      rw [mul_smul_comm]
    rw [hscale, LinearMap.map_smul_of_tower,
      Ψ_eq_smul_Tr1 _ (Algebra.FormallyUnramified.one_tmul_sub_tmul_one_mul_elem) (f0 b),
      Tr1map_elem_eq_one bB _ (Algebra.FormallyUnramified.one_tmul_sub_tmul_one_mul_elem) (Algebra.FormallyUnramified.lmul_elem)]
    rw [smul_eq_mul, mul_one]
    show p^n • f0 b = f0 (algebraMap A B (p^n) * b)
    rw [Algebra.smul_def, map_mul, ← IsScalarTower.algebraMap_apply A B Bp,
      IsScalarTower.algebraMap_apply A Ap Bp]
  have key2 : Tr1map Ap Bp ((f0 b ⊗ₜ[Ap] (1:Bp)) * (diagonalCompare p e))
      = f0 (Tr1map A B ((b ⊗ₜ[A] (1:B)) * e)) := by
    rw [show (f0 b ⊗ₜ[Ap] (1:Bp)) * (diagonalCompare p e)
        = diagonalCompare p ((b ⊗ₜ[A] (1:B)) * e) from by
      rw [map_mul, diagonalCompare_tmul, map_one]]
    exact Tr1map_diagonalCompare p _
  have := key1.symm.trans key2
  exact hf0inj this

/-! ## `Theorem 2.4(ii)` 第2文(`M^G/tr_G(M)`)、`M:=B` の場合を Ideal norm
で完成(2026-09-05、続々々々)。

`Theorem 2.4(ii)` の証明本文は「`m` が `A/tr_{B/A}(B)` を零化することを
示せば十分」で始まり、`remark_iii_trace_identity`(既に完成)の後、
「`N_{B/A}`(ノルム)を両辺に適用して結論を導く」と続く——この最後の
「ノルムを適用する」部分を、mathlibの`Ideal.relNorm`(`Mathlib.
RingTheory.Ideal.Norm.RelNorm`、`IsDedekindDomain`+`IsIntegrallyClosed`
な環同士の相対イデアルノルム)で実行した。 -/

/-- `Tr1map A B x` は常に `tr_{B/A}(B)` が生成する `A`-イデアルの
`B` への像の中にある(`Tr1map`の定義`x⊗y↦tr(x)•y`から直接)。 -/
theorem Tr1map_mem_traceIdeal_map {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (x : TensorProduct A B B) :
    Tr1map A B x ∈ (Ideal.span (Set.range (Algebra.trace A B)) : Ideal A).map (algebraMap A B) := by
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul a b =>
    rw [Tr1map_tmul, Algebra.smul_def]
    apply Ideal.mul_mem_right
    apply Ideal.mem_map_of_mem
    exact Ideal.subset_span ⟨a, rfl⟩
  | add x y hx hy => rw [map_add]; exact Ideal.add_mem _ hx hy

attribute [local instance] FractionRing.liftAlgebra

/-- **`Theorem 2.4(ii)` の証明が要求する「`m` は `A/tr_{B/A}(B)` を
零化する」の`Ideal norm`版**。`remark_iii_trace_identity`(`b:=1`)で
`p^n∈tr_{B/A}(B)·B`(`B`側のイデアル関係)を得た後、`Ideal.relNorm`
(`Ideal.relNorm_algebraMap`——`relNorm(I.map(algebraMap)) = I^finrank`)
を両辺に適用して`A`側へ引き戻し、`Ideal.pow_le_self`(高次のべきは
元のイデアルに含まれる)で閉じる。`n`(条件(iii)で任意に選べる)を
`finrank`倍した指数で、`p^{n・finrank}`が`tr_{B/A}(B)`に入ることを
示す——`n`を大きく取れば取るほど深い annihilation が言えるので、
Faltings の「`m` annihilates」の精神(`ε`を任意に小さく取れる)を
`finrank`個おきの指数列として実現する。 -/
theorem trace_ideal_pow_mem_traceIdeal {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Algebra.Etale A B] [Module.Finite A B] [Module.Free A B] (p : A)
    [IsDedekindDomain A] [IsDedekindDomain B] [Module.IsTorsionFree A B] :
    letI := awayAlgebra p (A := A) (B := B)
    haveI := isAlmostEtaleCovering_of_etale_general (A := A) (B := B) p |>.2.2.1
    haveI := (isAlmostEtaleCovering_of_etale_general (A := A) (B := B) p |>.2.1 : Module.Finite _ _)
    Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))) →
    ∀ (n : ℕ) (e : TensorProduct A B B),
    (diagonalCompare p e
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p))) →
    (Ideal.span {p ^ (n * Module.finrank (FractionRing A) (FractionRing B))} : Ideal A)
      ≤ Ideal.span (Set.range (Algebra.trace A B)) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  haveI hEtale := isAlmostEtaleCovering_of_etale_general (A := A) (B := B) p |>.2.2.1
  haveI hFinite := (isAlmostEtaleCovering_of_etale_general (A := A) (B := B) p |>.2.1 : Module.Finite _ _)
  intro hf0inj n e he
  set a := (Ideal.span (Set.range (Algebra.trace A B)) : Ideal A)
  have heq : (1:B) ⊗ₜ[A] (1:B) * e = e := by
    rw [show (1:B) ⊗ₜ[A] (1:B) = (1 : TensorProduct A B B) from rfl, one_mul]
  have hstep1 : algebraMap A B (p^n) ∈ a.map (algebraMap A B) := by
    have hthis := remark_iii_trace_identity (A := A) (B := B) p hf0inj n e 1 he
    rw [mul_one, heq] at hthis
    rw [hthis]
    exact Tr1map_mem_traceIdeal_map e
  have hstep2 : (Ideal.span {algebraMap A B (p^n)} : Ideal B) ≤ a.map (algebraMap A B) := by
    rw [Ideal.span_singleton_le_iff_mem]
    exact hstep1
  have hstep3 : (Ideal.span {algebraMap A B (p^n)} : Ideal B) = (Ideal.span {p^n} : Ideal A).map (algebraMap A B) := by
    rw [Ideal.map_span]
    congr 1
    simp
  rw [hstep3] at hstep2
  have hstep4 := Ideal.relNorm_mono A hstep2
  rw [Ideal.relNorm_algebraMap B (Ideal.span {p^n} : Ideal A)] at hstep4
  rw [Ideal.relNorm_algebraMap B a] at hstep4
  calc (Ideal.span {p ^ (n * Module.finrank (FractionRing A) (FractionRing B))} : Ideal A)
      = (Ideal.span {p^n} : Ideal A) ^ Module.finrank (FractionRing A) (FractionRing B) := by
        rw [Ideal.span_singleton_pow, pow_mul]
    _ ≤ a ^ Module.finrank (FractionRing A) (FractionRing B) := hstep4
    _ ≤ a := Ideal.pow_le_self (Module.finrank_pos.ne')

/-! ## remark 2.1(v)、Hochschild cohomology の消滅——**訂正**(2026-09-05、
続々々々々)。

これまで複数回、「remark 2.1(v)(`Theorem 2.2`-`2.4`が使う`m`が
Hochschild cohomology を零化するという事実)はFaltings自身が本文中で
証明せず外部参照に頼っている」と報告してきたが、**これは誤りだった**。
原文(物理p.6=印字p.259、260dpi精読)を再確認したところ、remark (v) は
「`e_{B/A}=Σxᵢ⊗yᵢ`が`B⊗AB`の元だったなら、`b₀⊗b₁⊗⋯⊗b_{n+1}↦Σxᵢ⊗
yᵢb₀⊗b₁⊗⋯⊗b_{n+1}`という null-homotopy が得られ、Hochschild
cohomology は消滅する」という**完全な、標準的な**構成を与えている——
分離代数論の古典的事実(バー分解の縮約ホモトピー)そのもので、
Faltings は省略していない。

この事実を、mathlib の `CategoryTheory.Abelian.Ext`(導来圏経由の
一般Ext理論、`Mathlib.Algebra.Homology.DerivedCategory.Ext.*`)を使って
**証明した**——`elem R S`(annihilation性質)から`S`が`S⊗_RS`-加群として
`μ:S⊗_RS→S`の**切断**を持つ(`hochSection`)ことを示し、これは
「`S`が`S⊗_RS`-加群として射影的である」ことを意味する(`hochModule_
projective`、`Module.Projective.of_split`経由)。射影加群からの`Ext`は
正の次数で消える(`Ext.eq_zero_of_projective`、標準的な一般論)ので、
`HH^n(S/R,M) := Ext^n_{S⊗_RS}(S,M)`(Theorem 15、`HHⁿ(A,M)≅Extⁿ_{Aᵉ}
(A,M)`という標準的な同値な定式化を直接採用)が`n>0`で消える。

★これは remark (v) の**honest な場合**(`S`が`R`上honestに formally
unramified、`p`が単元でなくても良い、`elem`の annihilation 性質が
exactに成り立つ場合)——`Theorem 2.2`-`2.4`が要求する**almost**な場合
(`B`が単に almost étale、`p^ε elem`のみ`B⊗AB`にある場合)への一般化
は、`remark_iii_trace_identity`と同型の「局所化を経由してinjectivity
で降ろす」議論を要し、次回への課題として残す——ただし今回の honest な
場合の完成により、remark (v) 自体が Faltings の言う通り**成立する
完全な定理**であることが実証され、以前の「未証明」という評価は
正式に撤回する。 -/

open CategoryTheory in
/-- annihilation 性質から、`elem R S` を使って`μ:S⊗_RS→S`の切断
`s(b):=(b⊗1)*elem`を構成する。`S⊗_RS`-線形性の検証が核心(swap
annihilation性質`one_tmul_mul_elem`を使う)。 -/
noncomputable def hochSection (R S : Type u) [CommRing R] [CommRing S] [Algebra R S]
    [Algebra.FormallyUnramified R S] [Algebra.EssFiniteType R S] :
    letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
    S →ₗ[TensorProduct R S S] TensorProduct R S S := by
  letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
  refine
  { toFun := fun b => (b ⊗ₜ[R] (1:S)) * Algebra.FormallyUnramified.elem R S
    map_add' := by intro x y; simp [TensorProduct.add_tmul, add_mul]
    map_smul' := ?_ }
  intro z b
  simp only [RingHom.id_apply, Algebra.smul_def]
  induction z using TensorProduct.induction_on with
  | zero => simp
  | tmul p q =>
    show ((Algebra.TensorProduct.lmul' R (p ⊗ₜ[R] q) * b) ⊗ₜ[R] (1:S)) * Algebra.FormallyUnramified.elem R S
      = (p ⊗ₜ[R] q) * ((b ⊗ₜ[R] (1:S)) * Algebra.FormallyUnramified.elem R S)
    rw [Algebra.TensorProduct.lmul'_apply_tmul]
    have hstep1 : ((p*q*b : S) ⊗ₜ[R] (1:S) : TensorProduct R S S)
        = ((p*b : S) ⊗ₜ[R] (1:S)) * (q ⊗ₜ[R] (1:S)) := by
      rw [Algebra.TensorProduct.tmul_mul_tmul]; ring_nf
    rw [hstep1, mul_assoc]
    have hstep2 : (q ⊗ₜ[R] (1:S) : TensorProduct R S S) * Algebra.FormallyUnramified.elem R S
        = (1 ⊗ₜ[R] q) * Algebra.FormallyUnramified.elem R S :=
      (Algebra.FormallyUnramified.one_tmul_mul_elem q).symm
    rw [hstep2, ← mul_assoc, Algebra.TensorProduct.tmul_mul_tmul, mul_one, one_mul]
    rw [← mul_assoc, Algebra.TensorProduct.tmul_mul_tmul, mul_one]
  | add x y hx hy =>
    simp only [map_add, add_mul, TensorProduct.add_tmul]
    rw [hx, hy]

theorem hochSection_augment (R S : Type u) [CommRing R] [CommRing S] [Algebra R S]
    [Algebra.FormallyUnramified R S] [Algebra.EssFiniteType R S] (b : S) :
    Algebra.TensorProduct.lmul' R (hochSection R S b) = b := by
  show Algebra.TensorProduct.lmul' R ((b ⊗ₜ[R] (1:S)) * Algebra.FormallyUnramified.elem R S) = b
  rw [map_mul, Algebra.TensorProduct.lmul'_apply_tmul, mul_one, Algebra.FormallyUnramified.lmul_elem, mul_one]

open CategoryTheory in
/-- `μ:S⊗_RS→S`自身を`S⊗_RS`-線形写像として。`hochSection`と合わせて
`S`が`S⊗_RS`-加群として`μ`の切断を持つことを示す。 -/
noncomputable def hochSectionLmul (R S : Type u) [CommRing R] [CommRing S] [Algebra R S] :
    letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
    TensorProduct R S S →ₗ[TensorProduct R S S] S := by
  letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
  exact { toFun := Algebra.TensorProduct.lmul' R,
          map_add' := map_add _,
          map_smul' := fun z x => by
            show Algebra.TensorProduct.lmul' R (z*x) = algebraMap (TensorProduct R S S) S z * Algebra.TensorProduct.lmul' R x
            rw [map_mul]
            rfl }

open CategoryTheory in
/-- **`S`は`S⊗_RS`-加群として射影的**(`hochSection`が`μ`の切断を
与えることから、`Module.Projective.of_split`経由)。remark(v)の
「`S`が`S⊗AS`のbimodule direct summand」という主張の honest な場合。 -/
theorem hochModule_projective (R S : Type u) [CommRing R] [CommRing S] [Algebra R S]
    [Algebra.FormallyUnramified R S] [Algebra.EssFiniteType R S] :
    letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
    Module.Projective (TensorProduct R S S) S := by
  letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
  apply Module.Projective.of_split (hochSection R S) (hochSectionLmul R S)
  ext b
  simp only [LinearMap.comp_apply, LinearMap.id_apply]
  show Algebra.TensorProduct.lmul' R (hochSection R S b) = b
  exact hochSection_augment R S b

open CategoryTheory in
theorem hochModule_catProjective (R S : Type u) [CommRing R] [CommRing S] [Algebra R S]
    [Algebra.FormallyUnramified R S] [Algebra.EssFiniteType R S] :
    letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
    Projective (ModuleCat.of (TensorProduct R S S) S) := by
  letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
  haveI := hochModule_projective R S
  exact ModuleCat.projective_of_categoryTheory_projective (ModuleCat.of (TensorProduct R S S) S)

set_option maxHeartbeats 1000000 in
open CategoryTheory in
/-- **remark 2.1(v)、honest な場合、完全に証明**: `S`が`R`上formally
unramified・ess finite typeなら(honestly「formally étale」、`p`が単元
である必要は無い)、`HH^n(S/R,M) := Ext^n_{S⊗_RS}(S,M)`は`n>0`で
消える——`S`が射影的な`S⊗_RS`-加群であることから、射影対象からの
`Ext`は正の次数で消えるという一般論(`Ext.eq_zero_of_projective`)
一発で閉じる。Faltings の remark(v)前半(「`e_{B/A}`が`B⊗AB`の元
だったならHochschild cohomologyは消滅する」)の完全な形式化。 -/
theorem hochschild_ext_eq_zero (R S : Type u) [CommRing R] [CommRing S] [Algebra R S]
    [Algebra.FormallyUnramified R S] [Algebra.EssFiniteType R S]
    (M : Type u) [AddCommGroup M] :
    letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
    haveI := CategoryTheory.HasExt.standard (ModuleCat (TensorProduct R S S))
    ∀ [Module (TensorProduct R S S) M] (n : ℕ) (e : CategoryTheory.Abelian.Ext
        (ModuleCat.of (TensorProduct R S S) S) (ModuleCat.of (TensorProduct R S S) M) (n+1)), e = 0 := by
  letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
  haveI := CategoryTheory.HasExt.standard (ModuleCat (TensorProduct R S S))
  intro _ n e
  haveI := hochModule_catProjective R S
  exact CategoryTheory.Abelian.Ext.eq_zero_of_projective e

/-! ## remark 2.1(v) の almost 化への足場(2026-09-05、新規)。

`hochschild_ext_eq_zero`(honest な場合)は`B`が`FormallyUnramified`
であること(`elem A B`自身がその annihilation性質を満たすこと)に
依っていた。`Theorem 2.2`-`2.4(i)`が実際に必要とするのは
`IsAlmostEtaleCovering`(`B/A`自体は honestly unramified である必要
は無く、局所化`Bp/Ap`だけが古典的に étale であればよい)の下での
almost 版——「`p^ε e_{B/A}`に対応する witness `w∈B⊗_AB`は、annihilation
性質を`B⊗_AB`の中で厳密に満たす」という事実である。これは`w`自体は
`diagonalCompare p w = p^n•elem(Ap,Bp)`としてしか特徴づけられておらず、
honest な annihilation の議論がそのままでは`B⊗_AB`に降りてこない
(`diagonalCompare`は環準同型であって同型ではないため)。

この障壁は**`diagonalCompare`の単射性**(`Module.Free A B`と
`algebraMap B Bp`の単射性さえあれば従う——`Module.Flat`の
`rTensor_preserves_injective_linearMap`/`lTensor_preserves_injective_linearMap`
を2回使い、`IsLocalization.moduleTensorEquiv`で「両成分が既に局所化
されたテンソル積からさらに局所化を消す」同型に帰着させる)で解消
できる。これが`diagonalCompare_injective`——`Bp⊗_ApBp`での厳密な
annihilation を`diagonalCompare`で引き戻し、単射性で`B⊗_AB`へ
そのまま降ろす。 -/

/-- **`diagonalCompare p`の単射性**。`B`が`A`上free、`algebraMap B Bp`が
単射(Faltings 自身が終始置いている「p-torsion-free」という標準仮定、
このファイルの他の場所(`remark_iii_trace_identity`)でも同じ形で
使っている)という2条件だけから従う。`diagonalCompare`を
`φ := (rTensor Bp f₀)∘(lTensor B f₀)`(`f₀:B→Bp`をテンソル積の両成分に
適用、`Module.Flat`の保存性で単射)と`e:=IsLocalization.moduleTensorEquiv`
(「両成分が既に局所化された`Ap⊗ApBpBp`から、さらなる局所化`Ap`を
消して`A⊗Bp⊗Bp`に戻す」という mathlib の同型)の合成`e.symm∘φ`として
分解することで示す。 -/
theorem diagonalCompare_injective {A B : Type u} [CommRing A] [CommRing B] [Algebra A B] [Module.Free A B]
    (p : A) (hf0inj : Function.Injective (algebraMap B (Localization.Away (algebraMap A B p)))) :
    letI := awayAlgebra p (A := A) (B := B)
    haveI := awayScalarTower p (A := A) (B := B)
    Function.Injective (diagonalCompare (A := A) (B := B) p) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  set Bp := Localization.Away (algebraMap A B p)
  set Ap := Localization.Away p
  set f0 : B →ₗ[A] Bp := (Algebra.linearMap B Bp).restrictScalars A with hf0def
  have hf0injL : Function.Injective f0 := hf0inj
  haveI hflatB : Module.Flat A B := inferInstance
  haveI hflatBp : Module.Flat A Bp := by
    haveI : Module.Flat B Bp := IsLocalization.flat Bp (Submonoid.powers (algebraMap A B p))
    exact Module.Flat.trans A B Bp
  have hlTinj : Function.Injective (LinearMap.lTensor B f0) :=
    Module.Flat.lTensor_preserves_injective_linearMap f0 hf0injL
  have hrTinj : Function.Injective (LinearMap.rTensor Bp f0) :=
    Module.Flat.rTensor_preserves_injective_linearMap f0 hf0injL
  set φ : TensorProduct A B B →ₗ[A] TensorProduct A Bp Bp :=
    LinearMap.rTensor Bp f0 ∘ₗ LinearMap.lTensor B f0
  have hφinj : Function.Injective φ := hrTinj.comp hlTinj
  set e := IsLocalization.moduleTensorEquiv (Submonoid.powers p) Ap Bp Bp
  have hcomp : ∀ z : TensorProduct A B B, diagonalCompare (A := A) (B := B) p z = e.symm (φ z) := by
    intro z
    induction z using TensorProduct.induction_on with
    | zero => simp
    | tmul b1 b2 =>
      rw [diagonalCompare_tmul]
      show f0 b1 ⊗ₜ[Ap] f0 b2 = e.symm (φ (b1 ⊗ₜ[A] b2))
      have hφtmul : φ (b1 ⊗ₜ[A] b2) = f0 b1 ⊗ₜ[A] f0 b2 := rfl
      rw [hφtmul]
      symm
      rw [LinearEquiv.symm_apply_eq]
      show e (f0 b1 ⊗ₜ[Ap] f0 b2) = f0 b1 ⊗ₜ[A] f0 b2
      rfl
    | add x y hx hy => rw [map_add, map_add, hx, hy, map_add]
  intro x y hxy
  have := hcomp x ▸ hcomp y ▸ hxy
  have h2 : φ x = φ y := e.symm.injective this
  exact hφinj h2

/-- **almost witness `w`の`B⊗_AB`での厳密な swap-annihilation**。
`IsAlmostEtaleCovering`条件(iii)の witness `w`(`diagonalCompare p w =
p^n•elem(Ap,Bp)`としてしか特徴づけられない)について、honest な
`elem`の annihilation性質(`one_tmul_sub_tmul_one_mul_elem`、`Ap,Bp`
レベル)を`diagonalCompare`で引き戻し、`diagonalCompare_injective`で
`B⊗_AB`へそのまま降ろす。これが前回セッションから持ち越されていた
「almost 化の技術的な壁」の解消そのもの。 -/
theorem almost_swap_annihilate {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Free A B] (p : A)
    (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p)))
    (q : B) :
    (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0 := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  haveI hAPFree := hAE.1
  haveI hAPFinite := (hAE.2.1 : Module.Finite _ _)
  haveI hAPEtale := hAE.2.2.1
  apply diagonalCompare_injective p hf0inj
  rw [map_zero, map_mul, map_sub, diagonalCompare_tmul, diagonalCompare_tmul, map_one, hw]
  have hkey := Algebra.FormallyUnramified.one_tmul_sub_tmul_one_mul_elem (R := Localization.Away p)
    (S := Localization.Away (algebraMap A B p)) (algebraMap B (Localization.Away (algebraMap A B p)) q)
  rw [mul_smul_comm]
  rw [show ((1:Localization.Away (algebraMap A B p)) ⊗ₜ[Localization.Away p] (algebraMap B (Localization.Away (algebraMap A B p)) q)
      - (algebraMap B (Localization.Away (algebraMap A B p)) q) ⊗ₜ[Localization.Away p] (1:Localization.Away (algebraMap A B p)))
      * Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p)) = 0 from hkey]
  simp

/-- `almost_swap_annihilate`の引き算形を等式形へ言い換えただけ
(`hochSectionOfWitness`が要求する形と揃える)。 -/
theorem almost_swap_mul_eq {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Free A B] (p : A)
    (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p)))
    (q : B) :
    (1:B) ⊗ₜ[A] q * w = q ⊗ₜ[A] (1:B) * w := by
  have h := almost_swap_annihilate p hAE hf0inj n w hw q
  rw [sub_mul] at h
  exact sub_eq_zero.mp h

/-- **augmentation の almost 版**: `μ(w) = p^n`(`B`の中で、局所化を
経ずに厳密に成り立つ)。`Ap`-値の等式`lmul' Ap (diagonalCompare p w) =
p^n•1 = algebraMap B Bp(lmul' A w)`(`lmul'_diagonalCompare`の自然性 +
`elem`の augmentation`lmul_elem`)を`algebraMap B Bp`の単射性で`B`へ
降ろす。`remark_iii_trace_identity`と並ぶ、Faltings remark(iii)/(v)
両方の土台になる「augmentation の厳密性」。 -/
theorem almost_swap_augment {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (p : A)
    (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p))) :
    Algebra.TensorProduct.lmul' A w = p ^ n • (1 : B) := by
  letI := awayAlgebra p (A := A) (B := B)
  haveI := awayScalarTower p (A := A) (B := B)
  haveI hEtale := hAE.2.2.1
  haveI hFinite := (hAE.2.1 : Module.Finite _ _)
  apply hf0inj
  have h1 := lmul'_diagonalCompare p w
  rw [hw] at h1
  have hpull := (Algebra.TensorProduct.lmul' (Localization.Away p) (S := Localization.Away (algebraMap A B p))).toLinearMap.map_smul_of_tower
    (p^n) (Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p)))
  rw [show (Algebra.TensorProduct.lmul' (Localization.Away p))
      (p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p)))
      = p ^ n • (Algebra.TensorProduct.lmul' (Localization.Away p)
          (Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p))))
      from hpull, Algebra.FormallyUnramified.lmul_elem] at h1
  have hLHS : (p ^ n • (1:Localization.Away (algebraMap A B p)))
      = algebraMap A (Localization.Away (algebraMap A B p)) (p ^ n) := by
    rw [Algebra.smul_def, mul_one]
  have hRHS : algebraMap B (Localization.Away (algebraMap A B p)) (p ^ n • (1:B))
      = algebraMap A (Localization.Away (algebraMap A B p)) (p ^ n) := by
    rw [Algebra.smul_def, mul_one, ← IsScalarTower.algebraMap_apply A B (Localization.Away (algebraMap A B p))]
  rw [← h1, hLHS, hRHS]

/-- **`hochSection`の一般化**: `elem`固定ではなく、任意の witness
`w`と、それが満たすべき swap 性質`hswap`だけから`S⊗_RS`-線形な
切断もどきを構成する。`hochSection R S = hochSectionOfWitness R S
(elem R S) one_tmul_mul_elem`として honest な場合を復元できる
(この事実自体は今回は使わないが、構成が同一であることの確認として
`hochSectionOfWitness_augment`で確かめられる)。 -/
noncomputable def hochSectionOfWitness (R S : Type u) [CommRing R] [CommRing S] [Algebra R S]
    (w : TensorProduct R S S) (hswap : ∀ q : S, (1:S) ⊗ₜ[R] q * w = q ⊗ₜ[R] (1:S) * w) :
    letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
    S →ₗ[TensorProduct R S S] TensorProduct R S S := by
  letI : Algebra (TensorProduct R S S) S := (Algebra.TensorProduct.lmul' R).toRingHom.toAlgebra
  refine
  { toFun := fun b => (b ⊗ₜ[R] (1:S)) * w
    map_add' := by intro x y; simp [TensorProduct.add_tmul, add_mul]
    map_smul' := ?_ }
  intro z b
  simp only [RingHom.id_apply, Algebra.smul_def]
  induction z using TensorProduct.induction_on with
  | zero => simp
  | tmul p q =>
    show ((Algebra.TensorProduct.lmul' R (p ⊗ₜ[R] q) * b) ⊗ₜ[R] (1:S)) * w
      = (p ⊗ₜ[R] q) * ((b ⊗ₜ[R] (1:S)) * w)
    rw [Algebra.TensorProduct.lmul'_apply_tmul]
    have hstep1 : ((p*q*b : S) ⊗ₜ[R] (1:S) : TensorProduct R S S)
        = ((p*b : S) ⊗ₜ[R] (1:S)) * (q ⊗ₜ[R] (1:S)) := by
      rw [Algebra.TensorProduct.tmul_mul_tmul]; ring_nf
    rw [hstep1, mul_assoc]
    have hstep2 : (q ⊗ₜ[R] (1:S) : TensorProduct R S S) * w
        = (1 ⊗ₜ[R] q) * w :=
      (hswap q).symm
    rw [hstep2, ← mul_assoc, Algebra.TensorProduct.tmul_mul_tmul, mul_one, one_mul]
    rw [← mul_assoc, Algebra.TensorProduct.tmul_mul_tmul, mul_one]
  | add x y hx hy =>
    simp only [map_add, add_mul, TensorProduct.add_tmul]
    rw [hx, hy]

theorem hochSectionOfWitness_augment (R S : Type u) [CommRing R] [CommRing S] [Algebra R S]
    (w : TensorProduct R S S) (hswap : ∀ q : S, (1:S) ⊗ₜ[R] q * w = q ⊗ₜ[R] (1:S) * w) (b : S) :
    Algebra.TensorProduct.lmul' R (hochSectionOfWitness R S w hswap b)
      = b * Algebra.TensorProduct.lmul' R w := by
  show Algebra.TensorProduct.lmul' R ((b ⊗ₜ[R] (1:S)) * w) = b * Algebra.TensorProduct.lmul' R w
  rw [map_mul, Algebra.TensorProduct.lmul'_apply_tmul, mul_one]

/-- **almost 版 augmentation、具体形**: `hochSectionOfWitness`を
almost witness `w`に適用したときの合成`μ∘s`が「`B`の元`p^n`との掛け算」
そのものになる。これが remark(v)の「`S`が`S⊗AS`の direct summand
"up to p^n"」という主張の、`Definition 2.1`一般の`IsAlmostEtaleCovering`
仮定だけからの完全な形式化——honest な`Algebra.FormallyUnramified A B`
を一切要求しない点が`hochModule_projective`(honest な場合)からの
真の一般化になっている。 -/
theorem hochSectionAlmost_augment {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Free A B] (p : A)
    (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p) (Localization.Away (algebraMap A B p)))
    (b : B) :
    Algebra.TensorProduct.lmul' A (hochSectionOfWitness A B w (almost_swap_mul_eq p hAE hf0inj n w hw) b)
      = algebraMap A B (p ^ n) * b := by
  rw [hochSectionOfWitness_augment, almost_swap_augment p hAE hf0inj n w hw]
  rw [Algebra.smul_def, mul_one, mul_comm]

/-! **現状のまとめ(2026-09-05)**: 上記4定理・2定義により、
「`IsAlmostEtaleCovering A B p`(honest な`Algebra.FormallyUnramified
A B`を要求しない、真に一般の almost 設定)の下で、`B`は`B⊗_AB`-線形
写像`hochSectionOfWitness`を経て`μ:B⊗_AB→B`の"p^n-almost section"を
持つ」ことが完全に証明された。remark(v)の almost 版を`Ext`の消滅
(`p^n•Ext^k(B,M)=0`、honest な場合の`hochschild_ext_eq_zero`の直接
一般化)まで持ち上げるには、`Ext^k_{T}(-,M)`(`T:=B⊗_AB`)の関手性
(`s:B→T`・`μ:T→B`からの`Ext^k(s)`・`Ext^k(μ)`と、その合成が
「`B`の元`algebraMap A B(p^n)`による掛け算」の pullback に一致する
という事実)を要する——`CategoryTheory.Abelian.Ext`のこの向きの
API(pre/post-composition による Hom 群上の作用)がどこまで揃って
いるかの調査が次回の入口。 -/

end ABC3.Found.Falt1
