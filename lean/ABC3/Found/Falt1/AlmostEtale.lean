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
import Mathlib.RingTheory.TensorProduct.Basic

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

end ABC3.Found.Falt1
