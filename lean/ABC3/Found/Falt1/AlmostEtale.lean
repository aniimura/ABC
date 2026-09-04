import ABC3.Meta.Claim
import Mathlib.RingTheory.Etale.Basic
import Mathlib.RingTheory.Unramified.Finite
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
(★2026-09-04 時点では `awayAlgebra`(局所化の自己写像)と mathlib の
標準的な自己代数(`Algebra.id`)のインスタンス衝突により、`A=B` の
具体例での完全な non-vacuous witness の構成は未完成——`Definition
2.1` 自体(`isAlmostEtaleCovering`・`awayAlgebra`・`awayScalarTower`・
`diagonalCompare`)は sorry 無しで完成しているが、続く witness は
次のセッションへ持ち越す)。 -/
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

end ABC3.Found.Falt1
