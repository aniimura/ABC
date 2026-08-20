/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex61Model
import ABC3.Found.Divisor.Ex63Model
import ABC3.Found.FrdI.Def24RlfCone

/-!
# `Φ^rlf` は `Example 6.1` / `Example 6.3` で `𝒟` 上の単系になる(鎖 `rlf` の `rlf-flat` の実例)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.103。

原文 (FrdI p.103):
> Then the realification Φrlf of Φ [cf. Definition 2.4, (i)] determines a monoid Φrlf on D

## ★★錐が sharp であることは「係数の総和」で出る

`Def24RlfCone.lean` は `Φ^rlf` を**錐模型**で立て、
残るのは**錐が sharp** であることだけにしてある(`phiRlfConeOnOfPos`)。
★その十分条件は「`Φ(A)` 上で狭義正の線型汎関数」であり、
★★**Cartier 因子でも算術因子でも「係数の総和」がそれである** ——
`Φ(A) = Γ ∩ ℤ≥0[S]` の元は成分がすべて `≥ 0` で、非零なら
どこか `> 0` だから、総和は `> 0` である。

## ★本ファイルで閉じること

| 定義 | 中身 |
|---|---|
| `sumHomR` | 係数の総和 `(S →₀ ℤ) →+ ℝ` |
| `effSubPosFunctional` | ★`Φ^gp → ℝ`(総和) |
| `phiRlfOfCartierDatum` | ★★★**`Φ^rlf` は `𝒟` 上の単系**(`CartierDatum` から) |
-/

namespace ABC3.Found.Divisor

open CategoryTheory ABC3.Found.FrdI

universe v u w

/-! ## ★1. 係数の総和 -/

/-- ★**係数の総和**。 -/
noncomputable def sumHomR (S : Type*) : (S →₀ ℤ) →+ ℝ :=
  Finsupp.liftAddHom (fun _ => Int.castAddHom ℝ)

@[simp] theorem sumHomR_apply {S : Type*} (d : S →₀ ℤ) :
    sumHomR S d = d.sum fun _ c => (c : ℝ) := rfl

/-- ★★**有効な非零因子の係数の総和は正**。 -/
theorem sumHomR_pos {S : Type*} {d : S →₀ ℤ} (hnn : ∀ s, 0 ≤ d s) (hne : d ≠ 0) :
    0 < sumHomR S d := by
  rw [sumHomR_apply, Finsupp.sum]
  obtain ⟨s, hs⟩ : ∃ s, d s ≠ 0 := by
    by_contra hc
    exact hne (Finsupp.ext (by simpa using hc))
  refine Finset.sum_pos' (fun i _ => by exact_mod_cast hnn i) ⟨s, ?_, ?_⟩
  · exact Finsupp.mem_support_iff.mpr hs
  · exact_mod_cast lt_of_le_of_ne (hnn s) (Ne.symm hs)

/-! ## ★2. `Φ^gp → ℝ` -/

variable {D : Type u} [Category.{v} D] (Δ : CartierDatum.{v, u, w} D) (hD : IsOfFSMType D)

/-- ★`Φ^gp → ℝ`(係数の総和)。 -/
noncomputable def effSubPosFunctional (A : D) : Gp ((Δ.phi hD).val A) →+ ℝ :=
  (sumHomR (Δ.primes A)).comp (((Δ.grp A).subtype).comp (effSubGpHom (Δ.grp A)))

theorem effSubPosFunctional_pos (A : D) (m : (Δ.phi hD).val A) (hm : m ≠ 0) :
    0 < effSubPosFunctional Δ hD A (toGp _ m) := by
  have hkey : ((effSubGpHom (Δ.grp A) (toGp ((Δ.phi hD).val A) m) : Δ.grp A)
      : Δ.primes A →₀ ℤ) = m.1 := phiGpHomC_toGp Δ hD m
  show 0 < sumHomR (Δ.primes A)
    ((effSubGpHom (Δ.grp A) (toGp ((Δ.phi hD).val A) m) : Δ.grp A) : Δ.primes A →₀ ℤ)
  rw [hkey]
  exact sumHomR_pos (fun s => (mem_effSub.mp m.2).2 s) (fun h => hm (Subtype.ext h))

/-- ★★★★**`Φ^rlf` は `𝒟` 上の単系**(`CartierDatum` から)。

★錐の sharp 性は「係数の総和」が狭義正の線型汎関数であることから出る。 -/
noncomputable def phiRlfOfCartierDatum : MonoidOn.{v, u, w} D :=
  phiRlfConeOnOfPos (Δ.phi hD) (Δ.phi_isDivisorialOn hD)
    (effSubPosFunctional Δ hD) (effSubPosFunctional_pos Δ hD)

/-! ## ★3. 算術の場合(`Example 6.3`) -/

/-- ★**係数の総和**(実係数)。 -/
noncomputable def sumHomRR (S : Type*) : (S →₀ ℝ) →+ ℝ :=
  Finsupp.liftAddHom (fun _ => AddMonoidHom.id ℝ)

@[simp] theorem sumHomRR_apply {S : Type*} (d : S →₀ ℝ) :
    sumHomRR S d = d.sum fun _ c => c := rfl

theorem sumHomRR_pos {S : Type*} {d : S →₀ ℝ} (hnn : ∀ s, 0 ≤ d s) (hne : d ≠ 0) :
    0 < sumHomRR S d := by
  rw [sumHomRR_apply, Finsupp.sum]
  obtain ⟨s, hs⟩ : ∃ s, d s ≠ 0 := by
    by_contra hc
    exact hne (Finsupp.ext (by simpa using hc))
  exact Finset.sum_pos' (fun i _ => hnn i)
    ⟨s, Finsupp.mem_support_iff.mpr hs, lt_of_le_of_ne (hnn s) (Ne.symm hs)⟩

section Arith

variable (Δ' : ArithDatum.{v, u, w} D) (hD' : IsOfFSMType D)

/-- ★`Φ^gp → ℝ`(算術の場合、係数の総和)。 -/
noncomputable def effRPosFunctional (A : D) : Gp ((Δ'.phi hD').val A) →+ ℝ :=
  (sumHomRR (Δ'.primes A)).comp (((Δ'.grp A).subtype).comp (effRGpHom (Δ'.grp A)))

theorem effRPosFunctional_pos (A : D) (m : (Δ'.phi hD').val A) (hm : m ≠ 0) :
    0 < effRPosFunctional Δ' hD' A (toGp _ m) := by
  have hkey : ((effRGpHom (Δ'.grp A) (toGp ((Δ'.phi hD').val A) m) : Δ'.grp A)
      : Δ'.primes A →₀ ℝ) = m.1 := phiGpHom_toGp Δ' hD' m
  show 0 < sumHomRR (Δ'.primes A)
    ((effRGpHom (Δ'.grp A) (toGp ((Δ'.phi hD').val A) m) : Δ'.grp A) : Δ'.primes A →₀ ℝ)
  rw [hkey]
  exact sumHomRR_pos (fun s => (mem_effR.mp m.2).2 s) (fun h => hm (Subtype.ext h))

/-- ★★★★**`Φ^rlf` は `𝒟` 上の単系**(`ArithDatum` から、`Example 6.3`)。 -/
noncomputable def phiRlfOfArithDatum : MonoidOn.{v, u, w} D :=
  phiRlfConeOnOfPos (Δ'.phi hD') (Δ'.phi_isDivisorialOn hD')
    (effRPosFunctional Δ' hD') (effRPosFunctional_pos Δ' hD')

end Arith

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Proposition 5.3` の `Φ^rlf`(`Example 6.1` の場合)。 -/
def phiRlfOfCartierDatum.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — Φ^rlf は monoid on 𝒟(Cartier 因子の場合)",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.Divisor
