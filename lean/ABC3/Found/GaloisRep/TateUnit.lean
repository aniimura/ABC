import ABC3.Found.GaloisRep.TateDelta
import Mathlib.RingTheory.PowerSeries.Inverse

/-!
# Galois (G6) 第 98 ブロック —— **★★★★局所高さ = 判別式の付値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★これが `Definition 3.3` の算術的な中身である

第 97 ブロックで `Δ(E_q) = q + O(q²)` を得た。★形式冪級数では

    Δ = X · u,      u は**単元**(定数項が 1)

★★特殊化すると `Δ(E_{q₀}) = q₀ · (R の単元)` であり、したがって

    v(Δ) = v(q₀)

★★★**局所高さ `v_K(q_E)` は判別式の付値である**——これが原文の
`Definition 3.3` を計算可能にする形である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `evalAdic_X` | ★`X` の値は `q` |
| `tateDelta_eq_X_mul_unit` | ★★★**`Δ = X · 単元`** |
| `tateCurveAt_Delta_eq_mul_unit` | ★★★★**`Δ(E_q) = q · 単元`** |
-/

namespace ABC3.Found.GaloisRep

open PowerSeries

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★`X` の値は `q` である。 -/
theorem evalAdic_X [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    evalAdic PowerSeries.X q hq = q := by
  refine evalAdic_unique _ _ _ _ (fun n => ?_)
  match n with
  | 0 =>
    have htop : (I ^ 0 • (⊤ : Submodule R R)) = ⊤ := by simp
    rw [SModEq.sub_mem, htop]
    exact Submodule.mem_top
  | 1 =>
    have hp : partialEval (PowerSeries.X : PowerSeries ℤ) q 1 = 0 := by
      rw [partialEval]
      simp
    rw [SModEq.sub_mem, hp]
    simpa using neg_mem (by simpa using hq : q ∈ I ^ 1)
  | (m + 2) =>
    have hp : partialEval (PowerSeries.X : PowerSeries ℤ) q (m + 2) = q := by
      rw [partialEval, Finset.sum_eq_single 1]
      · simp
      · intro b _ hb1
        rcases Nat.eq_zero_or_pos b with hb0 | hb0
        · subst hb0; simp
        · have : PowerSeries.coeff b (PowerSeries.X : PowerSeries ℤ) = 0 := by
            rw [PowerSeries.coeff_X]
            exact if_neg hb1
          rw [this]
          simp
      · intro h
        exact absurd (Finset.mem_range.2 (by omega)) h
    rw [hp]

/-- ★★★**`Δ = X · 単元`**。 -/
theorem tateDelta_eq_X_mul_unit :
    ∃ u : PowerSeries ℤ, IsUnit u ∧ tateCurve.Δ = PowerSeries.X * u := by
  have hdvd : PowerSeries.X ∣ tateCurve.Δ := by
    rw [PowerSeries.X_dvd_iff, ← PowerSeries.coeff_zero_eq_constantCoeff_apply]
    exact coeff_zero_tateDelta
  obtain ⟨u, hu⟩ := hdvd
  refine ⟨u, ?_, hu⟩
  rw [PowerSeries.isUnit_iff_constantCoeff]
  have h1 : PowerSeries.constantCoeff u = 1 := by
    have hc := coeff_one_tateDelta
    rw [hu, PowerSeries.coeff_succ_X_mul 0] at hc
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]
    exact hc
  rw [h1]
  exact isUnit_one

/-- ★★★★**`Δ(E_q) = q · 単元`**——局所高さが判別式の付値であることの中身。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateCurveAt_Delta_eq_mul_unit [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    ∃ v : R, IsUnit v ∧ (tateCurveAt q hq).Δ = q * v := by
  obtain ⟨u, hu, hΔ⟩ := tateDelta_eq_X_mul_unit
  refine ⟨evalAdicHom q hq u, hu.map (evalAdicHom q hq), ?_⟩
  rw [tateCurveAt_Delta, hΔ, map_mul]
  congr 1
  exact evalAdic_X q hq

/-! ## ★出典の紐付け(`.src`) -/

def tateCurveAt_Delta_eq_mul_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(局所高さが判別式の付値であること)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
