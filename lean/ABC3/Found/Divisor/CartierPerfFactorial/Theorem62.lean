/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.CartierMonoprime
import ABC3.Found.FrdI.Def24
import ABC3.Found.FrdI.Def24SuppElt
import ABC3.Found.Divisor.CartierPerfFactorial.Example61

/-!
# CartierPerfFactorial —— `[FrdI] Theorem 6.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.FrdI

open Finsupp
variable {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)}

/-! ## ★★★台(`SuppElt`)と係数の対応 —— `Theorem 6.2, (iii)` の rationally standard に要る -/

/-- ★★**`SuppElt` は係数が 0 でないことに他ならない**。 -/
theorem mem_suppElt_iff (hQ : IsQCartierSubgroup Γ) (a : effSub Γ) (p : Prime (effSub Γ)) :
    p ∈ SuppElt (iotaAt hQ) a ↔ ((a : effSub Γ) : S →₀ ℤ) (effSubPrimeEquiv hQ p) ≠ 0 := by
  have h1 : pfCoeffHom (Pf.of a) = intToRatFs ((a : effSub Γ) : S →₀ ℤ) := pfCoeff_of a
  have h2 := congrArg (fun f : S →₀ ℚ => f (effSubPrimeEquiv hQ p)) h1
  have h4 : (intToRatFs ((a : effSub Γ) : S →₀ ℤ)) (effSubPrimeEquiv hQ p)
      = ((((a : effSub Γ) : S →₀ ℤ) (effSubPrimeEquiv hQ p) : ℤ) : ℚ) := rfl
  have hkey : ratCoeffAt hQ p (Pf.mk a 1)
      = ((((a : effSub Γ) : S →₀ ℤ) (effSubPrimeEquiv hQ p) : ℤ) : ℚ) := by
    show pfCoeffHom (Pf.of a) (effSubPrimeEquiv hQ p) = _
    rw [h2, h4]
  show factorMap (iotaAt hQ) (Pf.mk a 1) p ≠ 0 ↔ _
  rw [factorMap_iotaAt hQ (Pf.mk a 1) p]
  constructor
  · intro hne hc
    refine hne ?_
    rw [iotaAt_eq_zero_iff, hkey, hc]
    rfl
  · intro hne hc
    rw [iotaAt_eq_zero_iff, hkey] at hc
    exact hne (by exact_mod_cast hc)

/-- ★★★★**`IsStrictlyRational` が要求する分割**(台の言葉で)。

★`Φ^birat` の元 `y` が素点 `p` の係数を持てば、`p` で正・`p` を含まない
`a, b ∈ Φ(L)` に割れる。★`Φ^birat` は部分群なので `±y` のどちらでもよい。 -/
theorem exists_split_suppElt_of_qc [DecidableEq S] (hQ : IsQCartierSubgroup Γ)
    {y : S →₀ ℤ} (hy : y ∈ Γ) (p : Prime (effSub Γ))
    (hyp : y (effSubPrimeEquiv hQ p) ≠ 0) :
    ∃ (a b : effSub Γ) (z : S →₀ ℤ), z ∈ Γ ∧ (z = y ∨ z = -y) ∧
      ((a : effSub Γ) : S →₀ ℤ) = ((b : effSub Γ) : S →₀ ℤ) + z ∧
      p ∈ SuppElt (iotaAt hQ) a ∧ p ∉ SuppElt (iotaAt hQ) b := by
  set s := effSubPrimeEquiv hQ p with hs
  rcases lt_or_gt_of_ne hyp with hneg | hpos
  · have hnegy : 0 < (-y) s := by
      show 0 < -(y s)
      linarith
    obtain ⟨a, b, ha, hb, hab, hap, hbp⟩ := exists_split_of_qc hQ (Γ.neg_mem hy) hnegy
    refine ⟨⟨a, ha⟩, ⟨b, hb⟩, -y, Γ.neg_mem hy, Or.inr rfl, hab, ?_, ?_⟩
    · exact (mem_suppElt_iff hQ ⟨a, ha⟩ p).mpr (ne_of_gt hap)
    · intro hc
      exact ((mem_suppElt_iff hQ ⟨b, hb⟩ p).mp hc) hbp
  · obtain ⟨a, b, ha, hb, hab, hap, hbp⟩ := exists_split_of_qc hQ hy hpos
    refine ⟨⟨a, ha⟩, ⟨b, hb⟩, y, hy, Or.inl rfl, hab, ?_, ?_⟩
    · exact (mem_suppElt_iff hQ ⟨a, ha⟩ p).mpr (ne_of_gt hap)
    · intro hc
      exact ((mem_suppElt_iff hQ ⟨b, hb⟩ p).mp hc) hbp

/-- ★★★locator —— `Theorem 6.2, (iii)` の rationally standard に要る分割(台の言葉)。 -/
def exists_split_suppElt_of_qc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — Φ^birat の元は p で正・p を含まない差に割れる",
    sectionId := "frdi-thm-6-2" }

end ABC3.Found.FrdI
