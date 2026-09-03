/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtFaltCovolume
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.IsogenyPeriodPair
import ABC3.Found.GenEll.LatticeScale
import ABC3.Found.GenEll.Uniformization
import ABC3.Meta.Claim
import ABC3.Found.GaloisRep.VeluNormalized.Proposition34

/-!
# VeluNormalized —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open NumberField IsDedekindDomain WeierstrassCurve ABC3.Found.GenEll
open scoped Classical
variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5` の高さ評価——`L` 上の位数 `l` の点から**

    `ht^Falt(E/⟨Q⟩) ≤ ht^Falt(E) + 2·log(l)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★★★★☆**これが `Lemma 3.5` の解析側の到達点である。**

★組み立て:

* `S ≔ ⟨Q⟩∖{O}` の座標（`L × L` の有限集合）、`E′ ≔ veluQuotientFull E S`
* 各 `σ` で `rhPoint σ` により `Q` を `E.map σ` の位数 `l` の点へ（第 696）、
  点集合は `σ(S)` になる（第 703）
* 第 702（`ℂ` 側最終形）で `latticeCurve P′ = C σ • veluQuotientFull (E.map σ) (σS)`
* `veluQuotientFull` は底変換と可換（第 679）なので
  `= C σ • (E′.map σ)`——**同じ変数変換 `C σ`、すなわち `α = 1`**
* 第 678 に渡す

☆残る仮定は `hmin`（`E` が極小）・`hint`（`E′` が整）だけである。 -/
theorem htFalt_veluQuotientFull_le (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ) (hl : 0 < l)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (P : (L →+* ℂ) → PeriodPair) (C : (L →+* ℂ) → VariableChange ℂ)
    (hΔ : ∀ σ, latticeDisc (P σ) ≠ 0)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hell1 : ∀ σ : L →+* ℂ, (E.map σ).IsElliptic)
    (hell2 : ∀ σ : L →+* ℂ, (C σ • (E.map σ)).IsElliptic)
    (hmin : ∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0)
    (hint : ∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  have key : ∀ σ : L →+* ℂ, ∃ (P' : PeriodPair) (A B Cc D : ℤ),
      (P σ).ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      (P σ).ω₂ = (Cc : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * Cc).natAbs = l ∧
      C σ • (E'.map σ) = latticeCurve P' := by
    intro σ
    haveI := hell1 σ
    haveI := hell2 σ
    have hQσ : addOrderOf (rhPoint σ E Q) = l := by
      rw [addOrderOf_rhPoint]; exact hQ
    obtain ⟨P', A, B, Cc, D, h1, h2, hdet, hEq⟩ :=
      exists_periodPair_veluQuotientFull (E.map σ) (C σ) (P σ) (hΔ σ) (hPC σ) hl hQσ
    refine ⟨P', A, B, Cc, D, h1, h2, hdet, ?_⟩
    rw [hEq, image_pointCoords_rhPoint_nsmul σ E hQ, hE', veluQuotientFull_map]
  choose P' A B Cc D h1 h2 hdet hPC' using key
  exact htFalt_isogeny_le_of_velu E E' l hl P P' C hPC hPC' A B Cc D h1 h2 hdet hmin hint

/-- ★★★★★★★★★★★★★★★★**`Lemma 3.5` の高さ評価——有限素点側を
1 本の不等式で受ける形**（第 1049）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1049）の測定**——第 704 は
`hmin`（`E` が大域極小）と `hint`（`E′` が整）を受けていたが、
それらは**ただ 1 本の不等式 `hfin` を出すためだけ**に使われていた
（`htFalt_isogeny_le_of_archDefect_minimal` の証明を実測）。
☆本定理はその `hfin` を直接受ける形である。 -/
theorem htFalt_veluQuotientFull_le_of_defect (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ) (hl : 0 < l)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (P : (L →+* ℂ) → PeriodPair) (C : (L →+* ℂ) → VariableChange ℂ)
    (hΔ : ∀ σ, latticeDisc (P σ) ≠ 0)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hell1 : ∀ σ : L →+* ℂ, (E.map σ).IsElliptic)
    (hell2 : ∀ σ : L →+* ℂ, (C σ • (E.map σ)).IsElliptic)
    (hfin : (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
          - (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        ≤ (3 / 2) * (Module.finrank ℚ L : ℝ) * Real.log l) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  have key : ∀ σ : L →+* ℂ, ∃ (P' : PeriodPair) (A B Cc D : ℤ),
      (P σ).ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      (P σ).ω₂ = (Cc : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * Cc).natAbs = l ∧
      C σ • (E'.map σ) = latticeCurve P' := by
    intro σ
    haveI := hell1 σ
    haveI := hell2 σ
    have hQσ : addOrderOf (rhPoint σ E Q) = l := by
      rw [addOrderOf_rhPoint]; exact hQ
    obtain ⟨P', A, B, Cc, D, h1, h2, hdet, hEq⟩ :=
      exists_periodPair_veluQuotientFull (E.map σ) (C σ) (P σ) (hΔ σ) (hPC σ) hl hQσ
    refine ⟨P', A, B, Cc, D, h1, h2, hdet, ?_⟩
    rw [hEq, image_pointCoords_rhPoint_nsmul σ E hQ, hE', veluQuotientFull_map]
  choose P' A B Cc D h1 h2 hdet hPC' using key
  refine htFalt_isogeny_le_of_archDefect E E' l ?_ hfin
  exact archDefect_isogeny E E' l hl P P' C C hPC hPC' (fun _ => 1) (fun _ => one_ne_zero)
    (fun _ => by ring) A B Cc D (fun σ => by rw [one_mul]; exact h1 σ)
    (fun σ => by rw [one_mul]; exact h2 σ) hdet

def htFalt_veluQuotientFull_le_of_defect.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(高さ評価——有限素点側を 1 本の不等式で受ける形)",
    sectionId := "genell-lemma-3-5" }

def htFalt_veluQuotientFull_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ht^Falt(E/⟨Q⟩) ≤ ht^Falt(E) + 2 log l——解析側の到達点)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
