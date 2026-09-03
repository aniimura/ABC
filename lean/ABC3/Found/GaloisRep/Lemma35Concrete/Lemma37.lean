/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Found.GaloisRep.HtFaltBounds
import ABC3.Found.GaloisRep.VeluNormalized
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Meta.Claim
import ABC3.Found.GaloisRep.Lemma35Concrete.Lemma35

/-!
# Lemma35Concrete —— `[GenEll] Lemma 3.7` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical
variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★★★★★★★★★`Lemma 3.7` が受けている `hdag` の形 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.7` が受けている `hdag`（`ε₀ = 1/6`、すなわち `12(1+ε₀) = 14`）**

    `(l/14)·deg_∞(E) ≤ ht^Falt(E) + 2·log(l) + C′`

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★★★★★★★☆`Found/GaloisRep/Lemma37C.lean` の `finite_j_of_condA` が
**唯一の入力**として受けている `hdag` はこの形である。
本補題はそれを `lemma_3_5_velu_bad`（第 707）から出す——
すなわち **`Lemma 3.7` の第 3 の主張も未証明の外部引用から自由になった**。

☆`C′ ≥ 0` は `max` で取り直しているだけである。 -/
theorem hdag_of_velu :
    ∃ C' : ℝ, 0 ≤ C' ∧ ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E') →
      ∀ S : Set (HeightOneSpectrum (𝓞 L)),
      (∀ p ∈ S, minDeltaExp p E' = l * minDeltaExp p E) →
      (∀ p ∉ S, minDeltaExp p E = 0) →
      (∀ p ∉ S, minDeltaExp p E' = 0) →
      ((l : ℝ) / 14) * degInfOf L E ≤ htFaltOf L E + 2 * Real.log l + C' := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu_bad (1 / 6) (by norm_num)
  refine ⟨max C 0, le_max_right _ _, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC
    hell1 hell2 hmin hint hss S hbad hgoodE hgoodE' => ?_⟩
  have h := hC L E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint hss S
    hbad hgoodE hgoodE'
  have hcoef : (1 / (12 * (1 + (1 / 6 : ℝ)))) * (l : ℝ) * degInfOf L E
      = ((l : ℝ) / 14) * degInfOf L E := by
    rw [show (12 * (1 + (1 / 6 : ℝ))) = 14 by norm_num]
    ring
  rw [hcoef] at h
  exact le_trans h (by have := le_max_left C 0; linarith)

def hdag_of_velu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 3 の主張が受けている hdag——Lemma 3.5 から出した形)",
    sectionId := "genell-lemma-3-7" }

def hdag_of_velu.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_velu_bad(Lemma 3.5、外部引用なし、§9-1149)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_velu_bad") 2 ]

end ABC3.Found.GaloisRep
