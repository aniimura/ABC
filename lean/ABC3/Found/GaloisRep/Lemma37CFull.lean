/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma37C
import ABC3.Found.GaloisRep.Lemma35Unconditional
import ABC3.Meta.Claim

/-!
# 第 1088 ブロック —— **`Lemma 3.7` の条件 (a) から `(†)` が消えた**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★これは何か

`Found/GaloisRep/Lemma37C.lean` は「受けているのは `(†)`（＝ `Lemma 3.5` の結論）だけ」
と書いていた。★第 1087 で `Lemma 3.5` が条なしになったので、**その `(†)` を埋める**。

## ★★★★`l ≠ 2` と `d + 1 < l` は条件 (a) から出る

`Lemma 3.5`（第 1085）は原典に無い仮説 `l ≠ 2`・`d + 1 < l` を置いている。
☆条件 (a) はそれらを**含意する**:

    `C ≔ C₀ + |B| + 1`（`B ≤ ht^Falt` は下界、`C₀ > 0`）と取ると
    `ht^Falt + C·d^ε ≥ B + C₀ + |B| + 1 ≥ 1`
    よって `l ≥ 100·d·1 = 100·d ≥ 100`。

★`l ≥ 100` から `l ≠ 2`、`l ≥ 100·d ≥ d + 99·d ≥ d + 99` から `d + 1 < l`。
☆すなわち **`Lemma 3.5` の逸脱は本ファイルで完全に吸収される**。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open ABC3.Found.GenEll ABC3.Meta
open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★
**条件 (a) だけで `ht^Falt ≤ C′`**（第 1088）——`(†)` は要らない。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★★★★**2026-09-01（第 1088）**——`Lemma 3.5`（第 1085）を代入した形である。
☆`Lemma 3.5` が置いている `l ≠ 2`・`d + 1 < l` は条件 (a) から従うので、
外から受ける仮説は**原典の Lemma 3.7 が置くものだけ**になる。 -/
theorem htFalt_le_of_condA_lcyclic (eps : ℝ) (heps : 0 < eps) :
    ∃ C C' : ℝ, 0 < C ∧
      ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
        [E.IsElliptic] [E'.IsElliptic] (l : ℕ) (p : HeightOneSpectrum (𝓞 L)),
        l.Prime →
        minDeltaExp p E ≠ 0 →
        100 * (Module.finrank ℚ L : ℝ)
            * (htFaltOf L E + C * (Module.finrank ℚ L : ℝ) ^ eps) ≤ (l : ℝ) →
        ∀ Q : E.toAffine.Point, addOrderOf Q = l →
        E' = veluQuotientFull E (((Finset.range l).erase 0).image
            (fun k : ℕ => pointCoords (k • Q))) →
        (∀ r, SemistableAt r E) →
        (∀ r, SemistableAt r E') →
        (∀ r : HeightOneSpectrum (𝓞 L), jExp r E < 0 → ¬ ((l : ℤ) ∣ jExp r E)) →
        htFaltOf L E ≤ C' := by
  obtain ⟨C₀, hC₀0, hC₀⟩ := htFalt_le_of_condA eps heps
  obtain ⟨B, hB⟩ := exists_htFalt_bddBelow
  obtain ⟨C₅, h₅⟩ := lemma_3_5_height_ineq (1/6) (by norm_num)
  refine ⟨C₀ + |B| + 1, max C₅ 0, by positivity, ?_⟩
  intro L _ _ E E' _ _ l p hl hne hcondA Q hQ hE' hssE hssE' hcop
  set d : ℝ := (Module.finrank ℚ L : ℝ) with hddef
  have hd1 : (1:ℝ) ≤ d := by rw [hddef]; exact_mod_cast Module.finrank_pos
  have hP1 : (1:ℝ) ≤ d ^ eps := Real.one_le_rpow hd1 heps.le
  have hBabs : (0:ℝ) ≤ B + |B| := by
    have := neg_abs_le B; linarith
  -- ☆条件 (a) から `100·d ≤ l`
  have hone : (1:ℝ) ≤ htFaltOf L E + (C₀ + |B| + 1) * d ^ eps := by
    have h1 : B ≤ htFaltOf L E := hB L E
    have h2 : C₀ + |B| + 1 ≤ (C₀ + |B| + 1) * d ^ eps := by nlinarith [abs_nonneg B]
    linarith [hBabs, hC₀0]
  have hbig : 100 * d ≤ (l : ℝ) := by nlinarith [hcondA, hd1, hone]
  have hbig' : (100 : ℝ) * (Module.finrank ℚ L : ℝ) ≤ (l : ℝ) := by
    rw [hddef] at hbig; exact hbig
  have hbigN : 100 * Module.finrank ℚ L ≤ l := by exact_mod_cast hbig'
  have hdN : 1 ≤ Module.finrank ℚ L := Module.finrank_pos
  -- ☆`Lemma 3.5` で `(†)` を埋める
  have hdag : ((l : ℝ) / 14) * degInfOf L E
      ≤ htFaltOf L E + 2 * Real.log l + max C₅ 0 := by
    have h := h₅ L E E' l hl Q hQ hE' hssE hssE' hcop
    have hcoef : (1 : ℝ) / (12 * (1 + 1/6)) = 1 / 14 := by norm_num
    rw [hcoef] at h
    have hrw : (1 : ℝ) / 14 * (l : ℝ) * degInfOf L E
        = ((l : ℝ) / 14) * degInfOf L E := by ring
    rw [hrw] at h
    exact le_trans h (by linarith [le_max_left C₅ (0:ℝ)])
  -- ☆`C₀` の側の条件 (a)
  have hcondA₀ : 100 * d * (htFaltOf L E + C₀ * d ^ eps) ≤ (l : ℝ) := by
    refine le_trans ?_ hcondA
    have hle : C₀ * d ^ eps ≤ (C₀ + |B| + 1) * d ^ eps := by nlinarith [abs_nonneg B]
    nlinarith [hd1]
  exact hC₀ (max C₅ 0) (le_max_right _ _) L E l p hl.one_lt.le hne hcondA₀ hdag

def htFalt_le_of_condA_lcyclic.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 3 の主張・条件 (a)——(†) を埋めた形。★外部入力なし)",
    sectionId := "genell-lemma-3-7" }

def htFalt_le_of_condA_lcyclic.needs : List ProofObligation :=
  [ .citation "[ABC3]" "htFalt_le_of_condA(数値の核、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_le_of_condA") 1,
    .citation "[ABC3]" "lemma_3_5_height_ineq((†)、第 1085、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_height_ineq") 1 ]

end ABC3.Found.GaloisRep
