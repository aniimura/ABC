/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Lemma37Full
import ABC3.Meta.Claim

/-!
# 第 1210 ブロック —— **`hdag` を「導ける形」に弱める**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★★★★★これは何か——仮説を実際に導ける形へ

第 1166 の `lemma_3_7_stableLine` は `hdag` を

    ∀ E l, l 素数 → HasLCyclicJ E l → (†)

の形で受け取っていた。☆しかし `Lemma 3.5` の不等式は
**`l` が局所高さと素であること**を要るので、この形の `hdag` は
そのままでは導けない。

★本ブロックは `hdag` に `PrimeToLocalHeights` を**追加の仮説として持たせる**:

    ∀ E l, l 素数 → HasLCyclicJ E l → PrimeToLocalHeights l → (†)

☆使う側では (a)(b) のどちらの枝でも `PrimeToLocalHeights` が手元にある
——(a) では `lemma_3_7_a_coprime`、(b) では条件の定義そのもの——ので、
**結論は第 1166 とまったく同じ**である。

★★★これで残るギャップ `hdag` は**実際に導ける形**になった。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve ABC3.Meta
open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.7 —— 安定直線の側、`hdag` は素性も受け取る**（第 1210）。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

☆第 1166 との違いは `hdag` が `PrimeToLocalHeights` も仮説に持つことだけ。
★使う側では (a)(b) のどちらの枝でもそれが手元にあるので、結論は同じである。

★★★これで残るギャップは**実際に導ける形**の `hdag` ただ 1 つになった。 -/
theorem lemma_3_7_stableLine_cop (KV : Set ℂ) (hKV : CompactlyBoundedJ KV)
    (eps : ℝ) (heps : 0 < eps) (C₅ : ℝ) (hC₅ : 0 ≤ C₅)
    (hdag : ∀ (E : DegCurve) (l : ℕ), Nat.Prime l → 5 ≤ l → HasLCyclicJ E.toSSCurve l →
      E.toSSCurve.PrimeToLocalHeights l →
      ((l : ℝ) / 14) * degInfOf E.toSSCurve.fld E.toSSCurve.W
        ≤ htFaltOf E.toSSCurve.fld E.toSSCurve.W + 2 * Real.log l + C₅) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set ℂ, GaloisFiniteJ Exc ∧
      ∀ (E : DegCurve) (l : ℕ), Nat.Prime l → 5 ≤ l →
        ∀ condA condB : Prop,
          (condA ↔ (100 * (E.deg : ℝ)
                      * (faltingsHeightJ E.j + C * (E.deg : ℝ) ^ eps) ≤ (l : ℝ)
                    ∧ E.toSSCurve.HasMultRed)) →
          (condB ↔ (E.j ∈ KV ∧ E.toSSCurve.PrimeToLocalHeights l)) →
          (condA → E.toSSCurve.PrimeToLocalHeights l)
        ∧ (condB → E.j ∉ Exc → E.toSSCurve.HasMultRed)
        ∧ ((condA ∨ condB) → HasLCyclicJ E.toSSCurve l → E.j ∈ Exc) := by
  obtain ⟨M, hM⟩ := hKV
  obtain ⟨Ca, hCa0, hCa⟩ := lemma_3_7_a_coprime eps heps
  obtain ⟨Cc, hCc0, hCc⟩ := htFalt_le_of_condA eps heps
  obtain ⟨C₂, hCb⟩ := htFalt_le_of_condB
  refine ⟨max Ca Cc, lt_of_lt_of_le hCa0 (le_max_left _ _),
    {x : ℂ | faltingsHeightJ x ≤ max C₅ (|M + C₂| / 5 + 28 / 5 + 1.4 * C₅)},
    galoisFiniteJ_htFalt_le _, ?_⟩
  intro E l hl hl5 condA condB hcA hcB
  have hFeq : faltingsHeightJ E.j = htFaltOf E.toSSCurve.fld E.toSSCurve.W :=
    faltingsHeightJ_eq E.toSSCurve
  have hdeq : (E.deg : ℝ) = (Module.finrank ℚ E.toSSCurve.fld : ℝ) := rfl
  have hd1 : (1 : ℝ) ≤ (E.deg : ℝ) := by exact_mod_cast E.deg_pos
  have hp0 : (0 : ℝ) ≤ (E.deg : ℝ) ^ eps := Real.rpow_nonneg (by linarith) eps
  have hweak : ∀ C₀ : ℝ, C₀ ≤ max Ca Cc →
      (100 * (E.deg : ℝ) * (faltingsHeightJ E.j + max Ca Cc * (E.deg : ℝ) ^ eps) ≤ (l : ℝ)) →
      100 * (Module.finrank ℚ E.toSSCurve.fld : ℝ)
        * (htFaltOf E.toSSCurve.fld E.toSSCurve.W
            + C₀ * (Module.finrank ℚ E.toSSCurve.fld : ℝ) ^ eps) ≤ (l : ℝ) := by
    intro C₀ hle hA
    rw [← hdeq, ← hFeq]
    refine le_trans ?_ hA
    have hmul : C₀ * (E.deg : ℝ) ^ eps ≤ max Ca Cc * (E.deg : ℝ) ^ eps :=
      mul_le_mul_of_nonneg_right hle hp0
    have h100 : (0 : ℝ) ≤ 100 * (E.deg : ℝ) := by linarith
    exact mul_le_mul_of_nonneg_left (by linarith) h100
  have hA1 : condA → E.toSSCurve.PrimeToLocalHeights l := by
    intro hc
    rw [hcA] at hc
    intro p hp
    exact hCa E.toSSCurve.fld E.toSSCurve.W l p hl (hweak Ca (le_max_left _ _) hc.1) hp
  refine ⟨hA1, fun _ _ => E.multRed, ?_⟩
  intro hAB hcyc
  show faltingsHeightJ E.j ≤ max C₅ (|M + C₂| / 5 + 28 / 5 + 1.4 * C₅)
  rw [hFeq]
  rcases hAB with hc | hc
  · -- ☆条件 (a)——素性は `lemma_3_7_a_coprime` から
    have hdagE := hdag E l hl hl5 hcyc (hA1 hc)
    have hcA' := hcA.1 hc
    obtain ⟨p, hp⟩ := hcA'.2
    have hmain := hCc C₅ hC₅ E.toSSCurve.fld E.toSSCurve.W l p hl.one_lt.le hp
      (hweak Cc (le_max_right _ _) hcA'.1) hdagE
    exact le_trans hmain (le_max_left _ _)
  · -- ☆条件 (b)——素性は条件の定義そのもの
    have hcB' := hcB.1 hc
    have hdagE := hdag E l hl hl5 hcyc hcB'.2
    have harch : htArchJ E.toSSCurve.fld E.toSSCurve.W ≤ M := hM E.toSSCurve hcB'.1
    have hmain := hCb M C₅ hC₅ E.toSSCurve.fld E.toSSCurve.W E.toSSCurve.ss harch l
      hl.two_le hdagE
    exact le_trans hmain (le_max_right _ _)

/-! ## ★出典の紐付け(`.src`) -/

def lemma_3_7_stableLine_cop.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(安定直線の側——hdag が PrimeToLocalHeights も受け取る形)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_stableLine_cop.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_7_a_coprime(第 1088、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_7_a_coprime") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1210）**——第 1166 の `hdag` は" ++
       "`Lemma 3.5` が要る**素性**を持っていなかったので、そのままでは導けなかった。" ++
       "☆使う側では (a)(b) のどちらの枝でも `PrimeToLocalHeights` が手元にある" ++
       "——(a) では `lemma_3_7_a_coprime`、(b) では条件の定義そのもの——ので、" ++
       "結論は第 1166 とまったく同じである。" ++
       "★★★これで残るギャップ `hdag` は**実際に導ける形**になった。") 2 ]

end ABC3.Found.GenEll
