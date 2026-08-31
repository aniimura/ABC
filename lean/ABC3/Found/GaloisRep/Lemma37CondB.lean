/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Meta.Claim

/-!
# `Lemma 3.7` の第 3 の主張・**条件 (b) の側**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★これは何か

`Found/GaloisRep/Lemma37C.lean` は**条件 (a)**（`l` が `ht^Falt` に比べて大きい）の側で
`ht^Falt ≤ C′` を出した。本ファイルは**条件 (b)**（`[E_L] ∈ K_V`、
すなわちアルキメデス側の高さが有界）の側で同じことをする。

## ★★★機構

`K_V` が compactly bounded なら `h_∞(j) ≤ M` である。半安定なら

    12·ht^Falt ≤ h(j) + C₂ = h_fin(j) + h_∞(j) + C₂ ≤ deg∞ + M + C₂

（`prop_3_4_chain_semistable` の 3 本目と `htFinJ_le_degInfOf`）。
★これと `Lemma 3.5` の

    (†)  (l/14)·deg∞ ≤ ht^Falt + 2·log(l) + C′

を合わせると、`m ≔ M + C₂` として

    ht^Falt·(12l/14 − 1) ≤ (l/14)·m + 2·log(l) + C′

★★`l ≥ 2` では `12l/14 − 1 ≥ 5l/14` なので、両辺を `5l/14` で割って

    ht^Falt ≤ |m|/5 + 28/5 + 1.4·C′

——**`l` に依らない上界**が出る。☆`log(l) ≤ l − 1` を使う。

## ★★条件 (a) との違い

条件 (a) では `l` が**大きい**ことから上界が出た。条件 (b) では `l` は小さくてよく、
代わりに**アルキメデス側が有界**であることから出る。★どちらの場合も
`ht^Falt` が有界になり、`northcott`（`§9-1179`）で類が有限個になる。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve

/-! ## ★★★★★★数値の核 -/

/-- ★★★★★★**条件 (b) の数値の核**。

`12h ≤ D + m` と `(lR/14)·D ≤ h + 2·Lg + C′`（`Lg ≤ lR − 1`、`lR ≥ 2`、`0 ≤ D`）から

    h ≤ |m|/5 + 28/5 + 1.4·C′

を出す。★`l` に依らない上界であることが要点である。 -/
theorem htFalt_le_condB_numeric (h D m Lg C' lR : ℝ)
    (hl2 : 2 ≤ lR) (hC'0 : 0 ≤ C') (hLg0 : 0 ≤ Lg) (hLg : Lg ≤ lR - 1)
    (hD0 : 0 ≤ D) (hchain : 12 * h ≤ D + m)
    (hdag : (lR / 14) * D ≤ h + 2 * Lg + C') :
    h ≤ |m| / 5 + 28 / 5 + 1.4 * C' := by
  have hmabs : m ≤ |m| := le_abs_self m
  have hmabs0 : (0:ℝ) ≤ |m| := abs_nonneg m
  have hlpos : (0:ℝ) < lR := by linarith
  rcases le_or_gt h 0 with hh | hhpos
  · nlinarith
  rcases le_or_gt (12 * h - m) 0 with hcase | hcase
  · -- `12h ≤ m` なら直ちに有界
    nlinarith
  · -- `12h > m`——(†) に `D ≥ 12h − m` を代入する
    have hD : 12 * h - m ≤ D := by linarith
    have hstep : (lR / 14) * (12 * h - m) ≤ h + 2 * Lg + C' := by
      have hmul : (lR / 14) * (12 * h - m) ≤ (lR / 14) * D :=
        mul_le_mul_of_nonneg_left hD (by positivity)
      linarith
    -- `h(12lR/14 − 1) ≤ (lR/14)m + 2Lg + C′`
    have hcoef : (5 * lR / 14) ≤ 12 * lR / 14 - 1 := by linarith
    have hcm : (5 * lR / 14) * h ≤ (12 * lR / 14 - 1) * h :=
      mul_le_mul_of_nonneg_right hcoef hhpos.le
    have hmain : (5 * lR / 14) * h ≤ (lR / 14) * m + 2 * Lg + C' := by nlinarith
    have hLglR : Lg ≤ lR := by linarith
    have h5 : (0:ℝ) < 5 * lR / 14 := by positivity
    rw [← le_div_iff₀' h5] at hmain
    have hbound : ((lR / 14) * m + 2 * Lg + C') / (5 * lR / 14)
        ≤ |m| / 5 + 28 / 5 + 1.4 * C' := by
      rw [div_le_iff₀ h5]
      have h1 : (lR / 14) * m ≤ (lR / 14) * |m| :=
        mul_le_mul_of_nonneg_left hmabs (by positivity)
      nlinarith
    linarith

/-! ## ★★★★★★★★曲線の水準 -/

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★★★★★★★**条件 (b) ＋ (†) なら `ht^Falt` は有界**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

☆受けているのは `hdag`（＝ `Lemma 3.5` の結論 (†)）と、
`h_∞(j) ≤ M`（`[E_L] ∈ K_V` の内容）だけである。 -/
theorem htFalt_le_of_condB :
    ∃ C₂ : ℝ, ∀ (M C' : ℝ), 0 ≤ C' →
      ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic],
        (∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E) →
        htArchJ L E ≤ M →
        ∀ l : ℕ, 2 ≤ l →
        ((l : ℝ) / 14) * degInfOf L E ≤ htFaltOf L E + 2 * Real.log l + C' →
        htFaltOf L E ≤ |M + C₂| / 5 + 28 / 5 + 1.4 * C' := by
  obtain ⟨C, hC⟩ := prop_3_4_chain_semistable 1 one_pos
  refine ⟨C / 2, fun M C' hC'0 L _ _ E _ hss harch l hl2 hdag => ?_⟩
  have hlR : (2:ℝ) ≤ (l : ℝ) := by exact_mod_cast hl2
  have hl1 : (1:ℝ) ≤ (l : ℝ) := by linarith
  have hLg0 : (0:ℝ) ≤ Real.log l := Real.log_nonneg hl1
  have hLg : Real.log l ≤ (l : ℝ) - 1 :=
    Real.log_le_sub_one_of_pos (by linarith)
  have hD0 : (0:ℝ) ≤ degInfOf L E := degInfOf_nonneg E
  -- `12·ht^Falt ≤ deg∞ + h_∞(j) + C/2`
  have hchain : 12 * htFaltOf L E ≤ degInfOf L E + (M + C / 2) := by
    have h3 := (hC L E hss).2.2
    have hfin : htFinJ L E ≤ degInfOf L E := htFinJ_le_degInfOf E
    have hjeq : htJ L E = htFinJ L E + htArchJ L E := rfl
    have := harch
    rw [hjeq] at h3
    linarith
  exact htFalt_le_condB_numeric (htFaltOf L E) (degInfOf L E) (M + C / 2)
    (Real.log l) C' (l : ℝ) hlR hC'0 hLg0 hLg hD0 hchain hdag

/-! ## ★出典の紐付け(`.src`) -/

def htFalt_le_condB_numeric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 3 の主張・条件 (b)——数値の核。★無条件)",
    sectionId := "genell-lemma-3-7" }

def htFalt_le_of_condB.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 3 の主張・条件 (b)——ht^Falt ≤ C″。★(†) を受ける)",
    sectionId := "genell-lemma-3-7" }

def htFalt_le_of_condB.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "Lemma 3.5 の結論 (†)((l/14)·deg∞ ≤ ht^Falt + 2log l + C′)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hdag_of_velu") 3 ]

end ABC3.Found.GaloisRep
