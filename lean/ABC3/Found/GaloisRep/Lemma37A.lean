/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtFaltBounds
import Mathlib.Analysis.Complex.ExponentialBounds
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★`Lemma 3.7` の**第 1 の主張**（`Found`、無条件）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

## ★★★★★★★★★★★★★★★★これは何か

`Lemma 3.7` は 3 つの主張を持つ:

| 主張 | 内容 | 状態 |
|---|---|---|
| **第 1** | 条件 (a) なら **`l` は乗法還元の全素点の局所高さより大きい** | ★★★**本ファイルで証明**（無条件） |
| 第 2 | 条件 (b) かつ `[E] ∉ Exc` なら乗法還元の素点を持つ | ☆例外集合が要る |
| 第 3 | (a) か (b) かつ l-cyclic なら `[E] ∈ Exc` | ☆`Lemma 3.5` を経由する |

★★**第 1 の主張は原文の証明の最初の 2 段だけで出る**——そしてその 2 段は
`§9-1010`（`Found/GaloisRep/HtFaltBounds.lean`）で**どちらも証明済み**である。

## ★原文の 2 段（どちらも済）

1. 「if `v` is any local height of `E_L`, then `d · deg∞([E_L]) ≥ v · log(2)`」
   ——★`minDeltaExp_log_two_le`（`§9-1010`）。
2. 「by Proposition 3.4 … `ht^Falt + C·d^ϵ ≥ ht^Falt + C ≥ (1/14)·deg∞`」
   ——★`exists_degInfOf_le_htFalt`（`§9-1010`）。★★こちらは **`14` ではなく `12`** で取れている。

★★★あとは原文どおり:

    `l ≥ 100d·(ht^Falt + C·d^ϵ) ≥ (100d/12)·deg∞ ≥ (100·log 2/12)·v > v`

（原文は `14`。`12` なら余裕はさらに大きい。）

## ★定数 `C` の取り方

★`ht^Falt` は**下に有界**（`exists_htFalt_bddBelow`、`§9-1010`）なので、
`ht^Falt < 0` の場合の超過分は `|B|` で吸収できる。
★★`C ≔ |A| + |B| + 1` で足りる（`A` は `deg∞ ≤ 12ht^Falt + A` の定数）。
★★★原文の「for an appropriate choice of `C`」がこれである。

## ☆残り

☆第 2・第 3 の主張には**例外集合**（`Exc`）が要り、それは `Lemma 3.5` を経由する。
`Lemma 3.5` に残る唯一の入力は Faltings 高さの同種写像評価である
（`Found/GaloisRep/Lemma35Concrete.lean` に型で固定してある）。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★★★★★★★★★第 1 の主張 -/

/-- ★★★★★★★★★★★★★★★★**[GenEll] Lemma 3.7 の第 1 の主張**——★**無条件**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

条件 (a) の不等式 `l ≥ 100d·(ht^Falt + C·d^ϵ)` から
**`l` は局所高さ `v_p(Δ_min)` より大きい**。

★`C` は `ϵ` にすら依らない普遍定数である（`A`・`B` は `§9-1010` の普遍定数）。 -/
theorem lemma_3_7_a (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
      [E.IsElliptic] (l : ℕ) (p : HeightOneSpectrum (𝓞 L)),
      100 * (Module.finrank ℚ L : ℝ)
          * (htFaltOf L E + C * (Module.finrank ℚ L : ℝ) ^ eps) ≤ (l : ℝ) →
      (minDeltaExp p E : ℝ) < (l : ℝ) := by
  obtain ⟨A, hA0, hA⟩ := exists_degInfOf_le_htFalt
  obtain ⟨B, hB⟩ := exists_htFalt_bddBelow
  have hLlo : (0.69 : ℝ) ≤ Real.log 2 := by linarith [Real.log_two_gt_d9]
  have hLhi : Real.log 2 ≤ (0.70 : ℝ) := by linarith [Real.log_two_lt_d9]
  refine ⟨|A| + |B| + 1, by positivity, ?_⟩
  intro L _ _ E _ l p hcond
  set C : ℝ := |A| + |B| + 1 with hC
  set d : ℝ := (Module.finrank ℚ L : ℝ) with hd
  set F : ℝ := htFaltOf L E with hF
  set D : ℝ := degInfOf L E with hD
  set v : ℝ := (minDeltaExp p E : ℝ) with hv
  have hd1 : (1:ℝ) ≤ d := by
    rw [hd]; exact_mod_cast Module.finrank_pos
  have hdpos : (0:ℝ) < d := by linarith
  have hP1 : (1:ℝ) ≤ d ^ eps := Real.one_le_rpow hd1 heps.le
  have hΔ : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  have hvL : v * Real.log 2 ≤ d * D := minDeltaExp_log_two_le E hΔ p
  have hDF : D ≤ 12 * F + A := hA L E
  have hBF : B ≤ F := hB L E
  have hvnn : (0:ℝ) ≤ v := by
    rw [hv]; exact_mod_cast minDeltaExp_nonneg p E
  have hCpos : (0:ℝ) < C := by rw [hC]; positivity
  have hAle : A ≤ |A| := le_abs_self A
  have hBge : -|B| ≤ B := neg_abs_le B
  have hAnn : (0:ℝ) ≤ |A| := abs_nonneg A
  have hBnn : (0:ℝ) ≤ |B| := abs_nonneg B
  have hlow : 100 * d * F + 100 * d * C ≤ (l : ℝ) := by
    refine le_trans ?_ hcond
    nlinarith [mul_le_mul_of_nonneg_left hP1 (by positivity : (0:ℝ) ≤ 100 * d * C)]
  have hvle : v * Real.log 2 ≤ d * (12 * F + A) := by nlinarith
  nlinarith [hvle, hlow, hvnn, hd1, hCpos, hAle, hBge, hLlo, hLhi, hBF]

/-- ★★★★★★★★★★**`l` は局所高さと素である**——第 1 の主張の系。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★「`l` is prime to the local heights」は原文が「`l` is > these local heights」から
出している。★★局所高さは**正**（悪い還元の素点）なので `hne` を置く。 -/
theorem lemma_3_7_a_coprime (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
      [E.IsElliptic] (l : ℕ) (p : HeightOneSpectrum (𝓞 L)), Nat.Prime l →
      100 * (Module.finrank ℚ L : ℝ)
          * (htFaltOf L E + C * (Module.finrank ℚ L : ℝ) ^ eps) ≤ (l : ℝ) →
      minDeltaExp p E ≠ 0 →
      Nat.Coprime l (minDeltaExp p E).toNat := by
  obtain ⟨C, hCpos, hC⟩ := lemma_3_7_a eps heps
  refine ⟨C, hCpos, fun L _ _ E _ l p hl hcond hne => ?_⟩
  have hlt : (minDeltaExp p E : ℝ) < (l : ℝ) := hC L E l p hcond
  have hltZ : minDeltaExp p E < (l : ℤ) := by exact_mod_cast hlt
  have hnn := minDeltaExp_nonneg p E
  have hpos : 0 < (minDeltaExp p E).toNat := by omega
  rw [Nat.Prime.coprime_iff_not_dvd hl]
  intro hdvd
  have := Nat.le_of_dvd hpos hdvd
  omega

/-! ## ★出典の紐付け(`.src`) -/

def lemma_3_7_a.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 1 の主張——条件 (a) なら l は局所高さより大きい。★無条件)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_a_coprime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 1 の主張の系——l は局所高さと素。★無条件)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_a.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_log_two_le(d·deg∞ ≥ v·log 2、§9-1010)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_log_two_le") 3,
    .citation "[ABC3]" "exists_degInfOf_le_htFalt(deg∞ ≤ 12·ht^Falt + A、§9-1010)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_degInfOf_le_htFalt") 3,
    .citation "[ABC3]" "exists_htFalt_bddBelow(ht^Falt は下に有界、§9-1010)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_htFalt_bddBelow") 2,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29、第 566): Lemma 3.7 の**第 1 の主張が無条件で取れた**。" ++
       "★原文の証明の最初の 2 段(d·deg∞ ≥ v·log 2 と deg∞ ≤ 14(ht^Falt + C))は" ++
       "どちらも §9-1010 で証明済みであり、しかも後者は **12** で取れている" ++
       "(原文は 14 なので余裕はさらに大きい)。" ++
       "★★定数 C = |A| + |B| + 1 が原文の「for an appropriate choice of C」である" ++
       "——ht^Falt < 0 の場合の超過分を |B| で吸収する") 8,
    .folklore
      ("☆第 2・第 3 の主張には**例外集合 Exc**(Galois-finite)が要り、" ++
       "それは Lemma 3.5 を経由する。★Lemma 3.5 に残る唯一の入力は " ++
       "Faltings 高さの同種写像評価であり、" ++
       "Found/GaloisRep/Lemma35Concrete.lean に型で固定してある。" ++
       "★★条件 (b) の側にはさらに compactly bounded subset の具体化が要る") 9 ]

end ABC3.Found.GaloisRep
