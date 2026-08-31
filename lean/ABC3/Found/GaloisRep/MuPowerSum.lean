/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuCharSum
import Mathlib.Tactic.LinearCombination

/-!
# Galois (G6) 第 847 ブロック —— **★★★★★★★★★★`∑_ζ ζ^a/(1−ζ)^k` の機械**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★★★★★これは何か

`sum_mu_dyterm`（第 846 の葉）は `240·∑_{ζ≠1} ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1` であり、
`p₄ = ∑_{ζ≠1} 1/(1−ζ)⁴` が要る。★これまでの道具（`sum_mu_frac` 等）は
`p₁, p₂, p₃` までで、`p₄` には新しい種が要る。

★★**種は要らない**——次の 3 つの性質だけで**すべての `k` が決まる**:

    `A^{(k)}_a := ∑_{ζ≠1} ζ^a/(1−ζ)^k`

| 性質 | 内容 | 根拠 |
|---|---|---|
| (i) | `A^{(k+1)}_a − A^{(k+1)}_{a+1} = A^{(k)}_a` | `ζ^a·(1−ζ)/(1−ζ)^{k+1} = ζ^a/(1−ζ)^k` |
| (ii) | `∑_{a<l} A^{(k)}_a = 0` | ★★**`∑_{a<l} ζ^a = 0`**（`ζ ≠ 1`） |
| (iii) | `A^{(0)}_a = l·[l∣a] − 1` | 指標和 |

(i) を繰り返すと `A^{(k+1)}_a = A^{(k+1)}_0 − ∑_{j<a} A^{(k)}_j`。これを (ii) に入れると

    ★★★★★★`l · A^{(k+1)}_0 = ∑_{j<l} (l − 1 − j) · A^{(k)}_j`

——★**下の層だけから上の層の定数項が決まる漸化式**である。

## ★★★★★★★★★★★★★★★★★★★★★★★★重み付きの和で登る（次のブロックの計画）

一般の `a` での閉じた式を作ると `∑ j^m` が増えていく。★そこで**重み付きの和**

    `M_m^{(k)} := ∑_{i<l} C(l−1−i, m) · A^{(k)}_i`

を使う。ホッケースティック `∑_{t<l} C(t,m) = C(l,m+1)` と (i)(ii) から

| 式 | 内容 |
|---|---|
| `M_0^{(k)} = 0` | (ii) そのもの |
| ★★`M_m^{(k+1)} = C(l,m+1)·A^{(k+1)}_0 − M_{m+1}^{(k)}` | (i) を代入してホッケースティック |
| `M_m^{(0)} = l·C(l−1,m) − C(l,m+1)` | (iii) から |
| ★`l·A^{(k+1)}_0 = M_1^{(k)}` | `M_0^{(k+1)} = 0` の別の書き方 |

★★★**一般の `a` の式は一切要らない**。`p₄` までの鎖は有限で短い:

```
M_4^{(0)} → M_3^{(1)} → M_2^{(2)} → M_1^{(3)} → A^{(4)}_0
M_3^{(0)} → M_2^{(1)} → M_1^{(2)} → A^{(3)}_0
M_2^{(0)} → M_1^{(1)} → A^{(2)}_0
M_1^{(0)} → A^{(1)}_0
```

☆検算: `A^{(1)}_0 = [l·C(l−1,1) − C(l,2)]/l = (l−1) − (l−1)/2 = (l−1)/2` ✓
☆`A^{(2)}_0 = [C(l,2)(l−1)/2 − l·C(l−1,2) + C(l,3)]/l = (l−1)(5−l)/12` ✓

## ★★★検算（`l = 2, 3, 5, 7, 11, 13, 17` で数値確認、2026-08-31）

    `p₁ = (l−1)/2`
    `p₂ = (l−1)(5−l)/12`
    `p₃ = (l−1)(3−l)/8`
    `p₄ = (l−1)(l³ + l² − 109l + 251)/720`

そして `3p₄ − 7p₃ + 5p₂ − p₁ = (l⁴−1)/240` が厳密に成り立つ
（`x = 1/(1−ζ)` と置くと `ζ²(2+ζ)/(1−ζ)⁴ = 3x⁴ − 7x³ + 5x² − x`）。

## ★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `muPow` | `A^{(k)}_a = ∑_{ζ≠1} ζ^a/(1−ζ)^k` |
| `muPow_succ_sub` | ★(i) |
| `muPow_sum_range` | ★★(ii)——`∑_{a<l} ζ^a = 0` から |
| `muPow_zero_level` | (iii) |
| `muPow_shift` | `A^{(k+1)}_a = A^{(k+1)}_0 − ∑_{j<a}A^{(k)}_j` |
| `muPow_zero_rec` | ★★★★★★**漸化式** |
| `sum_range_sum_range_lt` | `∑_{a<n}∑_{j<a}f j = ∑_{j<n}(n−1−j)f j` |
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {F : Type*} [Field F]

/-! ## ★★★★★機械の定義 -/

/-- ★★★★★★**`A^{(k)}_a = ∑_{ζ≠1} ζ^a/(1−ζ)^k`**。 -/
noncomputable def muPow (l : ℕ) (ζ : F) (k a : ℕ) : F :=
  ∑ i ∈ (range l).erase 0, (ζ ^ i) ^ a * ((1 - ζ ^ i)⁻¹) ^ k

theorem one_sub_zeta_ne_zero {l : ℕ} {ζ : F} (hζ : IsPrimitiveRoot ζ l)
    {i : ℕ} (hi : i ∈ (range l).erase 0) : (1 : F) - ζ ^ i ≠ 0 := by
  have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
  have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
  intro h
  exact hζ.pow_ne_one_of_pos_of_lt hi0 hil (sub_eq_zero.1 h).symm

/-- ★★★**(i) 一段下げる**——`ζ^a·(1−ζ)/(1−ζ)^{k+1} = ζ^a/(1−ζ)^k`。 -/
theorem muPow_succ_sub {l : ℕ} {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k a : ℕ) :
    muPow l ζ (k + 1) a - muPow l ζ (k + 1) (a + 1) = muPow l ζ k a := by
  simp only [muPow, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have hinv : (1 - ζ ^ i) * ((1 - ζ ^ i)⁻¹) = 1 :=
    mul_inv_cancel₀ (one_sub_zeta_ne_zero hζ hi)
  linear_combination ((ζ ^ i) ^ a * ((1 - ζ ^ i)⁻¹) ^ k) * hinv

/-- ★★★★★★**(ii) `a` について一周すると消える**——`∑_{a<l} ζ^a = 0` から。 -/
theorem muPow_sum_range {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k : ℕ) :
    ∑ a ∈ range l, muPow l ζ k a = 0 := by
  simp only [muPow]
  rw [Finset.sum_comm]
  refine Finset.sum_eq_zero (fun i hi => ?_)
  obtain ⟨-, hsum⟩ := zeta_pow_facts hl hζ hi
  rw [← Finset.sum_mul, hsum, zero_mul]

/-- ★★(iii) 一番下の層。 -/
theorem muPow_zero_level {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (a : ℕ) :
    muPow l ζ 0 a = (if l ∣ a then (l : F) else 0) - 1 := by
  simp only [muPow, pow_zero, mul_one]
  rw [← sum_mu_pow_erase_zero hl hζ a]
  exact Finset.sum_congr rfl (fun i _ => by rw [← pow_mul])

/-- ★★★★**(i) を繰り返した形**。 -/
theorem muPow_shift {l : ℕ} {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k a : ℕ) :
    muPow l ζ (k + 1) a = muPow l ζ (k + 1) 0 - ∑ j ∈ range a, muPow l ζ k j := by
  induction a with
  | zero => simp
  | succ n ih =>
      have h := muPow_succ_sub hζ k n
      rw [Finset.sum_range_succ]
      linear_combination -h + ih

/-- ★`∑_{a<n}∑_{j<a} f j = ∑_{j<n} (n − 1 − j)·f j`。 -/
theorem sum_range_sum_range_lt (n : ℕ) (f : ℕ → F) :
    ∑ a ∈ range n, ∑ j ∈ range a, f j = ∑ j ∈ range n, ((n : F) - 1 - (j : F)) * f j := by
  induction n with
  | zero => simp
  | succ m ih =>
      rw [Finset.sum_range_succ, ih, Finset.sum_range_succ]
      push_cast
      rw [show ((m : F) + 1 - 1 - (m : F)) = 0 by ring, zero_mul, add_zero,
        ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl (fun j _ => ?_)
      ring

/-- ★★★★★★★★★★**漸化式**——下の層だけから上の層の定数項が決まる。

    `l · A^{(k+1)}_0 = ∑_{j<l} (l − 1 − j) · A^{(k)}_j` -/
theorem muPow_zero_rec {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k : ℕ) :
    (l : F) * muPow l ζ (k + 1) 0
      = ∑ j ∈ range l, ((l : F) - 1 - (j : F)) * muPow l ζ k j := by
  have h0 := muPow_sum_range hl hζ (k + 1)
  rw [Finset.sum_congr rfl (fun a _ => muPow_shift hζ k a)] at h0
  rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range, nsmul_eq_mul,
    sum_range_sum_range_lt] at h0
  linear_combination h0

/-! ## ★出典の紐付け(`.src`) -/

def muPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(A^{(k)}_a = ∑_{ζ≠1} ζ^a/(1−ζ)^k。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muPow_succ_sub.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(一段下げる。★ζ が原始 l 乗根)",
    sectionId := "genell-lemma-3-2" }

def muPow_sum_range.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(a について一周すると消える。★l は素数)",
    sectionId := "genell-lemma-3-2" }

def muPow_zero_rec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(層をひとつ上げる漸化式。★l は素数)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
