import ABC3.Found.GaloisRep.CollinearGroupLaw

/-!
# Galois (G6) 第 258 ブロック —— **★★★★★★葉 (d) の主要項**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★葉 (d)(単射性)の設計

第 257 で「`x` 座標が相異なること」が群法則に要ると分かった。それは
**Tate 座標が単射であること**(葉 (d))から来る。

★★★単射性の骨格は**逐次近似**である:

    a − a' ∈ I^j  かつ  w − w' ∈ I^j
      ⟹ (X, Y が一致するなら) a − a' ∈ I^{j+1}

これを `j` について繰り返すと `a − a' ∈ ⋂ I^n = 0`(`IsHausdorff`)。
★本ブロックはその **1 段分の材料**を揃える。

## ★★★主要項

`f(t) = t/(1−t)² = t + 2t² + …`、`g(t) = t²/(1−t)³ = t² + …` なので

| 補題 | 内容 |
|---|---|
| `tateXterm_sub_self_mem` | `f(t) − t ∈ I^{2k}`(`t ∈ I^k`) |
| `tateYterm_mem_two_mul` | `g(t) ∈ I^{2k}` |
| `tateXpair_sub_add_mem` | ★★★★★★**`X(a,w) ≡ a + w`(2 位まで)** |
| `tateYpair_add_mem` | ★★★★★**`Y(a,w) ≡ −w`(2 位まで)** |

★`(1−t)²·f(t) = t` を使うと `f(t) − t = inv(1−t)²·t²(2−t)` と**明示的に**書けるので、
冪級数の係数を触らずに済む。

## ★★★★★差の形

逐次近似には「差」の形が要る:

    f(a) − f(b) − (a − b) ∈ I^{j+1}      (`a − b ∈ I^j`、`a, b ∈ I`)

★★これも明示式で出る。`D := (1−a)²(1−b)²` として

    D·(f(a) − f(b)) = (a−b)(1 − ab)
    ⟹ f(a) − f(b) − (a−b) = (a−b)·inv(D)·((1−ab) − D)

で、`(1−ab) − D = −ab + 2s − s²`(`s := a+b−ab ∈ I`)は `I` の元である。

## ★★配管 —— adic 和の冪

`adicSum_mem` は `∈ I` しか言わない。`∈ I^k` が要るので `adicSum_mem_pow` を足した
(部分和 `partialSum a k` が `I^k` に入り、残差も `I^k` に入る)。
`evalAdic_mem_pow` も同じ形である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `adicSum_mem_pow`・`evalAdic_mem_pow` | ★★★★adic 和の冪への所属 |
| `tateXtail_mem_pow`・`tateYtail_mem_pow` | ★★★尾の冪 |
| `tateXpair_sub_add_mem`・`tateYpair_add_mem` | ★★★★★★**主要項** |
| `tateXterm_diff_mem` | ★★★★★**差の形** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★adic 和の冪への所属 -/

/-- ★★★★**adic 和の冪への所属**——各項が `I^k` に入るなら和も入る。 -/
theorem adicSum_mem_pow [IsPrecomplete I R] (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) (k : ℕ)
    (hk : ∀ n, a n ∈ I ^ k) : adicSum a ha ∈ I ^ k := by
  have hspec := (smodEq_iff_sub_mem (I := I) (n := k) (partialSum a k) (adicSum a ha)).1
    (adicSum_spec a ha k)
  have hpart : partialSum a k ∈ I ^ k := by
    rw [partialSum]
    exact Ideal.sum_mem _ fun n _ => hk n
  have h2 : adicSum a ha - partialSum a k ∈ I ^ k := by
    have h3 := neg_mem hspec
    rwa [neg_sub] at h3
  have h4 := Ideal.add_mem (I ^ k) h2 hpart
  simpa using h4

theorem evalAdic_mem_pow [IsAdicComplete I R] (f : PowerSeries ℤ)
    (hf : PowerSeries.coeff 0 f = 0) (q : R) (hq : q ∈ I) (k : ℕ) (hqk : q ∈ I ^ k) :
    evalAdic f q hq ∈ I ^ k := by
  have h := evalAdic_spec (I := I) f q hq k
  have hpart : partialEval f q k ∈ I ^ k := by
    rw [partialEval]
    refine Ideal.sum_mem _ fun n _ => ?_
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp [hf]
    · refine Ideal.mul_mem_left _ _ ?_
      have hq1 : q ^ n = q ^ (n - 1) * q := by
        rw [← pow_succ]
        congr 1
        omega
      rw [hq1]
      exact Ideal.mul_mem_left _ _ hqk
  have hsub : evalAdic f q hq - partialEval f q k ∈ I ^ k := by
    have h2 := (smodEq_iff_sub_mem (I := I) (n := k) (partialEval f q k) (evalAdic f q hq)).1 h
    have h3 := neg_mem h2
    rwa [neg_sub] at h3
  have h4 := Ideal.add_mem (I ^ k) hsub hpart
  simpa using h4

theorem tateXtail_mem_pow [IsAdicComplete I R] (u q : R) (hq : q ∈ I) (k : ℕ)
    (hqk : q ∈ I ^ k) : tateXtail u q hq ∈ I ^ k := by
  refine adicSum_mem_pow _ _ k fun n => ?_
  refine tateXterm_mem_pow ?_
  have h : q ^ (n + 1) * u = q ^ n * (q * u) := by ring
  rw [h]
  exact Ideal.mul_mem_left _ _ (Ideal.mul_mem_right _ _ hqk)

theorem tateYtail_mem_pow [IsAdicComplete I R] (u q : R) (hq : q ∈ I) (k : ℕ)
    (hqk : q ∈ I ^ k) : tateYtail u q hq ∈ I ^ k := by
  refine adicSum_mem_pow _ _ k fun n => ?_
  refine tateYterm_mem_pow ?_
  have h : q ^ (n + 1) * u = q ^ n * (q * u) := by ring
  rw [h]
  exact Ideal.mul_mem_left _ _ (Ideal.mul_mem_right _ _ hqk)

/-! ## ★★★主要項 -/

theorem sq_mem_two_mul {k : ℕ} {t : R} (ht : t ∈ I ^ k) : t ^ 2 ∈ I ^ (2 * k) := by
  rw [two_mul, pow_add, sq]
  exact Ideal.mul_mem_mul ht ht

/-- ★★★**`f(t) − t` は `t²` の位**——`f(t) = t + 2t² + …`。

★`(1−t)²·f(t) = t` から `f(t) − t = inv(1−t)²·t²(2−t)` と明示的に書ける。 -/
theorem tateXterm_sub_self_mem [IsAdicComplete I R] {k : ℕ} {t : R} (ht : t ∈ I ^ k)
    (ht1 : t ∈ I) : tateXterm t - t ∈ I ^ (2 * k) := by
  have hu : IsUnit (1 - t) := isUnit_one_sub ht1
  have hkey : tateXterm t - t = Ring.inverse (1 - t) ^ 2 * (t ^ 2 * (2 - t)) := by
    have h1 : (1 - t) ^ 2 * (tateXterm t - t) = t ^ 2 * (2 - t) := by
      rw [mul_sub, mul_tateXterm ht1]
      ring
    calc tateXterm t - t
        = Ring.inverse (1 - t) ^ 2 * ((1 - t) ^ 2 * (tateXterm t - t)) := by
          rw [← mul_assoc, ← mul_pow, Ring.inverse_mul_cancel _ hu, one_pow, one_mul]
      _ = Ring.inverse (1 - t) ^ 2 * (t ^ 2 * (2 - t)) := by rw [h1]
  rw [hkey]
  exact Ideal.mul_mem_left _ _ (Ideal.mul_mem_right _ _ (sq_mem_two_mul ht))

/-- ★★★**`g(t)` は `t²` の位**。 -/
theorem tateYterm_mem_two_mul {k : ℕ} {t : R} (ht : t ∈ I ^ k) : tateYterm t ∈ I ^ (2 * k) := by
  rw [tateYterm]
  exact Ideal.mul_mem_right _ _ (sq_mem_two_mul ht)

/-- ★★★★★★**`X(a,w) ≡ a + w`(2 位まで)**——葉 (d) の主要項。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateXpair_sub_add_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (k : ℕ)
    (ha : a ∈ I ^ k) (hw : w ∈ I ^ k) (ha1 : a ∈ I) (hw1 : w ∈ I) (haw : a * w = q) :
    tateXpair a w q hq - (a + w) ∈ I ^ (2 * k) := by
  have hq2k : q ∈ I ^ (2 * k) := by
    rw [← haw, two_mul, pow_add]
    exact Ideal.mul_mem_mul ha hw
  have hA := tateXterm_sub_self_mem (I := I) ha ha1
  have hW := tateXterm_sub_self_mem (I := I) hw hw1
  have hTA := tateXtail_mem_pow (I := I) a q hq (2 * k) hq2k
  have hTW := tateXtail_mem_pow (I := I) w q hq (2 * k) hq2k
  have hS := evalAdic_mem_pow (I := I) (sigmaSeries 1) (by simp [sigmaSeries]) q hq (2 * k) hq2k
  have hkey : tateXpair a w q hq - (a + w)
      = (tateXterm a - a) + (tateXterm w - w) + tateXtail a q hq + tateXtail w q hq
        + (-2) * evalAdic (sigmaSeries 1) q hq := by
    rw [tateXpair]
    ring
  rw [hkey]
  exact Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _ hA hW) hTA) hTW)
    (Ideal.mul_mem_left _ _ hS)

/-- ★★★★★**`Y(a,w) ≡ −w`(2 位まで)**——`g(a)`・`g(w)` は 2 位である。 -/
theorem tateYpair_add_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (k : ℕ)
    (ha : a ∈ I ^ k) (hw : w ∈ I ^ k) (hw1 : w ∈ I) (haw : a * w = q) :
    tateYpair a w q hq + w ∈ I ^ (2 * k) := by
  have hq2k : q ∈ I ^ (2 * k) := by
    rw [← haw, two_mul, pow_add]
    exact Ideal.mul_mem_mul ha hw
  have hGA := tateYterm_mem_two_mul (I := I) ha
  have hGW := tateYterm_mem_two_mul (I := I) hw
  have hW := tateXterm_sub_self_mem (I := I) hw hw1
  have hTA := tateYtail_mem_pow (I := I) a q hq (2 * k) hq2k
  have hTX := tateXtail_mem_pow (I := I) w q hq (2 * k) hq2k
  have hTW := tateYtail_mem_pow (I := I) w q hq (2 * k) hq2k
  have hS := evalAdic_mem_pow (I := I) (sigmaSeries 1) (by simp [sigmaSeries]) q hq (2 * k) hq2k
  have hkey : tateYpair a w q hq + w
      = tateYterm a + tateYtail a q hq + (-1) * (tateXterm w - w) + (-1) * tateXtail w q hq
        + (-1) * tateYterm w + (-1) * tateYtail w q hq + evalAdic (sigmaSeries 1) q hq := by
    rw [tateYpair]
    ring
  rw [hkey]
  refine Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _
    (Ideal.add_mem _ (Ideal.add_mem _ hGA hTA) (Ideal.mul_mem_left _ _ hW))
    (Ideal.mul_mem_left _ _ hTX)) (Ideal.mul_mem_left _ _ hGW))
    (Ideal.mul_mem_left _ _ hTW)) hS

/-! ## ★★★★★差の形 -/

/-- ★★★★★**差の形**——`f(a) − f(b) − (a − b)` は 1 つ位が上がる。

★`D := (1−a)²(1−b)²` として `D·(f(a) − f(b)) = (a−b)(1 − ab)` なので

    f(a) − f(b) − (a−b) = (a−b)·inv(D)·((1−ab) − D)

であり、`(1−ab) − D = −ab + 2s − s²`(`s := a+b−ab`)は `I` の元である。 -/
theorem tateXterm_diff_mem [IsAdicComplete I R] {j : ℕ} {a b : R} (ha : a ∈ I) (hb : b ∈ I)
    (hab : a - b ∈ I ^ j) : tateXterm a - tateXterm b - (a - b) ∈ I ^ (j + 1) := by
  have hua : IsUnit (1 - a) := isUnit_one_sub ha
  have hub : IsUnit (1 - b) := isUnit_one_sub hb
  set D : R := (1 - a) ^ 2 * (1 - b) ^ 2 with hDdef
  have huD : IsUnit D := by
    simp only [hDdef]
    exact (hua.pow 2).mul (hub.pow 2)
  have hda : (1 - a) ^ 2 * tateXterm a = a := mul_tateXterm ha
  have hdb : (1 - b) ^ 2 * tateXterm b = b := mul_tateXterm hb
  have hD : D * (tateXterm a - tateXterm b) = (a - b) * (1 - a * b) := by
    simp only [hDdef]
    linear_combination (1 - b) ^ 2 * hda - (1 - a) ^ 2 * hdb
  have hE : (1 - a * b) - D ∈ I := by
    have hs : a + b - a * b ∈ I :=
      Ideal.sub_mem _ (Ideal.add_mem _ ha hb) (Ideal.mul_mem_right _ _ ha)
    have hEeq : (1 - a * b) - D
        = -(a * b) + 2 * (a + b - a * b) - (a + b - a * b) ^ 2 := by
      simp only [hDdef]; ring
    rw [hEeq]
    exact Ideal.sub_mem _ (Ideal.add_mem _ (neg_mem (Ideal.mul_mem_right _ _ ha))
      (Ideal.mul_mem_left _ _ hs)) (by rw [sq]; exact Ideal.mul_mem_right _ _ hs)
  have h1 : Ring.inverse D * D = 1 := Ring.inverse_mul_cancel _ huD
  have hkey : tateXterm a - tateXterm b - (a - b)
      = (a - b) * (Ring.inverse D * ((1 - a * b) - D)) := by
    calc tateXterm a - tateXterm b - (a - b)
        = Ring.inverse D * (D * (tateXterm a - tateXterm b)) - (a - b) := by
          rw [← mul_assoc, h1, one_mul]
      _ = Ring.inverse D * ((a - b) * (1 - a * b)) - (a - b) := by rw [hD]
      _ = (a - b) * (Ring.inverse D * ((1 - a * b) - D)) := by
          linear_combination (a - b) * h1
  rw [hkey, pow_succ]
  exact Ideal.mul_mem_mul hab (Ideal.mul_mem_left _ _ hE)

/-! ## ★出典の紐付け(`.src`) -/

def tateXpair_sub_add_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——X の主要項)",
    sectionId := "genell-def-3-3" }

def tateXterm_diff_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——f の差は位が上がる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
