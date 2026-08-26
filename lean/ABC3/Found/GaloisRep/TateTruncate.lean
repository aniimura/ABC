import ABC3.Found.GaloisRep.TateFunctorial

/-!
# Galois (G6) 第 222 ブロック —— **★★★★★★★★段 6 の骨格(葉 (b) を有限の主張に落とす)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★葉 (b) を有限の主張の族に落とす

葉 (b) は `tateDefect a w q = 0`(Weierstrass 方程式)である。`R` は `I` 進完備
(したがって `IsHausdorff`)なので、**`∀ n, tateDefect ∈ I^n` を示せば足りる**。

★`I^n` を法とすると `adicSum` は**部分和**に、`evalAdic` は `partialEval` に化ける。
つまり `tateDefect` は `n` 次の**切り詰め** `tateDefectTrunc n` と `I^n` を法として等しい:

    tateDefect a w q hq − tateDefectTrunc n a w q ∈ I^n        (`tateDefect_sub_trunc`)

★★`tateDefectTrunc n a w q` は `Ring.inverse` を除けば **`(a, w, q)` の多項式**である
——無限和も極限も入っていない。

## ★★★★★★★★そして特殊化の規準が立つ

`tateDefectTrunc` は**完備性を要求しない**ので、任意の環準同型で移る
(`map_tateDefectTrunc`)。したがって:

> 万有な環 `S` の元 `A, W` について `(AW)^n ∣ tateDefectTrunc n A W (AW)` が示せれば、
> `φ A = a`、`φ W = w` なる任意の特殊化で `tateDefect a w q = 0` が従う。

これが `tateDefect_eq_zero_of_specialize` である。★★**段 6 の骨格はこれで立った**——
残るのは「万有な環 `ℤ[A,W]`(`(1−A)(1−W)` 等で局所化)でその整除性を示す」ことだけになる。

## ★★方程式の差は対称である

`X(w,a) = X(a,w)`、`Y(w,a) = −X−Y` なので、差は `(a,w)` の交換で不変である
(`tateDefect_symm`)。★万有な環での議論で `A` と `W` の役割を入れ替えられる
——これは係数の消滅を両側から取るのに要る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateDefect_symm` | ★★★★差は `(a,w)` の交換で不変 |
| `eq_zero_of_mem_pow` | ★すべての `I^n` に属する元は 0 |
| `tateXtrunc` / `tateYtrunc` / `tateDefectTrunc` | ★★`n` 次の切り詰め |
| `tateXpair_sub_trunc` / `tateYpair_sub_trunc` | ★★★級数と切り詰めの差は `I^n` |
| `tateDefect_sub_trunc` | ★★★★★**差と切り詰めの差は `I^n`** |
| `tateDefect_eq_zero_of_trunc` | ★★★★★★★**葉 (b) は有限の主張の族に落ちる** |
| `map_tateDefectTrunc` | ★★★★★★切り詰めは任意の環準同型で移る |
| `tateDefect_eq_zero_of_specialize` | ★★★★★★★★**特殊化の規準** |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★差の対称性 -/

/-- ★★★★**方程式の差は `(a, w)` の交換で不変**——`X` は不変、`Y ↦ −X−Y` だから。 -/
theorem tateDefect_symm [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) :
    tateDefect w a q hq = tateDefect a w q hq := by
  rw [tateDefect, tateDefect, tateXpair_symm a w q hq, tateYpair_swap a w q hq]
  ring

/-! ## ★★`n` 次の切り詰め -/

/-- ★`I` 進的に 0 に近い元は 0(`IsHausdorff`)。 -/
theorem eq_zero_of_mem_pow [IsAdicComplete I R] {x : R} (h : ∀ n, x ∈ I ^ n) : x = 0 := by
  refine (IsHausdorff.eq_iff_smodEq (I := I)).2 fun n => ?_
  rw [smodEq_iff_sub_mem, sub_zero]
  exact h n

/-- ★★**`X(u,q)` の `n` 次の切り詰め**——`adicSum` を部分和に、`evalAdic` を `partialEval` に。 -/
noncomputable def tateXtrunc (n : ℕ) (a w q : R) : R :=
  (tateXterm a + partialSum (fun m => tateXterm (q ^ (m + 1) * a)) n)
    + (tateXterm w + partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
    - 2 * partialEval (sigmaSeries 1) q n

/-- ★★**`Y(u,q)` の `n` 次の切り詰め**。 -/
noncomputable def tateYtrunc (n : ℕ) (a w q : R) : R :=
  (tateYterm a + partialSum (fun m => tateYterm (q ^ (m + 1) * a)) n)
    - (tateXterm w + partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
    - (tateYterm w + partialSum (fun m => tateYterm (q ^ (m + 1) * w)) n)
    + partialEval (sigmaSeries 1) q n

theorem tateXtail_sub_partialSum [IsAdicComplete I R] (n : ℕ) (u q : R) (hq : q ∈ I) :
    tateXtail u q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * u)) n ∈ I ^ n := by
  have h := (smodEq_iff_sub_mem (I := I) _ _).1
    (adicSum_spec (fun m => tateXterm (q ^ (m + 1) * u)) (tateXtail_aux hq) n)
  have h2 := neg_mem h
  rw [neg_sub] at h2
  exact h2

theorem tateYtail_sub_partialSum [IsAdicComplete I R] (n : ℕ) (u q : R) (hq : q ∈ I) :
    tateYtail u q hq - partialSum (fun m => tateYterm (q ^ (m + 1) * u)) n ∈ I ^ n := by
  have h := (smodEq_iff_sub_mem (I := I) _ _).1
    (adicSum_spec (fun m => tateYterm (q ^ (m + 1) * u)) (tateYtail_aux hq) n)
  have h2 := neg_mem h
  rw [neg_sub] at h2
  exact h2

theorem evalAdic_sub_partialEval [IsAdicComplete I R] (n : ℕ) (f : PowerSeries ℤ) (q : R)
    (hq : q ∈ I) : evalAdic f q hq - partialEval f q n ∈ I ^ n := by
  have h := (smodEq_iff_sub_mem (I := I) _ _).1 (evalAdic_spec f q hq n)
  have h2 := neg_mem h
  rw [neg_sub] at h2
  exact h2

theorem tateXpair_sub_trunc [IsAdicComplete I R] (n : ℕ) (a w q : R) (hq : q ∈ I) :
    tateXpair a w q hq - tateXtrunc n a w q ∈ I ^ n := by
  have h1 := tateXtail_sub_partialSum (I := I) n a q hq
  have h2 := tateXtail_sub_partialSum (I := I) n w q hq
  have h3 := evalAdic_sub_partialEval (I := I) n (sigmaSeries 1) q hq
  have hkey : tateXpair a w q hq - tateXtrunc n a w q
      = (tateXtail a q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * a)) n)
        + (tateXtail w q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
        + (-2) * (evalAdic (sigmaSeries 1) q hq - partialEval (sigmaSeries 1) q n) := by
    rw [tateXpair, tateXtrunc]
    ring
  rw [hkey]
  exact Ideal.add_mem _ (Ideal.add_mem _ h1 h2) (Ideal.mul_mem_left _ _ h3)

theorem tateYpair_sub_trunc [IsAdicComplete I R] (n : ℕ) (a w q : R) (hq : q ∈ I) :
    tateYpair a w q hq - tateYtrunc n a w q ∈ I ^ n := by
  have h1 := tateYtail_sub_partialSum (I := I) n a q hq
  have h2 := tateYtail_sub_partialSum (I := I) n w q hq
  have h2x := tateXtail_sub_partialSum (I := I) n w q hq
  have h3 := evalAdic_sub_partialEval (I := I) n (sigmaSeries 1) q hq
  have hkey : tateYpair a w q hq - tateYtrunc n a w q
      = (tateYtail a q hq - partialSum (fun m => tateYterm (q ^ (m + 1) * a)) n)
        + (-1) * (tateXtail w q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
        + (-1) * (tateYtail w q hq - partialSum (fun m => tateYterm (q ^ (m + 1) * w)) n)
        + (evalAdic (sigmaSeries 1) q hq - partialEval (sigmaSeries 1) q n) := by
    rw [tateYpair, tateYtrunc]
    ring
  rw [hkey]
  exact Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _ h1
    (Ideal.mul_mem_left _ _ h2x)) (Ideal.mul_mem_left _ _ h2)) h3

/-! ## ★★★★★★★葉 (b) を有限の主張の族に落とす -/

/-- ★★★**方程式の差の `n` 次の切り詰め**——`Ring.inverse` を除けば `(a,w,q)` の多項式。 -/
noncomputable def tateDefectTrunc (n : ℕ) (a w q : R) : R :=
  tateYtrunc n a w q ^ 2 + tateXtrunc n a w q * tateYtrunc n a w q
    - (tateXtrunc n a w q ^ 3 + partialEval tateA4 q n * tateXtrunc n a w q
      + partialEval tateA6 q n)

/-- ★★★★★**差と切り詰めの差は `I^n` の元**。

★ここは `Ring.inverse` を含まない多項式の関係なので、`Ideal.Quotient` を使ってよい
(第 212・221 ブロックで `Ring.inverse` を避けたのとは事情が違う)。 -/
theorem tateDefect_sub_trunc [IsAdicComplete I R] (n : ℕ) (a w q : R) (hq : q ∈ I) :
    tateDefect a w q hq - tateDefectTrunc n a w q ∈ I ^ n := by
  refine Ideal.Quotient.eq.1 ?_
  have hX : (Ideal.Quotient.mk (I ^ n)) (tateXpair a w q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateXtrunc n a w q) :=
    Ideal.Quotient.eq.2 (tateXpair_sub_trunc n a w q hq)
  have hY : (Ideal.Quotient.mk (I ^ n)) (tateYpair a w q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateYtrunc n a w q) :=
    Ideal.Quotient.eq.2 (tateYpair_sub_trunc n a w q hq)
  have h4 : (Ideal.Quotient.mk (I ^ n)) ((tateCurveAt q hq).a₄)
      = (Ideal.Quotient.mk (I ^ n)) (partialEval tateA4 q n) :=
    Ideal.Quotient.eq.2 (evalAdic_sub_partialEval n tateA4 q hq)
  have h6 : (Ideal.Quotient.mk (I ^ n)) ((tateCurveAt q hq).a₆)
      = (Ideal.Quotient.mk (I ^ n)) (partialEval tateA6 q n) :=
    Ideal.Quotient.eq.2 (evalAdic_sub_partialEval n tateA6 q hq)
  simp only [tateDefect, tateDefectTrunc, map_sub, map_add, map_mul, map_pow, hX, hY, h4, h6]

/-- ★★★★★★★**葉 (b) は有限の代数的主張の族に落ちる**。

`tateDefectTrunc n` は無限和も極限も含まない。★これがすべての `n` で `I^n` に入れば、
`IsHausdorff` により差そのものが 0 になる。 -/
theorem tateDefect_eq_zero_of_trunc [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (h : ∀ n, tateDefectTrunc n a w q ∈ I ^ n) : tateDefect a w q hq = 0 := by
  refine eq_zero_of_mem_pow (I := I) fun n => ?_
  have hd := tateDefect_sub_trunc (I := I) n a w q hq
  have hsum := Ideal.add_mem (I ^ n) hd (h n)
  simpa using hsum

end ABC3.Found.GaloisRep

/-! ## ★★★★★★切り詰めの函手性(完備性は要らない) -/

namespace ABC3.Found.GaloisRep

variable {R R' : Type} [CommRing R] [CommRing R']

theorem map_tateXtrunc (φ : R →+* R') (n : ℕ) (a w q : R)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hqa : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * a))
    (hqw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * w)) :
    φ (tateXtrunc n a w q) = tateXtrunc n (φ a) (φ w) (φ q) := by
  have hsa : φ (partialSum (fun m => tateXterm (q ^ (m + 1) * a)) n)
      = partialSum (fun m => tateXterm (φ q ^ (m + 1) * φ a)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateXterm φ (hqa m (Finset.mem_range.1 hm)), map_mul, map_pow]
  have hsw : φ (partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
      = partialSum (fun m => tateXterm (φ q ^ (m + 1) * φ w)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateXterm φ (hqw m (Finset.mem_range.1 hm)), map_mul, map_pow]
  rw [tateXtrunc, tateXtrunc, map_sub, map_add, map_add, map_add, map_mul,
    map_tateXterm φ ha, map_tateXterm φ hw, hsa, hsw, map_partialEval, map_ofNat φ 2]

theorem map_tateYtrunc (φ : R →+* R') (n : ℕ) (a w q : R)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hqa : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * a))
    (hqw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * w)) :
    φ (tateYtrunc n a w q) = tateYtrunc n (φ a) (φ w) (φ q) := by
  have hsa : φ (partialSum (fun m => tateYterm (q ^ (m + 1) * a)) n)
      = partialSum (fun m => tateYterm (φ q ^ (m + 1) * φ a)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateYterm φ (hqa m (Finset.mem_range.1 hm)), map_mul, map_pow]
  have hsw : φ (partialSum (fun m => tateYterm (q ^ (m + 1) * w)) n)
      = partialSum (fun m => tateYterm (φ q ^ (m + 1) * φ w)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateYterm φ (hqw m (Finset.mem_range.1 hm)), map_mul, map_pow]
  have hsxw : φ (partialSum (fun m => tateXterm (q ^ (m + 1) * w)) n)
      = partialSum (fun m => tateXterm (φ q ^ (m + 1) * φ w)) n := by
    rw [partialSum, partialSum, map_sum]
    refine Finset.sum_congr rfl fun m hm => ?_
    rw [map_tateXterm φ (hqw m (Finset.mem_range.1 hm)), map_mul, map_pow]
  rw [tateYtrunc, tateYtrunc, map_add, map_sub, map_sub, map_add, map_add, map_add,
    map_tateYterm φ ha, map_tateXterm φ hw, map_tateYterm φ hw, hsa, hsw, hsxw,
    map_partialEval]

/-- ★★★★★★**切り詰めた差は任意の環準同型で移る**——完備性は要らない。

★これが万有な環からの特殊化を可能にする。 -/
theorem map_tateDefectTrunc (φ : R →+* R') (n : ℕ) (a w q : R)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hqa : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * a))
    (hqw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * w)) :
    φ (tateDefectTrunc n a w q) = tateDefectTrunc n (φ a) (φ w) (φ q) := by
  simp only [tateDefectTrunc, map_sub, map_add, map_mul, map_pow,
    map_tateXtrunc φ n a w q ha hw hqa hqw, map_tateYtrunc φ n a w q ha hw hqa hqw,
    map_partialEval]

end ABC3.Found.GaloisRep

/-! ## ★★★★★★★★特殊化の規準 -/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★**特殊化の規準**——万有な環 `S` の上で切り詰めが `(AW)^n` の倍数なら、
任意の完備環で方程式の差が 0 になる。

★★これが段 6 の骨格である。残るのは「万有な環でその整除性を示す」ことだけになる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateDefect_eq_zero_of_specialize [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q)
    {S : Type} [CommRing S] (A W : S) (φ : S →+* R)
    (hA : φ A = a) (hW : φ W = w)
    (hAu : IsUnit (1 - A)) (hWu : IsUnit (1 - W))
    (hAq : ∀ m : ℕ, IsUnit (1 - (A * W) ^ (m + 1) * A))
    (hWq : ∀ m : ℕ, IsUnit (1 - (A * W) ^ (m + 1) * W))
    (H : ∀ n : ℕ, ((A * W) ^ n) ∣ tateDefectTrunc n A W (A * W)) :
    tateDefect a w q hq = 0 := by
  refine tateDefect_eq_zero_of_trunc (I := I) a w q hq fun n => ?_
  have hQ : φ (A * W) = q := by rw [map_mul, hA, hW, haw]
  have hmap := map_tateDefectTrunc φ n A W (A * W) hAu hWu
    (fun m _ => hAq m) (fun m _ => hWq m)
  rw [hA, hW, hQ] at hmap
  obtain ⟨c, hc⟩ := H n
  rw [← hmap, hc, map_mul, map_pow, hQ]
  exact Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)

/-! ## ★出典の紐付け(`.src`) -/

def tateDefect_eq_zero_of_trunc.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——葉 (b) を有限の主張の族に落とす)",
    sectionId := "genell-def-3-3" }

def tateDefect_eq_zero_of_specialize.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——万有な環からの特殊化の規準)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
