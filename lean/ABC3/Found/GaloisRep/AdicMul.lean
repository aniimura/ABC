import ABC3.Found.GaloisRep.TateSigma

/-!
# Galois (G6) 第 112 ブロック —— **★★★★★`I` 進級数の Cauchy 積**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★葉 (b) は**積**を要求する

葉 (b) は `Y² + XY = X³ + a₄X + a₆` の verify である。
★`X`・`Y` は第 110–111 ブロックで `q` の冪級数として書けたが、
**`Y²`・`XY`・`X³` を係数の言葉に直すには積の公式が要る**:

    (∑ aₙ)(∑ bₙ) = ∑ₙ (∑_{i≤n} aᵢ·b_{n−i})

★★`aₙ ∈ Iⁿ`・`bₙ ∈ Iⁿ` という**次数付き**の仮定があるので、これは成り立つ。

## ★★★★証明——三角形と正方形の差

    A_k·B_k = ∑_{i<k} aᵢ·B_k            (正方形)
    C_k     = ∑_{i<k} aᵢ·(∑_{j<k−i} bⱼ) (三角形、`cauchy_partial`)

★差は `∑_{i<k} aᵢ·(∑_{k−i ≤ j < k} bⱼ)` であり、
各項は `i + j ≥ k` なので `I^k` に入る。★★あとは
`A·B − A_k·B_k = (A−A_k)·B + A_k·(B−B_k) ∈ I^k` と合わせて `IsHausdorff`。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `cauchy_partial` | ★★★有限版(三角形 = 正方形の一部) |
| `adicSum_mul` | ★★★★★**Cauchy 積** |
| `adicSum_neg` / `adicSum_sub` | ★符号 |
| `adicSum_add_const` | ★★定数を第 0 項に押し込む |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★有限版 -/

/-- ★★★**三角形の和は「係数×部分和」の形になる**。 -/
theorem cauchy_partial (a b : ℕ → R) (k : ℕ) :
    ∑ n ∈ Finset.range k, (∑ i ∈ Finset.range (n + 1), a i * b (n - i))
      = ∑ i ∈ Finset.range k, a i * (∑ j ∈ Finset.range (k - i), b j) := by
  induction k with
  | zero => simp
  | succ m ih =>
    have hL : ∑ n ∈ Finset.range (m + 1), (∑ i ∈ Finset.range (n + 1), a i * b (n - i))
        = (∑ i ∈ Finset.range m, a i * (∑ j ∈ Finset.range (m - i), b j))
          + ∑ i ∈ Finset.range (m + 1), a i * b (m - i) := by
      rw [Finset.sum_range_succ, ih]
    have hR : ∑ i ∈ Finset.range (m + 1), a i * (∑ j ∈ Finset.range (m + 1 - i), b j)
        = (∑ i ∈ Finset.range m, a i * (∑ j ∈ Finset.range (m + 1 - i), b j))
          + a m * (∑ j ∈ Finset.range (m + 1 - m), b j) := Finset.sum_range_succ _ _
    rw [hL, hR]
    have hstep : ∀ i ∈ Finset.range m, a i * (∑ j ∈ Finset.range (m + 1 - i), b j)
        = a i * (∑ j ∈ Finset.range (m - i), b j) + a i * b (m - i) := by
      intro i hi
      have him : i < m := Finset.mem_range.1 hi
      have h1 : m + 1 - i = (m - i) + 1 := by omega
      rw [h1, Finset.sum_range_succ, mul_add]
    rw [Finset.sum_congr rfl hstep, Finset.sum_add_distrib, Finset.sum_range_succ]
    have h2 : m + 1 - m = 1 := by omega
    rw [h2, Finset.sum_range_one, Nat.sub_self]
    ring

theorem cauchy_mem {a b : ℕ → R} (ha : ∀ n, a n ∈ I ^ n) (hb : ∀ n, b n ∈ I ^ n) (n : ℕ) :
    (∑ i ∈ Finset.range (n + 1), a i * b (n - i)) ∈ I ^ n := by
  refine Submodule.sum_mem _ (fun i hi => ?_)
  have hin : i ≤ n := by
    have := Finset.mem_range.1 hi
    omega
  have hmem : a i * b (n - i) ∈ I ^ (i + (n - i)) := by
    rw [pow_add]
    exact Ideal.mul_mem_mul (ha i) (hb (n - i))
  have he : i + (n - i) = n := by omega
  rwa [he] at hmem

/-! ## ★★★★★Cauchy 積 -/

/-- ★★★★★**`I` 進級数の Cauchy 積**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem adicSum_mul [IsAdicComplete I R] (a b : ℕ → R)
    (ha : ∀ n, a n ∈ I ^ n) (hb : ∀ n, b n ∈ I ^ n) :
    adicSum a ha * adicSum b hb
      = adicSum (fun n => ∑ i ∈ Finset.range (n + 1), a i * b (n - i)) (cauchy_mem ha hb) := by
  refine (IsHausdorff.eq_iff_smodEq (I := I)).2 (fun k => ?_)
  rw [SModEq.sub_mem]
  set A := adicSum a ha with hA
  set B := adicSum b hb with hB
  set C := adicSum (fun n => ∑ i ∈ Finset.range (n + 1), a i * b (n - i)) (cauchy_mem ha hb) with hC
  set Ak := partialSum a k with hAk
  set Bk := partialSum b k with hBk
  set Ck := partialSum (fun n => ∑ i ∈ Finset.range (n + 1), a i * b (n - i)) k with hCk
  have h1 : A - Ak ∈ I ^ k := by
    have := adicSum_spec a ha k
    rw [SModEq.sub_mem] at this
    simpa using neg_mem (show Ak - A ∈ I ^ k by simpa using this)
  have h2 : B - Bk ∈ I ^ k := by
    have := adicSum_spec b hb k
    rw [SModEq.sub_mem] at this
    simpa using neg_mem (show Bk - B ∈ I ^ k by simpa using this)
  have h3 : C - Ck ∈ I ^ k := by
    have := adicSum_spec (fun n => ∑ i ∈ Finset.range (n + 1), a i * b (n - i))
      (cauchy_mem ha hb) k
    rw [SModEq.sub_mem] at this
    simpa using neg_mem (show Ck - C ∈ I ^ k by simpa using this)
  have hABmem : A * B - Ak * Bk ∈ I ^ k := by
    have he : A * B - Ak * Bk = (A - Ak) * B + Ak * (B - Bk) := by ring
    rw [he]
    exact Submodule.add_mem _ (Ideal.mul_mem_right _ _ h1) (Ideal.mul_mem_left _ _ h2)
  have hkey : Ak * Bk - Ck ∈ I ^ k := by
    have e1 : Ak * Bk = ∑ i ∈ Finset.range k, a i * Bk := by
      rw [hAk, partialSum, Finset.sum_mul]
    have e2 : Ck = ∑ i ∈ Finset.range k, a i * partialSum b (k - i) := by
      rw [hCk, partialSum]
      exact cauchy_partial a b k
    rw [e1, e2, ← Finset.sum_sub_distrib]
    refine Submodule.sum_mem _ (fun i hi => ?_)
    have hik : i < k := Finset.mem_range.1 hi
    have e3 : a i * Bk - a i * partialSum b (k - i) = a i * (Bk - partialSum b (k - i)) := by ring
    rw [e3, hBk, partialSum_sub_partialSum (show k - i ≤ k by omega), Finset.mul_sum]
    refine Submodule.sum_mem _ (fun j hj => ?_)
    have hj1 : k - i ≤ j := (Finset.mem_Ico.1 hj).1
    have hmem : a i * b j ∈ I ^ (i + j) := by
      rw [pow_add]
      exact Ideal.mul_mem_mul (ha i) (hb j)
    exact Ideal.pow_le_pow_right (by omega) hmem
  have hfin : A * B - C = (A * B - Ak * Bk) - (C - Ck) + (Ak * Bk - Ck) := by ring
  have hgoal : A * B - C ∈ I ^ k := by
    rw [hfin]
    exact Submodule.add_mem _ (Submodule.sub_mem _ hABmem h3) hkey
  simpa using hgoal

/-! ## ★符号と定数 -/

theorem adicSum_neg [IsAdicComplete I R] (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) :
    adicSum (fun n => -(a n)) (fun n => neg_mem (ha n)) = - adicSum a ha := by
  refine adicSum_unique _ _ _ (fun n => ?_)
  have hps : partialSum (fun n => -(a n)) n = - partialSum a n := by
    simp only [partialSum, Finset.sum_neg_distrib]
  rw [hps, SModEq.sub_mem]
  have h := adicSum_spec a ha n
  rw [SModEq.sub_mem] at h
  have h2 : partialSum a n - adicSum a ha ∈ I ^ n := by simpa using h
  have he : - partialSum a n - -adicSum a ha = -(partialSum a n - adicSum a ha) := by ring
  rw [he]
  simpa using neg_mem h2

theorem adicSum_sub [IsAdicComplete I R] (a b : ℕ → R)
    (ha : ∀ n, a n ∈ I ^ n) (hb : ∀ n, b n ∈ I ^ n) :
    adicSum (fun n => a n - b n) (fun n => Submodule.sub_mem _ (ha n) (hb n))
      = adicSum a ha - adicSum b hb := by
  refine adicSum_unique _ _ _ (fun n => ?_)
  have hps : partialSum (fun n => a n - b n) n = partialSum a n - partialSum b n := by
    simp only [partialSum]
    exact Finset.sum_sub_distrib (s := Finset.range n) (f := a) (g := b)
  rw [hps]
  exact SModEq.sub (adicSum_spec a ha n) (adicSum_spec b hb n)

theorem const_mem (c : R) (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) (n : ℕ) :
    (if n = 0 then c + a 0 else a n) ∈ I ^ n := by
  by_cases h : n = 0
  · subst h
    simp
  · rw [if_neg h]
    exact ha n

/-- ★★定数を第 0 項に押し込む。 -/
theorem adicSum_add_const [IsAdicComplete I R] (c : R) (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) :
    c + adicSum a ha
      = adicSum (fun n => if n = 0 then c + a 0 else a n) (const_mem c a ha) := by
  refine (adicSum_unique _ _ _ (fun n => ?_)).symm
  rw [SModEq.sub_mem]
  rcases Nat.eq_zero_or_pos n with h0 | h0
  · subst h0
    simp
  · have hps : partialSum (fun n => if n = 0 then c + a 0 else a n) n = c + partialSum a n := by
      rw [partialSum, partialSum]
      rcases Nat.exists_eq_succ_of_ne_zero (by omega : n ≠ 0) with ⟨m, rfl⟩
      rw [Finset.sum_range_succ', Finset.sum_range_succ']
      have hcong : ∀ i ∈ Finset.range m,
          (if i + 1 = 0 then c + a 0 else a (i + 1)) = a (i + 1) := by
        intro i _
        rw [if_neg (by omega)]
      rw [Finset.sum_congr rfl hcong]
      simp
      ring
    rw [hps]
    have h := adicSum_spec a ha n
    rw [SModEq.sub_mem] at h
    have h2 : partialSum a n - adicSum a ha ∈ I ^ n := by simpa using h
    have he : c + partialSum a n - (c + adicSum a ha) = partialSum a n - adicSum a ha := by ring
    rw [he]
    simpa using h2

/-! ## ★出典の紐付け(`.src`) -/

def adicSum_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 級数の積——葉 (b) の道具)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
