import ABC3.Found.GaloisRep.TatePair

/-!
# Galois (G6) 第 108 ブロック —— **★★★★Tate 級数の項の q 展開**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★葉 (b) の入口——閉じた式から級数へ

葉 (b)(`Y² + XY = X³ + a₄X + a₆`)を verify するには、
`X`・`Y` を **`q` の冪級数として展開**しなければならない。★その第一歩が

    f(t) = t/(1−t)²  =  ∑_{n≥0} n·tⁿ
    g(t) = t²/(1−t)³ =  ∑_{n≥0} C(n,2)·tⁿ

である。★★これを入れると古典的な q 展開

    X = u/(1−u)² + ∑_{N≥1} (∑_{d|N} d(u^d + u^{−d} − 2)) q^N

が出る(その組み立ては次の層)。

## ★★★有限版から `I` 進極限へ

証明は**有限和の恒等式**を作って `I` 進極限を取る:

    (1−t)²·∑_{n<N} n tⁿ    = t − N tᴺ + (N−1) tᴺ⁺¹
    (1−t)³·∑_{n<N} C(n,2)tⁿ = t² − tᴺ·(C(N,2)(1−t)² + N t(1−t) + t²)

★誤差項はどちらも `I^N` に入るので、`IsHausdorff` で `(1−t)²·S = t` が出る。
★★あとは第 105 ブロックの `(1−t)²·f(t) = t` と突き合わせ、
`1−t` が単元(`isUnit_one_sub`)であることから割って結論する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `geomDeriv_partial` | ★★有限版(`f`) |
| `geomDeriv2_partial` | ★★有限版(`g`) |
| `tateXterm_eq_adicSum` | ★★★★**`t/(1−t)² = ∑ n tⁿ`** |
| `tateYterm_eq_adicSum` | ★★★★**`t²/(1−t)³ = ∑ C(n,2) tⁿ`** |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★有限版の恒等式 -/

/-- ★★`(1−t)²·∑_{n<N} n tⁿ = t − N tᴺ + (N−1) tᴺ⁺¹`。 -/
theorem geomDeriv_partial (t : R) (N : ℕ) :
    (1 - t) ^ 2 * (∑ n ∈ Finset.range N, (n : R) * t ^ n)
      = t - (N : R) * t ^ N + ((N : R) - 1) * t ^ (N + 1) := by
  induction N with
  | zero => simp
  | succ m ih =>
    rw [Finset.sum_range_succ, mul_add, ih]
    push_cast
    ring

/-- ★★`(1−t)³·∑_{n<N} C(n,2) tⁿ = t² − tᴺ·(C(N,2)(1−t)² + N t(1−t) + t²)`。 -/
theorem geomDeriv2_partial (t : R) (N : ℕ) :
    (1 - t) ^ 3 * (∑ n ∈ Finset.range N, (n.choose 2 : R) * t ^ n)
      = t ^ 2 - t ^ N * ((N.choose 2 : R) * (1 - t) ^ 2 + (N : R) * t * (1 - t) + t ^ 2) := by
  induction N with
  | zero => simp
  | succ m ih =>
    rw [Finset.sum_range_succ, mul_add, ih]
    have hc : ((m + 1).choose 2 : R) = (m.choose 2 : R) + (m : R) := by
      rw [Nat.choose_succ_succ, Nat.choose_one_right]
      push_cast
      ring
    rw [hc]
    push_cast
    ring

/-! ## ★★★★q 展開 -/

/-- ★★★★**`t/(1−t)² = ∑_{n≥0} n tⁿ`** —— `X` の項の q 展開。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateXterm_eq_adicSum [IsAdicComplete I R] {t : R} (ht : t ∈ I) :
    tateXterm t
      = adicSum (fun n => (n : R) * t ^ n)
          (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) := by
  set S := adicSum (fun n => (n : R) * t ^ n)
    (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) with hS
  have hkey : (1 - t) ^ 2 * S = t := by
    refine (IsHausdorff.eq_iff_smodEq (I := I)).2 (fun N => ?_)
    rw [SModEq.sub_mem]
    set P := partialSum (fun n => (n : R) * t ^ n) N with hP
    have hspec := adicSum_spec (fun n => (n : R) * t ^ n)
      (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) N
    rw [SModEq.sub_mem] at hspec
    have h1 : P - S ∈ I ^ N := by simpa using hspec
    have h1' : S - P ∈ I ^ N := by simpa using neg_mem h1
    have hgd : (1 - t) ^ 2 * P = t - (N : R) * t ^ N + ((N : R) - 1) * t ^ (N + 1) := by
      rw [hP, partialSum]
      exact geomDeriv_partial t N
    have h2 : (1 - t) ^ 2 * S - t
        = (1 - t) ^ 2 * (S - P) + (-(N : R) * t ^ N + ((N : R) - 1) * t ^ (N + 1)) := by
      rw [mul_sub, hgd]
      ring
    have h3 : (1 - t) ^ 2 * S - t ∈ I ^ N := by
      rw [h2]
      refine Submodule.add_mem _ (Ideal.mul_mem_left _ _ h1') (Submodule.add_mem _ ?_ ?_)
      · exact Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht N)
      · exact Ideal.pow_le_pow_right (Nat.le_succ N)
          (Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht (N + 1)))
    simpa using h3
  have hu : IsUnit ((1 - t) ^ 2) := (isUnit_one_sub ht).pow 2
  exact ((IsUnit.mul_right_inj hu).mp ((mul_tateXterm ht).trans hkey.symm))

/-- ★★★★**`t²/(1−t)³ = ∑_{n≥0} C(n,2) tⁿ`** —— `Y` の項の q 展開。 -/
theorem tateYterm_eq_adicSum [IsAdicComplete I R] {t : R} (ht : t ∈ I) :
    tateYterm t
      = adicSum (fun n => (n.choose 2 : R) * t ^ n)
          (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) := by
  set S := adicSum (fun n => (n.choose 2 : R) * t ^ n)
    (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) with hS
  have hkey : (1 - t) ^ 3 * S = t ^ 2 := by
    refine (IsHausdorff.eq_iff_smodEq (I := I)).2 (fun N => ?_)
    rw [SModEq.sub_mem]
    set P := partialSum (fun n => (n.choose 2 : R) * t ^ n) N with hP
    have hspec := adicSum_spec (fun n => (n.choose 2 : R) * t ^ n)
      (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) N
    rw [SModEq.sub_mem] at hspec
    have h1 : P - S ∈ I ^ N := by simpa using hspec
    have h1' : S - P ∈ I ^ N := by simpa using neg_mem h1
    have hgd : (1 - t) ^ 3 * P
        = t ^ 2 - t ^ N * ((N.choose 2 : R) * (1 - t) ^ 2 + (N : R) * t * (1 - t) + t ^ 2) := by
      rw [hP, partialSum]
      exact geomDeriv2_partial t N
    have h2 : (1 - t) ^ 3 * S - t ^ 2
        = (1 - t) ^ 3 * (S - P)
          + (- (t ^ N * ((N.choose 2 : R) * (1 - t) ^ 2 + (N : R) * t * (1 - t) + t ^ 2))) := by
      rw [mul_sub, hgd]
      ring
    have h3 : (1 - t) ^ 3 * S - t ^ 2 ∈ I ^ N := by
      rw [h2]
      refine Submodule.add_mem _ (Ideal.mul_mem_left _ _ h1') (neg_mem ?_)
      exact Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow ht N)
    simpa using h3
  have hu : IsUnit ((1 - t) ^ 3) := (isUnit_one_sub ht).pow 3
  exact ((IsUnit.mul_right_inj hu).mp ((mul_tateYterm ht).trans hkey.symm))

/-! ## ★出典の紐付け(`.src`) -/

def tateXterm_eq_adicSum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 級数の項の q 展開——葉 (b) の入口)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
