import ABC3.Found.GaloisRep.AdicContraction

/-!
# Galois (G6) 第 104 ブロック —— **★★★★`I` 進級数の和**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★葉 (a) の道具を先に作る

`Skeleton/GaloisRep/TateUniformization.lean` の葉 (a) は

    X(u,q) = ∑_{n∈ℤ} qⁿu/(1−qⁿu)² − 2s₁(q)

の**収束**である。★これは第 95 ブロックの `evalAdic` と同じ構造だが、
項が `c_n·qⁿ` の形をしていない——**一般の級数**である。

★★そこで `evalAdic` を一般化する:

    a : ℕ → R,  a n ∈ I^n  ⟹  ∑ a n は I 進収束する

★★★`evalAdic` は `a n = c_n·qⁿ` の場合として**この定理の系になる**
(`evalAdic_eq_adicSum`)。

## ★★★★★位相を使わない

第 95 ブロックと同じく `IsPrecomplete` で極限、`IsHausdorff` で一意性を出す。
★**Tate 一意化の級数もこの枠内で定義できる**——`qᵐu` と `qᵐu⁻¹` はどちらも
`I^m` に入るからである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `partialSum` / `partialSum_cauchy` | ★部分和は `I` 進 Cauchy |
| `adicSum` | ★★★**`I` 進級数の和** |
| `adicSum_spec` / `adicSum_unique` | ★★極限であること・一意性 |
| `adicSum_congr` / `adicSum_add` / `adicSum_smul` | ★係数計算 |
| `adicSum_shift` | ★★項を 1 つずらす `∑ a = a 0 + ∑ a(·+1)` |
| `adicSum_mem` | ★`a 0 ∈ I` なら和も `I` に入る |
| `evalAdic_eq_adicSum` | ★★★第 95 ブロックが系になる |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★部分和 -/

/-- ★部分和 `∑_{i<N} a i`。 -/
def partialSum (a : ℕ → R) (N : ℕ) : R := ∑ i ∈ Finset.range N, a i

/-- ★★**項が `I^n` に入るなら部分和は `I` 進 Cauchy である**。 -/
theorem partialSum_cauchy {a : ℕ → R} (ha : ∀ n, a n ∈ I ^ n) {m n : ℕ} (hmn : m ≤ n) :
    partialSum a m ≡ partialSum a n [SMOD (I ^ m • ⊤ : Submodule R R)] := by
  rw [SModEq.sub_mem]
  have hsub : partialSum a m - partialSum a n = -∑ k ∈ Finset.Ico m n, a k := by
    rw [partialSum, partialSum, Finset.range_eq_Ico, Finset.range_eq_Ico,
      ← Finset.sum_Ico_consecutive a (Nat.zero_le m) hmn]
    ring
  rw [hsub]
  have hmem : ∀ k ∈ Finset.Ico m n, a k ∈ I ^ m := by
    intro k hk
    exact Ideal.pow_le_pow_right (Finset.mem_Ico.1 hk).1 (ha k)
  simpa using neg_mem (Submodule.sum_mem _ hmem)

/-! ## ★★★和 -/

/-- ★★★**`I` 進級数の和**——項が `I^n` に入るなら部分和は収束する。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def adicSum [IsPrecomplete I R] (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) : R :=
  Classical.choose (IsPrecomplete.prec' (partialSum a) (partialSum_cauchy ha))

theorem adicSum_spec [IsPrecomplete I R] (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) (n : ℕ) :
    partialSum a n ≡ adicSum a ha [SMOD (I ^ n • ⊤ : Submodule R R)] :=
  Classical.choose_spec (IsPrecomplete.prec' (partialSum a) (partialSum_cauchy ha)) n

/-- ★★極限は一意である(`IsHausdorff`)。 -/
theorem adicSum_unique [IsAdicComplete I R] (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) (L : R)
    (hL : ∀ n, partialSum a n ≡ L [SMOD (I ^ n • ⊤ : Submodule R R)]) :
    adicSum a ha = L :=
  (IsHausdorff.eq_iff_smodEq (I := I)).2 (fun n => (adicSum_spec a ha n).symm.trans (hL n))

/-! ## ★係数計算 -/

theorem adicSum_congr [IsAdicComplete I R] {a b : ℕ → R} (ha : ∀ n, a n ∈ I ^ n)
    (hb : ∀ n, b n ∈ I ^ n) (hab : ∀ n, a n = b n) : adicSum a ha = adicSum b hb := by
  refine adicSum_unique _ _ _ (fun n => ?_)
  have hps : partialSum a n = partialSum b n := by
    simp only [partialSum]
    exact Finset.sum_congr rfl (fun i _ => hab i)
  rw [hps]
  exact adicSum_spec b hb n

theorem adicSum_add [IsAdicComplete I R] (a b : ℕ → R) (ha : ∀ n, a n ∈ I ^ n)
    (hb : ∀ n, b n ∈ I ^ n) :
    adicSum (fun n => a n + b n) (fun n => Submodule.add_mem _ (ha n) (hb n))
      = adicSum a ha + adicSum b hb := by
  refine adicSum_unique _ _ _ (fun n => ?_)
  have hps : partialSum (fun n => a n + b n) n = partialSum a n + partialSum b n := by
    simp only [partialSum]
    exact Finset.sum_add_distrib
  rw [hps]
  exact SModEq.add (adicSum_spec a ha n) (adicSum_spec b hb n)

theorem adicSum_smul [IsAdicComplete I R] (c : R) (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) :
    adicSum (fun n => c * a n) (fun n => Ideal.mul_mem_left _ _ (ha n)) = c * adicSum a ha := by
  refine adicSum_unique _ _ _ (fun n => ?_)
  have hps : partialSum (fun n => c * a n) n = c * partialSum a n := by
    simp only [partialSum, Finset.mul_sum]
  rw [hps, SModEq.sub_mem, ← mul_sub]
  have h := adicSum_spec a ha n
  rw [SModEq.sub_mem] at h
  have h2 : partialSum a n - adicSum a ha ∈ I ^ n := by simpa using h
  simpa using Ideal.mul_mem_left _ c h2

theorem adicSum_shift [IsAdicComplete I R] (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n) :
    adicSum a ha
      = a 0 + adicSum (fun n => a (n + 1))
          (fun n => Ideal.pow_le_pow_right (Nat.le_succ n) (ha (n + 1))) := by
  set hb : ∀ n, a (n + 1) ∈ I ^ n :=
    fun n => Ideal.pow_le_pow_right (Nat.le_succ n) (ha (n + 1)) with hbdef
  set T := adicSum (fun n => a (n + 1)) hb with hTdef
  refine adicSum_unique _ _ _ (fun n => ?_)
  rw [SModEq.sub_mem]
  have hps : partialSum a n = partialSum (fun n => a (n + 1)) n + a 0 - a n := by
    induction n with
    | zero => simp [partialSum]
    | succ m ih =>
      simp only [partialSum, Finset.sum_range_succ] at ih ⊢
      rw [ih]
      ring
  have hspec := adicSum_spec (fun n => a (n + 1)) hb n
  rw [SModEq.sub_mem] at hspec
  have h1 : partialSum (fun n => a (n + 1)) n - T ∈ I ^ n := by simpa using hspec
  have h2 : partialSum a n - (a 0 + T)
      = (partialSum (fun n => a (n + 1)) n - T) - a n := by
    rw [hps]
    ring
  have h3 : partialSum a n - (a 0 + T) ∈ I ^ n := by
    rw [h2]
    exact Submodule.sub_mem _ h1 (ha n)
  simpa using h3

/-- ★第 0 項が `I` に入るなら和も `I` に入る。 -/
theorem adicSum_mem [IsPrecomplete I R] (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n)
    (ha0 : a 0 ∈ I) : adicSum a ha ∈ I := by
  have h := adicSum_spec a ha 1
  rw [SModEq.sub_mem] at h
  have h1 : partialSum a 1 - adicSum a ha ∈ I := by simpa using h
  have h2 : partialSum a 1 = a 0 := by simp [partialSum]
  have h3 : adicSum a ha = a 0 - (partialSum a 1 - adicSum a ha) := by rw [h2]; ring
  rw [h3]
  exact Submodule.sub_mem _ ha0 h1

/-! ## ★★★第 95 ブロックが系になる -/

/-- ★★★**冪級数の値は `I` 進級数の和の特別な場合である**。 -/
theorem evalAdic_eq_adicSum [IsAdicComplete I R] (f : PowerSeries ℤ) (q : R) (hq : q ∈ I) :
    evalAdic f q hq
      = adicSum (fun n => ((PowerSeries.coeff n f : ℤ) : R) * q ^ n)
          (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n)) :=
  (adicSum_unique _ _ _ (fun n => evalAdic_spec f q hq n)).symm

/-! ## ★出典の紐付け(`.src`) -/

def adicSum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化の級数が収束すること——葉 (a) の道具)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
