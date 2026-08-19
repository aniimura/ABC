import ABC3.Found.GaloisRep.TateSeries
import Mathlib.RingTheory.AdicCompletion.Basic

/-!
# Galois (G6) 第 95 ブロック —— **★★★完備環での冪級数の値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★位相を使わずに特殊化する

第 94 ブロックは Tate 曲線を **`ℤ⟦q⟧` 上**に作った。
★これを完備離散付値環 `R` の元 `q ∈ 𝔪` で特殊化したい。

★★mathlib の `PowerSeries.eval₂Hom` は**位相**(`CompleteSpace`・`IsLinearTopology`)を要求するが、
(G6) の設定にあるのは **`IsAdicComplete`(代数的)**である。

★★★そこで**部分和の Cauchy 列**を直接使う:

    ∑_{n<N} c_n q^n     は  I 進 Cauchy(差は I^N に入る)

`IsPrecomplete` が極限を与え、`IsHausdorff` が一意性を与える。
★★★★**位相を経由せずに評価が定義できる。**

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `partialEval` | ★部分和 |
| `partialEval_cauchy` | ★★部分和は `I` 進 Cauchy |
| `evalAdic` | ★★★**完備環での値** |
| `evalAdic_spec` / `evalAdic_unique` | ★★極限であること・一意性 |
| `evalAdic_add` | ★加法性 |
-/

namespace ABC3.Found.GaloisRep

open PowerSeries

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★部分和 `∑_{n<N} c_n q^n`。 -/
noncomputable def partialEval (f : PowerSeries ℤ) (q : R) (N : ℕ) : R :=
  ∑ n ∈ Finset.range N, ((PowerSeries.coeff n f : ℤ) : R) * q ^ n

/-- ★★**部分和は `I` 進 Cauchy である**。 -/
theorem partialEval_cauchy (f : PowerSeries ℤ) {q : R} (hq : q ∈ I) {m n : ℕ} (hmn : m ≤ n) :
    partialEval f q m ≡ partialEval f q n [SMOD (I ^ m • ⊤ : Submodule R R)] := by
  rw [SModEq.sub_mem]
  have hsub : partialEval f q m - partialEval f q n
      = -∑ k ∈ Finset.Ico m n, ((PowerSeries.coeff k f : ℤ) : R) * q ^ k := by
    rw [partialEval, partialEval, Finset.range_eq_Ico, Finset.range_eq_Ico,
      ← Finset.sum_Ico_consecutive (fun k => ((PowerSeries.coeff k f : ℤ) : R) * q ^ k)
        (Nat.zero_le m) hmn]
    ring
  rw [hsub]
  have hmem : ∀ k ∈ Finset.Ico m n, ((PowerSeries.coeff k f : ℤ) : R) * q ^ k ∈ I ^ m := by
    intro k hk
    have hkm : m ≤ k := (Finset.mem_Ico.1 hk).1
    have hpow : q ^ k = q ^ m * q ^ (k - m) := by
      rw [← pow_add]
      congr 1
      omega
    rw [hpow, ← mul_assoc]
    exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq m))
  have hs : (∑ k ∈ Finset.Ico m n, ((PowerSeries.coeff k f : ℤ) : R) * q ^ k) ∈ I ^ m :=
    Submodule.sum_mem _ hmem
  simpa using neg_mem hs

/-- ★★★**完備環での冪級数の値**——部分和の `I` 進極限。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def evalAdic [IsPrecomplete I R] (f : PowerSeries ℤ) (q : R) (hq : q ∈ I) : R :=
  Classical.choose (IsPrecomplete.prec' (partialEval f q) (partialEval_cauchy f hq))

theorem evalAdic_spec [IsPrecomplete I R] (f : PowerSeries ℤ) (q : R) (hq : q ∈ I) (n : ℕ) :
    partialEval f q n ≡ evalAdic f q hq [SMOD (I ^ n • ⊤ : Submodule R R)] :=
  Classical.choose_spec (IsPrecomplete.prec' (partialEval f q) (partialEval_cauchy f hq)) n

/-- ★★極限は一意である(`IsHausdorff`)。 -/
theorem evalAdic_unique [IsAdicComplete I R] (f : PowerSeries ℤ) (q : R) (hq : q ∈ I) (L : R)
    (hL : ∀ n, partialEval f q n ≡ L [SMOD (I ^ n • ⊤ : Submodule R R)]) :
    evalAdic f q hq = L :=
  (IsHausdorff.eq_iff_smodEq (I := I)).2 (fun n => (evalAdic_spec f q hq n).symm.trans (hL n))

theorem partialEval_add (f g : PowerSeries ℤ) (q : R) (n : ℕ) :
    partialEval (f + g) q n = partialEval f q n + partialEval g q n := by
  simp only [partialEval, map_add, Int.cast_add, add_mul]
  rw [Finset.sum_add_distrib]

/-- ★加法性。 -/
theorem evalAdic_add [IsAdicComplete I R] (f g : PowerSeries ℤ) (q : R) (hq : q ∈ I) :
    evalAdic (f + g) q hq = evalAdic f q hq + evalAdic g q hq := by
  refine evalAdic_unique _ _ _ _ (fun n => ?_)
  rw [partialEval_add]
  exact SModEq.add (evalAdic_spec f q hq n) (evalAdic_spec g q hq n)

/-! ## ★出典の紐付け(`.src`) -/

def evalAdic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線——完備環への特殊化)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
