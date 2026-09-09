import ABC3.Found.PGC.PadicValueGroup
import Mathlib.RingTheory.Polynomial.Cyclotomic.Basic

/-!
# [pGC] 最後の「数学」が定理になった —— `Φ_N(X+1)` から関係式を出す

## 持ち場（前波の表で唯一「数学」と残っていた行）

「次数 `n` のモニック関係式」。

## ★結果

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `pow_eq_sum_of_monic_aeval_zero` | ★抽象核。モニック多項式が消せば `x^n = Σ_{j<n} (−coeff j) x^j` |
| §2 | `relation_of_primitive_root` | ★円分の場合。多項式は `Φ_N(X+1)` |
| §3 | `finrank_padic_of_zeta` | ★★★**絶対次数が円分の材料だけから出る** |

★§2 で使う多項式は `ZetaSubOnePrime` が Eisenstein 性で使ったのと**同じもの**だが、
★ここでは Eisenstein 性は要らず、**モニック性・次数・根であること**の 3 つだけを使う。
⇒ `cyclotomic.monic` / `natDegree_comp` ＋ `natDegree_cyclotomic` /
`IsPrimitiveRoot.isRoot_cyclotomic` ＋ `map_cyclotomic` の 4 本で閉じた。

## ★★★残っているものを**数えた**（本体が「何本か測っていない」と書いた点）

`MainPartCoeffs.loss_le_of_main_coeffs`（主定理）まで繋ぐのに要るものを全部並べる。
★**数学の入力はもう無い。** 残りは**模型の構成**である。

**インスタンス（12 本）**

1. `[NormedField L]` `[IsUltrametricDist L]` `[CharZero L]` `[NormOneClass L]`
2. `[Field E₁]`
3. `[NormedAlgebra ℚ_p L]` `[Algebra ℚ_p E₁]` `[Algebra E₁ L]` `[IsScalarTower ℚ_p E₁ L]`
4. `[FiniteDimensional ℚ_p L]` `[FiniteDimensional E₁ L]` `[FiniteDimensional ℚ_p E₁]`

**定義的な事実（6 本）**

5. `IsPrimitiveRoot ζ (p^{m+2})`、`IsPrimitiveRoot ξ (p^{m+1})`（`ξ = ζ^p`）
6. `Algebra.adjoin ℚ_p {ζ−1} = ⊤`（`L = ℚ_p(ζ)`）
7. `q ≠ 1`（`L ≠ E₁`）
8. `‖(p : L)‖ < 1`、`(p : L) ≠ 0`

★これだけである。★**`p ≥ 3` の `loss ≤ 2p−2` の鎖は、模型を作れば閉じる。**

## ★配管 #358 が**同じ波で 2 度**効いた

§3 で `omit [Fact p.Prime] in` と書いたら

```
error: cannot omit referenced section variable `inst✝`
```

（`ℚ_[p]` の記法自体が `[Fact p.Prime]` を参照するので外せない）。★そのとき
`#print axioms` に `sorryAx` が出たが、★#358 のとおり**上のエラーを先に読んだ**。
修正後に 4 件すべて `[propext, Classical.choice, Quot.sound]` を確認した。

## 逸脱の記録

- §2 は `K` を一般の体のままにした（`ℚ_p` に限らない）。★`IsDomain M` と
  `IsPrimitiveRoot` があれば成り立つ。
- §3 の `hadj`（`adjoin = ⊤`）は**塔の定義**として受ける。★これを「証明する」のは
  `L` を `AdjoinRoot` などで**構成**するときに自動的に付いてくる性質である。
-/

namespace ABC3.Found.PGC

namespace MonicRelation

open Polynomial

/-! ## §1 抽象核 —— モニック多項式が消せば冪の関係式が出る -/

section Kernel

theorem pow_eq_sum_of_monic_aeval_zero {K M : Type*} [CommRing K] [CommRing M] [Algebra K M]
    {P : Polynomial K} {n : ℕ} (hmonic : P.Monic) (hdeg : P.natDegree = n)
    {x : M} (hx : Polynomial.aeval x P = 0) :
    x ^ n = ∑ j ∈ Finset.range n, algebraMap K M (-P.coeff j) * x ^ j := by
  have hsum := Polynomial.aeval_eq_sum_range (p := P) x
  rw [hx, hdeg, Finset.sum_range_succ] at hsum
  have hcn : P.coeff n = 1 := by
    rw [← hdeg]
    exact hmonic.coeff_natDegree
  rw [hcn, one_smul] at hsum
  have hxn : x ^ n = -∑ i ∈ Finset.range n, P.coeff i • x ^ i := by
    linear_combination -hsum
  have hfinal : ∑ j ∈ Finset.range n, algebraMap K M (-P.coeff j) * x ^ j
      = -∑ i ∈ Finset.range n, P.coeff i • x ^ i := by
    rw [← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl fun i _ => by rw [Algebra.smul_def, map_neg]; ring
  rw [hfinal, hxn]

end Kernel

/-! ## §2 円分の場合 —— `π = ζ − 1` の関係式 -/

section Cyclotomic

/-- ★★**残った最後の「数学」**。`ζ` が原始 `N` 乗根なら、`π = ζ − 1` は
次数 `φ(N)` のモニックな関係式を満たす。

★多項式は `Φ_N(X+1)`（`ZetaSubOnePrime` が Eisenstein 性で使ったのと**同じもの**）。
★ここでは Eisenstein 性は要らず、★**モニック性と次数と根であること**だけを使う。 -/
theorem relation_of_primitive_root {K M : Type*} [Field K] [CommRing M] [IsDomain M]
    [Algebra K M] {N : ℕ} (hN : 0 < N) {ζ : M} (hζ : IsPrimitiveRoot ζ N) :
    ∃ c : ℕ → K, (ζ - 1) ^ N.totient
      = ∑ j ∈ Finset.range N.totient, algebraMap K M (c j) * (ζ - 1) ^ j := by
  classical
  set P : Polynomial K := (cyclotomic N K).comp (X + C 1) with hP
  have hmonic : P.Monic := (cyclotomic.monic N K).comp_X_add_C 1
  have hdeg : P.natDegree = N.totient := by
    rw [hP, natDegree_comp, natDegree_cyclotomic, natDegree_X_add_C, mul_one]
  have hroot : Polynomial.aeval (ζ - 1) P = 0 := by
    rw [hP, aeval_comp]
    have h1 : Polynomial.aeval (ζ - 1) (X + C (1 : K)) = ζ := by simp
    rw [h1, Polynomial.aeval_def, ← Polynomial.eval_map, Polynomial.map_cyclotomic]
    exact hζ.isRoot_cyclotomic hN
  exact ⟨fun j => -P.coeff j, pow_eq_sum_of_monic_aeval_zero hmonic hdeg hroot⟩

end Cyclotomic

/-! ## §3 ★組み上げ —— 絶対次数が円分の材料だけから出る -/

section Assemble

variable {p : ℕ} [Fact p.Prime]

/-- `algebraMap ℚ_p M (p : ℚ_p) = (p : M)`（自然数の像）。 -/
theorem norm_algebraMap_natCast {M : Type*} [NormedField M] [Algebra ℚ_[p] M] :
    ‖algebraMap ℚ_[p] M (p : ℚ_[p])‖ = ‖(p : M)‖ := by
  rw [map_natCast]

/-- ★★★**絶対次数 `finrank ℚ_p L = φ(p^{m+1})` が円分の材料だけから出る**。

入力は `IsPrimitiveRoot ζ (p^{m+1})`、`‖p‖ < 1`、`(p : M) ≠ 0`、
`Algebra.adjoin ℚ_p {ζ−1} = ⊤`（塔の定義）だけ。
★`‖ζ−1‖^{φ} = ‖p‖` は `ZetaSubOnePrime`（自分）、
★`‖ζ−1‖ < 1` は `CyclotomicSubstitution.norm_sub_one_lt_one`（自分）。 -/
theorem finrank_padic_of_zeta {M : Type*} [NormedField M] [IsUltrametricDist M] [CharZero M]
    [NormedAlgebra ℚ_[p] M] [NormOneClass M] [FiniteDimensional ℚ_[p] M]
    {m : ℕ} {ζ : M} (hζ : IsPrimitiveRoot ζ (p ^ (m + 1)))
    (hlt : ‖(p : M)‖ < 1) (hne : (p : M) ≠ 0)
    (hadj : Algebra.adjoin ℚ_[p] ({ζ - 1} : Set M) = ⊤) :
    Module.finrank ℚ_[p] M = Nat.totient (p ^ (m + 1)) := by
  have hp : Nat.Prime p := Fact.out
  have hpow : ‖ζ - 1‖ ^ (Nat.totient (p ^ (m + 1))) = ‖(p : M)‖ :=
    ZetaSubOnePrime.norm_zeta_sub_one_pow hζ hlt hne
  have hsmall : ‖ζ - 1‖ < 1 := CyclotomicSubstitution.norm_sub_one_lt_one hζ hlt hne
  have hπ0 : 0 < ‖ζ - 1‖ := by
    rcases eq_or_lt_of_le (norm_nonneg (ζ - 1)) with h | h
    · exfalso
      rw [← h, zero_pow (Nat.totient_pos.mpr (pow_pos hp.pos _)).ne'] at hpow
      exact hne (norm_eq_zero.mp hpow.symm)
    · exact h
  obtain ⟨c, hc⟩ :=
    relation_of_primitive_root (K := ℚ_[p]) (pow_pos hp.pos (m + 1)) hζ
  exact PadicValueGroup.finrank_padic_eq hπ0 (ne_of_lt hsmall)
    (by rw [norm_algebraMap_natCast]; exact hpow) c hc hadj

end Assemble

/-! ## §4 使っている公理の一覧 -/

#print axioms pow_eq_sum_of_monic_aeval_zero
#print axioms relation_of_primitive_root
#print axioms norm_algebraMap_natCast
#print axioms finrank_padic_of_zeta

end MonicRelation

end ABC3.Found.PGC
