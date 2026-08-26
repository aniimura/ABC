/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.SplitInfinite
import Mathlib.RingTheory.Localization.Integral
import Mathlib.FieldTheory.Minpoly.IsIntegrallyClosed

/-!
# ほとんどすべての `p` で mod `p` は分離的(鎖 `cheb` の `cheb-spl-infinite` を閉じる)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★残っていた 1 段

`SplitInfinite.lean` は「`minpoly ℤ θ` が mod `p` で**単純な**根を持てば `p` は完全分解する」
まで取った。★残るのは「**単純**であること」——すなわち
「ほとんどすべての `p` で `minpoly ℤ θ` mod `p` が分離的」である。

★★中身は **Bézout の等式を整数係数で書く**ことに尽きる:

  `u·f + v·f′ = c`  (`u, v ∈ ℤ[X]`、`c ∈ ℤ`、`c ≠ 0`)

★`f` は標数 `0` で既約だから分離的、すなわち `ℚ[X]` で `f` と `f′` は互いに素。
その Bézout 係数の**分母を払う**と上の形になる(`IsLocalization.integerNormalization`)。
★★あとは `p ∤ c` を課せば、mod `p` でも `c` が可逆なので互いに素性が残る。

## ★★★これで迂回路が閉じる

| 段 | 内容 | 場所 |
|---|---|---|
| 1 | `f` の値を割る素数は無限個(Schur) | `PrimeDivisorsOfValues.lean` |
| 2 | 根があれば剰余次数 1 の素がある(Dedekind–Kummer) | `SplitInfinite.lean` |
| 3 | ★**ほとんどすべての `p` で mod `p` は分離的** | ★**本ファイル** |
| 4 | Galois なら完全分解 | `SplitInfinite.lean` |

★★**Chebotarev の密度定理は 1 度も使わない。**

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `exists_bezout_int` | Bézout の等式を整数係数で |
| `separable_map_of_not_dvd` | `p ∤ c` なら mod `p` は分離的 |
| `multiplicity_X_sub_C_eq_one` | 分離的なら根の重複度は 1 |
| `infinite_splitsCompletely` | ★★★★**完全分解する素数は無限個** |
-/

namespace ABC3.Found.NF

open _root_.NumberField Polynomial Ideal
open scoped _root_.NumberField

/-! ## ★1. Bézout の等式を整数係数で -/

/-- ★★★**分離的な整数係数多項式には整数係数の Bézout 等式がある**。

  `u·f + v·f′ = C c`  (`c ≠ 0`)

★`ℚ[X]` の Bézout 係数の分母を `IsLocalization.integerNormalization` で払うだけ。 -/
theorem exists_bezout_int (f : ℤ[X])
    (hsep : (f.map (Int.castRingHom ℚ)).Separable) :
    ∃ (u v : ℤ[X]) (c : ℤ), c ≠ 0 ∧ u * f + v * derivative f = C c := by
  obtain ⟨A, B, hAB⟩ := hsep
  obtain ⟨b₁, hb₁, hA⟩ :=
    IsLocalization.integerNormalization_spec (nonZeroDivisors ℤ) (S := ℚ) A
  obtain ⟨b₂, hb₂, hB⟩ :=
    IsLocalization.integerNormalization_spec (nonZeroDivisors ℤ) (S := ℚ) B
  set A' := IsLocalization.integerNormalization (nonZeroDivisors ℤ) (S := ℚ) A with hA'def
  set B' := IsLocalization.integerNormalization (nonZeroDivisors ℤ) (S := ℚ) B with hB'def
  have hcast : (algebraMap ℤ ℚ : ℤ →+* ℚ) = Int.castRingHom ℚ := Subsingleton.elim _ _
  rw [hcast] at hA hB
  refine ⟨C b₂ * A', C b₁ * B', b₁ * b₂,
    mul_ne_zero (nonZeroDivisors.ne_zero hb₁) (nonZeroDivisors.ne_zero hb₂), ?_⟩
  refine Polynomial.map_injective (Int.castRingHom ℚ) Int.cast_injective ?_
  have hd : Polynomial.map (Int.castRingHom ℚ) (derivative f)
      = derivative (Polynomial.map (Int.castRingHom ℚ) f) := (derivative_map f _).symm
  simp only [Polynomial.map_add, Polynomial.map_mul, hd, hA, hB,
    eq_intCast, Polynomial.map_intCast, Int.cast_mul]
  linear_combination ((b₁ : ℚ[X]) * (b₂ : ℚ[X])) * hAB

/-! ## ★2. `p ∤ c` なら mod `p` は分離的 -/

/-- ★★**整数係数の Bézout 等式は、`p ∤ c` のとき mod `p` で分離性を与える**。 -/
theorem separable_map_of_not_dvd {f u v : ℤ[X]} {c : ℤ}
    (h : u * f + v * derivative f = C c) {p : ℕ} [Fact p.Prime]
    (hc : ¬ (p : ℤ) ∣ c) : (f.map (Int.castRingHom (ZMod p))).Separable := by
  have hcz : ((c : ZMod p)) ≠ 0 := by
    simpa [ZMod.intCast_zmod_eq_zero_iff_dvd] using hc
  have hmap := congrArg (Polynomial.map (Int.castRingHom (ZMod p))) h
  simp only [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C,
    ← derivative_map] at hmap
  refine ⟨C ((c : ZMod p))⁻¹ * u.map (Int.castRingHom (ZMod p)),
    C ((c : ZMod p))⁻¹ * v.map (Int.castRingHom (ZMod p)), ?_⟩
  rw [mul_assoc, mul_assoc, ← mul_add, hmap, ← C_mul,
    eq_intCast (Int.castRingHom (ZMod p)) c, inv_mul_cancel₀ hcz, map_one]

/-! ## ★3. 分離的なら根の重複度は 1 -/

/-- ★**squarefree な多項式の根の重複度は 1**。 -/
theorem multiplicity_X_sub_C_eq_one {R : Type*} [CommRing R] [IsDomain R]
    {g : R[X]} (hsq : Squarefree g) {a : R} (ha : g.IsRoot a) :
    multiplicity (X - C a) g = 1 := by
  refine multiplicity_eq_of_dvd_of_not_dvd ?_ ?_
  · simpa using dvd_iff_isRoot.mpr ha
  · intro hdvd
    have : IsUnit (X - C a) := by
      refine hsq (X - C a) ?_
      rw [← sq]
      exact hdvd
    exact not_isUnit_X_sub_C a this

/-! ## ★4. 完全分解する素数は無限個 -/

variable {K : Type*} [Field K] [NumberField K]

/-- ★素数の有限性 —— 零でない整数を割る自然数は有限個。 -/
theorem finite_primes_dvd {c : ℤ} (hc : c ≠ 0) : {p : ℕ | (p : ℤ) ∣ c}.Finite := by
  classical
  refine Set.Finite.subset (Set.finite_Icc 0 c.natAbs) ?_
  intro p hp
  have hp' : p ∣ c.natAbs := by
    have := Int.natAbs_dvd_natAbs.mpr hp
    simpa using this
  refine ⟨Nat.zero_le _, Nat.le_of_dvd ?_ hp'⟩
  simpa [Int.natAbs_pos] using hc

/-- ★★★★★★**Galois 拡大で完全分解する素数は無限個ある**(Chebotarev の迂回路)。

★★**Chebotarev の密度定理は 1 度も使わない** ——
Schur の補題(値を割る素数は無限個)＋ Dedekind–Kummer ＋ Galois の基本等式だけ。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) = -/
theorem infinite_splitsCompletely [IsGalois ℚ K] (θ : 𝓞 K)
    (hdeg : 0 < (minpoly ℤ θ).natDegree)
    (hsep : ((minpoly ℤ θ).map (Int.castRingHom ℚ)).Separable)
    (hexp : _root_.RingOfIntegers.exponent θ ≠ 0) :
    {p : ℕ | p.Prime ∧
      (primesOver (span {(p : ℤ)}) (𝓞 K)).ncard = Module.finrank ℚ K}.Infinite := by
  classical
  obtain ⟨u, v, c, hc0, hbez⟩ := exists_bezout_int (minpoly ℤ θ) hsep
  have hSchur := _root_.ABC3.Found.NumberField.infinite_primes_dvd_eval (minpoly ℤ θ) hdeg
  have hbad : ({p : ℕ | (p : ℤ) ∣ c} ∪ {p : ℕ | p ∣ _root_.RingOfIntegers.exponent θ}).Finite := by
    refine Set.Finite.union (finite_primes_dvd hc0) ?_
    refine Set.Finite.subset (Set.finite_Icc 0 (_root_.RingOfIntegers.exponent θ)) ?_
    intro q hq
    exact ⟨Nat.zero_le _, Nat.le_of_dvd (Nat.pos_of_ne_zero hexp) hq⟩
  refine ((hSchur.sdiff hbad).mono ?_)
  rintro q ⟨⟨hqp, n, hn⟩, hqbad⟩
  simp only [Set.mem_union, Set.mem_setOf_eq, not_or] at hqbad
  haveI : Fact q.Prime := ⟨hqp⟩
  refine ⟨hqp, ?_⟩
  have hroot : ((minpoly ℤ θ).map (Int.castRingHom (ZMod q))).IsRoot ((n : ZMod q)) := by
    have h0 : ((((minpoly ℤ θ).eval n : ℤ) : ZMod q)) = 0 := by
      simpa [ZMod.intCast_zmod_eq_zero_iff_dvd] using hn
    have hev : eval ((Int.castRingHom (ZMod q)) n)
        ((minpoly ℤ θ).map (Int.castRingHom (ZMod q)))
        = (Int.castRingHom (ZMod q)) ((minpoly ℤ θ).eval n) := by
      rw [eval_map, eval₂_at_apply]
    show eval ((n : ZMod q)) ((minpoly ℤ θ).map (Int.castRingHom (ZMod q))) = 0
    rw [show ((n : ZMod q)) = (Int.castRingHom (ZMod q)) n from rfl, hev]
    simpa using h0
  have hsq : Squarefree ((minpoly ℤ θ).map (Int.castRingHom (ZMod q))) :=
    (separable_map_of_not_dvd hbez hqbad.1).squarefree
  exact splitsCompletely_of_isRoot θ hqbad.2 hroot
    (multiplicity_X_sub_C_eq_one hsq hroot)

/-! ## ★5. 仮定の回収 —— 分離性と次数は自動で出る -/

/-- ★**代数的整数の最小多項式は `ℚ` 上で分離的**（標数 `0`）。 -/
theorem minpoly_int_map_separable (θ : 𝓞 K) :
    ((minpoly ℤ θ).map (Int.castRingHom ℚ)).Separable := by
  have h1 : minpoly ℚ ((θ : K)) = (minpoly ℤ θ).map (algebraMap ℤ ℚ) :=
    minpoly.isIntegrallyClosed_eq_field_fractions ℚ K θ.isIntegral
  have h2 : (minpoly ℚ ((θ : K))).Separable := Algebra.IsSeparable.isSeparable ℚ (θ : K)
  rw [h1] at h2
  have hcast : (algebraMap ℤ ℚ : ℤ →+* ℚ) = Int.castRingHom ℚ := Subsingleton.elim _ _
  rwa [hcast] at h2

/-- ★★★★★★**完全分解する素数は無限個**（仮定を回収した形）。

★残る仮定は `exponent θ ≠ 0`（すなわち `ℤ[θ] ⊆ 𝓞 K` が有限指数）だけで、
★これは **`θ` が原始元である**ことに相当する。 -/
theorem infinite_splitsCompletely' [IsGalois ℚ K] (θ : 𝓞 K)
    (hexp : _root_.RingOfIntegers.exponent θ ≠ 0) :
    {p : ℕ | p.Prime ∧
      (primesOver (span {(p : ℤ)}) (𝓞 K)).ncard = Module.finrank ℚ K}.Infinite :=
  infinite_splitsCompletely θ (minpoly.natDegree_pos θ.isIntegral)
    (minpoly_int_map_separable θ) hexp

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Theorem 6.4, (iv)` の (a)(完全分解する素数の存在)。 -/
def infinite_splitsCompletely.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素数は無限個(Chebotarev を使わない迂回路)",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
