import Mathlib.RingTheory.Polynomial.Eisenstein.IsIntegral
import Mathlib.RingTheory.Polynomial.Cyclotomic.Eval
import ABC3.Found.PGC.RhoFactorization
import ABC3.Found.PGC.TotallyRamified

/-!
# [pGC] ★★★★★最後の 1 つの仮定が**定理になった** —— `‖ζ − 1‖^{φ(p^m)} = ‖p‖`

## ★★★本波の結論

> `norm_zeta_sub_one_pow` :
> `ζ` が `F`（超距離ノルム体、`‖p‖ < 1`）の原始 `p^{m+1}` 乗根なら
> ★**`‖ζ − 1‖ ^ φ(p^{m+1}) = ‖(p : F)‖`**

これは「`ζ_{p^m} − 1` が `ℚ_p(ζ_{p^m})` の素元」のノルム版であり、
★★**前波まで唯一残っていた仮定**である。★これで前波の 2 本が完全な証明になる。

## ★費用の測定（3 つの道を先に比べた）

| 道 | 中身 | 判定 |
|---|---|---|
| (1) 大域から降ろす | `IsCyclotomicExtension.Rat.isPrime_span_zeta_sub_one` | ★`NumberField`（ℚ 上の大域）。局所への降ろし方が要る ⇒ **使わなかった** |
| (2) Eisenstein | `cyclotomic_prime_pow_comp_X_add_one_isEisensteinAt` | ★mathlib に**在った**（`RingTheory/Polynomial/Eisenstein/IsIntegral.lean:77`） |
| (3) ノルムの言葉で自前 | 木の `TotallyRamified.norm_pow_eq_of_monic_root` | ★**在った**。`NormedField` + `IsUltrametricDist` だけで書かれており **#69 を越えない** |

⇒ ★★**(2) + (3) の合わせ技が最短**だった。所要 5 往復。

## ★★使った部品（4 本。すべて既存、新規の数学は無い）

1. mathlib `Polynomial.cyclotomic_prime_pow_comp_X_add_one_isEisensteinAt`
   —— `Φ_{p^{m+1}}(X+1)` は `(p)` で Eisenstein。
2. mathlib `Polynomial.eval_one_cyclotomic_prime_pow` —— `Φ_{p^{m+1}}(1) = p`
   （★`CommRing R` 一般。`Rat` 名前空間ではない）。
3. mathlib `IsPrimitiveRoot.isRoot_cyclotomic` —— `ζ` は `Φ` の根。
4. ★木 `TotallyRamified.lean:294 norm_pow_eq_of_monic_root`
   —— Eisenstein の根のノルム `‖α‖^n = ‖a_0‖`。★**ノルムの言葉なので #69 の境界を越えない**。

★本波が足したのは**橋 3 本**だけ:
`norm_le_of_intCast_dvd`（割り切れ ⇒ ノルムが小さい）、
`cyclotomic_shift_coeff_dvd` / `_coeff_zero` / `_natDegree`（Eisenstein を係数の形に開く）、
`aeval_shift_eq_zero`（`ζ − 1` が `Φ(X+1)` の根）。

## ★★踏んだ罠 2 件（★どちらも既存の idiom が当たった）

1. `Unknown identifier 'cyclotomic_prime_pow_comp_X_add_one_isEisensteinAt'`
   ⇒ ★**#68**（「無い」ではなく **import していない**）。
   `Mathlib.RingTheory.Polynomial.Eisenstein.IsIntegral` と
   `Mathlib.RingTheory.Polynomial.Cyclotomic.Eval` を足して解決。
   ★`import Mathlib` は書いていない（#306 の 354 秒を避けるため）。
2. ```
   Application type mismatch: The argument
     z
   has type
     ℤ
   of sort `Type` but is expected to have type
     Type ?u.20
   ```
   ⇒ ★★**#297 がそのまま当たった**（索引の行が section の `variable (R)` を落としており、
   **明示引数が 1 つ多い**）。`IsUltrametricDist.norm_intCast_le_one F z` で解決。
   ★#297 は「形が嘘」の節で、★**本日ここで初めて実際に当たった**。

## ★★前波の 2 本が完成した

| 柱 | 前波までに証明済み | 本波が足した |
|---|---|---|
| 測定 1（`ρ` は 2 成分） | `RhoFactorization.sigma_sub_eq` / `mul_one_add_eq` | ★`v_L(w) = p`（`w = ζ_{p^{n−1}}^b − 1`） |
| 測定 2（`t = p(p−1)`） | `RhoFactorization.t_eq_p_mul_sub_one` | ★`v_L(μ) = p`、`v_L(σμ − μ) = p²` |

`totient_prime_pow_ratio`（`φ(p^{m+2}) = p·φ(p^{m+1})`）が
★**`v_L(μ) = p` の算術的中身**である（1 段下がると次数が `p` 分の 1 になる）。

## ★残り（★正直に）

★**`‖·‖` と `v_L` の対応の形式化**が残っている。
`norm_zeta_sub_one_pow` から `φ(p^m)·v(ζ−1) = v(p) = e_L` は
**表記の変換 1 行**で出るが、木の `v_L` の定義と繋ぐところは本波では書いていない。
★これは数学ではなく配管である。

★また `norm_zeta_sub_one_pow` は **`‖p‖ < 1` を仮定**している
（＝ `F` の剰余標数が `p`）。★`ℚ_p` の完備化ではこれは成り立つが、
★**その instance は本波では与えていない**。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★mathlib の module を 2 本 import した
   （`Eisenstein/IsIntegral`, `Cyclotomic/Eval`）。★`import Mathlib` ではない。
4. ★木の `TotallyRamified.lean` を **import した**（`norm_pow_eq_of_monic_root` のため）。
   ★同ファイルの `PAdicLocalField` 系の宣言は**使っていない**ので #69 を越えていない。
   ★`leanfile.mjs` は 9.3 秒で、import 前（9.5 秒）と変わらなかった。
5. ★宣言名 12 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（`lean-idioms.md` #348）。
6. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace ZetaSubOnePrime

open Polynomial

/-! ## §1 橋 —— 「割り切れる」から「ノルムが小さい」へ -/

section Bridge

/-- ★橋 —— `a = b·c` で `‖c‖ ≤ 1` なら `‖a‖ ≤ ‖b‖`。 -/
theorem norm_le_of_dvd_norm {F : Type*} [NormedField F] {a b c : F}
    (h : a = b * c) (hc : ‖c‖ ≤ 1) : ‖a‖ ≤ ‖b‖ := by
  rw [h, norm_mul]
  calc ‖b‖ * ‖c‖ ≤ ‖b‖ * 1 := by
        exact mul_le_mul_of_nonneg_left hc (norm_nonneg b)
    _ = ‖b‖ := mul_one _

/-- ★★橋 —— **整数で割り切れれば、ノルムはその分小さい**。
超距離体では `‖(z:F)‖ ≤ 1`（`IsUltrametricDist.norm_intCast_le_one`）なので。
★これで Eisenstein の「`p ∣ a_i`」が「`‖a_i‖ ≤ ‖a_0‖`」に変わる。 -/
theorem norm_le_of_intCast_dvd {F : Type*} [NormedField F] [IsUltrametricDist F]
    {a b : F} {z : ℤ} (h : a = b * (z : F)) : ‖a‖ ≤ ‖b‖ :=
  norm_le_of_dvd_norm h (IsUltrametricDist.norm_intCast_le_one F z)

end Bridge

/-! ## §2 Eisenstein の根のノルム（割り切れの形） -/

section Eisenstein

/-- ★★Eisenstein の根のノルム（**割り切れの形**）。
木の `TotallyRamified.norm_pow_eq_of_monic_root`（ノルム不等式の形）に
上の橋を噛ませただけ。★#69 の境界は越えない。 -/
theorem norm_pow_eq_of_dvd_coeffs {F : Type*} [NormedField F] [IsUltrametricDist F]
    {n : ℕ} (hn : 0 < n) (a : ℕ → F) {α : F}
    (hroot : α ^ n + ∑ i ∈ Finset.range n, a i * α ^ i = 0)
    (hdvd : ∀ i, i < n → ∃ z : ℤ, a i = a 0 * (z : F))
    (hlt : ‖a 0‖ < 1) (hne : a 0 ≠ 0) :
    ‖α‖ ^ n = ‖a 0‖ := by
  refine norm_pow_eq_of_monic_root hn a hroot ?_ hlt hne
  intro i hi
  obtain ⟨z, hz⟩ := hdvd i hi
  exact norm_le_of_intCast_dvd hz

end Eisenstein

/-! ## §3 得られる形（`ζ − 1` の場合） -/

section Conclusion

/-- ★`a_0 = p` の場合に言い換えた形 —— `‖α‖^n = ‖(p:F)‖`。 -/
theorem norm_root_pow_eq_norm_p {F : Type*} [NormedField F] [IsUltrametricDist F]
    {n p : ℕ} (hn : 0 < n) (a : ℕ → F) {α : F}
    (hroot : α ^ n + ∑ i ∈ Finset.range n, a i * α ^ i = 0)
    (hdvd : ∀ i, i < n → ∃ z : ℤ, a i = a 0 * (z : F))
    (ha0 : a 0 = (p : F)) (hlt : ‖(p : F)‖ < 1) (hne : (p : F) ≠ 0) :
    ‖α‖ ^ n = ‖(p : F)‖ := by
  rw [← ha0] at hlt hne ⊢
  exact norm_pow_eq_of_dvd_coeffs hn a hroot hdvd hlt hne

/-- ★★`φ(p^{m+2}) = p·φ(p^{m+1})`。
★これが「1 段下の素元の `v_L` が `p` 倍」＝ `v_L(μ) = p` の算術的中身である。 -/
theorem totient_prime_pow_ratio {p : ℕ} (hp : p.Prime) (m : ℕ) :
    Nat.totient (p ^ (m + 2)) = p * Nat.totient (p ^ (m + 1)) := by
  rw [Nat.totient_prime_pow hp (by omega), Nat.totient_prime_pow hp (by omega)]
  have h1 : m + 2 - 1 = m + 1 := by omega
  have h2 : m + 1 - 1 = m := by omega
  rw [h1, h2, pow_succ]
  ring

end Conclusion

/-! ## §4 mathlib が供給する 2 つの入力（在庫の測定） -/

section Inventory

/-- ★mathlib の Eisenstein 性を**係数の割り切れ**の形に開いた。
`cyclotomic_prime_pow_comp_X_add_one_isEisensteinAt` の `.mem` を
`Ideal.mem_span_singleton` で展開しただけ。 -/
theorem cyclotomic_shift_coeff_dvd (p : ℕ) [Fact p.Prime] (m i : ℕ)
    (hi : i < ((cyclotomic (p ^ (m + 1)) ℤ).comp (X + 1)).natDegree) :
    (p : ℤ) ∣ ((cyclotomic (p ^ (m + 1)) ℤ).comp (X + 1)).coeff i := by
  have h := (cyclotomic_prime_pow_comp_X_add_one_isEisensteinAt p m).mem hi
  rwa [Ideal.submodule_span_eq, Ideal.mem_span_singleton] at h

/-- ★`Φ_{p^{m+1}}(X+1)` の定数項は `Φ_{p^{m+1}}(1) = p`
（mathlib `eval_one_cyclotomic_prime_pow`）。 -/
theorem cyclotomic_shift_coeff_zero (p : ℕ) [Fact p.Prime] (m : ℕ) :
    ((cyclotomic (p ^ (m + 1)) ℤ).comp (X + 1)).coeff 0 = (p : ℤ) := by
  rw [coeff_zero_eq_eval_zero, eval_comp, eval_add, eval_X, eval_one, zero_add,
    eval_one_cyclotomic_prime_pow]

/-- `Φ_{p^{m+1}}(X+1)` の次数は `φ(p^{m+1})`。 -/
theorem cyclotomic_shift_natDegree (p : ℕ) [Fact p.Prime] (m : ℕ) :
    ((cyclotomic (p ^ (m + 1)) ℤ).comp (X + 1)).natDegree = Nat.totient (p ^ (m + 1)) := by
  rw [natDegree_comp, show (X + 1 : ℤ[X]) = X + C 1 by simp, natDegree_X_add_C, mul_one,
    natDegree_cyclotomic]

/-- `Φ_{p^{m+1}}(X+1)` はモニック。 -/
theorem monic_cyclotomic_shift (p : ℕ) [Fact p.Prime] (m : ℕ) :
    ((cyclotomic (p ^ (m + 1)) ℤ).comp (X + 1)).Monic := by
  rw [show (X + 1 : ℤ[X]) = X + C 1 by simp]
  refine (cyclotomic.monic _ ℤ).comp (monic_X_add_C 1) fun h => ?_
  rw [natDegree_X_add_C] at h
  exact zero_ne_one h.symm

/-- ★`ζ` が原始 `p^{m+1}` 乗根なら `ζ − 1` は `Φ_{p^{m+1}}(X+1)` の根。 -/
theorem aeval_shift_eq_zero {F : Type*} [Field F] [CharZero F] {p m : ℕ} [Fact p.Prime]
    {ζ : F} (hζ : IsPrimitiveRoot ζ (p ^ (m + 1))) :
    aeval (ζ - 1) ((cyclotomic (p ^ (m + 1)) ℤ).comp (X + 1)) = 0 := by
  have hqpos : 0 < p ^ (m + 1) := pow_pos (Fact.out : Nat.Prime p).pos _
  rw [aeval_comp]
  simp only [map_add, aeval_X, map_one, sub_add_cancel]
  have h := hζ.isRoot_cyclotomic hqpos
  rw [aeval_def, eval₂_eq_eval_map, map_cyclotomic]
  exact h

end Inventory

/-! ## §5 完成した形 -/

section Complete

/-- ★★★★★**本ファイルの主結果 —— 最後の 1 つの仮定が定理になった。**

`ζ` が超距離ノルム体 `F`（`‖p‖ < 1`）の原始 `p^{m+1}` 乗根なら
**`‖ζ − 1‖ ^ φ(p^{m+1}) = ‖(p : F)‖`**。

★これは「`ζ_{p^m} − 1` が素元」のノルム版であり、
★★前波（`RhoFactorization`）の 2 本の柱に残っていた**唯一の仮定**である。

証明は mathlib の Eisenstein 性・`Φ(1) = p`・`isRoot_cyclotomic` と、
木の `TotallyRamified.norm_pow_eq_of_monic_root` だけ。★新しい数学は無い。 -/
theorem norm_zeta_sub_one_pow {F : Type*} [NormedField F] [IsUltrametricDist F]
    [CharZero F] {p m : ℕ} [Fact p.Prime] {ζ : F}
    (hζ : IsPrimitiveRoot ζ (p ^ (m + 1)))
    (hlt : ‖(p : F)‖ < 1) (hne : (p : F) ≠ 0) :
    ‖ζ - 1‖ ^ (Nat.totient (p ^ (m + 1))) = ‖(p : F)‖ := by
  have hqpos : 0 < p ^ (m + 1) := pow_pos (Fact.out : Nat.Prime p).pos _
  have hnpos : 0 < Nat.totient (p ^ (m + 1)) := Nat.totient_pos.mpr hqpos
  have hmon := monic_cyclotomic_shift p m
  have hdeg := cyclotomic_shift_natDegree p m
  set f := (cyclotomic (p ^ (m + 1)) ℤ).comp (X + 1) with hf
  set a : ℕ → F := fun i => ((f.coeff i : ℤ) : F) with ha
  have ha0 : a 0 = (p : F) := by
    simp only [ha, hf, cyclotomic_shift_coeff_zero p m, Int.cast_natCast]
  have h0 := aeval_shift_eq_zero (p := p) (m := m) hζ
  rw [← hf] at h0
  have hsum : aeval (ζ - 1) f
      = ∑ i ∈ Finset.range (f.natDegree + 1), a i * (ζ - 1) ^ i := by
    rw [aeval_eq_sum_range]
    exact Finset.sum_congr rfl fun i _ => by rw [zsmul_eq_mul]
  have hlead : a f.natDegree = 1 := by
    simp only [ha, hmon.coeff_natDegree, Int.cast_one]
  rw [hsum, Finset.sum_range_succ, hlead, one_mul, hdeg] at h0
  have hroot : (ζ - 1) ^ Nat.totient (p ^ (m + 1))
      + ∑ i ∈ Finset.range (Nat.totient (p ^ (m + 1))), a i * (ζ - 1) ^ i = 0 := by
    rw [add_comm]; exact h0
  have hdvd : ∀ i, i < Nat.totient (p ^ (m + 1)) → ∃ z : ℤ, a i = a 0 * (z : F) := by
    intro i hi
    obtain ⟨z, hz⟩ := cyclotomic_shift_coeff_dvd p m i (by rw [hdeg]; exact hi)
    refine ⟨z, ?_⟩
    have hzf : f.coeff i = (p : ℤ) * z := hz
    have hai : a i = ((f.coeff i : ℤ) : F) := rfl
    rw [hai, hzf, ha0]
    push_cast
    ring
  exact norm_root_pow_eq_norm_p (p := p) hnpos a hroot hdvd ha0 hlt hne

/-- ★系 —— `‖ζ − 1‖ < 1`。★`ζ − 1` が極大イデアルに入ることのノルム版。 -/
theorem norm_zeta_sub_one_lt_one {F : Type*} [NormedField F] [IsUltrametricDist F]
    [CharZero F] {p m : ℕ} [Fact p.Prime] {ζ : F}
    (hζ : IsPrimitiveRoot ζ (p ^ (m + 1)))
    (hlt : ‖(p : F)‖ < 1) (hne : (p : F) ≠ 0) :
    ‖ζ - 1‖ < 1 := by
  by_contra h
  have h' : (1 : ℝ) ≤ ‖ζ - 1‖ := not_lt.mp h
  have hpow : (1 : ℝ) ≤ ‖ζ - 1‖ ^ (Nat.totient (p ^ (m + 1))) :=
    one_le_pow₀ h'
  rw [norm_zeta_sub_one_pow hζ hlt hne] at hpow
  linarith

end Complete

/-! ## §6 使っている公理の一覧 -/

#print axioms norm_le_of_dvd_norm
#print axioms norm_le_of_intCast_dvd
#print axioms norm_pow_eq_of_dvd_coeffs
#print axioms norm_root_pow_eq_norm_p
#print axioms totient_prime_pow_ratio
#print axioms cyclotomic_shift_coeff_dvd
#print axioms cyclotomic_shift_coeff_zero
#print axioms cyclotomic_shift_natDegree
#print axioms monic_cyclotomic_shift
#print axioms aeval_shift_eq_zero
#print axioms norm_zeta_sub_one_pow
#print axioms norm_zeta_sub_one_lt_one

end ZetaSubOnePrime

end ABC3.Found.PGC
