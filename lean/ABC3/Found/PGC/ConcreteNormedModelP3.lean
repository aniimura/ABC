import ABC3.Found.PGC.ConcreteNormedModel

/-!
# [pGC] ★★★`p = 3` の具体模型 —— `ℚ₃(ζ₃)((√−3)^{1/3}) / ℚ₃(ζ₃)`（3 次巡回・全分岐）

`ConcreteNormedModel`（`p = 2`）に続く 2 本目。★**奇素数**で、しかも
★**底が `ℚ_p` ではなく `ℚ₃(ζ₃)`**（`e = 2`）という、持ち場が求めていた形である。
`JumpFromValueGroup.exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree`
に**仮説を 1 つも残さず**代入した（`model3_exists`）。

## 構成（2 段の `AdjoinRoot` ＋ 2 段の `spectralNorm`）

```
ℚ_[3]  --(X² + 3 が既約: v(−3)=1 は奇数)-->  F3 = ℚ₃(√−3) = ℚ₃(ζ₃)
       spectralNorm ℚ_[3] F3 で NormedField / IsUltrametricDist
F3     --(X³ − π が既約: π は F3 の 3 乗でない)-->  M3 = F3(π^{1/3})
       spectralNorm F3 M3 で NormedField / IsUltrametricDist
```

* `ζ₃ = (−1 + √−3)/2 ∈ F3`（`zeta3_prim : IsPrimitiveRoot zeta3 3`）
* `g = kummerAut zeta3_prim π`（`α ↦ ζ₃·α`、`orderOf g = 3`）
* `π = √−3`, `α = π^{1/3}`, `v(α) = 1`, `v(π) = 3`, `v(3) = 6`

## ★★★★これが可能になった理由 —— 抽象核 1 本

`p = 2` の版で「残り」として測定した

> `F = ℚ_p(π)`（`π² = c`、`v(c)` 奇数）の値群は `‖π‖^ℤ`

は、`ConcreteNormedModel.exists_zpow_norm_of_quadratic`
（★分岐・付値・Galois の語彙が 1 語も出ない抽象核）＋ `quadratic_decomp`
で埋まった。本ファイルはそれを `val_F3` として使うだけである。
★値群が要るのは 2 か所:(a) `hvalK`、(b)「`π` が `F3` の 3 乗でない」(`piF_not_cube`)。

## ★★原典との差（記録すべき逸脱）

持ち場が名指ししたのは `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`（9 次、`k = 2`）だが、
本ファイルは `k = 0`（3 次）で、しかも上の体は **`ℚ₃(ζ₉)` ではない**。
★2 つは**跳びが違う**（どちらも `F3` 上 3 次巡回・全分岐だが別の拡大）:

| 模型 | 素元 | `‖gπ − π‖` | 跳び `u₀` |
|---|---|---|---|
| `ℚ₃(ζ₉)/ℚ₃(ζ₃)` | `λ = ζ₉ − 1` | `v(ζ₉) + v(ζ₃−1) = 0 + 3` | ★`2` |
| ★本ファイル `F3(π^{1/3})` | `α` | `v(ζ₃−1) + v(α) = 3 + 1` | ★`3` |

★**本ファイルの方が `harith` の上界にぴったり乗る**:
`(p−1)·u₀ = 2·3 = 6` が `p^{k+1}·e = 3·2 = 6` と**等号**になる
（`ℚ₃(ζ₉)` 側は `2·2 = 4 ≤ 6` で余裕がある）。
⇒ ★出口定理の定数勘定を**境界で 1 度検算した**ことになる。

`ℚ₃(ζ₉)` そのものを作らなかったのは、`Irreducible (cyclotomic 9 ℚ_[3])` が
mathlib に無い（`NumberTheory/Padics/` 12 ファイルに `cyclotom` は 0 件）ため。
★Kummer 側（`FieldTheory/KummerPolynomial.lean` / `KummerExtension.lean`）は
既約性・自己同型・位数がすべて揃っているので、そちらを通した。

## 在庫の測定（コマンドを残す）

```
grep -n "IsPrimitiveRoot.iff_orderOf" .cache/mathlib-index.txt
  → IsPrimitiveRoot ζ k ↔ orderOf ζ = k。★これで orderOf_eq_prime に落とせた。
grep -n "IsPrimitiveRoot.neg_one" .cache/mathlib-index.txt  （p=2 側で使用）
grep -n "spectralNorm.completeSpace" .cache/mathlib-index.txt
  → FiniteDimensional なら CompleteSpace。★2 段目の spectralNorm（底が F3）に要る。
  本ファイルは `example : CompleteSpace F3 := by infer_instance` で実測した
  （ProperSpace F3 ← FiniteDimensional.proper ℚ_[3] F3 経由で通る）。
```

★**測ったが要らなかったもの**: `‖(2 : ℚ_[3])‖ = 1`。
`ζ₃ − 1 = (−3 + π)/2` から直接ノルムを出すと `‖2‖` が要るが、
`(ζ₃−1)² = −3·ζ₃`（`zeta3_sub_one_sq`）を使うと **2 が現れない**。
★2 乗の等式に直すと分母が消える、という手（本セッションで 1 回効いた）。
-/

open Polynomial

namespace ABC3.Found.PGC
namespace ConcreteNormedModel

section Base

variable {p : ℕ} [Fact p.Prime]

/-- 付値 1 の元は ℚ_p の平方でない（一般形）。 -/
theorem not_sq_of_norm_eq_inv_p (c : ℚ_[p]) (hc : ‖c‖ = ((p : ℕ) : ℝ) ^ (-1 : ℤ)) (b : ℚ_[p]) :
    b ^ 2 ≠ c := by
  have hp1 : (1 : ℝ) < ((p : ℕ) : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  intro hb
  have hb0 : b ≠ 0 := by
    intro h
    rw [h] at hb
    rw [← hb] at hc
    simp at hc
    nlinarith [hc, zpow_pos (by linarith : (0:ℝ) < ((p:ℕ):ℝ)) (-1 : ℤ)]
  have h1 : ‖b‖ ^ 2 = ‖c‖ := by rw [← norm_pow, hb]
  rw [Padic.norm_eq_zpow_neg_valuation hb0, hc, ← zpow_natCast _ 2, ← zpow_mul] at h1
  have h2 := zpow_right_injective₀ (a := ((p:ℕ):ℝ)) (by linarith) (by linarith) h1
  push_cast at h2
  omega

end Base

section F3

instance factPrimeThree : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩

theorem norm_neg_three : ‖(-3 : ℚ_[3])‖ = ((3 : ℕ) : ℝ) ^ (-1 : ℤ) := by
  rw [norm_neg, show (3 : ℚ_[3]) = ((3 : ℕ) : ℚ_[3]) by norm_num, Padic.norm_p]
  norm_num

instance factIrrF3 : Fact (Irreducible (X ^ 2 - C (-3 : ℚ_[3]))) :=
  ⟨X_pow_sub_C_irreducible_of_prime Nat.prime_two
    (not_sq_of_norm_eq_inv_p (-3 : ℚ_[3]) norm_neg_three)⟩

/-- 底の体 `F3 = ℚ₃(√−3) = ℚ₃(ζ₃)`。 -/
abbrev F3 : Type := AdjoinRoot (X ^ 2 - C (-3 : ℚ_[3]))

instance : FiniteDimensional ℚ_[3] F3 :=
  (AdjoinRoot.powerBasis (poly_ne_zero 2 (-3 : ℚ_[3]))).finite

attribute [local instance] Algebra.IsAlgebraic.of_finite

@[implicit_reducible] noncomputable instance normedFieldF3 : NormedField F3 :=
  spectralNorm.normedField ℚ_[3] F3

noncomputable instance ultraF3 : IsUltrametricDist F3 :=
  IsUltrametricDist.isUltrametricDist_of_forall_norm_add_le_max_norm
    (isNonarchimedean_spectralNorm (K := ℚ_[3]) (L := F3))

theorem norm_algebraMap_F3 (x : ℚ_[3]) : ‖algebraMap ℚ_[3] F3 x‖ = ‖x‖ :=
  spectralNorm_extends x

/-- `F3` の素元 `π = √−3`。 -/
noncomputable def piF : F3 := AdjoinRoot.root (X ^ 2 - C (-3 : ℚ_[3]))

theorem piF_sq : piF ^ 2 = algebraMap ℚ_[3] F3 (-3) := root_pow_eq (n := 2) (-3 : ℚ_[3])

theorem piF_ne_zero : piF ≠ 0 := root_X_pow_sub_C_ne_zero (by norm_num) _

theorem norm_piF_sq : ‖piF‖ ^ 2 = 1 / 3 := by
  rw [← norm_pow, piF_sq, norm_algebraMap_F3, norm_neg_three]
  norm_num

theorem norm_piF_pos : 0 < ‖piF‖ := norm_pos_iff.mpr piF_ne_zero

theorem norm_piF_lt_one : ‖piF‖ < 1 := by nlinarith [norm_piF_sq, norm_piF_pos]

theorem norm_piF_ne_one : ‖piF‖ ≠ 1 := ne_of_lt norm_piF_lt_one

end F3

section ValueGroup

/-- 底 ℚ₃ の値群は ‖π‖^(2ℤ)。 -/
theorem base_val_F3 (x : ℚ_[3]) (hx : x ≠ 0) :
    ∃ m : ℤ, ‖algebraMap ℚ_[3] F3 x‖ = ‖piF‖ ^ ((2 : ℤ) * m) := by
  refine ⟨Padic.valuation x, ?_⟩
  rw [norm_algebraMap_F3, Padic.norm_eq_zpow_neg_valuation hx, zpow_mul]
  have h2 : ‖piF‖ ^ (2 : ℤ) = ((3 : ℕ) : ℝ) ^ (-1 : ℤ) := by
    rw [show (2 : ℤ) = ((2 : ℕ) : ℤ) by norm_num, zpow_natCast, norm_piF_sq]
    norm_num
  rw [h2, ← zpow_mul]
  ring_nf

theorem decomp_F3 (z : F3) :
    ∃ x y : ℚ_[3], z = algebraMap ℚ_[3] F3 x + algebraMap ℚ_[3] F3 y * piF :=
  quadratic_decomp (-3 : ℚ_[3]) z

/-- ★F3 の値群は ‖π‖^ℤ（＝ F3/ℚ₃ は全分岐）。抽象核 `exists_zpow_norm_of_quadratic` の適用。 -/
theorem val_F3 {z : F3} (hz : z ≠ 0) : ∃ m : ℤ, ‖z‖ = ‖piF‖ ^ m :=
  exists_zpow_norm_of_quadratic norm_piF_pos norm_piF_ne_one base_val_F3 decomp_F3 hz

/-- ★π は F3 の 3 乗でない（値群の指数が 3 で割れない）。 -/
theorem piF_not_cube (b : F3) : b ^ 3 ≠ piF := by
  intro hb
  have hb0 : b ≠ 0 := by
    intro h
    rw [h] at hb
    exact piF_ne_zero (by simpa using hb.symm)
  obtain ⟨m, hm⟩ := val_F3 hb0
  have h1 : ‖piF‖ ^ (m * 3) = ‖piF‖ ^ (1 : ℤ) := by
    rw [zpow_mul, ← hm, show (3 : ℤ) = ((3 : ℕ) : ℤ) by norm_num, zpow_natCast, ← norm_pow, hb,
      zpow_one]
  have h2 := zpow_right_injective₀ norm_piF_pos norm_piF_ne_one h1
  omega

end ValueGroup

section Zeta

theorem algebraMap_neg_three : algebraMap ℚ_[3] F3 (-3) = -3 := by
  rw [map_neg, map_ofNat]

theorem piF_sq_eq : piF ^ 2 = -3 := by rw [piF_sq, algebraMap_neg_three]

theorem two_ne_zero_F3 : (2 : F3) ≠ 0 := by
  intro h
  have h1 : algebraMap ℚ_[3] F3 (2 : ℚ_[3]) = algebraMap ℚ_[3] F3 0 := by
    rw [map_ofNat, map_zero]; exact h
  have h2 := (algebraMap ℚ_[3] F3).injective h1
  norm_num at h2

theorem three_ne_zero_F3 : (3 : F3) ≠ 0 := by
  intro h
  have h1 : algebraMap ℚ_[3] F3 (3 : ℚ_[3]) = algebraMap ℚ_[3] F3 0 := by
    rw [map_ofNat, map_zero]; exact h
  have h2 := (algebraMap ℚ_[3] F3).injective h1
  norm_num at h2

/-- ζ₃ = (−1 + √−3)/2 ∈ F3。 -/
noncomputable def zeta3 : F3 := (-1 + piF) / 2

theorem zeta3_sq_add_one : zeta3 ^ 2 + zeta3 + 1 = 0 := by
  rw [zeta3]
  have h2 := two_ne_zero_F3
  field_simp
  linear_combination piF_sq_eq

theorem zeta3_cube : zeta3 ^ 3 = 1 := by linear_combination (zeta3 - 1) * zeta3_sq_add_one

theorem zeta3_ne_one : zeta3 ≠ 1 := by
  intro h
  apply three_ne_zero_F3
  have h2 := zeta3_sq_add_one
  rw [h] at h2
  linear_combination h2

/-- ★ζ₃ は F3 の原始 3 乗根。 -/
theorem zeta3_prim : IsPrimitiveRoot zeta3 3 :=
  IsPrimitiveRoot.iff_orderOf.mpr (orderOf_eq_prime zeta3_cube zeta3_ne_one)

end Zeta

section TopField

attribute [local instance] Algebra.IsAlgebraic.of_finite

@[implicit_reducible] noncomputable instance normedAlgebraF3 : NormedAlgebra ℚ_[3] F3 :=
  spectralNorm.normedAlgebra ℚ_[3] F3

@[implicit_reducible] noncomputable instance nontriviallyF3 : NontriviallyNormedField F3 where
  toNormedField := normedFieldF3
  non_trivial := by
    refine ⟨algebraMap ℚ_[3] F3 ((3 : ℚ_[3])⁻¹), ?_⟩
    rw [norm_algebraMap_F3, norm_inv, show (3 : ℚ_[3]) = ((3 : ℕ) : ℚ_[3]) by norm_num,
      Padic.norm_p, inv_inv]
    norm_num

noncomputable instance properF3 : ProperSpace F3 := FiniteDimensional.proper ℚ_[3] F3

example : CompleteSpace F3 := by infer_instance

instance factIrrM3 : Fact (Irreducible (X ^ 3 - C piF)) :=
  ⟨X_pow_sub_C_irreducible_of_prime Nat.prime_three piF_not_cube⟩

/-- ★模型の台 `M3 = F3(π^{1/3}) = ℚ₃(ζ₃)(√−3^{1/3})`（F3 上 3 次巡回・全分岐）。 -/
abbrev M3 : Type := AdjoinRoot (X ^ 3 - C piF)

instance : FiniteDimensional F3 M3 := (AdjoinRoot.powerBasis (poly_ne_zero 3 piF)).finite

@[implicit_reducible] noncomputable instance normedFieldM3 : NormedField M3 :=
  spectralNorm.normedField F3 M3

noncomputable instance ultraM3 : IsUltrametricDist M3 :=
  IsUltrametricDist.isUltrametricDist_of_forall_norm_add_le_max_norm
    (isNonarchimedean_spectralNorm (K := F3) (L := M3))

theorem norm_algebraMap_M3 (x : F3) : ‖algebraMap F3 M3 x‖ = ‖x‖ := spectralNorm_extends x

end TopField

section Alpha

/-- `M3` の素元 `α = π^{1/3}`。 -/
noncomputable def alpha : M3 := AdjoinRoot.root (X ^ 3 - C piF)

theorem alpha_cube : alpha ^ 3 = algebraMap F3 M3 piF := root_pow_eq (n := 3) piF

theorem alpha_ne_zero : alpha ≠ 0 := root_X_pow_sub_C_ne_zero (by norm_num) _

theorem norm_alpha_cube : ‖alpha‖ ^ 3 = ‖piF‖ := by
  rw [← norm_pow, alpha_cube, norm_algebraMap_M3]

theorem norm_alpha_pos : 0 < ‖alpha‖ := norm_pos_iff.mpr alpha_ne_zero

theorem norm_alpha_lt_one : ‖alpha‖ < 1 := by
  rcases lt_or_ge ‖alpha‖ 1 with h | h
  · exact h
  · exfalso
    have h3 : (1 : ℝ) ≤ ‖alpha‖ ^ 3 := one_le_pow₀ h
    rw [norm_alpha_cube] at h3
    linarith [norm_piF_lt_one]

theorem norm_alpha_ne_one : ‖alpha‖ ≠ 1 := ne_of_lt norm_alpha_lt_one

/-- 正規化: ‖3‖ = 1/3。 -/
theorem norm_three_F3 : ‖((3 : ℕ) : F3)‖ = ((3 : ℕ) : ℝ)⁻¹ := by
  rw [← map_natCast (algebraMap ℚ_[3] F3) 3, norm_algebraMap_F3]
  exact Padic.norm_p

theorem norm_three_M3 : ‖((3 : ℕ) : M3)‖ = ((3 : ℕ) : ℝ)⁻¹ := by
  rw [← map_natCast (algebraMap F3 M3) 3, norm_algebraMap_M3]
  exact norm_three_F3

/-- `heM`: ‖3‖ = ‖α‖^{p^{k+1}·e} = ‖α‖^6。 -/
theorem heM3 : ‖((3 : ℕ) : M3)‖ = ‖alpha‖ ^ (3 ^ (0 + 1) * 2) := by
  rw [norm_three_M3, show (3 : ℕ) ^ (0 + 1) * 2 = 3 * 2 from rfl, pow_mul, norm_alpha_cube]
  have h := norm_piF_sq
  norm_num at h ⊢
  linarith [h]

/-- `hvalK`: 底 F3 の値群は ‖α‖^{3ℤ}。 -/
theorem valK3 (z : F3) (hz : z ≠ 0) :
    ∃ m : ℤ, ‖algebraMap F3 M3 z‖ = ‖alpha‖ ^ ((3 : ℤ) * m) := by
  obtain ⟨m, hm⟩ := val_F3 hz
  refine ⟨m, ?_⟩
  rw [norm_algebraMap_M3, hm, zpow_mul,
    show (3 : ℤ) = ((3 : ℕ) : ℤ) by norm_num, zpow_natCast, norm_alpha_cube]

theorem finrank_M3 : Module.finrank F3 M3 = 3 ^ (0 + 1) := by
  rw [finrank_kummer (n := 3) piF]
  norm_num

end Alpha

section Aut3

open GainedTowerModel.GaloisTower JumpFromValueGroup

/-- 模型の生成元 `g : α ↦ ζ₃·α`。 -/
noncomputable def g3 : M3 ≃ₐ[F3] M3 := kummerAut zeta3_prim piF

theorem g3_alpha : g3 alpha = algebraMap F3 M3 zeta3 * alpha := kummerAut_root zeta3_prim piF

theorem orderOf_g3 : orderOf g3 = 3 := orderOf_kummerAut Nat.prime_three zeta3_prim

theorem norm_zeta3 : ‖zeta3‖ = 1 := by
  have h : ‖zeta3‖ ^ 3 = 1 := by rw [← norm_pow, zeta3_cube, norm_one]
  have hnn : (0 : ℝ) ≤ ‖zeta3‖ := norm_nonneg _
  rcases lt_trichotomy ‖zeta3‖ 1 with hlt | heq | hgt
  · exact absurd h (by have := pow_lt_one₀ hnn hlt (by norm_num : (3 : ℕ) ≠ 0); linarith)
  · exact heq
  · exact absurd h (by have := one_lt_pow₀ hgt (by norm_num : (3 : ℕ) ≠ 0); linarith)

theorem zeta3_sub_one_sq : (zeta3 - 1) ^ 2 = -3 * zeta3 := by linear_combination zeta3_sq_add_one

/-- ★`ζ₃ − 1` は `π = √−3` と同じノルムを持つ（`(ζ−1)² = −3ζ`）。 -/
theorem norm_zeta3_sub_one : ‖zeta3 - 1‖ = ‖piF‖ := by
  have h3 : ‖(-3 : F3)‖ = 1 / 3 := by
    rw [show (-3 : F3) = -((3 : ℕ) : F3) by push_cast; ring, norm_neg, norm_three_F3]
    norm_num
  have h1 : ‖zeta3 - 1‖ ^ 2 = ‖piF‖ ^ 2 := by
    rw [← norm_pow, zeta3_sub_one_sq, norm_mul, h3, norm_zeta3, mul_one, norm_piF_sq]
  nlinarith [h1, norm_nonneg (zeta3 - 1), norm_nonneg piF]

/-- 跳び: ‖gα − α‖ = ‖α‖^4。 -/
theorem g3_alpha_sub : ‖g3 alpha - alpha‖ = ‖alpha‖ ^ (4 : ℕ) := by
  have h : g3 alpha - alpha = algebraMap F3 M3 (zeta3 - 1) * alpha := by
    rw [map_sub, map_one, g3_alpha]; ring
  rw [h, norm_mul, norm_algebraMap_M3, norm_zeta3_sub_one, ← norm_alpha_cube]
  ring

theorem u0_eq_three (u : ℕ → ℤ) (hu : ‖g3 alpha - alpha‖ = ‖alpha‖ ^ (u 0 + 1)) : u 0 = 3 := by
  have h4 : ‖alpha‖ ^ (u 0 + 1) = ‖alpha‖ ^ ((4 : ℕ) : ℤ) := by
    rw [← hu, g3_alpha_sub, zpow_natCast]
  have h5 := zpow_right_injective₀ norm_alpha_pos norm_alpha_ne_one h4
  omega

noncomputable def tau3 : M3 →+ M3 := AddMonoidHom.mk' (fun z => g3 z) (fun x y => map_add g3 x y)

noncomputable def s3 (j : ℕ) : M3 →+* M3 := ((g3 ^ 3 ^ j).toAlgHom : M3 →ₐ[F3] M3).toRingHom

/-- ★★★`p = 3` の具体模型の出口 —— `ℚ₃(ζ₃)((√−3)^{1/3})/ℚ₃(ζ₃)`（3 次巡回・全分岐、`e = 2`）。 -/
theorem model3_exists (x : M3) :
    ∃ y : twr g3 3 0, ‖x - algebraMap (twr g3 3 0) M3 y‖
      ≤ (∏ j ∈ Finset.Icc 1 (0 + 1), axDecay 3 j) * ‖tau3 x - x‖ :=
  exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree
    (p := 3) (e := 2) (k := 0) (π := alpha) g3 (by simpa using orderOf_g3)
    tau3 (fun _ => rfl) s3 (fun _ _ => rfl) (by norm_num)
    (by
      intro u hu
      have h0 : ‖g3 alpha - alpha‖ = ‖alpha‖ ^ (u 0 + 1) := by
        simpa [s3] using hu 0 le_rfl
      have hval := u0_eq_three u h0
      refine ⟨by omega, fun m hm => absurd hm (Nat.not_lt_zero m),
        fun m hm => absurd hm (Nat.not_lt_zero m), ?_⟩
      rw [hval]
      norm_num)
    finrank_M3 (by simpa using valK3) norm_three_M3 heM3 (adjoin_root_top (n := 3) piF) x

end Aut3

/-! ## 使っている公理の一覧 -/

#print axioms not_sq_of_norm_eq_inv_p
#print axioms val_F3
#print axioms piF_not_cube
#print axioms zeta3_prim
#print axioms norm_zeta3_sub_one
#print axioms u0_eq_three
#print axioms model3_exists

end ConcreteNormedModel
end ABC3.Found.PGC
