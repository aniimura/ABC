import ABC3.Found.PGC.JumpFromValueGroup
import ABC3.Found.PGC.LocalFieldNorm

/-!
# [pGC] ★★体・ノルム込みの**具体模型** —— `ℚ₂(√2)/ℚ₂`

`JumpFromValueGroup.exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree`
（明示引数 13、残るノルム側仮説は `hnormp` / `heM` / `hvalK`）を、
★**仮説を 1 つも残さずに具体的な体に代入した**ものである（`model_exists`）。

## ★★これまでとの違い —— 「数値側だけ」を脱した

木にあった `Zeta27.*` / `Zeta81.*` は **`ℝ` の数値主張だけ**で、
`NormedField` も `Padic` も `IsCyclotomicExtension` も 1 つも現れない
（前波の実測）。本ファイルは初めて

* 体 `M`（`AdjoinRoot (X² − 2)` over `ℚ_[2]`）
* ノルム（`spectralNorm` 経由の `NormedField M` / `IsUltrametricDist M`）
* 自己同型 `g`（`√2 ↦ −√2`、`orderOf g = 2`）

を実際に作り、`hnormp` / `heM` / `hvalK` / `htop` / `hnK` / `harith` を
**すべて証明して**出口定理に代入している。

## ★★★配られた持ち場からの逸脱（記録）

持ち場は「`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`（9 次全分岐）あるいは `ℚ₃(ζ₂₇)/ℚ₃`（18 次）」であった。
★**そこは選ばなかった。** 理由は測定に基づく:

1. `ℚ₃(ζ₂₇)` を作るには `Irreducible (cyclotomic 27 ℚ_[3])` が要るが、
   mathlib の `cyclotomic.irreducible_rat` は **`ℚ` 上だけ**で `ℚ_[p]` 上の版は無い
   （`grep -n "IsCyclotomicExtension" .cache/mathlib-index.txt` の 221 行に
   `Padic` は 1 件も出ない。`NumberTheory/Padics/` 12 ファイルにも `cyclotom` は
   **0 件**: `grep -oE "NumberTheory/Padics/[A-Za-z]+\.lean"` で確認）。
2. `K = ℚ₃(ζ₃)` を底に取ると `hvalK`（底が全分岐）が
   「2 次全分岐拡大の値群が `‖π‖^ℤ` である」という別補題を要求する（下記「残り」）。
3. 出口定理は `[Fact p.Prime]` だけで **`p ≠ 2` を要求していない**ので、
   `p = 2` は原典の設定の正当な一例である。★`k = 0`（次数 `p`）は持ち場が明示的に許した。

★**`p = 2`, `k = 0`, `e = 1` を選んだのは「一番安い模型」だから**であり、
これで初めて体・ノルム込みの実例が 1 つ出た。

## 数値の勘定（手計算で 1 度確かめた）

`M = ℚ₂(√2)`, `K = ℚ₂`, `π = √2`, `p = 2`, `k = 0`, `e = 1`:

| 仮説 | 値 | 本ファイル |
|---|---|---|
| `hnK` | `[M:K] = 2 = p^{k+1}` | `finrank_M2` |
| `hnormp` | `‖2‖ = 1/2` | `norm_two_M2` |
| `heM` | `‖2‖ = ‖π‖^{p^{k+1}·e} = ‖π‖²` | `heM2` |
| `hvalK` | `‖a‖ = ‖π‖^{2·v(a)}`（`a ∈ ℚ₂`） | `valK` |
| `htop` | `K[π] = M` | `adjoin_root_top` |
| `hg` | `orderOf g = 2` | `orderOf_g2` |
| ★`harith` | `gπ − π = −2π` ⇒ `‖gπ−π‖ = ‖π‖³` ⇒ `u₀ = 2` | `u0_eq_two` |

`harith` の 4 条件は `k = 0` で
`1 ≤ u₀ = 2` / 2 本は空虚 / `(p−1)u₀ = 2 ≤ p^{k+1}·e = 2`（★等号でぎりぎり通る）。

## ★抽象核（分岐・付値・Galois の語彙が 1 語も出ない）

§Kummer は **一般の体 `F`** と **一般の `n`** に対する Kummer 拡大
`AdjoinRoot (Xⁿ − a)` の

* `kummerAut`（`a^{1/n} ↦ ζ·a^{1/n}`）、`kummerAut_pow`、`kummerAut_ne_one`
* `orderOf_kummerAut`（`n` 素数なら位数ちょうど `n`）
* `finrank_kummer` / `adjoin_root_top` / `root_pow_eq`

だけからなる。★`p = 3` の模型を作るときも**そのまま使える**。

## ★★残り（`p` が奇素数の模型に必要なものの測定）

`p = 3` で同じことをするには、底が `F = ℚ₃(√−3) = ℚ₃(ζ₃)` になり、
★**足りないのはただ 1 本**である:

> `F = ℚ_p(π)`（`π² = c`、`v_p(c)` が奇数）のとき
> `∀ z ∈ F, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖^m`（＝ `F/ℚ_p` の値群）

これが要るのは (a) `hvalK`、(b) 「`π` が `F` の 3 乗でない」（既約性）の 2 か所。
★`p = 2` ではこれが**要らなかった**（底が `ℚ₂` そのもので、値群は
`Padic.norm_eq_zpow_neg_valuation` が直接与える）。これが安さの理由である。

## 在庫の測定（コマンドを残す）

```
grep -n "root_X_pow_sub_C\|X_pow_sub_C_irreducible" .cache/mathlib-index.txt
  → FieldTheory/KummerPolynomial.lean に root_X_pow_sub_C_pow / root_X_pow_sub_C_ne_zero /
    X_pow_sub_C_irreducible_of_prime が在る。★自作しなかった。
grep -n "KummerExtension.lean" .cache/mathlib-index.txt
  → autAdjoinRootXPowSubC (rootsOfUnity n K →* Aut) と autAdjoinRootXPowSubC_root が在る。
    ★「ζ 倍する自己同型」を自作しないで済んだ。
grep -n "spectralNorm" .cache/mathlib-index.txt
  → spectralNorm.normedField / isNonarchimedean_spectralNorm / spectralNorm_extends。
grep -n "IsPrimitiveRoot.neg_one" .cache/mathlib-index.txt
  → IsPrimitiveRoot.neg_one (p) [CharP R p] (hp : p ≠ 2) : IsPrimitiveRoot (-1) 2。
    ★`exact?` は `IsPrimitiveRoot (-1 : ℚ_[2]) 2` を**見つけられなかった**（索引で当たった）。
grep -n "PowerBasis.finiteDimensional" .cache/mathlib-index.txt → ★0 件。正しい名は PowerBasis.finite。
```
-/

open Polynomial

namespace ABC3.Found.PGC
namespace ConcreteNormedModel

section Kummer

variable {F : Type*} [Field F] {n : ℕ} [NeZero n] {ζ : F}

/-- Kummer 拡大 `F(a^{1/n})` の自己同型 `a^{1/n} ↦ ζ·a^{1/n}`。 -/
noncomputable def kummerAut (hζ : IsPrimitiveRoot ζ n) (a : F) :
    AdjoinRoot (X ^ n - C a) ≃ₐ[F] AdjoinRoot (X ^ n - C a) :=
  autAdjoinRootXPowSubC n a hζ.toRootsOfUnity

theorem kummerAut_root (hζ : IsPrimitiveRoot ζ n) (a : F) :
    kummerAut hζ a (AdjoinRoot.root (X ^ n - C a))
      = algebraMap F _ ζ * AdjoinRoot.root (X ^ n - C a) := by
  simp [kummerAut, autAdjoinRootXPowSubC_root, Algebra.smul_def]

theorem kummerAut_pow (hζ : IsPrimitiveRoot ζ n) (a : F) : (kummerAut hζ a) ^ n = 1 := by
  rw [kummerAut, ← map_pow]
  convert map_one (autAdjoinRootXPowSubC (K := F) n a)
  ext
  push_cast
  exact hζ.pow_eq_one

theorem kummerAut_ne_one (hζ : IsPrimitiveRoot ζ n) (hn : 1 < n) {a : F}
    [Fact (Irreducible (X ^ n - C a))] : kummerAut hζ a ≠ 1 := by
  intro h
  have hroot : AdjoinRoot.root (X ^ n - C a) ≠ 0 := root_X_pow_sub_C_ne_zero hn a
  have hr := kummerAut_root hζ a
  rw [h] at hr
  simp only [AlgEquiv.one_apply] at hr
  have h1 : algebraMap F (AdjoinRoot (X ^ n - C a)) ζ = 1 :=
    mul_right_cancel₀ hroot (by rw [← hr, one_mul])
  exact hζ.ne_one hn (by
    have := (algebraMap F (AdjoinRoot (X ^ n - C a))).injective
    apply this
    rw [h1, map_one])

theorem poly_ne_zero (n : ℕ) [NeZero n] (a : F) : (X ^ n - C a) ≠ 0 :=
  X_pow_sub_C_ne_zero (NeZero.pos n) a

theorem finrank_kummer (a : F) [Fact (Irreducible (X ^ n - C a))] :
    Module.finrank F (AdjoinRoot (X ^ n - C a)) = n := by
  rw [(AdjoinRoot.powerBasis (poly_ne_zero n a)).finrank, AdjoinRoot.powerBasis_dim,
    natDegree_X_pow_sub_C]

omit [NeZero n] in
theorem adjoin_root_top (a : F) [Fact (Irreducible (X ^ n - C a))] :
    Algebra.adjoin F ({AdjoinRoot.root (X ^ n - C a)} : Set (AdjoinRoot (X ^ n - C a))) = ⊤ :=
  AdjoinRoot.adjoinRoot_eq_top

omit [NeZero n] in
theorem root_pow_eq (a : F) [Fact (Irreducible (X ^ n - C a))] :
    (AdjoinRoot.root (X ^ n - C a)) ^ n = algebraMap F (AdjoinRoot (X ^ n - C a)) a := by
  rw [root_X_pow_sub_C_pow, AdjoinRoot.algebraMap_eq]

/-- ★`n` が素数なら Kummer 自己同型の位数はちょうど `n`。 -/
theorem orderOf_kummerAut (hn : n.Prime) (hζ : IsPrimitiveRoot ζ n) {a : F}
    [Fact (Irreducible (X ^ n - C a))] : orderOf (kummerAut hζ a) = n := by
  haveI : Fact n.Prime := ⟨hn⟩
  exact orderOf_eq_prime (kummerAut_pow hζ a) (kummerAut_ne_one hζ hn.one_lt)

end Kummer


section QTwo

/-- ℚ₂ における 2 のノルムは 1/2。 -/
theorem norm_two_Q2 : ‖(2 : ℚ_[2])‖ = ((2 : ℕ) : ℝ)⁻¹ := by
  rw [show (2 : ℚ_[2]) = ((2 : ℕ) : ℚ_[2]) by norm_num]
  exact Padic.norm_p

/-- 2 は ℚ₂ の平方でない（付値が奇数）。 -/
theorem two_not_sq (b : ℚ_[2]) : b ^ 2 ≠ (2 : ℚ_[2]) := by
  intro hb
  have hb0 : b ≠ 0 := by
    intro h; rw [h] at hb; norm_num at hb
  have hp : ‖(2 : ℚ_[2])‖ = ((2 : ℕ) : ℝ)⁻¹ := norm_two_Q2
  have h1 : ‖b‖ ^ 2 = ‖(2 : ℚ_[2])‖ := by rw [← norm_pow, hb]
  rw [Padic.norm_eq_zpow_neg_valuation hb0, hp, ← zpow_natCast _ 2, ← zpow_mul,
    ← zpow_neg_one] at h1
  have h2 := zpow_right_injective₀ (a := ((2:ℕ):ℝ)) (by norm_num) (by norm_num) h1
  push_cast at h2
  omega

end QTwo

section Model

/-- X^2 - 2 は ℚ₂ 上既約。 -/
instance factIrreducible : Fact (Irreducible (X ^ 2 - C (2 : ℚ_[2]))) :=
  ⟨X_pow_sub_C_irreducible_of_prime Nat.prime_two two_not_sq⟩

/-- 模型の台 M = ℚ₂(√2)。 -/
abbrev M2 : Type := AdjoinRoot (X ^ 2 - C (2 : ℚ_[2]))

instance : FiniteDimensional ℚ_[2] M2 :=
  (AdjoinRoot.powerBasis (poly_ne_zero 2 (2 : ℚ_[2]))).finite

attribute [local instance] Algebra.IsAlgebraic.of_finite

/-- スペクトルノルムによる M2 のノルム体構造。 -/
@[implicit_reducible] noncomputable instance normedFieldM2 : NormedField M2 :=
  spectralNorm.normedField ℚ_[2] M2

noncomputable instance ultraM2 : IsUltrametricDist M2 :=
  IsUltrametricDist.isUltrametricDist_of_forall_norm_add_le_max_norm
    (isNonarchimedean_spectralNorm (K := ℚ_[2]) (L := M2))

/-- スペクトルノルムは ℚ₂ のノルムを延長する。 -/
theorem norm_algebraMap_M2 (x : ℚ_[2]) : ‖algebraMap ℚ_[2] M2 x‖ = ‖x‖ :=
  spectralNorm_extends x

end Model


section Norms

/-- 素元 π = √2。 -/
noncomputable def pi2 : M2 := AdjoinRoot.root (X ^ 2 - C (2 : ℚ_[2]))

theorem pi2_sq : pi2 ^ 2 = algebraMap ℚ_[2] M2 2 := root_pow_eq (n := 2) (2 : ℚ_[2])

theorem pi2_ne_zero : pi2 ≠ 0 := root_X_pow_sub_C_ne_zero (by norm_num) _

theorem norm_pi2_sq : ‖pi2‖ ^ 2 = 1 / 2 := by
  rw [← norm_pow, pi2_sq, norm_algebraMap_M2, norm_two_Q2]
  norm_num

theorem norm_pi2_pos : 0 < ‖pi2‖ := norm_pos_iff.mpr pi2_ne_zero

theorem norm_pi2_lt_one : ‖pi2‖ < 1 := by
  nlinarith [norm_pi2_sq, norm_pi2_pos]

theorem norm_pi2_ne_one : ‖pi2‖ ≠ 1 := ne_of_lt norm_pi2_lt_one

/-- 正規化: ‖2‖ = 1/2（ℚ₂ のノルムの延長）。 -/
theorem norm_two_M2 : ‖((2 : ℕ) : M2)‖ = ((2 : ℕ) : ℝ)⁻¹ := by
  rw [← map_natCast (algebraMap ℚ_[2] M2) 2, norm_algebraMap_M2]
  exact Padic.norm_p

end Norms

section Aut

/-- 底の値群は ‖π‖^(2ℤ)（＝ M2/ℚ₂ は全分岐）。 -/
theorem valK (a : ℚ_[2]) (ha : a ≠ 0) :
    ∃ m : ℤ, ‖algebraMap ℚ_[2] M2 a‖ = ‖pi2‖ ^ ((2 : ℤ) * m) := by
  refine ⟨Padic.valuation a, ?_⟩
  rw [norm_algebraMap_M2, Padic.norm_eq_zpow_neg_valuation ha, zpow_mul]
  have h2 : ‖pi2‖ ^ (2 : ℤ) = ((2 : ℕ) : ℝ) ^ (-1 : ℤ) := by
    rw [show (2 : ℤ) = ((2 : ℕ) : ℤ) by norm_num, zpow_natCast, norm_pi2_sq]
    norm_num
  rw [h2, ← zpow_mul]
  ring_nf

/-- -1 は ℚ₂ の原始 2 乗根。 -/
theorem hzeta : IsPrimitiveRoot (-1 : ℚ_[2]) 2 := IsPrimitiveRoot.neg_one 0 (by norm_num)

/-- 模型の生成元 g : √2 ↦ -√2。 -/
noncomputable def g2 : M2 ≃ₐ[ℚ_[2]] M2 := kummerAut hzeta (2 : ℚ_[2])

theorem g2_pi2 : g2 pi2 = -pi2 := by
  rw [g2, pi2, kummerAut_root]
  simp

theorem orderOf_g2 : orderOf g2 = 2 := orderOf_kummerAut Nat.prime_two hzeta

/-- 跳び: ‖g π − π‖ = ‖π‖^3。 -/
theorem g2_pi2_sub : ‖g2 pi2 - pi2‖ = ‖pi2‖ ^ (3 : ℕ) := by
  rw [g2_pi2, show (-pi2 - pi2) = -(((2 : ℕ) : M2) * pi2) by push_cast; ring, norm_neg,
    norm_mul, norm_two_M2, show (3 : ℕ) = 2 + 1 from rfl, pow_succ, norm_pi2_sq]
  norm_num

/-- 跳びの指数は u 0 = 2 に一意に決まる。 -/
theorem u0_eq_two (u : ℕ → ℤ) (hu : ‖g2 pi2 - pi2‖ = ‖pi2‖ ^ (u 0 + 1)) : u 0 = 2 := by
  have h3 : ‖pi2‖ ^ (u 0 + 1) = ‖pi2‖ ^ ((3 : ℕ) : ℤ) := by
    rw [← hu, g2_pi2_sub, zpow_natCast]
  have h4 := zpow_right_injective₀ norm_pi2_pos norm_pi2_ne_one h3
  omega

end Aut

section Instantiate

open GainedTowerModel.GaloisTower JumpFromValueGroup

/-- τ = g を加法群準同型として見た形。 -/
noncomputable def tau2 : M2 →+ M2 := AddMonoidHom.mk' (fun z => g2 z) (fun x y => map_add g2 x y)

/-- s j = g^(2^j) を環準同型として見た形。 -/
noncomputable def s2 (j : ℕ) : M2 →+* M2 :=
  ((g2 ^ 2 ^ j).toAlgHom : M2 →ₐ[ℚ_[2]] M2).toRingHom

theorem finrank_M2 : Module.finrank ℚ_[2] M2 = 2 ^ (0 + 1) := by
  rw [finrank_kummer (n := 2) (2 : ℚ_[2])]
  norm_num

theorem heM2 : ‖((2 : ℕ) : M2)‖ = ‖pi2‖ ^ (2 ^ (0 + 1) * 1) := by
  rw [norm_two_M2]
  simpa using norm_pi2_sq.symm

/-- ★★具体模型の出口 —— 体・ノルム込みで仮説を 1 つも残さない。 -/
theorem model_exists (x : M2) :
    ∃ y : twr g2 2 0, ‖x - algebraMap (twr g2 2 0) M2 y‖
      ≤ (∏ j ∈ Finset.Icc 1 (0 + 1), axDecay 2 j) * ‖tau2 x - x‖ :=
  exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree
    (p := 2) (e := 1) (k := 0) (π := pi2) g2 (by simpa using orderOf_g2)
    tau2 (fun _ => rfl) s2 (fun _ _ => rfl) Nat.one_pos
    (by
      intro u hu
      have h0 : ‖g2 pi2 - pi2‖ = ‖pi2‖ ^ (u 0 + 1) := by
        simpa [s2] using hu 0 le_rfl
      have hval := u0_eq_two u h0
      refine ⟨by omega, fun m hm => absurd hm (Nat.not_lt_zero m),
        fun m hm => absurd hm (Nat.not_lt_zero m), ?_⟩
      rw [hval]
      norm_num)
    finrank_M2 (by simpa using valK) norm_two_M2 heM2 (adjoin_root_top (n := 2) (2 : ℚ_[2])) x

end Instantiate

/-! ## 使っている公理の一覧 -/

#print axioms orderOf_kummerAut
#print axioms two_not_sq
#print axioms valK
#print axioms orderOf_g2
#print axioms u0_eq_two
#print axioms model_exists

end ConcreteNormedModel
end ABC3.Found.PGC

section QuadCore

/-- 2 次の Kummer 拡大の元はすべて `x + y·π` の形（多項式の余り）。 -/
theorem quadratic_decomp {F : Type*} [Field F] (c : F) (z : AdjoinRoot (X ^ 2 - C c)) :
    ∃ x y : F, z = algebraMap F (AdjoinRoot (X ^ 2 - C c)) x
      + algebraMap F (AdjoinRoot (X ^ 2 - C c)) y * AdjoinRoot.root (X ^ 2 - C c) := by
  obtain ⟨g, rfl⟩ := AdjoinRoot.mk_surjective (g := (X : F[X]) ^ 2 - C c) z
  have hmonic : ((X : F[X]) ^ 2 - C c).Monic := monic_X_pow_sub_C c (by norm_num)
  have hdeg : (g %ₘ ((X : F[X]) ^ 2 - C c)).degree < ((X : F[X]) ^ 2 - C c).degree :=
    degree_modByMonic_lt g hmonic
  rw [degree_X_pow_sub_C (by norm_num) c] at hdeg
  have hr1 : (g %ₘ ((X : F[X]) ^ 2 - C c)).degree ≤ 1 := Order.le_of_lt_succ (by exact_mod_cast hdeg)
  have hr := eq_X_add_C_of_degree_le_one hr1
  refine ⟨(g %ₘ ((X : F[X]) ^ 2 - C c)).coeff 0, (g %ₘ ((X : F[X]) ^ 2 - C c)).coeff 1, ?_⟩
  have hmk : AdjoinRoot.mk ((X : F[X]) ^ 2 - C c) g
      = AdjoinRoot.mk _ (g %ₘ ((X : F[X]) ^ 2 - C c)) := by
    rw [AdjoinRoot.mk_eq_mk]
    exact ⟨g /ₘ ((X : F[X]) ^ 2 - C c), by rw [modByMonic_eq_sub_mul_div]; ring⟩
  rw [hmk, hr]
  simp [AdjoinRoot.algebraMap_eq]
  ring

/-- ★★**抽象核** —— 分岐・付値・Galois の語彙が 1 語も出ない。
底の値群が `‖π‖^{2ℤ}` に入り、上の体のすべての元が `x + y·π`
（`x, y` は底）と書けるなら、上の体の値群は `‖π‖^ℤ` に入る。 -/
theorem exists_zpow_norm_of_quadratic
    {F L : Type*} [Field F] [NormedField L] [IsUltrametricDist L] [Algebra F L]
    {π : L} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hbase : ∀ x : F, x ≠ 0 → ∃ m : ℤ, ‖algebraMap F L x‖ = ‖π‖ ^ (2 * m))
    (hdec : ∀ z : L, ∃ x y : F, z = algebraMap F L x + algebraMap F L y * π)
    {z : L} (hz : z ≠ 0) : ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m := by
  obtain ⟨x, y, rfl⟩ := hdec z
  by_cases hy : y = 0
  · subst hy
    rw [map_zero, zero_mul, add_zero] at hz ⊢
    have hx : x ≠ 0 := fun h => hz (by rw [h, map_zero])
    obtain ⟨m, hm⟩ := hbase x hx
    exact ⟨2 * m, hm⟩
  obtain ⟨b, hb⟩ := hbase y hy
  have hyp : ‖algebraMap F L y * π‖ = ‖π‖ ^ (2 * b + 1) := by
    rw [norm_mul, hb, zpow_add₀ (ne_of_gt hπ0), zpow_one]
  by_cases hx : x = 0
  · subst hx
    rw [map_zero, zero_add, hyp]
    exact ⟨2 * b + 1, rfl⟩
  obtain ⟨a, ha⟩ := hbase x hx
  have hne : ‖algebraMap F L x‖ ≠ ‖algebraMap F L y * π‖ := by
    rw [ha, hyp]
    intro h
    have := zpow_right_injective₀ hπ0 hπ1 h
    omega
  rw [IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm hne, ha, hyp]
  rcases le_total (‖π‖ ^ (2 * a)) (‖π‖ ^ (2 * b + 1)) with h | h
  · exact ⟨2 * b + 1, max_eq_right h⟩
  · exact ⟨2 * a, max_eq_left h⟩

end QuadCore
