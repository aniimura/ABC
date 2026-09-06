import ABC3.Found.PGC.AbelianDecomposition
import ABC3.Found.PGC.UnramifiedSubextension

/-!
# Milne Lemma 4.11 の第 1 段——Frobenius の持ち上げと「巡回部分群 × 惰性」の格子計算

経路 Λ7(`K^ab = K_π · K^ur`)の起点である Milne, *Class Field Theory*,
Lemma 4.11(p.57)の**群論側**をここに置く。体の側(Lemma 4.11 本体)は
`Found/PGC/AbelianSplitUnramified.lean`。

## 原典(Milne, Lemma 4.11)

原文は「L を K の指数 m の有限アーベル拡大、K_m を K の次数 m の不分岐拡大と
するとき、K 上完全分岐なアーベル拡大 L_t が在って L · K_m = L_t · K_m」。
原典の証明は「`Gal(LK_m/K)` は指数 `m` のアーベル群で、`σ|_{K_m}` が Frobenius に
なる `σ` は位数 `m`、ゆえに `Gal(LK_m/K) = ⟨σ⟩ × H`(直積)」と進む。

## ★見立てとの差(位数勘定すら要らなかった)

段取りでは第 3 段を「`IsCompl ⟨σ⟩ N` を**位数の積**で出す」としていたが、
実際には**位数を一度も数えずに済む**。絶対 Galois 群 `Γ_K` の中で

* `S := ⟨σ⟩ ⊔ H_M`(`H_M := (L ⊔ K_m).fixingSubgroup`)
* `S ⊓ H_m = H_M` ⟸ 「`σ^k ∈ H_m ⟹ m ∣ k`」+「`σ^m ∈ H_M`」
* `S ⊔ H_m = ⊤`  ⟸ 「`Γ_K/H_m` は `σ` の像で生成される」

の 2 本だけを示せばよく、どちらも**部分群の元を追うだけ**で終わる
(`zpowers_sup_inf_eq` / `zpowers_sup_eq_top`)。直積分解も
`Nat.card` の積も、入射加群も `Module.Baer` も使っていない。

## 本ファイルの内容

* `zpowers_sup_inf_eq` / `zpowers_sup_eq_top`——上の 2 本(純群論)
* `mul_comm_of_commutator_eq_one` / `commutator_eq_one_of_mul_comm`——交換子の言い換え
* `mul_comm_of_forall_mem_zpowers`——1 元生成なら可換
* `isCyclic_gal_of_isUnramifiedAdjoin`——不分岐なら `Gal(K(z)/K)` は巡回
  (`exists_isCyclic_gal` は「或る生成元について」しか言わないので、一意性
  `adjoin_eq_of_isUnramified` で任意の生成元へ運ぶ)
* `exists_unramified_frobenius_lift`——次数 `m` の不分岐拡大 `K_m` と、
  `Γ_K` の元 `σ` で「`σ^k ∈ Gal(K̄/K_m) ⟺ m ∣ k`」「`⟨σ⟩ ⊔ Gal(K̄/K_m) = ⊤`」
  「`∀ g, g^m ∈ Gal(K̄/K_m)`」「`∀ a b, [a,b] ∈ Gal(K̄/K_m)`」を満たすもの
* `finrank_dvd_of_isUnramified_of_pow_mem`——`Gal(K̄/K(z))` が `m` 乗をすべて
  含み `K(z)/K` が不分岐なら `[K(z):K] ∣ m`(巡回群の指数勘定)

## 逸脱(記録)

* 原典は `Gal(LK_m/K)` という**有限商**の中で議論するが、本ファイルは
  `Γ_K` の中の**開部分群の束**で議論している(商群を作らないという
  `AbelianDecomposition.lean` の設計上の制約に合わせた)。主張は同値。
* 原典の `σ|_{K_m} = Frobenius` は、ここでは「`Gal(K_m/K)`(位数 `m` の巡回群)の
  **生成元**の持ち上げ」として実現している。剰余体の `x ↦ x^q` との一致は
  主張していない(`Found/PGC/UnramifiedFrobenius.lean` の担当)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued
open Pointwise

/-! ## 0. 交換子と 1 元生成(純群論) -/

/-- `[a,b] = 1` なら `a` と `b` は可換。 -/
theorem mul_comm_of_commutator_eq_one {G : Type*} [Group G] {a b : G}
    (h : a * b * a⁻¹ * b⁻¹ = 1) : a * b = b * a := by
  calc a * b = (a * b * a⁻¹ * b⁻¹) * (b * a) := by group
    _ = 1 * (b * a) := by rw [h]
    _ = b * a := one_mul _

/-- 可換なら `[a,b] = 1`。 -/
theorem commutator_eq_one_of_mul_comm {G : Type*} [Group G] {a b : G}
    (h : a * b = b * a) : a * b * a⁻¹ * b⁻¹ = 1 := by
  rw [h]; group

/-- 1 元 `g` の冪で尽くされる群は可換。 -/
theorem mul_comm_of_forall_mem_zpowers {G : Type*} [Group G] {g : G}
    (h : ∀ x : G, x ∈ Subgroup.zpowers g) (a b : G) : a * b = b * a := by
  obtain ⟨i, hi⟩ := h a
  obtain ⟨j, hj⟩ := h b
  have hi' : g ^ i = a := hi
  have hj' : g ^ j = b := hj
  rw [← hi', ← hj', ← zpow_add, ← zpow_add, add_comm]

/-! ## 1. 「巡回部分群 × 正規部分群」の格子計算(位数を数えない)

★ここが段取りとの最大の差。`IsCompl ⟨σ⟩ N` を位数の積で出す代わりに、
`fixedField` を当てる直前の形——`⊓` と `⊔` の値——を直接計算する。 -/

/-- ★**`(⟨σ⟩ ⊔ N) ⊓ M = N`**——`N ◁ Γ`、`N ≤ M`、そして
「`σ` の冪が `M` に入るなら `N` に入る」ならば。

`Subgroup.mul_normal` で `↑(⟨σ⟩ ⊔ N) = ↑⟨σ⟩ * ↑N` と分解し、
`x = σ^k n` から `σ^k = x n⁻¹ ∈ M` を経由するだけ。

退化の自己検査:仮定 `h` を落とすと**偽**(`σ` の位数が `M` の側で
`m` より大きいと `⟨σ⟩ ⊓ M ⊄ N`)。これが Milne Example 4.13 の反例の芯。 -/
theorem zpowers_sup_inf_eq {Γ : Type*} [Group Γ] {N M : Subgroup Γ} [N.Normal] (hNM : N ≤ M)
    {σ : Γ} (h : ∀ k : ℤ, σ ^ k ∈ M → σ ^ k ∈ N) :
    (Subgroup.zpowers σ ⊔ N) ⊓ M = N := by
  refine le_antisymm ?_ (le_inf le_sup_right hNM)
  intro x hx
  obtain ⟨hxS, hxM⟩ := hx
  have hmul : x ∈ ((Subgroup.zpowers σ : Set Γ) * (N : Set Γ)) := by
    rw [← Subgroup.mul_normal]; exact hxS
  obtain ⟨a, ha, n, hn, rfl⟩ := hmul
  obtain ⟨k, rfl⟩ := ha
  have han : σ ^ k ∈ M := by
    have hrw : σ ^ k = (σ ^ k * n) * n⁻¹ := by group
    rw [hrw]
    exact Subgroup.mul_mem M hxM (Subgroup.inv_mem M (hNM hn))
  exact Subgroup.mul_mem N (h k han) hn

/-- ★**`⟨σ⟩ ⊔ M = ⊤`**——`Γ/M` が `σ` の像で生成されるならば。 -/
theorem zpowers_sup_eq_top {Γ : Type*} [Group Γ] {M : Subgroup Γ} {σ : Γ}
    (h : ∀ x : Γ, ∃ k : ℤ, x * (σ ^ k)⁻¹ ∈ M) : Subgroup.zpowers σ ⊔ M = ⊤ := by
  rw [Subgroup.eq_top_iff']
  intro x
  obtain ⟨k, hk⟩ := h x
  have hx : x = (x * (σ ^ k)⁻¹) * σ ^ k := (inv_mul_cancel_right x (σ ^ k)).symm
  rw [hx]
  exact Subgroup.mul_mem _ (le_sup_right (a := Subgroup.zpowers σ) (b := M) hk)
    (le_sup_left (a := Subgroup.zpowers σ) (b := M) (Subgroup.mem_zpowers_iff.mpr ⟨k, rfl⟩))

/-! ## 2. 不分岐拡大の Galois 群は巡回(任意の生成元について) -/

variable {p : ℕ} [Fact p.Prime]

/-- **`K(z)/K` が不分岐なら `Gal(K(z)/K)` は巡回群**。

`exists_isCyclic_gal` は「次数 `n` の不分岐拡大が**或る**生成元 `x` で取れて
`Gal(K(x)/K)` が巡回」しか言わないので、一意性(`adjoin_eq_of_isUnramified`)で
任意の `z` に運ぶ。 -/
theorem isCyclic_gal_of_isUnramifiedAdjoin (K : PAdicLocalField p) (z : K.closure)
    (hz : IsUnramifiedAdjoin K z) :
    IsCyclic ((IntermediateField.adjoin K.carrier ({z} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({z} : Set K.closure))) := by
  obtain ⟨x, hrank, hux, hcyc, -⟩ := exists_isCyclic_gal K
    (Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({z} : Set K.closure)))
    Module.finrank_pos.ne'
  have heq : IntermediateField.adjoin K.carrier ({z} : Set K.closure)
      = IntermediateField.adjoin K.carrier ({x} : Set K.closure) :=
    adjoin_eq_of_isUnramified K z x hz hux hrank.symm
  exact (MulEquiv.isCyclic (AlgEquiv.autCongr (IntermediateField.equivOfEq heq))).mpr hcyc

/-! ## 3. Frobenius の持ち上げ -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**次数 `m` の不分岐拡大 `K_m = K(w)` と、`Γ_K` の中の Frobenius の持ち上げ `σ`**。

`σ` は `Gal(K_m/K)`(位数 `m` の巡回群)の生成元の持ち上げであり、

* `σ^k ∈ Gal(K̄/K_m) ⟺ m ∣ k`(`σ` の `Gal(K_m/K)` での位数がちょうど `m`)
* `⟨σ⟩ ⊔ Gal(K̄/K_m) = ⊤`(`σ` の像が `Gal(K_m/K)` を生成する)
* `∀ g, g^m ∈ Gal(K̄/K_m)`(`Gal(K_m/K)` の位数が `m`)
* `∀ a b, [a,b] ∈ Gal(K̄/K_m)`(`Gal(K_m/K)` は巡回ゆえ可換)

を満たす。★`w` と `σ` は `∃` の内側に閉じ込めてある。

退化の自己検査:`m = 0` を許すと「次数 `0` の不分岐拡大」が無いので**偽**。 -/
theorem exists_unramified_frobenius_lift (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0) :
    ∃ (w : K.closure) (σ : K.absGal),
      IsUnramifiedAdjoin K w ∧
      Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({w} : Set K.closure)) = m ∧
      (∀ k : ℤ, σ ^ k ∈ (IntermediateField.adjoin K.carrier
          ({w} : Set K.closure)).fixingSubgroup ↔ (m : ℤ) ∣ k) ∧
      Subgroup.zpowers σ ⊔ (IntermediateField.adjoin K.carrier
        ({w} : Set K.closure)).fixingSubgroup = ⊤ ∧
      (∀ g : K.absGal, g ^ m ∈ (IntermediateField.adjoin K.carrier
        ({w} : Set K.closure)).fixingSubgroup) ∧
      (∀ a b : K.absGal, a * b * a⁻¹ * b⁻¹ ∈ (IntermediateField.adjoin K.carrier
        ({w} : Set K.closure)).fixingSubgroup) := by
  haveI := isGalois_closure K
  haveI := IsAlgClosure.normal K.carrier K.closure
  obtain ⟨w, hrank, huw, hcyc, hcard⟩ := exists_isCyclic_gal K m hm
  haveI := hcyc
  haveI := normal_of_isUnramifiedAdjoin K w huw
  set Km := IntermediateField.adjoin K.carrier ({w} : Set K.closure) with hKmdef
  set φ := AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (Km : Type _) with hφ
  have hker : φ.ker = Km.fixingSubgroup := ker_restrictNormalHom_eq_fixingSubgroup _
  have hsurj : Function.Surjective φ :=
    AlgEquiv.restrictNormalHom_surjective (F := K.carrier) (K₁ := (Km : Type _)) K.closure
  obtain ⟨g0, hg0⟩ := IsCyclic.exists_generator (α := (Km ≃ₐ[K.carrier] Km))
  have htop0 : Subgroup.zpowers g0 = ⊤ := by rw [Subgroup.eq_top_iff']; exact hg0
  have hord : orderOf g0 = m := by
    rw [← Nat.card_zpowers, htop0, ← hcard]
    exact Nat.card_congr (Subgroup.topEquiv).toEquiv
  obtain ⟨σ, hσ⟩ := hsurj g0
  refine ⟨w, σ, huw, hrank, ?_, ?_, ?_, ?_⟩
  · intro k
    rw [← hker, MonoidHom.mem_ker, map_zpow, hσ, ← orderOf_dvd_iff_zpow_eq_one, hord]
  · refine zpowers_sup_eq_top ?_
    intro x
    obtain ⟨k, hk⟩ := hg0 (φ x)
    have hk' : φ x = g0 ^ k := hk.symm
    exact ⟨k, by rw [← hker, MonoidHom.mem_ker, map_mul, map_inv, map_zpow, hσ, hk',
      mul_inv_cancel]⟩
  · intro g
    rw [← hker, MonoidHom.mem_ker, map_pow, ← hcard]
    exact pow_card_eq_one'
  · intro a b
    have h1 : φ (a * b * a⁻¹ * b⁻¹) = 1 := by
      simp only [map_mul, map_inv]
      exact commutator_eq_one_of_mul_comm (mul_comm_of_forall_mem_zpowers hg0 (φ a) (φ b))
    rw [← hker, MonoidHom.mem_ker]
    exact h1

/-! ## 4. 「`m` 乗が固定部分群に入る不分岐拡大の次数は `m` を割る」 -/

/-- ★**`Gal(K̄/K(z))` が `Γ_K` の `m` 乗をすべて含み、`K(z)/K` が不分岐なら
`[K(z):K] ∣ m`**。

`Gal(K(z)/K)` は巡回(`isCyclic_gal_of_isUnramifiedAdjoin`)で、その生成元 `g₀` は
`Γ_K` の元 `g` の像だから `g₀^m = 1`。ゆえに `[K(z):K] = |Gal| = orderOf g₀ ∣ m`。

退化の自己検査:不分岐性(⟹ 巡回)を落とすと**偽**——指数 `m` の非巡回群
(例えば `(ℤ/2)^2`、`m = 2`)は位数が `m` を割らない。 -/
theorem finrank_dvd_of_isUnramified_of_pow_mem (K : PAdicLocalField p) {z : K.closure}
    (huz : IsUnramifiedAdjoin K z) {m : ℕ}
    (hpow : ∀ g : K.absGal, g ^ m ∈ (IntermediateField.adjoin K.carrier
      ({z} : Set K.closure)).fixingSubgroup) :
    Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({z} : Set K.closure)) ∣ m := by
  haveI := isGalois_closure K
  haveI := IsAlgClosure.normal K.carrier K.closure
  haveI := normal_of_isUnramifiedAdjoin K z huz
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({z} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier _
  haveI : IsGalois K.carrier (IntermediateField.adjoin K.carrier ({z} : Set K.closure)) := ⟨⟩
  haveI := isCyclic_gal_of_isUnramifiedAdjoin K z huz
  set Kz := IntermediateField.adjoin K.carrier ({z} : Set K.closure) with hKz
  set φ := AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (Kz : Type _) with hφ
  have hker : φ.ker = Kz.fixingSubgroup := ker_restrictNormalHom_eq_fixingSubgroup _
  have hsurj : Function.Surjective φ :=
    AlgEquiv.restrictNormalHom_surjective (F := K.carrier) (K₁ := (Kz : Type _)) K.closure
  obtain ⟨g0, hg0⟩ := IsCyclic.exists_generator (α := (Kz ≃ₐ[K.carrier] Kz))
  have htop0 : Subgroup.zpowers g0 = ⊤ := by rw [Subgroup.eq_top_iff']; exact hg0
  obtain ⟨g, hg⟩ := hsurj g0
  have hg0m : g0 ^ m = 1 := by
    rw [← hg, ← map_pow, ← MonoidHom.mem_ker, hker]; exact hpow g
  have hcard : Nat.card (Kz ≃ₐ[K.carrier] Kz) = orderOf g0 := by
    rw [← Nat.card_zpowers, htop0]
    exact (Nat.card_congr (Subgroup.topEquiv).toEquiv).symm
  rw [← IsGalois.card_aut_eq_finrank, hcard]
  exact orderOf_dvd_of_pow_eq_one hg0m

end ABC3.Found.PGC
