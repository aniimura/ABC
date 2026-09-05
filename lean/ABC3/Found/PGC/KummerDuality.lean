import ABC3.Found.PGC.ContinuousHomCount
import Mathlib.FieldTheory.KummerExtension
import Mathlib.FieldTheory.Galois.Infinite
import Mathlib.FieldTheory.KrullTopology
import Mathlib.GroupTheory.SpecificGroups.Cyclic
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import Mathlib.RingTheory.RootsOfUnity.Basic

/-!
# Kummer 双対 —— `#Hom_cont(Γ_F, ℤ/n) = [F^× : (F^×)^n]`(`μ_n ⊆ F` のとき)

[pGC] Proposition 1.2 への経路 C(`ResearchPaper/pgc-goal.md`)の**要石**。
`μ_n ⊆ F` なる標数 0 の体 `F` について

  `card_contHom_eq_index_powRange :
     contHomCard Γ_F n = ((powMonoidHom n : Fˣ →* Fˣ).range).index`

を証明する。ここで `Γ_F = Gal(F̄/F)`、`contHomCard`(`ContinuousHomCount.lean`)は
「核が開な準同型 `Γ_F →* ℤ/n` の個数」である。

## 設計(mathlib の在庫に合わせて選んだ道。逸脱は下に記録)

### 1. Hilbert 90 を使わない

mathlib の Hilbert 90(`groupCohomology.H1ofAutOnUnitsUnique` など)は
`[FiniteDimensional K L]` を要求する。`Γ_F` は無限次なので当たらない
(mathlib 自身が無限次への拡張を TODO にしている)。

代わりに **各指標をそれが切る有限次拡大に落とす**:
`f : Γ_F →* μ_n` の核 `H` が開 ⟹ 固定体 `M = F̄^H` は `F` 上**有限次**で
`Gal(M/F)` は巡回 ⟹ `exists_root_adjoin_eq_top_of_isCyclic`
(`Mathlib/FieldTheory/KummerExtension.lean`)で `M = F(α)`, `α^d ∈ F^×`。
`Γ_F` は無限次のままでよい。

### 2. `X^n - a` の既約性に触らない

mathlib の既約判定は**奇数 `n` 限定**である
(`X_pow_sub_C_irreducible_iff_of_odd` 等、ファイルに `TODO: criteria for even n`)。
既約性を使う設計は `p = 2` で詰まるので、Kummer 指標を**素手で**定義する:

  `κ_a(σ) := σ(α)/α`   (`α` は `X^n - a` の**任意の**根)

* well-defined(`kummerChar_congr`)—— 根の取り替えは `μ_n ⊆ F` だけで済む
* 準同型(`kummerHom` の `map_mul'`)—— 同上
* 核 `= (F^×)^n`(`ker_kummerMap`)—— `κ_a = 1 ⟺ α ∈ F` だけ
* 全射(`range_kummerMap`)—— 上の `exists_root_adjoin_eq_top_of_isCyclic`

`autEquivZmod` / `autEquivRootsOfUnity` は使わない(偶数 `n` の壁に当たるため)。

### 3. 係数群は `μ_n(F̄)`、最後に `ℤ/n` へ移す

本体は `A = rootsOfUnity n F̄`(`Ωˣ` の部分群)で走らせる。`ζ` が原始 `n` 乗根
なので `A ≃* Multiplicative (ZMod n)`(`zmodMulEquivRootsOfUnity`)であり、
`contHomEquivOfMulEquiv` で個数を移す。

## 証明の骨格

```
Fˣ  --κ-->  Hom_cont(Γ_F, μ_n)
ker κ = (F^×)^n            (ker_kummerMap)
range κ = Hom_cont         (range_kummerMap)
  ⊆ : ker κ_a ⊇ Gal(F̄/F(α)) は開           (isOpen_ker_kummerChar)
  ⊇ : f ↦ 固定体 M ↦ M = F(α), α^d ∈ F      (exists_kummerMap_ker_eq)
      同じ核をもつ 2 指標は互いの冪           (exists_zpow_eq_of_ker_eq)
⟹ Fˣ ⧸ (F^×)^n ≃ Hom_cont(Γ_F, μ_n)         (第一同型定理)
```

`Nat.card` の等式を**全単射経由**で出しているので、両辺が無限でも
(`Nat.card = 0` の規約のもとで)正しい。
-/

namespace ABC3.Found.PGC

namespace KummerDual

open IntermediateField

/-! ### 群論の補題(体とは無関係) -/

/-- 有限巡回群 `A` で `d ∣ |A|` なら、`d` 乗写像の核の位数はちょうど `d`。 -/
lemma card_ker_powMonoidHom {A : Type*} [CommGroup A] [Finite A] [IsCyclic A] {d : ℕ}
    (hd : d ∣ Nat.card A) : Nat.card ↥((powMonoidHom d : A →* A).ker) = d := by
  have h1 : Nat.card ↥((powMonoidHom d : A →* A).ker)
      * ((powMonoidHom d : A →* A).ker).index = Nat.card A := Subgroup.card_mul_index _
  have h2 : ((powMonoidHom d : A →* A).ker).index
      = Nat.card ↥((powMonoidHom d : A →* A).range) :=
    Nat.card_congr (QuotientGroup.quotientKerEquivRange _).toEquiv
  have h3 : Nat.card ↥((powMonoidHom d : A →* A).range) * d = Nat.card A := by
    have h := Subgroup.card_mul_index ((powMonoidHom d : A →* A).range)
    rwa [IsCyclic.index_powMonoidHom_range A d, Nat.gcd_eq_right hd] at h
  rw [h2] at h1
  refine Nat.eq_of_mul_eq_mul_right (Nat.card_pos (α := ↥((powMonoidHom d : A →* A).range))) ?_
  rw [h1, ← h3, mul_comm]

/-- 巡回群への準同型の像は巡回、したがって商も巡回。 -/
lemma isCyclic_quotient_ker {G : Type*} [Group G] {A : Type*} [CommGroup A] [IsCyclic A]
    (f : G →* A) : IsCyclic (G ⧸ f.ker) := by
  haveI := Subgroup.isCyclic f.range
  exact isCyclic_of_surjective (QuotientGroup.quotientKerEquivRange f).symm
    (QuotientGroup.quotientKerEquivRange f).symm.surjective

/-- **有限巡回群 `A` への 2 つの準同型が同じ核をもてば、一方は他方の冪**。

`A` は `|A|` の各約数に対して位数がその約数である部分群をちょうど 1 つもつので、
2 つの像は一致する。あとは商 `G/ker` の生成元での値を比べればよい。 -/
lemma exists_zpow_eq_of_ker_eq {G : Type*} [Group G] {A : Type*} [CommGroup A] [Finite A]
    [IsCyclic A] (f g : G →* A) (h : f.ker = g.ker) : ∃ k : ℤ, f = g ^ k := by
  haveI : Finite (G ⧸ g.ker) :=
    Finite.of_injective (QuotientGroup.kerLift g) (QuotientGroup.kerLift_injective g)
  haveI : IsCyclic (G ⧸ g.ker) := isCyclic_quotient_ker g
  obtain ⟨x, hx⟩ := IsCyclic.exists_generator (α := G ⧸ g.ker)
  obtain ⟨σ₀, rfl⟩ := QuotientGroup.mk_surjective x
  set d := Nat.card (G ⧸ g.ker) with hd
  have hrep : ∀ σ : G, ∃ m : ℤ, f σ = f σ₀ ^ m ∧ g σ = g σ₀ ^ m := by
    intro σ
    obtain ⟨m, hm⟩ := hx (σ : G ⧸ g.ker)
    have hm' : ((σ₀ : G ⧸ g.ker)) ^ m = (σ : G ⧸ g.ker) := hm
    rw [← QuotientGroup.mk_zpow, QuotientGroup.eq] at hm'
    refine ⟨m, ?_, ?_⟩
    · have h2 : f ((σ₀ ^ m)⁻¹ * σ) = 1 := by rw [← MonoidHom.mem_ker, h]; exact hm'
      rw [map_mul, map_inv, map_zpow, inv_mul_eq_one] at h2
      exact h2.symm
    · have h2 : g ((σ₀ ^ m)⁻¹ * σ) = 1 := hm'
      rw [map_mul, map_inv, map_zpow, inv_mul_eq_one] at h2
      exact h2.symm
  have hpowd : ∀ σ : G, σ ^ d ∈ g.ker := by
    intro σ
    rw [← QuotientGroup.eq_one_iff, QuotientGroup.mk_pow]
    exact pow_card_eq_one'
  have hdvd : d ∣ Nat.card A := by
    have hc : Nat.card (G ⧸ g.ker) = Nat.card ↥g.range :=
      Nat.card_congr (QuotientGroup.quotientKerEquivRange g).toEquiv
    rw [hd, hc]
    exact Subgroup.card_subgroup_dvd_card _
  have hzp : g.range = Subgroup.zpowers (g σ₀) := by
    apply le_antisymm
    · rintro y ⟨σ, rfl⟩
      obtain ⟨m, -, hm⟩ := hrep σ
      exact ⟨m, hm.symm⟩
    · rw [Subgroup.zpowers_le]
      exact ⟨σ₀, rfl⟩
  have hrange : g.range = (powMonoidHom d : A →* A).ker := by
    refine Subgroup.eq_of_le_of_card_ge ?_ ?_
    · rintro y ⟨σ, rfl⟩
      rw [MonoidHom.mem_ker]
      show (g σ) ^ d = 1
      rw [← map_pow]
      exact hpowd σ
    · rw [card_ker_powMonoidHom hdvd,
        ← Nat.card_congr (QuotientGroup.quotientKerEquivRange g).toEquiv]
  have hfmem : f σ₀ ∈ Subgroup.zpowers (g σ₀) := by
    rw [← hzp, hrange, MonoidHom.mem_ker]
    show (f σ₀) ^ d = 1
    rw [← map_pow, ← MonoidHom.mem_ker, h]
    exact hpowd σ₀
  obtain ⟨k, hk⟩ := hfmem
  refine ⟨k, MonoidHom.ext fun σ => ?_⟩
  obtain ⟨m, hfm, hgm⟩ := hrep σ
  rw [MonoidHom.zpow_apply, hfm, hgm, ← hk, ← zpow_mul, ← zpow_mul, mul_comm]

/-- 係数群を同型で取り替えても「核が開」は変わらない。 -/
def contHomEquivOfMulEquiv {G : Type*} [Group G] [TopologicalSpace G]
    [SeparatelyContinuousMul G] {A B : Type*} [CommGroup A] [CommGroup B] (e : A ≃* B) :
    contHom G A ≃ contHom G B where
  toFun f := ⟨e.toMonoidHom.comp (f : G →* A), by
    rw [mem_contHom, MonoidHom.ker_comp_of_injective _ _ e.injective]; exact f.2⟩
  invFun g := ⟨e.symm.toMonoidHom.comp (g : G →* B), by
    rw [mem_contHom, MonoidHom.ker_comp_of_injective _ _ e.symm.injective]; exact g.2⟩
  left_inv f := Subtype.ext (MonoidHom.ext fun x => by simp)
  right_inv g := Subtype.ext (MonoidHom.ext fun x => by simp)

section Base

variable {F : Type} [Field F] {n : ℕ} {ζ : F}

/-! ### `μ_n(F̄) ⊆ F` —— 根の取り替えを支える唯一の事実 -/

/-- `Ω = F̄` の `n` 乗根 1 はすべて `ζ` の冪、したがって `F` の元の像。 -/
lemma exists_algebraMap_eq_of_pow_eq_one (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0)
    {x : AlgebraicClosure F} (hx : x ^ n = 1) :
    ∃ c : F, algebraMap F (AlgebraicClosure F) c = x := by
  haveI : NeZero n := ⟨hn⟩
  have hζ' : IsPrimitiveRoot (algebraMap F (AlgebraicClosure F) ζ) n :=
    hζ.map_of_injective (algebraMap F (AlgebraicClosure F)).injective
  obtain ⟨i, -, hi⟩ := hζ'.eq_pow_of_pow_eq_one hx
  exact ⟨ζ ^ i, by rw [map_pow, hi]⟩

/-- `n` 乗根 1 は `Γ_F` の作用で固定される(`μ_n ⊆ F` だから)。 -/
lemma apply_eq_self_of_pow_eq_one (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0)
    {x : AlgebraicClosure F} (hx : x ^ n = 1)
    (σ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) : σ x = x := by
  obtain ⟨c, rfl⟩ := exists_algebraMap_eq_of_pow_eq_one hζ hn hx
  exact σ.commutes c

/-! ### Kummer 指標 -/

lemma ne_zero_of_pow_eq {a : Fˣ} {α : AlgebraicClosure F} (hn : n ≠ 0)
    (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F)) : α ≠ 0 := by
  intro h0
  rw [h0, zero_pow hn] at h
  exact a.ne_zero ((algebraMap F (AlgebraicClosure F)).injective (by simpa using h.symm))

lemma fix_pow_of_pow_eq {a : Fˣ} {α : AlgebraicClosure F}
    (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F))
    (σ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) : σ (α ^ n) = α ^ n := by
  rw [h, AlgEquiv.commutes]

lemma div_pow_eq_one {a : Fˣ} {α : AlgebraicClosure F} (hn : n ≠ 0)
    (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F))
    (σ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) : (σ α / α) ^ n = 1 := by
  rw [div_pow, ← map_pow, fix_pow_of_pow_eq h σ, div_self (pow_ne_zero n (ne_zero_of_pow_eq hn h))]

/-- **Kummer 指標**(素手の定義)。`α^n = a ∈ F^×` なる `α` に対し `κ_a(σ) = σα/α`。
`X^n - a` の既約性は一切使わない —— `map_mul'` は「`σα/α ∈ μ_n ⊆ F` は `τ` で
固定される」だけから出る。 -/
noncomputable def kummerHom (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) (a : Fˣ)
    {α : AlgebraicClosure F} (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F)) :
    (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) →* (AlgebraicClosure F)ˣ where
  toFun σ := Units.mk0 (σ α / α)
    (div_ne_zero (by simpa using ne_zero_of_pow_eq hn h) (ne_zero_of_pow_eq hn h))
  map_one' := by apply Units.ext; simp [div_self (ne_zero_of_pow_eq hn h)]
  map_mul' σ τ := by
    have hα := ne_zero_of_pow_eq hn h
    apply Units.ext
    have h1 : σ (τ α / α) = τ α / α :=
      apply_eq_self_of_pow_eq_one hζ hn (div_pow_eq_one hn h τ) σ
    have h2 : (σ * τ) α = (τ α / α) * σ α := by
      rw [AlgEquiv.mul_apply]
      nth_rewrite 1 [show τ α = (τ α / α) * α from (div_mul_cancel₀ _ hα).symm]
      rw [map_mul, h1]
    simp only [Units.val_mk0, Units.val_mul, h2]
    field_simp

@[simp] lemma kummerHom_val (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) (a : Fˣ)
    {α : AlgebraicClosure F} (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F))
    (σ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) :
    ((kummerHom hζ hn a h σ : (AlgebraicClosure F)ˣ) : AlgebraicClosure F) = σ α / α := rfl

/-- Kummer 指標は `μ_n(Ω)` に値をとる。 -/
noncomputable def kummerChar (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) (a : Fˣ)
    {α : AlgebraicClosure F} (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F)) :
    (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) →* ↥(rootsOfUnity n (AlgebraicClosure F)) :=
  MonoidHom.codRestrict (kummerHom hζ hn a h) _ (fun σ => by
    rw [mem_rootsOfUnity]
    exact Units.ext (by simpa using div_pow_eq_one hn h σ))

@[simp] lemma kummerChar_val (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) (a : Fˣ)
    {α : AlgebraicClosure F} (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F))
    (σ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) :
    (((kummerChar hζ hn a h σ : ↥(rootsOfUnity n (AlgebraicClosure F))) :
        (AlgebraicClosure F)ˣ) : AlgebraicClosure F) = σ α / α := rfl

/-- `μ_n(Ω)` に値をとる準同型の外延性(値を `Ω` の元として比べればよい)。 -/
lemma hom_ext {f g : (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) →*
      ↥(rootsOfUnity n (AlgebraicClosure F))}
    (h : ∀ σ, ((f σ : (AlgebraicClosure F)ˣ) : AlgebraicClosure F)
        = ((g σ : (AlgebraicClosure F)ˣ) : AlgebraicClosure F)) : f = g :=
  MonoidHom.ext fun σ => Subtype.ext (Units.ext (h σ))

/-- 根 `α` の取り方によらない(既約性不要:`β/α ∈ μ_n ⊆ F` だけ)。 -/
lemma kummerChar_congr (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) (a : Fˣ)
    {α β : AlgebraicClosure F} (hα : α ^ n = algebraMap F (AlgebraicClosure F) (a : F))
    (hβ : β ^ n = algebraMap F (AlgebraicClosure F) (a : F)) :
    kummerChar hζ hn a hα = kummerChar hζ hn a hβ := by
  have hα0 := ne_zero_of_pow_eq hn hα
  have hβ0 := ne_zero_of_pow_eq hn hβ
  refine hom_ext fun σ => ?_
  have hr : (β / α) ^ n = 1 := by
    rw [div_pow, hα, hβ, div_self (by rw [← hα]; exact pow_ne_zero n hα0)]
  have hfr : σ (β / α) = β / α := apply_eq_self_of_pow_eq_one hζ hn hr σ
  simp only [kummerChar_val]
  have hb : σ β = (β / α) * σ α := by
    nth_rewrite 1 [show β = (β / α) * α from (div_mul_cancel₀ _ hα0).symm]
    rw [map_mul, hfr]
  rw [hb]
  field_simp

/-- `κ_a` の核はちょうど `α` の固定群。 -/
lemma mem_ker_kummerChar (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) (a : Fˣ)
    {α : AlgebraicClosure F} (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F))
    (σ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) :
    σ ∈ (kummerChar hζ hn a h).ker ↔ σ α = α := by
  rw [MonoidHom.mem_ker, ← Subtype.coe_inj]
  change kummerHom hζ hn a h σ = 1 ↔ _
  rw [← Units.val_inj]
  change σ α / α = 1 ↔ _
  rw [div_eq_one_iff_eq (ne_zero_of_pow_eq hn h)]

/-- `κ_a` の核は `F(α)` の固定群にちょうど一致する。 -/
lemma ker_kummerChar_eq (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) (a : Fˣ)
    {α : AlgebraicClosure F} (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F)) :
    (kummerChar hζ hn a h).ker = F⟮α⟯.fixingSubgroup := by
  apply le_antisymm
  · rw [← IntermediateField.le_iff_le, IntermediateField.adjoin_le_iff]
    intro y hy
    rw [Set.mem_singleton_iff] at hy
    subst hy
    rw [SetLike.mem_coe, IntermediateField.mem_fixedField_iff]
    intro f hf
    exact (mem_ker_kummerChar hζ hn a h f).1 hf
  · intro σ hσ
    rw [mem_ker_kummerChar hζ hn a h]
    exact (IntermediateField.mem_fixingSubgroup_iff _ σ).1 hσ _
      (IntermediateField.mem_adjoin_simple_self F α)

/-- `a ∈ F^×` の `n` 乗根(代数閉包で選ぶ)。 -/
noncomputable def nthRoot (hn : n ≠ 0) (a : Fˣ) : AlgebraicClosure F :=
  (IsAlgClosed.exists_pow_nat_eq (algebraMap F (AlgebraicClosure F) (a : F))
    (Nat.pos_of_ne_zero hn)).choose

lemma nthRoot_pow (hn : n ≠ 0) (a : Fˣ) :
    (nthRoot hn a) ^ n = algebraMap F (AlgebraicClosure F) (a : F) :=
  (IsAlgClosed.exists_pow_nat_eq (algebraMap F (AlgebraicClosure F) (a : F))
    (Nat.pos_of_ne_zero hn)).choose_spec

/-- **Kummer 写像** `κ : F^× →* Hom(Γ_F, μ_n)`, `a ↦ κ_a`。 -/
noncomputable def kummerMap (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    Fˣ →* ((AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) →*
      ↥(rootsOfUnity n (AlgebraicClosure F))) where
  toFun a := kummerChar hζ hn a (nthRoot_pow hn a)
  map_one' := by
    refine hom_ext fun σ => ?_
    have h1 : (nthRoot hn (1 : Fˣ)) ^ n = 1 := by rw [nthRoot_pow]; simp
    have h2 := apply_eq_self_of_pow_eq_one hζ hn h1 σ
    rw [kummerChar_val, h2, div_self (ne_zero_of_pow_eq hn (nthRoot_pow hn (1 : Fˣ)))]
    simp
  map_mul' a b := by
    have hab : (nthRoot hn a * nthRoot hn b) ^ n
        = algebraMap F (AlgebraicClosure F) (((a * b : Fˣ) : F)) := by
      rw [mul_pow, nthRoot_pow, nthRoot_pow, Units.val_mul, map_mul]
    rw [kummerChar_congr hζ hn (a * b) (nthRoot_pow hn (a * b)) hab]
    refine hom_ext fun σ => ?_
    have ha0 := ne_zero_of_pow_eq hn (nthRoot_pow hn a)
    have hb0 := ne_zero_of_pow_eq hn (nthRoot_pow hn b)
    rw [kummerChar_val, MonoidHom.mul_apply]
    push_cast
    rw [kummerChar_val, kummerChar_val, map_mul]
    field_simp

/-- `κ_a` は `X^n - a` の**任意の**根で計算してよい。 -/
lemma kummerMap_eq (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) (a : Fˣ)
    {α : AlgebraicClosure F} (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F)) :
    kummerMap hζ hn a = kummerChar hζ hn a h :=
  kummerChar_congr hζ hn a (nthRoot_pow hn a) h

/-! ### `μ_n(Ω) ≃* ℤ/n` -/

/-- `ζ` を `Ω^×` の元として見たもの。 -/
noncomputable def zetaUnit (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) : (AlgebraicClosure F)ˣ :=
  Units.mk0 (algebraMap F (AlgebraicClosure F) ζ)
    ((map_ne_zero_iff _ (algebraMap F (AlgebraicClosure F)).injective).2 (hζ.ne_zero hn))

lemma zetaUnit_isPrimitiveRoot (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    IsPrimitiveRoot (zetaUnit hζ hn) n :=
  IsPrimitiveRoot.coe_units_iff.1
    (hζ.map_of_injective (algebraMap F (AlgebraicClosure F)).injective)

/-- `μ_n(Ω) = ⟨ζ⟩`。 -/
lemma zpowers_zetaUnit (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    Subgroup.zpowers (zetaUnit hζ hn) = rootsOfUnity n (AlgebraicClosure F) := by
  haveI : NeZero n := ⟨hn⟩
  exact (zetaUnit_isPrimitiveRoot hζ hn).zpowers_eq

lemma card_rootsOfUnity (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    Nat.card ↥(rootsOfUnity n (AlgebraicClosure F)) = n := by
  rw [← zpowers_zetaUnit hζ hn, Nat.card_zpowers, ← (zetaUnit_isPrimitiveRoot hζ hn).eq_orderOf]

lemma isCyclic_rootsOfUnity (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    IsCyclic ↥(rootsOfUnity n (AlgebraicClosure F)) := by
  rw [← zpowers_zetaUnit hζ hn]
  infer_instance

lemma finite_rootsOfUnity (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    Finite ↥(rootsOfUnity n (AlgebraicClosure F)) :=
  Nat.finite_of_card_ne_zero (by rw [card_rootsOfUnity hζ hn]; exact hn)

/-- `μ_n(Ω) ≃* ℤ/n`(乗法的に書いたもの)。 -/
noncomputable def zmodMulEquivRootsOfUnity (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    Multiplicative (ZMod n) ≃* ↥(rootsOfUnity n (AlgebraicClosure F)) :=
  (AddEquiv.toMultiplicativeLeft (zetaUnit_isPrimitiveRoot hζ hn).zmodEquivZPowers).trans
    (MulEquiv.subgroupCongr (zpowers_zetaUnit hζ hn))

end Base

section CharZero

variable {F : Type} [Field F] [CharZero F] {n : ℕ} {ζ : F}

/-! ### 核が `(F^×)^n` であること -/

/-- `Γ_F` の全元で固定される `Ω` の元は `F` の元(`fixedField ⊤ = ⊥`)。 -/
lemma exists_algebraMap_eq_of_forall_fixed {x : AlgebraicClosure F}
    (h : ∀ σ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F, σ x = x) :
    ∃ c : F, algebraMap F (AlgebraicClosure F) c = x := by
  have key : IntermediateField.fixedField
      (⊤ : Subgroup (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F))
      = (⊥ : IntermediateField F (AlgebraicClosure F)) := by
    rw [← IntermediateField.fixingSubgroup_bot]
    exact InfiniteGalois.fixedField_fixingSubgroup ⊥
  have hx : x ∈ IntermediateField.fixedField
      (⊤ : Subgroup (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)) := by
    rw [IntermediateField.mem_fixedField_iff]
    intro f _
    exact h f
  rw [key, IntermediateField.mem_bot] at hx
  exact hx

/-- **Kummer 写像の核は `(F^×)^n`**。既約性を使わない:`κ_a = 1 ⟺ α ∈ F` だけ。 -/
lemma ker_kummerMap (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    (kummerMap hζ hn).ker = (powMonoidHom n : Fˣ →* Fˣ).range := by
  ext a
  simp only [MonoidHom.mem_ker, MonoidHom.mem_range]
  constructor
  · intro h
    have hfix : ∀ σ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F,
        σ (nthRoot hn a) = nthRoot hn a := by
      intro σ
      rw [← mem_ker_kummerChar hζ hn a (nthRoot_pow hn a) σ, MonoidHom.mem_ker]
      have hone : kummerChar hζ hn a (nthRoot_pow hn a) = 1 := h
      rw [hone]
      rfl
    obtain ⟨c, hc⟩ := exists_algebraMap_eq_of_forall_fixed hfix
    have hcn : c ^ n = (a : F) := by
      apply (algebraMap F (AlgebraicClosure F)).injective
      rw [map_pow, hc, nthRoot_pow]
    have hc0 : c ≠ 0 := by
      intro h0
      rw [h0, zero_pow hn] at hcn
      exact a.ne_zero hcn.symm
    exact ⟨Units.mk0 c hc0, by ext; simpa using hcn⟩
  · rintro ⟨c, rfl⟩
    have hpow : (algebraMap F (AlgebraicClosure F) (c : F)) ^ n
        = algebraMap F (AlgebraicClosure F) (((powMonoidHom n c : Fˣ) : F)) := by
      rw [← map_pow]
      rfl
    rw [kummerMap_eq hζ hn _ hpow]
    refine hom_ext fun σ => ?_
    rw [kummerChar_val, AlgEquiv.commutes, div_self (ne_zero_of_pow_eq hn hpow)]
    rfl

/-! ### 連続性(核が開)と全射性 -/

/-- `κ_a` は連続(核が開)。`F(α)/F` が有限次だから。 -/
lemma isOpen_ker_kummerChar (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) (a : Fˣ)
    {α : AlgebraicClosure F} (h : α ^ n = algebraMap F (AlgebraicClosure F) (a : F)) :
    IsOpen (((kummerChar hζ hn a h).ker :
      Subgroup (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)) :
      Set (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)) := by
  haveI : FiniteDimensional F F⟮α⟯ :=
    IntermediateField.adjoin.finiteDimensional (Algebra.IsIntegral.isIntegral α)
  rw [ker_kummerChar_eq hζ hn a h]
  exact (InfiniteGalois.isOpen_iff_finite F⟮α⟯).2 inferInstance

/-- 有限次巡回拡大 `M/F`(次数が `n` を割る)の固定群は、ある Kummer 指標の核として
実現される。ここだけが `Mathlib/FieldTheory/KummerExtension.lean` の
`exists_root_adjoin_eq_top_of_isCyclic` を使う —— この補題は `X^n - a` の
既約性を要求しないので、偶数 `n` でも成り立つ。 -/
lemma exists_kummerMap_ker_eq_fixingSubgroup (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0)
    (M : IntermediateField F (AlgebraicClosure F)) [FiniteDimensional F ↥M] [IsGalois F ↥M]
    [IsCyclic (↥M ≃ₐ[F] ↥M)] (hdvd : Module.finrank F ↥M ∣ n) :
    ∃ c : Fˣ, (kummerMap hζ hn c).ker = M.fixingSubgroup := by
  have hdpos : 0 < Module.finrank F ↥M := Module.finrank_pos
  have hprim : (primitiveRoots (Module.finrank F ↥M) F).Nonempty := by
    obtain ⟨k, hk⟩ := hdvd
    exact ⟨ζ ^ k, (mem_primitiveRoots hdpos).2
      (IsPrimitiveRoot.pow (Nat.pos_of_ne_zero hn) hζ (hk.trans (mul_comm _ _)))⟩
  by_cases hd1 : Module.finrank F ↥M = 1
  · refine ⟨1, ?_⟩
    rw [map_one, MonoidHom.ker_one]
    symm
    rw [← Subgroup.index_eq_one, ← IntermediateField.finrank_eq_fixingSubgroup_index, hd1]
  · obtain ⟨α, hα1, hα2⟩ := exists_root_adjoin_eq_top_of_isCyclic F ↥M hprim
    have hα0 : α ≠ 0 := by
      rintro rfl
      apply hd1
      have hbt : (⊥ : IntermediateField F ↥M) = ⊤ := by
        rw [← hα2]
        exact (IntermediateField.adjoin_simple_eq_bot_iff.2 (zero_mem _)).symm
      rw [← IntermediateField.finrank_top' (F := F) (E := ↥M), ← hbt,
        IntermediateField.finrank_bot]
    obtain ⟨b, hb⟩ := hα1
    have hb0 : b ≠ 0 := by
      rintro rfl
      rw [map_zero] at hb
      exact pow_ne_zero _ hα0 hb.symm
    have hαd : ((α : AlgebraicClosure F)) ^ (Module.finrank F ↥M)
        = algebraMap F (AlgebraicClosure F) b := by
      have h : ((α ^ (Module.finrank F ↥M) : ↥M) : AlgebraicClosure F)
          = ((algebraMap F ↥M b : ↥M) : AlgebraicClosure F) := by rw [hb]
      rw [IsScalarTower.algebraMap_apply F ↥M (AlgebraicClosure F) b,
        IntermediateField.algebraMap_apply, ← h]
      push_cast
      ring
    have hαΩ : ((α : AlgebraicClosure F)) ^ n
        = algebraMap F (AlgebraicClosure F)
          (((Units.mk0 b hb0) ^ (n / Module.finrank F ↥M) : Fˣ) : F) := by
      have hnn : Module.finrank F ↥M * (n / Module.finrank F ↥M) = n :=
        Nat.mul_div_cancel' hdvd
      calc ((α : AlgebraicClosure F)) ^ n
          = (((α : AlgebraicClosure F)) ^ (Module.finrank F ↥M)) ^ (n / Module.finrank F ↥M) := by
            rw [← pow_mul, hnn]
        _ = (algebraMap F (AlgebraicClosure F) b) ^ (n / Module.finrank F ↥M) := by rw [hαd]
        _ = algebraMap F (AlgebraicClosure F) (b ^ (n / Module.finrank F ↥M)) := by rw [map_pow]
        _ = _ := by norm_cast
    have hFα : F⟮(α : AlgebraicClosure F)⟯ = M := by
      have h := IntermediateField.lift_adjoin F M ({α} : Set ↥M)
      rw [Set.image_singleton, show IntermediateField.adjoin F ({α} : Set ↥M) = ⊤ from hα2,
        IntermediateField.lift_top] at h
      exact h.symm
    exact ⟨_, by rw [kummerMap_eq hζ hn _ hαΩ, ker_kummerChar_eq hζ hn _ hαΩ, hFα]⟩

/-- **全射性の核心**:核が開な指標 `f` に対し、`f` と同じ核をもつ Kummer 指標がある。

`ker f` は開 ⟹ 固定体 `M` は有限次・Galois・`Gal(M/F)` は巡回(`Γ/ker f` の商だから)
⟹ 上の補題。**`Γ_F` は無限次のままでよい**。 -/
lemma exists_kummerMap_ker_eq (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0)
    (f : (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) →*
      ↥(rootsOfUnity n (AlgebraicClosure F)))
    (hf : IsOpen ((f.ker : Subgroup (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)) :
      Set (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F))) :
    ∃ c : Fˣ, (kummerMap hζ hn c).ker = f.ker := by
  haveI := isCyclic_rootsOfUnity hζ hn
  haveI := finite_rootsOfUnity hζ hn
  have hHclosed : IsClosed ((f.ker : Subgroup (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)) :
      Set (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)) := Subgroup.isClosed_of_isOpen _ hf
  have hMfix : (IntermediateField.fixedField f.ker).fixingSubgroup = f.ker :=
    InfiniteGalois.fixingSubgroup_fixedField
      (⟨f.ker, hHclosed⟩ : ClosedSubgroup (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F))
  haveI hfd : FiniteDimensional F ↥(IntermediateField.fixedField f.ker) :=
    (InfiniteGalois.isOpen_iff_finite _).1 (by rw [hMfix]; exact hf)
  haveI hgal : IsGalois F ↥(IntermediateField.fixedField f.ker) :=
    (InfiniteGalois.normal_iff_isGalois _).1 (by rw [hMfix]; infer_instance)
  have hkerle : f.ker ≤
      (AlgEquiv.restrictNormalHom (↥(IntermediateField.fixedField f.ker))).ker := by
    intro σ hσ
    rw [MonoidHom.mem_ker]
    refine AlgEquiv.ext fun y => Subtype.ext ?_
    rw [AlgEquiv.restrictNormalHom_apply]
    exact (IntermediateField.mem_fixedField_iff _ _).1 y.2 σ hσ
  haveI : IsCyclic ((AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) ⧸ f.ker) :=
    isCyclic_quotient_ker f
  haveI : IsCyclic (↥(IntermediateField.fixedField f.ker) ≃ₐ[F]
      ↥(IntermediateField.fixedField f.ker)) := by
    refine isCyclic_of_surjective (QuotientGroup.lift f.ker _ hkerle) ?_
    intro y
    obtain ⟨σ, hσ⟩ := AlgEquiv.restrictNormalHom_surjective
      (K₁ := ↥(IntermediateField.fixedField f.ker)) (E := AlgebraicClosure F) y
    exact ⟨(σ : (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) ⧸ f.ker), hσ⟩
  have hdindex : Module.finrank F ↥(IntermediateField.fixedField f.ker) = f.ker.index := by
    have h := IntermediateField.finrank_eq_fixingSubgroup_index
      (IntermediateField.fixedField f.ker)
    rwa [hMfix] at h
  have hdvd : Module.finrank F ↥(IntermediateField.fixedField f.ker) ∣ n := by
    have h1 : f.ker.index = Nat.card ↥f.range :=
      Nat.card_congr (QuotientGroup.quotientKerEquivRange f).toEquiv
    have h2 : Nat.card ↥f.range ∣ Nat.card ↥(rootsOfUnity n (AlgebraicClosure F)) :=
      Subgroup.card_subgroup_dvd_card _
    rw [card_rootsOfUnity hζ hn] at h2
    rw [hdindex, h1]
    exact h2
  obtain ⟨c, hc⟩ := exists_kummerMap_ker_eq_fixingSubgroup hζ hn
    (IntermediateField.fixedField f.ker) hdvd
  exact ⟨c, by rw [hc, hMfix]⟩

/-! ### 双対 -/

/-- **Kummer 双対**(部分群の形):`κ` の像はちょうど連続指標全体。 -/
theorem range_kummerMap (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    (kummerMap hζ hn).range
      = contHom (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)
          ↥(rootsOfUnity n (AlgebraicClosure F)) := by
  haveI := isCyclic_rootsOfUnity hζ hn
  haveI := finite_rootsOfUnity hζ hn
  apply le_antisymm
  · rintro g ⟨a, rfl⟩
    exact isOpen_ker_kummerChar hζ hn a (nthRoot_pow hn a)
  · intro f hf
    obtain ⟨c, hc⟩ := exists_kummerMap_ker_eq hζ hn f hf
    obtain ⟨k, hk⟩ := exists_zpow_eq_of_ker_eq f (kummerMap hζ hn c) hc.symm
    exact ⟨c ^ k, by rw [map_zpow, ← hk]⟩

/-- **Kummer 双対**(`μ_n` 版の計数)。 -/
theorem card_contHom_rootsOfUnity (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0) :
    Nat.card ↥(contHom (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)
        ↥(rootsOfUnity n (AlgebraicClosure F)))
      = ((powMonoidHom n : Fˣ →* Fˣ).range).index := by
  rw [← range_kummerMap hζ hn,
    ← Nat.card_congr (QuotientGroup.quotientKerEquivRange (kummerMap hζ hn)).toEquiv,
    ← ker_kummerMap hζ hn]
  rfl

end CharZero

end KummerDual

open KummerDual in
/-- **Kummer 双対**(経路 C の要石)。

`μ_n ⊆ F`(標数 0)のとき、`Γ_F = Gal(F̄/F)` からの連続指標の個数は
`F^×` の `n` 乗部分群の指数に一致する:

  `#Hom_cont(Γ_F, ℤ/n) = [F^× : (F^×)^n]`.

証明は `F^× ⧸ (F^×)^n ≃ Hom_cont(Γ_F, μ_n)` という**全単射**を作るので、
両辺が無限でも(`Nat.card = 0` の規約のもとで)正しい。 -/
theorem card_contHom_eq_index_powRange (F : Type) [Field F] [CharZero F] {n : ℕ} (hn : n ≠ 0)
    {ζ : F} (hζ : IsPrimitiveRoot ζ n) :
    contHomCard (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) n
      = ((powMonoidHom n : Fˣ →* Fˣ).range).index := by
  have h : contHomCard (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) n
      = Nat.card ↥(contHom (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)
          (Multiplicative (ZMod n))) := rfl
  rw [h, Nat.card_congr (contHomEquivOfMulEquiv (zmodMulEquivRootsOfUnity hζ hn))]
  exact card_contHom_rootsOfUnity hζ hn

end ABC3.Found.PGC
