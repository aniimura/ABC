import ABC3.Found.PGC.UnramifiedCriterion

/-!
# 不分岐拡大の部分拡大は不分岐——`K(x) ≤ K^ur ⟺ K(x)/K が不分岐`

`UnramifiedExtension.lean` は「不分岐拡大は各次数にちょうど一つ」
(存在 `exists_isUnramifiedAdjoin`・一意 `adjoin_eq_of_isUnramified`)まで
出していたが、**部分拡大の不分岐性**——`K(x) ⊆ K(y)` で `K(y)/K` が不分岐
なら `K(x)/K` も不分岐——は無かった。分岐指数の塔の乗法性を経由するのが
教科書の道だが、`Ideal.ramificationIdx_algebra_tower'` を使うには
`adjoinIntegers K x → adjoinIntegers K y` の代数構造と `LiesOver` を
組む配管が要る。

**本ファイルは塔を使わない**——不分岐拡大の Galois 群が**巡回群**である
ことを使う:

1. `m := [K(x):K]` は `n := [K(y):K]` を割る(`finrank_dvd_of_adjoin_le`)。
2. 次数 `n` の不分岐拡大 `K(w)` は `Gal` が巡回(`exists_isCyclic_gal`)で、
   一意性から `K(y) = K(w)`。
3. 次数 `m` の不分岐拡大 `K(z)` を取ると `K(z) ⊆ K(w)`(`adjoin_le_of_dvd`)。
4. `K(x)` と `K(z)` はどちらも `K(w)` の中の次数 `m` の部分拡大。
   `Γ_K/Gal(K̄/K(w)) ≅ Gal(K(w)/K)` は**位数 `n` の巡回群**なので、
   指数の等しい部分群は一致する ⟹ `K(x) = K(z)` ⟹ 不分岐。

## ★巡回群では位数の等しい部分群は一致する

mathlib に直接の補題が無かったので `IsCyclic.card_pow_eq_one_le`
(`x^d = 1` の解は高々 `d` 個)から作った:
`|H| = d` なら Lagrange で `H ⊆ {a | a^d = 1}`、後者の元は高々 `d` 個
なので一致。`H` は `d` だけで決まる。

## 結果

* **`isUnramifiedAdjoin_of_le`**——部分拡大の不分岐性
* **`adjoin_le_unramifiedClosure_iff`**——`K(x) ≤ K^ur ⟺ 不分岐`
* `mem_unramifiedClosure_iff_isUnramified`——`x ∈ K^ur ⟺ 不分岐`
* `isUnramifiedAt_iff_fixedField_le`——原文の判定条件は「固定体が `K^ur`
  に入る」ことと**同値**(`UnramifiedCriterion.lean` では片方向だった)
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NormedField Valued

/-! ## 巡回群の部分群は位数で決まる -/

/-- **有限巡回群では、位数の等しい部分群は一致する**。
mathlib に直接の補題が無いので `IsCyclic.card_pow_eq_one_le` から作る:
`|H| = d` なら Lagrange で `H ⊆ {a | a^d = 1}`、後者は高々 `d` 元。 -/
theorem eq_of_natCard_eq_of_isCyclic {G : Type*} [Group G] [Finite G] [IsCyclic G]
    {H H' : Subgroup G} (h : Nat.card H = Nat.card H') : H = H' := by
  classical
  haveI : Fintype G := Fintype.ofFinite G
  set d := Nat.card H with hd
  have hdpos : 0 < d := Nat.card_pos
  have hle : ({a : G | a ^ d = 1}).ncard ≤ d := by
    have h1 := IsCyclic.card_pow_eq_one_le (α := G) (n := d) hdpos
    rw [Set.ncard_eq_toFinset_card']
    convert h1 using 2
    ext a
    simp
  have key : ∀ J : Subgroup G, Nat.card J = d → (J : Set G) = {a : G | a ^ d = 1} := by
    intro J hJ
    have hsub : (J : Set G) ⊆ {a : G | a ^ d = 1} := by
      intro a ha
      have hp : (⟨a, ha⟩ : J) ^ d = 1 := by rw [← hJ]; exact pow_card_eq_one'
      have hv := congrArg (Subtype.val : J → G) hp
      simpa using hv
    have hcard : (J : Set G).ncard = Nat.card J := rfl
    refine Set.eq_of_subset_of_ncard_le hsub ?_ (Set.toFinite _)
    rw [hcard, hJ]; exact hle
  exact SetLike.ext' ((key H rfl).trans (key H' h.symm).symm)

/-- 商が有限巡回群なら、その核を含む部分群は**指数**で決まる。 -/
theorem subgroup_eq_of_index_eq_of_isCyclic_quotient {G : Type*} [Group G]
    {N : Subgroup G} [N.Normal] [Finite (G ⧸ N)] [IsCyclic (G ⧸ N)]
    {H H' : Subgroup G} (hH : N ≤ H) (hH' : N ≤ H') (h : H.index = H'.index) : H = H' := by
  classical
  set q := QuotientGroup.mk' N with hq
  have hker : q.ker = N := QuotientGroup.ker_mk' N
  have hkH : q.ker ≤ H := by rw [hker]; exact hH
  have hkH' : q.ker ≤ H' := by rw [hker]; exact hH'
  have hrange : q.range = ⊤ := MonoidHom.range_eq_top.mpr (QuotientGroup.mk'_surjective N)
  have hidx : ∀ J : Subgroup G, N ≤ J → (Subgroup.map q J).index = J.index := by
    intro J hJ
    rw [Subgroup.index_map, hker, sup_eq_left.mpr hJ, hrange, Subgroup.index_top, mul_one]
  have e1 := Subgroup.card_mul_index (Subgroup.map q H)
  have e2 := Subgroup.card_mul_index (Subgroup.map q H')
  rw [hidx H hH] at e1
  rw [hidx H' hH'] at e2
  have hpos : 0 < H.index := by
    rw [← hidx H hH]
    exact Nat.pos_of_ne_zero (Subgroup.index_ne_zero_iff_finite.mpr inferInstance)
  have hcard : Nat.card (Subgroup.map q H) = Nat.card (Subgroup.map q H') := by
    rw [← h] at e2
    exact Nat.eq_of_mul_eq_mul_right hpos (e1.trans e2.symm)
  have hmap := eq_of_natCard_eq_of_isCyclic hcard
  calc H = Subgroup.comap q (Subgroup.map q H) := (Subgroup.comap_map_eq_self hkH).symm
    _ = Subgroup.comap q (Subgroup.map q H') := by rw [hmap]
    _ = H' := Subgroup.comap_map_eq_self hkH'

/-! ## 部分拡大の不分岐性 -/

variable {p : ℕ} [Fact p.Prime]

/-- 不分岐性は**生成体だけ**で決まる(`residueDegree` と `finrank` の
両方が生成体だけで決まるから)。`intermediateLocalField_congr` 経由。 -/
theorem isUnramifiedAdjoin_of_adjoin_eq (K : PAdicLocalField p) {x z : K.closure}
    (h : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      = IntermediateField.adjoin K.carrier ({z} : Set K.closure))
    (hz : IsUnramifiedAdjoin K z) : IsUnramifiedAdjoin K x := by
  rw [isUnramifiedAdjoin_iff_residueDegree] at hz ⊢
  have hcard : residueDegree K x = residueDegree K z := by
    rw [← card_residueField_adjoinField, ← card_residueField_adjoinField,
      adjoinField_eq K x, adjoinField_eq K z]
    exact congrArg (fun L : PAdicLocalField p => Nat.card 𝓀[L.carrier])
      (intermediateLocalField_congr K _ _ h)
  rw [hcard, h]
  exact hz

/-- `[K(x):K] = [Γ_K : Gal(K̄/K(x))]`——無限次 Galois 対応の言い換え。 -/
theorem index_fixingSubgroup_adjoin (K : PAdicLocalField p) (x : K.closure) :
    ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup).index
      = Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  haveI := isGalois_closure K
  have hop : IsOpen
      (((IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup)
        : Set K.absGal) := IntermediateField.fixingSubgroup_isOpen _
  have h := finrank_fixedField_eq_index K _ hop
  rw [InfiniteGalois.fixedField_fixingSubgroup] at h
  exact h.symm

/-- **巡回 Galois 拡大の中では、次数の等しい部分拡大は一致する**。 -/
theorem adjoin_eq_of_le_of_finrank_eq (K : PAdicLocalField p) {x z w : K.closure}
    (hxw : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({w} : Set K.closure))
    (hzw : IntermediateField.adjoin K.carrier ({z} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({w} : Set K.closure))
    (hnw : Normal K.carrier (IntermediateField.adjoin K.carrier ({w} : Set K.closure)))
    (hcyc : IsCyclic ((IntermediateField.adjoin K.carrier ({w} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({w} : Set K.closure))))
    (hrank : Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({z} : Set K.closure))) :
    IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      = IntermediateField.adjoin K.carrier ({z} : Set K.closure) := by
  haveI := isGalois_closure K
  haveI := hnw
  haveI := IsAlgClosure.normal K.carrier K.closure
  set N := (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
    (IntermediateField.adjoin K.carrier ({w} : Set K.closure))).ker with hN
  have hNw : N = (IntermediateField.adjoin K.carrier ({w} : Set K.closure)).fixingSubgroup :=
    IntermediateField.restrictNormalHom_ker _
  have hNx : N ≤ (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup := by
    rw [hNw]; exact IntermediateField.fixingSubgroup_le hxw
  have hNz : N ≤ (IntermediateField.adjoin K.carrier ({z} : Set K.closure)).fixingSubgroup := by
    rw [hNw]; exact IntermediateField.fixingSubgroup_le hzw
  have hquot : (K.absGal ⧸ N) ≃* ((IntermediateField.adjoin K.carrier ({w} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({w} : Set K.closure))) :=
    QuotientGroup.quotientKerEquivOfSurjective _
      (AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
        (K₁ := (IntermediateField.adjoin K.carrier ({w} : Set K.closure))) (E := K.closure))
  haveI : Finite (K.absGal ⧸ N) := Finite.of_equiv _ hquot.symm.toEquiv
  haveI : IsCyclic (K.absGal ⧸ N) := (MulEquiv.isCyclic hquot).mpr hcyc
  have hsub := subgroup_eq_of_index_eq_of_isCyclic_quotient hNx hNz
    (by rw [index_fixingSubgroup_adjoin, index_fixingSubgroup_adjoin, hrank])
  have hff := congrArg IntermediateField.fixedField hsub
  rwa [InfiniteGalois.fixedField_fixingSubgroup, InfiniteGalois.fixedField_fixingSubgroup] at hff

/-- **★★★★★不分岐拡大の部分拡大は不分岐**。
分岐指数の塔の乗法性ではなく、不分岐拡大の Galois 群が巡回であることを使う。 -/
theorem isUnramifiedAdjoin_of_le (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (huy : IsUnramifiedAdjoin K y) : IsUnramifiedAdjoin K x := by
  have hmn : Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ∣ Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :=
    finrank_dvd_of_adjoin_le K hle
  have hmpos : Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ≠ 0 := Module.finrank_pos.ne'
  have hnpos : Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) ≠ 0 := Module.finrank_pos.ne'
  obtain ⟨w, hwrank, hwu, hwcyc, -⟩ := exists_isCyclic_gal K
    (Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))) hnpos
  have hyw : IntermediateField.adjoin K.carrier ({y} : Set K.closure)
      = IntermediateField.adjoin K.carrier ({w} : Set K.closure) :=
    adjoin_eq_of_isUnramified K y w huy hwu hwrank.symm
  obtain ⟨z, hzrank, hzu, -, -⟩ := exists_isUnramifiedAdjoin K
    (Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) hmpos
  have hzw : IntermediateField.adjoin K.carrier ({z} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({w} : Set K.closure) :=
    adjoin_le_of_dvd K z w hzu hwu (by rw [hzrank, hwrank]; exact hmn)
  have hxw : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({w} : Set K.closure) := hyw ▸ hle
  have hxz : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      = IntermediateField.adjoin K.carrier ({z} : Set K.closure) :=
    adjoin_eq_of_le_of_finrank_eq K hxw hzw (normal_of_isUnramifiedAdjoin K w hwu) hwcyc
      hzrank.symm
  exact isUnramifiedAdjoin_of_adjoin_eq K hxz hzu

/-! ## `K^ur` の完全な特徴づけ -/

/-- **★★★★★★`K(x) ≤ K^ur ⟺ K(x)/K が不分岐`**。 -/
theorem adjoin_le_unramifiedClosure_iff (K : PAdicLocalField p) (x : K.closure) :
    IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K
      ↔ IsUnramifiedAdjoin K x := by
  constructor
  · intro h
    have hx : x ∈ unramifiedClosure K :=
      h (IntermediateField.mem_adjoin_simple_self K.carrier x)
    obtain ⟨y, huy, hxy⟩ := (mem_unramifiedClosure_iff K x).mp hx
    refine isUnramifiedAdjoin_of_le K ?_ huy
    rw [IntermediateField.adjoin_simple_le_iff]
    exact hxy
  · exact adjoin_le_unramifiedClosure K

/-- **`x ∈ K^ur ⟺ K(x)/K が不分岐`**。 -/
theorem mem_unramifiedClosure_iff_isUnramified (K : PAdicLocalField p) (x : K.closure) :
    x ∈ unramifiedClosure K ↔ IsUnramifiedAdjoin K x := by
  rw [← adjoin_le_unramifiedClosure_iff, IntermediateField.adjoin_simple_le_iff]

/-- **原文の判定条件は「固定体が `K^ur` に入る」ことと同値**——
`UnramifiedCriterion.lean` では片方向だけだったところが、部分拡大の
不分岐性で閉じた。 -/
theorem isUnramifiedAt_iff_fixedField_le (K : PAdicLocalField p) {H : Subgroup K.absGal}
    (hH : IsOpen (H : Set K.absGal)) (hne : H ≠ ⊤) :
    IsUnramifiedAt (residueCardinality p) (subgroupCorrespondence p) K H hH
      ↔ IntermediateField.fixedField H ≤ unramifiedClosure K := by
  obtain ⟨x, hx⟩ := exists_adjoin_eq_fixedField K H hH
  rw [isUnramifiedAt_iff_isUnramifiedAdjoin K hH hne x hx, ← adjoin_le_unramifiedClosure_iff, hx]

end ABC3.Found.PGC
