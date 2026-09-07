import ABC3.Found.PGC.StageNaturality

/-!
# 有限段の輸送 —— `ᾱ` は `K(x)` を `K'(ᾱ x)` に運ぶ

`Found/PGC/StageNaturality.lean` は
`IsNaturalFiltration (ramificationFiltration p)` を **`IsNaturalStage p`**
(開正規部分群 1 つぶんの主張)に落とした。本ファイルはその `IsNaturalStage` を
さらに **有限段だけの主張** に落とす。

## 何が入ったか

* §0 **抽象核(純群論)** `map_comap_eq_comap_map` ——
  `f' ∘ Φ = ρ ∘ f` なる同型 `Φ`, `ρ` と準同型 `f`, `f'` があれば
  `map Φ (comap f H) = comap f' (map ρ H)`。★全射性も要らない。
* §1 **中間体の輸送** `map_adjoin_toSubfield` / `mem_adjoin_map_iff` ——
  `ᾱ` は `K(x)` を `K'(ᾱ x)` の上へ写す。中身は `Subfield.closure` の
  Galois 接続だけ(`IntermediateField.adjoin_toSubfield`)。
* §2 **有限段の体同型と Galois 群同型** `stageFieldEquiv` / `stageGalEquiv`
  および `stageGalEquiv_restrictNormalHom`
  ——★`Gal(K̄/K) → Gal(L/K)` の制限射と可換になる。
* §3 **段生成元の輸送** `fixingSubgroup_adjoin_map` / `stageGeneratorMap`
  —— `g : StageGenerator K N` から `StageGenerator K' (Φ N)` を作る。
  ★`Normal K'(ᾱx)/K'` は無限 Galois 対応(`InfiniteGalois.normal_iff_isGalois`)で出る。
* §4 ★★★**到達点** `isNaturalFiltration_of_isNaturalStageUpper` ——

```
IsNaturalStageUpper p → IsNaturalFiltration (ramificationFiltration p)
```

## ★★★残った穴 —— ★**同じシフトの `Found/PGC/StageUpperNaturality.lean` で埋まった**

★★`IsNaturalStageUpper p` は `isNaturalStageUpper_holds` として**無条件に証明済み**である。
以下は本ファイルを書いた時点での見通しで、実際その 4 手順どおりになった。

`IsNaturalStageUpper p` とは

> `stageGalEquiv α ᾱ hfwd x : Gal(K(x)/K) ≃ Gal(K'(ᾱx)/K')` が
> `stageUpperRamification K x v` を `stageUpperRamification K' (ᾱ x) v` に写す

という **有限次拡大 1 つぶんの主張**であり、絶対 Galois 群も逆極限も
開正規部分群も**もう出てこない**。★`ᾱ` は `extendToClosure α` に固定してあるので、
`Found/PGC/Prop22FixedForm.lean::norm_extendToClosure`(延長はスペクトルノルムを保つ)が
そのまま使える。

★次のノードの手順(測っていない、**着手していない**):

1. `norm_extendToClosure` から `ᾱ` が `𝒪_{K(x)} ≃+* 𝒪_{K'(ᾱx)}` を誘導すること
   (`adjoinIntegers` は `{y | ‖y‖ ≤ 1}` なので、ノルム保存がそのまま効く)。
2. その環同型が `stageGalEquiv` と同変であること。
3. `Found/PGC/UnramifiedBaseChangeInvariance.lean::ramIndex_ringHom_eq`
   (同変かつ付値を保つ環準同型は `ramIndex` を保つ)と
   `map_upperRamificationGroup_eq`(`ramIndex` が対応すれば `G^m` も対応する)を当てる。
4. 一意化元の取り方は
   `Found/PGC/RamificationImageStage.lean::stageUpperRamification_eq_map_upperRamificationGroup`
   で自由になっているので、`ᾱ` で写した素元をそのまま使えばよい。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC

/-! ## §0 ★★抽象核(純群論) —— 分岐・付値・Galois・体が 1 語も出てこない -/

/-- ★★**抽象核(純群論)** —— 可換図式 `f' ∘ Φ = ρ ∘ f`(`Φ`, `ρ` は同型、
`f`, `f'` はただの準同型)があれば、**引き戻しは像と交換する**。

★`f`・`f'` の全射性は要らない。 -/
theorem map_comap_eq_comap_map {Γ Γ' G G' : Type*} [Group Γ] [Group Γ'] [Group G] [Group G']
    (Φ : Γ ≃* Γ') (ρ : G ≃* G') (f : Γ →* G) (f' : Γ' →* G')
    (hcomm : ∀ σ : Γ, f' (Φ σ) = ρ (f σ)) (H : Subgroup G) :
    Subgroup.map (Φ : Γ →* Γ') (Subgroup.comap f H)
      = Subgroup.comap f' (Subgroup.map (ρ : G →* G') H) := by
  apply le_antisymm
  · rintro _ ⟨σ, hσ, rfl⟩
    show f' (Φ σ) ∈ Subgroup.map (ρ : G →* G') H
    rw [hcomm σ]
    exact ⟨f σ, hσ, rfl⟩
  · intro σ' hσ'
    refine ⟨Φ.symm σ', ?_, by simp⟩
    obtain ⟨h, hh, hheq⟩ := hσ'
    show f (Φ.symm σ') ∈ H
    have h1 : ρ (f (Φ.symm σ')) = f' σ' := by
      rw [← hcomm (Φ.symm σ'), Φ.apply_symm_apply]
    have h2 : ρ (f (Φ.symm σ')) = ρ h := by rw [h1, ← hheq]; rfl
    rwa [ρ.injective h2]

variable {p : ℕ} [Fact p.Prime] {K K' : PAdicLocalField p}
  (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
  (hfwd : ∀ y : K.carrier, ᾱ (algebraMap K.carrier K.closure y)
    = algebraMap K'.carrier K'.closure (α y))

/-! ## §1 中間体の輸送 -/

include hfwd in
/-- ★**`ᾱ` は `K(x)` を `K'(ᾱ x)` の上へ写す**。

中身は `IntermediateField.adjoin_toSubfield`(`adjoin` は
`range (algebraMap) ∪ S` の生成する部分体)と `Subfield` の Galois 接続だけ。 -/
theorem map_adjoin_toSubfield (x : K.closure) :
    Subfield.map (ᾱ : K.closure →+* K'.closure)
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).toSubfield
      = (IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)).toSubfield := by
  rw [IntermediateField.adjoin_toSubfield, IntermediateField.adjoin_toSubfield]
  apply le_antisymm
  · refine (Subfield.gc_map_comap _ _ _).mpr ?_
    rw [Subfield.closure_le]
    rintro z (⟨c, rfl⟩ | hz)
    · simp only [SetLike.mem_coe, Subfield.mem_comap]
      have hc : (ᾱ : K.closure →+* K'.closure) (algebraMap K.carrier K.closure c)
          = algebraMap K'.carrier K'.closure (α c) := hfwd c
      rw [hc]
      exact Subfield.subset_closure (Or.inl ⟨α c, rfl⟩)
    · rw [Set.mem_singleton_iff] at hz
      subst hz
      simp only [SetLike.mem_coe, Subfield.mem_comap]
      exact Subfield.subset_closure (Or.inr rfl)
  · rw [Subfield.closure_le]
    rintro z (⟨c', rfl⟩ | hz)
    · simp only [SetLike.mem_coe, Subfield.mem_map]
      refine ⟨algebraMap K.carrier K.closure (α.symm c'),
        Subfield.subset_closure (Or.inl ⟨α.symm c', rfl⟩), ?_⟩
      have hc : (ᾱ : K.closure →+* K'.closure) (algebraMap K.carrier K.closure (α.symm c'))
          = algebraMap K'.carrier K'.closure (α (α.symm c')) := hfwd _
      rw [hc, AlgEquiv.apply_symm_apply]
    · rw [Set.mem_singleton_iff] at hz
      subst hz
      simp only [SetLike.mem_coe, Subfield.mem_map]
      exact ⟨x, Subfield.subset_closure (Or.inr rfl), rfl⟩

include hfwd in
/-- ★`ᾱ y ∈ K'(ᾱ x) ↔ y ∈ K(x)`。 -/
theorem mem_adjoin_map_iff (x y : K.closure) :
    ᾱ y ∈ IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)
      ↔ y ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  have h := map_adjoin_toSubfield α ᾱ hfwd x
  constructor
  · intro hy
    have hm : ᾱ y ∈ Subfield.map (ᾱ : K.closure →+* K'.closure)
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).toSubfield := by
      rw [h]; exact hy
    obtain ⟨z, hz, hzy⟩ := hm
    have hzy' : z = y := ᾱ.injective hzy
    exact hzy' ▸ hz
  · intro hy
    have hm : ᾱ y ∈ Subfield.map (ᾱ : K.closure →+* K'.closure)
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).toSubfield := ⟨y, hy, rfl⟩
    rw [h] at hm
    exact hm

/-! ## §2 有限段の体同型と Galois 群同型 -/

include hfwd in
/-- ★**`ᾱ` の有限段への制限** `K(x) ≃+* K'(ᾱ x)`(`α` 上半線形)。 -/
noncomputable def stageFieldEquiv (x : K.closure) :
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃+* (IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)) where
  toFun z := ⟨ᾱ (z : K.closure), (mem_adjoin_map_iff α ᾱ hfwd x _).mpr z.2⟩
  invFun w := ⟨ᾱ.symm (w : K'.closure), by
    refine (mem_adjoin_map_iff α ᾱ hfwd x _).mp ?_
    rw [RingEquiv.apply_symm_apply]
    exact w.2⟩
  left_inv _ := Subtype.ext (by simp)
  right_inv _ := Subtype.ext (by simp)
  map_mul' _ _ := Subtype.ext (by simp)
  map_add' _ _ := Subtype.ext (by simp)

@[simp] theorem stageFieldEquiv_coe (x : K.closure)
    (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ((stageFieldEquiv α ᾱ hfwd x z : IntermediateField.adjoin K'.carrier
        ({ᾱ x} : Set K'.closure)) : K'.closure) = ᾱ (z : K.closure) := rfl

@[simp] theorem stageFieldEquiv_symm_coe (x : K.closure)
    (w : IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)) :
    (((stageFieldEquiv α ᾱ hfwd x).symm w : IntermediateField.adjoin K.carrier
        ({x} : Set K.closure)) : K.closure) = ᾱ.symm (w : K'.closure) := rfl

include hfwd in
/-- 制限は基礎体の同型 `α` と両立する。 -/
theorem stageFieldEquiv_algebraMap (x : K.closure) (c : K.carrier) :
    stageFieldEquiv α ᾱ hfwd x
        (algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) c)
      = algebraMap K'.carrier
          (IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)) (α c) := by
  apply Subtype.ext
  rw [stageFieldEquiv_coe]
  show ᾱ (algebraMap K.carrier K.closure c) = algebraMap K'.carrier K'.closure (α c)
  exact hfwd c

include hfwd in
theorem stageFieldEquiv_symm_algebraMap (x : K.closure) (c' : K'.carrier) :
    (stageFieldEquiv α ᾱ hfwd x).symm
        (algebraMap K'.carrier (IntermediateField.adjoin K'.carrier
          ({ᾱ x} : Set K'.closure)) c')
      = algebraMap K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (α.symm c') := by
  apply (stageFieldEquiv α ᾱ hfwd x).injective
  rw [RingEquiv.apply_symm_apply, stageFieldEquiv_algebraMap, AlgEquiv.apply_symm_apply]

include hfwd in
/-- `β τ β⁻¹` は `K'.carrier` を固定する。 -/
theorem stageConj_fixes (x : K.closure)
    (τ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (c' : K'.carrier) :
    ((stageFieldEquiv α ᾱ hfwd x).symm.trans
        (τ.toRingEquiv.trans (stageFieldEquiv α ᾱ hfwd x)))
        (algebraMap K'.carrier (IntermediateField.adjoin K'.carrier
          ({ᾱ x} : Set K'.closure)) c')
      = algebraMap K'.carrier (IntermediateField.adjoin K'.carrier
          ({ᾱ x} : Set K'.closure)) c' := by
  show stageFieldEquiv α ᾱ hfwd x (τ ((stageFieldEquiv α ᾱ hfwd x).symm _)) = _
  rw [stageFieldEquiv_symm_algebraMap α ᾱ hfwd x c']
  show stageFieldEquiv α ᾱ hfwd x (τ (algebraMap K.carrier _ (α.symm c'))) = _
  rw [τ.commutes, stageFieldEquiv_algebraMap α ᾱ hfwd x, AlgEquiv.apply_symm_apply]

include hfwd in
/-- `β⁻¹ τ' β` は `K.carrier` を固定する。 -/
theorem stageConj_symm_fixes (x : K.closure)
    (τ' : (IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure))
      ≃ₐ[K'.carrier] (IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)))
    (c : K.carrier) :
    ((stageFieldEquiv α ᾱ hfwd x).trans
        (τ'.toRingEquiv.trans (stageFieldEquiv α ᾱ hfwd x).symm))
        (algebraMap K.carrier (IntermediateField.adjoin K.carrier
          ({x} : Set K.closure)) c)
      = algebraMap K.carrier (IntermediateField.adjoin K.carrier
          ({x} : Set K.closure)) c := by
  show (stageFieldEquiv α ᾱ hfwd x).symm (τ' (stageFieldEquiv α ᾱ hfwd x _)) = _
  rw [stageFieldEquiv_algebraMap α ᾱ hfwd x c, τ'.commutes,
    stageFieldEquiv_symm_algebraMap α ᾱ hfwd x, AlgEquiv.symm_apply_apply]

include hfwd in
/-- ★★**有限段の Galois 群の同型** `Gal(K(x)/K) ≃* Gal(K'(ᾱx)/K')`。 -/
noncomputable def stageGalEquiv (x : K.closure) :
    ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      ≃* ((IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure))
        ≃ₐ[K'.carrier] (IntermediateField.adjoin K'.carrier
          ({ᾱ x} : Set K'.closure))) where
  toFun τ := AlgEquiv.ofRingEquiv (f := (stageFieldEquiv α ᾱ hfwd x).symm.trans
    (τ.toRingEquiv.trans (stageFieldEquiv α ᾱ hfwd x))) (stageConj_fixes α ᾱ hfwd x τ)
  invFun τ' := AlgEquiv.ofRingEquiv (f := (stageFieldEquiv α ᾱ hfwd x).trans
    (τ'.toRingEquiv.trans (stageFieldEquiv α ᾱ hfwd x).symm))
    (stageConj_symm_fixes α ᾱ hfwd x τ')
  left_inv τ := by
    apply AlgEquiv.ext
    intro z
    show (stageFieldEquiv α ᾱ hfwd x).symm
      (stageFieldEquiv α ᾱ hfwd x (τ ((stageFieldEquiv α ᾱ hfwd x).symm
        (stageFieldEquiv α ᾱ hfwd x z)))) = τ z
    rw [RingEquiv.symm_apply_apply, RingEquiv.symm_apply_apply]
  right_inv τ' := by
    apply AlgEquiv.ext
    intro w
    show stageFieldEquiv α ᾱ hfwd x ((stageFieldEquiv α ᾱ hfwd x).symm
      (τ' (stageFieldEquiv α ᾱ hfwd x ((stageFieldEquiv α ᾱ hfwd x).symm w)))) = τ' w
    rw [RingEquiv.apply_symm_apply, RingEquiv.apply_symm_apply]
  map_mul' τ₁ τ₂ := by
    apply AlgEquiv.ext
    intro w
    show stageFieldEquiv α ᾱ hfwd x ((τ₁ * τ₂) ((stageFieldEquiv α ᾱ hfwd x).symm w))
      = stageFieldEquiv α ᾱ hfwd x (τ₁ ((stageFieldEquiv α ᾱ hfwd x).symm
          (stageFieldEquiv α ᾱ hfwd x (τ₂ ((stageFieldEquiv α ᾱ hfwd x).symm w)))))
    rw [RingEquiv.symm_apply_apply]
    rfl

@[simp] theorem stageGalEquiv_apply_coe (x : K.closure)
    (τ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (w : IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)) :
    ((stageGalEquiv α ᾱ hfwd x τ w : IntermediateField.adjoin K'.carrier
        ({ᾱ x} : Set K'.closure)) : K'.closure)
      = ᾱ ((τ ((stageFieldEquiv α ᾱ hfwd x).symm w) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) := rfl

include hfwd in
/-- ★★★**制限射との可換性** —— `Gal(K̄/K) → Gal(L/K)` と `Gal(K̄'/K') → Gal(L'/K')`
は `galMulEquivOf` / `stageGalEquiv` で繋がる。

★これが「有限段だけの主張に落とす」ための鍵である。 -/
theorem stageGalEquiv_restrictNormalHom (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K'.carrier (IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure))]
    (σ : K.absGal) :
    stageGalEquiv α ᾱ hfwd x
        (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
          ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : Type _) σ)
      = AlgEquiv.restrictNormalHom (F := K'.carrier) (K₁ := K'.closure)
          ((IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)) : Type _)
          (galMulEquivOf α ᾱ hfwd σ) := by
  apply AlgEquiv.ext
  intro w
  apply Subtype.ext
  have h1 := AlgEquiv.restrictNormal_commutes σ
    ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : Type _)
    ((stageFieldEquiv α ᾱ hfwd x).symm w)
  have h2 := AlgEquiv.restrictNormal_commutes (galMulEquivOf α ᾱ hfwd σ)
    ((IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)) : Type _) w
  simp only [IntermediateField.algebraMap_apply] at h1 h2
  have hh1 : (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : Type _) σ)
      = AlgEquiv.restrictNormal σ
        ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : Type _) := rfl
  have hh2 : (AlgEquiv.restrictNormalHom (F := K'.carrier) (K₁ := K'.closure)
      ((IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)) : Type _)
      (galMulEquivOf α ᾱ hfwd σ))
      = AlgEquiv.restrictNormal (galMulEquivOf α ᾱ hfwd σ)
        ((IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)) : Type _) := rfl
  rw [stageGalEquiv_apply_coe, hh1, hh2, h1, h2, stageFieldEquiv_symm_coe]
  rfl

/-! ## §3 段生成元の輸送 -/

include hfwd in
/-- ★★**固定化部分群の輸送** —— `K'(ᾱ x)` の固定化部分群は `K(x)` のそれの像。 -/
theorem fixingSubgroup_adjoin_map (x : K.closure) :
    (IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure)).fixingSubgroup
      = Subgroup.map ((galMulEquivOf α ᾱ hfwd : K.absGal ≃* K'.absGal) : K.absGal →* K'.absGal)
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup := by
  ext σ'
  constructor
  · intro hσ'
    simp only [IntermediateField.mem_fixingSubgroup_iff] at hσ'
    refine ⟨(galMulEquivOf α ᾱ hfwd).symm σ', ?_, by simp⟩
    simp only [SetLike.mem_coe, IntermediateField.mem_fixingSubgroup_iff]
    intro y hy
    have hy' : ᾱ y ∈ IntermediateField.adjoin K'.carrier ({ᾱ x} : Set K'.closure) :=
      (mem_adjoin_map_iff α ᾱ hfwd x y).mpr hy
    have hval : ((galMulEquivOf α ᾱ hfwd).symm σ') y = ᾱ.symm (σ' (ᾱ y)) := by
      show ᾱ.symm (σ' (ᾱ.symm.symm y)) = ᾱ.symm (σ' (ᾱ y))
      rw [RingEquiv.symm_symm]
    rw [hval, hσ' _ hy', RingEquiv.symm_apply_apply]
  · rintro ⟨σ, hσ, rfl⟩
    simp only [SetLike.mem_coe, IntermediateField.mem_fixingSubgroup_iff] at hσ
    simp only [IntermediateField.mem_fixingSubgroup_iff]
    intro y' hy'
    have hy : ᾱ.symm y' ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
      refine (mem_adjoin_map_iff α ᾱ hfwd x _).mp ?_
      rw [RingEquiv.apply_symm_apply]
      exact hy'
    show ᾱ (σ (ᾱ.symm y')) = y'
    rw [hσ _ hy, RingEquiv.apply_symm_apply]

include hfwd in
/-- ★★★**段生成元の輸送** —— `x` の段生成元から `ᾱ x` の段生成元を作る。

`finiteDimensional` は `ᾱ x` が `K'` 上整であることから、`normal` は
無限 Galois 対応(`InfiniteGalois.normal_iff_isGalois`)から出る。 -/
noncomputable def stageGeneratorMap {N : Subgroup K.absGal} (hN : N.Normal)
    (g : StageGenerator K N) :
    StageGenerator K' (Subgroup.map ((galMulEquivOf α ᾱ hfwd : K.absGal ≃* K'.absGal) :
      K.absGal →* K'.absGal) N) where
  gen := ᾱ g.gen
  finiteDimensional :=
    IntermediateField.adjoin.finiteDimensional (Algebra.IsIntegral.isIntegral (ᾱ g.gen))
  normal := by
    haveI := isGalois_closure K'
    haveI hfix : ((IntermediateField.adjoin K'.carrier
        ({ᾱ g.gen} : Set K'.closure)).fixingSubgroup).Normal := by
      rw [fixingSubgroup_adjoin_map α ᾱ hfwd g.gen, g.fixingSubgroup_eq]
      exact Subgroup.Normal.map hN _ (galMulEquivOf α ᾱ hfwd).surjective
    haveI : IsGalois K'.carrier
        (IntermediateField.adjoin K'.carrier ({ᾱ g.gen} : Set K'.closure)) :=
      (InfiniteGalois.normal_iff_isGalois _).mp hfix
    infer_instance
  fixingSubgroup_eq := by
    rw [fixingSubgroup_adjoin_map α ᾱ hfwd g.gen, g.fixingSubgroup_eq]

@[simp] theorem stageGeneratorMap_gen {N : Subgroup K.absGal} (hN : N.Normal)
    (g : StageGenerator K N) : (stageGeneratorMap α ᾱ hfwd hN g).gen = ᾱ g.gen := rfl

/-! ## §4 ★★★残った穴を有限段だけの主張にする -/

/-- ★★★★**残った唯一の穴** —— 有限段の**上付き分岐群**の自然性。

`Gal(K(x)/K) ≃ Gal(K'(ᾱ x)/K')`(`stageGalEquiv`)が
`Gal(K(x)/K(x)_0)^v` を `Gal(K'(ᾱx)/K'(ᾱx)_0)^v` に写す、という主張。

★★分岐の言葉しか出てこない —— 絶対 Galois 群も逆極限も開正規部分群も消えている。 -/
def IsNaturalStageUpper (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K'.carrier (IntermediateField.adjoin K'.carrier
      ({extendToClosure α x} : Set K'.closure))]
    (v : ℝ),
      Subgroup.map (((stageGalEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x :
          ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
              ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
            ≃* ((IntermediateField.adjoin K'.carrier
                ({extendToClosure α x} : Set K'.closure))
              ≃ₐ[K'.carrier] (IntermediateField.adjoin K'.carrier
                ({extendToClosure α x} : Set K'.closure)))) :
          ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
              ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
            →* ((IntermediateField.adjoin K'.carrier
                ({extendToClosure α x} : Set K'.closure))
              ≃ₐ[K'.carrier] (IntermediateField.adjoin K'.carrier
                ({extendToClosure α x} : Set K'.closure)))))
        (stageUpperRamification K x v)
      = stageUpperRamification K' (extendToClosure α x) v

def IsNaturalStageUpper.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

theorem galContinuousMulEquiv_toMulEquiv_eq (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    ((galContinuousMulEquiv α).toMulEquiv : K.absGal ≃* K'.absGal)
      = galMulEquivOf α (extendToClosure α) (extendToClosure_algebraMap α) := rfl

/-- ★★★**段の輸送** —— 有限段の上付き分岐群が対応すれば、段データも対応する。

中身は §0 の抽象核 `map_comap_eq_comap_map` に
`stageGalEquiv_restrictNormalHom` を差し込むだけ。 -/
theorem map_stage_eq (h : IsNaturalStageUpper p) {N : Subgroup K.absGal} (hN : N.Normal)
    (g : StageGenerator K N) (v : ℝ) :
    Subgroup.map ((galMulEquivOf α (extendToClosure α) (extendToClosure_algebraMap α) :
        K.absGal ≃* K'.absGal) : K.absGal →* K'.absGal) (g.stage v)
      = (stageGeneratorMap α (extendToClosure α) (extendToClosure_algebraMap α) hN g).stage v := by
  letI := g.finiteDimensional
  letI := g.normal
  letI : FiniteDimensional K'.carrier (IntermediateField.adjoin K'.carrier
      ({extendToClosure α g.gen} : Set K'.closure)) :=
    (stageGeneratorMap α (extendToClosure α) (extendToClosure_algebraMap α) hN g).finiteDimensional
  letI : Normal K'.carrier (IntermediateField.adjoin K'.carrier
      ({extendToClosure α g.gen} : Set K'.closure)) :=
    (stageGeneratorMap α (extendToClosure α) (extendToClosure_algebraMap α) hN g).normal
  show Subgroup.map ((galMulEquivOf α (extendToClosure α) (extendToClosure_algebraMap α) :
        K.absGal ≃* K'.absGal) : K.absGal →* K'.absGal)
      (Subgroup.comap (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        ((IntermediateField.adjoin K.carrier ({g.gen} : Set K.closure)) : Type _))
        (stageUpperRamification K g.gen v))
    = Subgroup.comap (AlgEquiv.restrictNormalHom (F := K'.carrier) (K₁ := K'.closure)
        ((IntermediateField.adjoin K'.carrier
          ({extendToClosure α g.gen} : Set K'.closure)) : Type _))
        (stageUpperRamification K' (extendToClosure α g.gen) v)
  rw [map_comap_eq_comap_map (galMulEquivOf α (extendToClosure α) (extendToClosure_algebraMap α))
    (stageGalEquiv α (extendToClosure α) (extendToClosure_algebraMap α) g.gen) _ _
    (fun σ => (stageGalEquiv_restrictNormalHom α (extendToClosure α)
      (extendToClosure_algebraMap α) g.gen σ).symm),
    h α g.gen v]

/-- ★★★★**有限段の自然性 ⟹ 段データの自然性**。

`Found/PGC/StageNaturality.lean::isNaturalFiltration_of_isNaturalStage` と繋ぐと
`IsNaturalStageUpper p → IsNaturalFiltration (ramificationFiltration p)`。 -/
theorem isNaturalStage_of_isNaturalStageUpper (h : IsNaturalStageUpper p) :
    IsNaturalStage p := by
  intro K K' α N hN v _
  obtain ⟨g⟩ := nonempty_stageGenerator K hN
  rw [galContinuousMulEquiv_toMulEquiv_eq α]
  rw [absGalStage_eq_stage_any K g v,
    map_stage_eq α h hN.2 g v,
    absGalStage_eq_stage_any K' (stageGeneratorMap α (extendToClosure α)
      (extendToClosure_algebraMap α) hN.2 g) v]

/-- ★★★★★**到達点** —— 有限段の上付き分岐群の自然性から
`IsNaturalFiltration (ramificationFiltration p)` が出る。 -/
theorem isNaturalFiltration_of_isNaturalStageUpper (h : IsNaturalStageUpper p) :
    IsNaturalFiltration (ramificationFiltration p) :=
  isNaturalFiltration_of_isNaturalStage (isNaturalStage_of_isNaturalStageUpper h)

def isNaturalFiltration_of_isNaturalStageUpper.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

#print axioms map_comap_eq_comap_map
#print axioms map_adjoin_toSubfield
#print axioms mem_adjoin_map_iff
#print axioms stageFieldEquiv
#print axioms stageGalEquiv
#print axioms stageGalEquiv_restrictNormalHom
#print axioms fixingSubgroup_adjoin_map
#print axioms stageGeneratorMap
#print axioms map_stage_eq
#print axioms isNaturalStage_of_isNaturalStageUpper
#print axioms isNaturalFiltration_of_isNaturalStageUpper

end ABC3.Found.PGC
