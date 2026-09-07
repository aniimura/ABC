import ABC3.Found.PGC.StageTransport
import ABC3.Found.PGC.LubinTateQuotientDescent

/-!
# ★★★★★★分岐フィルトレーションの**自然性**が閉じた

到達点は

```
theorem isNaturalFiltration_ramificationFiltration (p : ℕ) [Fact p.Prime] :
    IsNaturalFiltration (ramificationFiltration p)
```

であり、**仮定は 1 つも無い**。`RF` は実物
`ABC3.Found.PGC.ramificationFiltration p`
(`Found/PGC/UnramifiedBaseChangeInvariance.lean`、仮定ゼロ)である。

⇒ ★★`Skeleton/PGC/Section4.lean::theorem_4_2` の仮説
`IsNaturalFiltration` が **1 つ落ちる**。
`Found/PGC/Prop22FilteredHypothesis.lean::filteredIsoOfAlgEquiv` /
`intKbarTransportFiltered_algEquiv` の `hnat` も落ちる。

## 経路(3 ファイル)

| ファイル | 落とした量 |
|---|---|
| `Found/PGC/Prop22FilteredHypothesis.lean`(直前の波) | `v ≤ 0` を落とし `v > 0` だけにした |
| `Found/PGC/StageNaturality.lean` | 逆極限を落とし**段データ 1 段**(`IsNaturalStage`)にした |
| `Found/PGC/StageTransport.lean` | 絶対 Galois 群を落とし**有限段**(`IsNaturalStageUpper`)にした |
| ★本ファイル | **有限段を埋めた** |

## ★★本ファイルの中身(すべて既存の在庫の組み合わせ)

1. `norm_stageFieldEquiv` —— `extendToClosure` はスペクトルノルムを保つ
   (`Found/PGC/Prop22FixedForm.lean::norm_extendToClosure`)。
   ★中間体のノルムは `K̄` のノルムの制限なので **1 層の `rfl`** で降りる。
2. `stageIntegersEquiv` —— したがって `𝒪_{K(x)} ≃+* 𝒪_{K'(ᾱx)}`
   (`adjoinIntegers` は `{y | ‖y‖ ≤ 1}` なので、ノルム保存がそのまま効く)。
3. `map_inertiaGalAdjoin` —— 惰性群が対応する。★これは
   **Corollary 1.3**(`inertia_recoverable_real`、無条件)と
   `stageGalEquiv_restrictNormalHom` の合成**だけ**である
   (`inertiaGal K L = I_K の像`)。
4. `map_stageUpperRamification` —— `ramIndex_ringHom_eq`(同変かつ付値を保つ環準同型は
   `i(σ)` を保つ)＋ `map_upperRamificationGroup_eq`(`i` が対応すれば `G^m` も対応する)。
   ★一意化元は `ᾱ` で写したものをそのまま使う
   (`stageUpperRamification_eq_map_upperRamificationGroup` が素元の選択を自由にしている)。

## ★★★原典より短い道(名指し)

★**原典 §4 は「自然な射」を局所類体論(Artin 写像)の自然性から述べるが、
本ファイルは局所類体論を 1 度も使っていない。** 使ったのは

* 無限 Galois 対応(`InfiniteGalois.normal_iff_isGalois`)、
* スペクトルノルムの保存、
* Herbrand 関数の全射準同型に沿った不変性(`map_upperRamificationGroup_eq`)

の 3 つだけである。★`Found/PGC/LubinTate*.lean`(57 本)は
`addVal_ringEquiv` 1 本を借りただけで、Lubin–Tate 塔そのものは使っていない。

## 逸脱の記録

1. `Skeleton/**` は 1 行も書き換えていない。
2. `IsNaturalStageUpper` は延長を `extendToClosure α` に固定している。
   `Found/PGC/StageNaturality.lean::map_gv_eq_of_ext` により
   **どの延長を選んでも `Γ_K^v` の像は同じ**なので、これは逸脱ではない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC IsLocalRing
open scoped NormedField Valued

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

variable {p : ℕ} [Fact p.Prime] {K K' : PAdicLocalField p}
  (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)

/-! ## §1 ノルムの保存 —— 中間体へ 1 層降りるだけ -/

/-- 中間体のノルムは `K̄` のノルムの制限(**1 層**の射影なので `rfl`。
2 層をまたぐと kernel が止まる —— `lean-idioms.md` #59)。 -/
theorem norm_coe_adjoin (x : K.closure)
    (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ‖z‖ = ‖(z : K.closure)‖ := rfl

/-- ★**有限段への制限もスペクトルノルムを保つ**。 -/
theorem norm_stageFieldEquiv (x : K.closure)
    (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ‖stageFieldEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x z‖ = ‖z‖ := by
  rw [norm_coe_adjoin, norm_coe_adjoin, stageFieldEquiv_coe]
  exact norm_extendToClosure α _

theorem norm_stageFieldEquiv_symm (x : K.closure)
    (w : IntermediateField.adjoin K'.carrier
      ({extendToClosure α x} : Set K'.closure)) :
    ‖(stageFieldEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x).symm w‖ = ‖w‖ := by
  have h := norm_stageFieldEquiv α x
    ((stageFieldEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x).symm w)
  rw [RingEquiv.apply_symm_apply] at h
  exact h.symm

/-! ## §2 整数環の輸送 -/

/-- ★★**`ᾱ` は整数環の同型 `𝒪_{K(x)} ≃+* 𝒪_{K'(ᾱx)}` を誘導する**
——`adjoinIntegers` が `{y | ‖y‖ ≤ 1}` だから、ノルム保存がそのまま効く。 -/
noncomputable def stageIntegersEquiv (x : K.closure) :
    adjoinIntegers K x ≃+* adjoinIntegers K' (extendToClosure α x) where
  toFun z := ⟨stageFieldEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x
      (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)), by
    show ‖_‖ ≤ 1
    rw [norm_stageFieldEquiv]
    exact z.2⟩
  invFun w := ⟨(stageFieldEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x).symm
      (w : IntermediateField.adjoin K'.carrier
        ({extendToClosure α x} : Set K'.closure)), by
    show ‖_‖ ≤ 1
    rw [norm_stageFieldEquiv_symm]
    exact w.2⟩
  left_inv _ := Subtype.ext (by simp)
  right_inv _ := Subtype.ext (by simp)
  map_mul' _ _ := Subtype.ext (by simp)
  map_add' _ _ := Subtype.ext (by simp)

@[simp] theorem stageIntegersEquiv_coe (x : K.closure) (z : adjoinIntegers K x) :
    ((stageIntegersEquiv α x z : adjoinIntegers K' (extendToClosure α x)) :
        IntermediateField.adjoin K'.carrier ({extendToClosure α x} : Set K'.closure))
      = stageFieldEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x
          (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := rfl

/-- ★整数環の同型は `Gal(L/K) ≃ Gal(L'/K')` と**同変**。 -/
theorem stageIntegersEquiv_smul (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [FiniteDimensional K'.carrier (IntermediateField.adjoin K'.carrier
      ({extendToClosure α x} : Set K'.closure))]
    (τ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (a : adjoinIntegers K x) :
    (stageGalEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x τ) •
        (stageIntegersEquiv α x a)
      = stageIntegersEquiv α x (τ • a) := by
  apply Subtype.ext
  rw [coe_smul_adjoinIntegers, stageIntegersEquiv_coe, stageIntegersEquiv_coe,
    coe_smul_adjoinIntegers]
  show stageFieldEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x
      (τ ((stageFieldEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x).symm
        (stageFieldEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x
          (a : IntermediateField.adjoin K.carrier ({x} : Set K.closure))))) = _
  rw [RingEquiv.symm_apply_apply]

/-! ## §3 惰性群の輸送 —— Corollary 1.3 の押し出しだけ -/

/-- ★★**惰性群が対応する**。

`inertiaGal K L = I_K の像`(`InertiaReduction.lean`)なので、
Corollary 1.3(`inertia_recoverable_real`、無条件)を
`stageGalEquiv_restrictNormalHom` で押し出すだけである。 -/
theorem map_inertiaGalAdjoin (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K'.carrier (IntermediateField.adjoin K'.carrier
      ({extendToClosure α x} : Set K'.closure))] :
    Subgroup.map ((stageGalEquiv α (extendToClosure α)
        (extendToClosure_algebraMap α) x : _ ≃* _) : _ →* _) (inertiaGalAdjoin K x)
      = inertiaGalAdjoin K' (extendToClosure α x) := by
  show Subgroup.map _ (Subgroup.map (AlgEquiv.restrictNormalHom (F := K.carrier)
      (K₁ := K.closure) ((IntermediateField.adjoin K.carrier
        ({x} : Set K.closure)) : Type _)) (absInertia K)) = _
  rw [Subgroup.map_map]
  have hcomp : (((stageGalEquiv α (extendToClosure α)
        (extendToClosure_algebraMap α) x : _ ≃* _) : _ →* _).comp
      (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : Type _)))
      = (AlgEquiv.restrictNormalHom (F := K'.carrier) (K₁ := K'.closure)
        ((IntermediateField.adjoin K'.carrier
          ({extendToClosure α x} : Set K'.closure)) : Type _)).comp
        ((galMulEquivOf α (extendToClosure α) (extendToClosure_algebraMap α) :
          K.absGal ≃* K'.absGal) : K.absGal →* K'.absGal) := by
    apply MonoidHom.ext
    intro σ
    exact stageGalEquiv_restrictNormalHom α (extendToClosure α)
      (extendToClosure_algebraMap α) x σ
  rw [hcomp, ← Subgroup.map_map]
  have hin : Subgroup.map ((galMulEquivOf α (extendToClosure α)
      (extendToClosure_algebraMap α) : K.absGal ≃* K'.absGal) : K.absGal →* K'.absGal)
      (absInertia K) = absInertia K' := by
    rw [← galContinuousMulEquiv_toMulEquiv_eq α, absInertia_eq_inertia, absInertia_eq_inertia]
    exact inertia_recoverable_real (galContinuousMulEquiv α)
  rw [hin]
  rfl

def map_inertiaGalAdjoin.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Corollary 1.3", sectionId := "cor-1-3" }

/-- 惰性群の間の同型。 -/
noncomputable def stageInertiaEquiv (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K'.carrier (IntermediateField.adjoin K'.carrier
      ({extendToClosure α x} : Set K'.closure))] :
    ↥(inertiaGalAdjoin K x) ≃* ↥(inertiaGalAdjoin K' (extendToClosure α x)) :=
  (MulEquiv.subgroupMap (stageGalEquiv α (extendToClosure α)
      (extendToClosure_algebraMap α) x) (inertiaGalAdjoin K x)).trans
    (MulEquiv.subgroupCongr (map_inertiaGalAdjoin α x))

@[simp] theorem stageInertiaEquiv_coe (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K'.carrier (IntermediateField.adjoin K'.carrier
      ({extendToClosure α x} : Set K'.closure))]
    (σ : ↥(inertiaGalAdjoin K x)) :
    ((stageInertiaEquiv α x σ : ↥(inertiaGalAdjoin K' (extendToClosure α x))) :
        (IntermediateField.adjoin K'.carrier ({extendToClosure α x} : Set K'.closure))
          ≃ₐ[K'.carrier] (IntermediateField.adjoin K'.carrier
            ({extendToClosure α x} : Set K'.closure)))
      = stageGalEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x
          (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) := rfl

theorem stageGal_comp_subtype (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K'.carrier (IntermediateField.adjoin K'.carrier
      ({extendToClosure α x} : Set K'.closure))] :
    ((stageGalEquiv α (extendToClosure α)
        (extendToClosure_algebraMap α) x : _ ≃* _) : _ →* _).comp
        (inertiaGalAdjoin K x).subtype
      = (inertiaGalAdjoin K' (extendToClosure α x)).subtype.comp
          ((stageInertiaEquiv α x : _ ≃* _) : _ →* _) :=
  MonoidHom.ext fun _ => rfl

/-! ## §4 ★★★★有限段の上付き分岐群の自然性 -/

/-- ★★★★**有限段の上付き分岐群の自然性**。

`ramIndex_ringHom_eq`(同変かつ付値を保つ環準同型は `i(σ)` を保つ)＋
`map_upperRamificationGroup_eq`(`i` が対応すれば `G^m` も対応する)。 -/
theorem map_stageUpperRamification (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K'.carrier (IntermediateField.adjoin K'.carrier
      ({extendToClosure α x} : Set K'.closure))] (v : ℝ) :
    Subgroup.map ((stageGalEquiv α (extendToClosure α)
        (extendToClosure_algebraMap α) x : _ ≃* _) : _ →* _) (stageUpperRamification K x v)
      = stageUpperRamification K' (extendToClosure α x) v := by
  have hirr : Irreducible (stageIntegersEquiv α x (stageUniformizer K x)) := by
    refine irreducible_of_addVal_eq_one ?_
    rw [addVal_ringEquiv (stageIntegersEquiv α x) (stageUniformizer K x)]
    exact IsDiscreteValuationRing.addVal_uniformizer (irreducible_stageUniformizer K)
  have hspan : maximalIdeal (adjoinIntegers K' (extendToClosure α x))
      = Ideal.span {stageIntegersEquiv α x (stageUniformizer K x)} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer _).mp hirr
  rw [stageUpperRamification_eq_map_upperRamificationGroup K' (extendToClosure α x) hspan v]
  show Subgroup.map _ (Subgroup.map (inertiaGalAdjoin K x).subtype
      (upperRamificationGroup ↥(inertiaGalAdjoin K x) (stageUniformizer K x) v)) = _
  rw [Subgroup.map_map, stageGal_comp_subtype, ← Subgroup.map_map]
  congr 1
  refine map_upperRamificationGroup_eq ?_ ?_ v
  · exact (stageInertiaEquiv α x).surjective
  intro σ
  have hval := addVal_ringEquiv (stageIntegersEquiv α x)
  have hequiv : ∀ (σ' : ↥(inertiaGalAdjoin K' (extendToClosure α x)))
      (a : adjoinIntegers K x),
      σ' • (stageIntegersEquiv α x a)
        = stageIntegersEquiv α x (((((stageInertiaEquiv α x).symm :
            ↥(inertiaGalAdjoin K' (extendToClosure α x)) ≃* ↥(inertiaGalAdjoin K x)) :
            ↥(inertiaGalAdjoin K' (extendToClosure α x)) →* ↥(inertiaGalAdjoin K x))
            σ') • a) := by
    intro σ' a
    have h := stageIntegersEquiv_smul α x
      ((((stageInertiaEquiv α x).symm σ' : ↥(inertiaGalAdjoin K x)) :
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
          ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) a
    have h2 : stageGalEquiv α (extendToClosure α) (extendToClosure_algebraMap α) x
        (((stageInertiaEquiv α x).symm σ' : ↥(inertiaGalAdjoin K x)) :
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
        = ((σ' : ↥(inertiaGalAdjoin K' (extendToClosure α x))) :
          (IntermediateField.adjoin K'.carrier
              ({extendToClosure α x} : Set K'.closure))
            ≃ₐ[K'.carrier] (IntermediateField.adjoin K'.carrier
              ({extendToClosure α x} : Set K'.closure))) := by
      rw [← stageInertiaEquiv_coe α x ((stageInertiaEquiv α x).symm σ'),
        MulEquiv.apply_symm_apply]
    rw [h2] at h
    exact h
  have key := ramIndex_ringHom_eq
    ((stageIntegersEquiv α x : adjoinIntegers K x ≃+* adjoinIntegers K' (extendToClosure α x)) :
      adjoinIntegers K x →+* adjoinIntegers K' (extendToClosure α x))
    (((stageInertiaEquiv α x).symm :
        ↥(inertiaGalAdjoin K' (extendToClosure α x)) ≃* ↥(inertiaGalAdjoin K x)) :
      ↥(inertiaGalAdjoin K' (extendToClosure α x)) →* ↥(inertiaGalAdjoin K x))
    hval hequiv (stageUniformizer K x) (stageInertiaEquiv α x σ)
  have h3 : (((stageInertiaEquiv α x).symm :
        ↥(inertiaGalAdjoin K' (extendToClosure α x)) ≃* ↥(inertiaGalAdjoin K x)) :
      ↥(inertiaGalAdjoin K' (extendToClosure α x)) →* ↥(inertiaGalAdjoin K x))
      (stageInertiaEquiv α x σ) = σ := by
    show (stageInertiaEquiv α x).symm (stageInertiaEquiv α x σ) = σ
    rw [MulEquiv.symm_apply_apply]
  rw [h3] at key
  exact key.symm

def map_stageUpperRamification.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-! ## §5 ★★★★★★到達点 -/

variable (p) in
/-- ★★★**`IsNaturalStageUpper p` は無条件に成り立つ**。 -/
theorem isNaturalStageUpper_holds : IsNaturalStageUpper p :=
  by
  intro K K' β x i1 i2 v
  haveI := i1
  haveI := i2
  exact map_stageUpperRamification β x v

variable (p) in
/-- ★★★★★★**[pGC] 分岐フィルトレーションの自然性** ——
`ℚ_p`-代数同型 `α : K ≃ K'` から誘導される `Γ_K ≃ₜ* Γ_{K'}` は、
**実物の**高次分岐フィルトレーション `Γ_K^v` を `Γ_{K'}^v` に写す。

★★**仮定は 1 つも無い。**
★これで `Skeleton/PGC/Section4.lean::theorem_4_2` の仮説
`IsNaturalFiltration` が落ちる。 -/
theorem isNaturalFiltration_ramificationFiltration :
    IsNaturalFiltration (ramificationFiltration p) :=
  isNaturalFiltration_of_isNaturalStageUpper (isNaturalStageUpper_holds p)

def isNaturalFiltration_ramificationFiltration.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-- ★★**自然な射(filtered group の同型)が無条件に作れる**。 -/
noncomputable def naturalFilteredIsoOfAlgEquiv {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    FilteredGroup.Iso (filtOf (ramificationFiltration p) K)
      (filtOf (ramificationFiltration p) K') :=
  naturalFilteredIso (ramificationFiltration p) (isNaturalFiltration_ramificationFiltration p) β

/-- ★★**Theorem 4.2 の「自然な射」が無条件に作れる**。 -/
noncomputable def naturalOuterIsoOfAlgEquiv {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    FilteredGroup.OuterIso (filtOf (ramificationFiltration p) K)
      (filtOf (ramificationFiltration p) K') :=
  naturalOuterIso (ramificationFiltration p) (isNaturalFiltration_ramificationFiltration p) β

def naturalOuterIsoOfAlgEquiv.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-- ★★**濾過つきの同型が体の同型から無条件に作れる**
(`Found/PGC/Prop22FilteredHypothesis.lean::filteredIsoOfAlgEquiv` の仮定なし版)。 -/
noncomputable def filteredIsoOfAlgEquiv' {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K') :=
  filteredIsoOfAlgEquiv (isNaturalFiltration_ramificationFiltration p) β

/-- ★★体の同型から来る `α` では、濾過つきでも `𝒪_{K̄}` が移送される(仮定なし)。 -/
theorem intKbarTransportFiltered_algEquiv' {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    IntKbarTransportFiltered (filteredIsoOfAlgEquiv' β) :=
  intKbarTransportFiltered_algEquiv (isNaturalFiltration_ramificationFiltration p) β

#print axioms norm_stageFieldEquiv
#print axioms stageIntegersEquiv
#print axioms stageIntegersEquiv_smul
#print axioms map_inertiaGalAdjoin
#print axioms stageInertiaEquiv
#print axioms stageGal_comp_subtype
#print axioms map_stageUpperRamification
#print axioms isNaturalStageUpper_holds
#print axioms isNaturalFiltration_ramificationFiltration
#print axioms naturalOuterIsoOfAlgEquiv
#print axioms filteredIsoOfAlgEquiv'
#print axioms intKbarTransportFiltered_algEquiv'

end ABC3.Found.PGC
