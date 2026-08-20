import ABC3.Found.Arakelov.UltraCoord

/-!
# Arakelov (C2) 第 89–90 ブロック —— **★★★★★★★★固有 ⟹ `X^arc` コンパクト**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★Chow の補題を使わずに (C2) の中身が出た

    X が ℤ 上固有 ⟹ X^arc = Hom(Spec ℂ, X) はコンパクト

★教科書の道は **Chow の補題**(固有 → 射影的な変更)で、mathlib に無い。
★★**超フィルターと付値判定法**で迂回した:

| 段 | ブロック |
|---|---|
| 超積の付値環 `O` と標準部分 | 80 |
| 超積体の点と付値判定法による持ち上げ | 81–82 |
| chart 間の翻訳の係数環自然性 | 83 |
| 転送(超積の点が `V` なら a.e. の点も `V`) | 84 |
| 基本開集合・閉点 | 85–86 |
| **超積の座標 = 座標の超積** | 87–88 |
| **収束**(位相の段)と結論 | ★本ブロック |

## ★★★収束の中身

`V` chart の座標 `b` について:

    各点の座標の超積 = 持ち上げ `l` の座標(`O` の元)
      ⟹ 有界 ⟹ 標準部分へ収束
      ⟹ その標準部分は極限点 `q` の座標

★★したがって `𝒰` は `q` に収束する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `le_nhds_affine` / `le_nhds_chart` / `le_nhds_X` | ★★座標の収束から近傍への収束 |
| `coordFun` ほか | ★点の座標(全域版) |
| `compactSpace_arc` | ★★★★★★★★**固有 ⟹ `X^arc` コンパクト** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Filter

variable {X : Scheme.{0}}

/-! ## ★★位相の段 -/

/-- ★★アフィンでは、座標の収束が近傍への収束を与える。 -/
theorem le_nhds_affine {A : CommRingCat.{0}}
    (𝒱 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ Spec A))
    (r : Spec (CommRingCat.of ℂ) ⟶ Spec A)
    (h : ∀ a : A, Tendsto (fun s => evalAffine A s a) (𝒱 : Filter _) (nhds (evalAffine A r a))) :
    (𝒱 : Filter (Spec (CommRingCat.of ℂ) ⟶ Spec A)) ≤ @nhds _ (arcTopologyAffine A) r := by
  letI := arcTopologyAffine A
  rw [nhds_induced, ← Filter.map_le_iff_le_comap]
  exact tendsto_pi_nhds.2 h

/-- ★★chart 上でも、座標の収束が近傍への収束を与える。 -/
theorem le_nhds_chart (V : X.affineOpens)
    (𝒱 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme))
    (r : Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme)
    (h : ∀ b : Γ(X, V.1), Tendsto (fun s => evalAffine Γ(X, V.1) (s ≫ V.2.isoSpec.hom) b)
        (𝒱 : Filter (Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme))
        (nhds (evalAffine Γ(X, V.1) (r ≫ V.2.isoSpec.hom) b))) :
    (𝒱 : Filter (Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme)) ≤ @nhds _ (arcTopologyOpen V) r := by
  letI := arcTopologyAffine (Γ(X, V.1))
  have hmap : (Ultrafilter.map (fun s => s ≫ V.2.isoSpec.hom) 𝒱 : Filter _)
      ≤ nhds (r ≫ V.2.isoSpec.hom) := by
    refine le_nhds_affine _ _ ?_
    intro a
    exact (tendsto_map'_iff).2 (h a)
  show (𝒱 : Filter _)
    ≤ @nhds _ (TopologicalSpace.induced (fun p => p ≫ V.2.isoSpec.hom) inferInstance) r
  rw [nhds_induced, ← Filter.map_le_iff_le_comap]
  exact hmap

/-- ★★★chart からの押し出しで `X` の近傍へ。 -/
theorem le_nhds_X (V : X.affineOpens)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X))
    (𝒱 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme))
    (hmap : Ultrafilter.map (fun s => s ≫ V.1.ι) 𝒱 = 𝒰)
    (r : Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme)
    (h : (𝒱 : Filter (Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme)) ≤ @nhds _ (arcTopologyOpen V) r) :
    (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ≤ @nhds _ (arcTopology X) (r ≫ V.1.ι) := by
  letI := arcTopology V.1.toScheme
  letI := arcTopology X
  rw [← hmap]
  rw [arcTopologyOpen_eq_arcTopology] at h
  have hcont := continuous_comp_openImmersion V.1.ι
  calc (Ultrafilter.map (fun s => s ≫ V.1.ι) 𝒱 : Filter (Spec (CommRingCat.of ℂ) ⟶ X))
      = Filter.map (fun s => s ≫ V.1.ι)
          (𝒱 : Filter (Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme)) := rfl
    _ ≤ Filter.map (fun s => s ≫ V.1.ι) (nhds r) := Filter.map_mono h
    _ ≤ nhds (r ≫ V.1.ι) := hcont.tendsto r

/-! ## ★点の座標(全域版) -/

theorem range_comp_ι_subset (V : X.Opens) (s : Spec (CommRingCat.of ℂ) ⟶ V.toScheme) :
    Set.range (s ≫ V.ι).base ⊆ (V : Set X) := by
  rintro _ ⟨z, rfl⟩
  have hz : (s ≫ V.ι).base z = V.ι.base (s.base z) := rfl
  rw [hz, ← Scheme.Opens.range_ι V]
  exact ⟨_, rfl⟩

theorem isUnit_pi_of_ne_zero {ι : Type} (g : ι → ℂ) (h : ∀ i, g i ≠ 0) : IsUnit g := by
  refine ⟨⟨g, fun i => (g i)⁻¹, ?_, ?_⟩, rfl⟩
  · funext i; exact mul_inv_cancel₀ (h i)
  · funext i; exact inv_mul_cancel₀ (h i)

open scoped Classical in
/-- ★点 `p` の `V` 座標(定義できないところでは 0)。 -/
noncomputable def coordFun (U V : X.affineOpens)
    (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme) (b : Γ(X, V.1))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : ℂ :=
  if h : Set.range (chartPoint U (CommRingCat.ofHom (arcHomOf U d p))).base ⊆ (V.1 : Set X)
  then (chartHom U V (CommRingCat.ofHom (arcHomOf U d p)) h).hom b else 0

open scoped Classical in
theorem coordFun_eq (U V : X.affineOpens) (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (b : Γ(X, V.1)) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (h : Set.range (chartPoint U (CommRingCat.ofHom (arcHomOf U d p))).base ⊆ (V.1 : Set X)) :
    coordFun U V d b p = (chartHom U V (CommRingCat.ofHom (arcHomOf U d p)) h).hom b := by
  rw [coordFun, dif_pos h]

theorem chartHom_congr {R : CommRingCat.{0}} (U V : X.affineOpens) (φ ψ : Γ(X, U.1) ⟶ R)
    (hφψ : φ = ψ) (hφ : Set.range (chartPoint U φ).base ⊆ (V.1 : Set X))
    (hψ : Set.range (chartPoint U ψ).base ⊆ (V.1 : Set X)) :
    chartHom U V φ hφ = chartHom U V ψ hψ := by
  subst hφψ
  rfl

theorem coordHom_congr {R : CommRingCat.{0}} (V : X.affineOpens) (g g' : Spec R ⟶ X)
    (hgg : g = g') (h : Set.range g.base ⊆ (V.1 : Set X))
    (h' : Set.range g'.base ⊆ (V.1 : Set X)) :
    coordHom V g h = coordHom V g' h' := by
  subst hgg
  rfl

theorem germ_ne_zero_eventually {α : Type} (𝒰 : Ultrafilter α) (g : α → ℂ)
    (h : ((g : α → ℂ) : Germ (𝒰 : Filter α) ℂ) ≠ 0) : ∀ᶠ a in (𝒰 : Filter α), g a ≠ 0 := by
  refine Ultrafilter.eventually_not.2 ?_
  intro hz
  apply h
  show ((g : α → ℂ) : Germ (𝒰 : Filter α) ℂ) = ((0 : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
  rw [Germ.coe_eq]
  filter_upwards [hz] with a ha
  exact ha

open scoped Classical in
theorem germRestrict_apply {α : Type} (𝒰 : Ultrafilter α) (T : Set α) (hT : T ∈ 𝒰)
    (g : ↥T → ℂ) :
    germRestrict 𝒰 T hT g
      = (((fun a => if h : a ∈ T then g ⟨a, h⟩ else 0) : α → ℂ) :
          Germ (𝒰 : Filter α) ℂ) := rfl

theorem evalAffine_eq_coordHom (V : X.affineOpens) (s : Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme)
    (h : Set.range (s ≫ V.1.ι).base ⊆ (V.1 : Set X)) (b : Γ(X, V.1)) :
    evalAffine Γ(X, V.1) (s ≫ V.2.isoSpec.hom) b = (coordHom V (s ≫ V.1.ι) h).hom b := by
  unfold coordHom evalAffine evalHom
  have hs : liftV V.1 (s ≫ V.1.ι) h = s := by
    apply (cancel_mono V.1.ι).1
    rw [liftV_fac]
  rw [hs]

theorem coordFun_eq_coordHom (U V : X.affineOpens) (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hpU : arcPt p ∈ U.1) (b : Γ(X, V.1))
    (h : Set.range p.base ⊆ (V.1 : Set X)) :
    coordFun U V d b p = (coordHom V p h).hom b := by
  have hcp := chartPoint_arcHomOf X U d p hpU
  have hcond : Set.range (chartPoint U (CommRingCat.ofHom (arcHomOf U d p))).base
      ⊆ (V.1 : Set X) := by
    rw [hcp]; exact h
  rw [coordFun_eq U V d b p hcond, chartHom_eq_coordHom]
  exact congrArg (fun (m : Γ(X, V.1) ⟶ CommRingCat.of ℂ) => m.hom b)
    (coordHom_congr V _ _ hcp hcond h)

/-! ## ★★★超積の座標 -/

open scoped Classical in
/-- ★★★**座標の超積は超積の座標**。 -/
theorem germ_coordFun (U V : X.affineOpens) (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X))
    (hstarV : Set.range (chartPoint U (CommRingCat.ofHom (starHom U d 𝒰))).base ⊆ (V.1 : Set X))
    (b : Γ(X, V.1)) :
    ((coordFun U V d b : (Spec (CommRingCat.of ℂ) ⟶ X) → ℂ) :
        Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ)
      = (chartHom U V (CommRingCat.ofHom (starHom U d 𝒰)) hstarV).hom b := by
  have hopen : IsOpen ((U.2.isoSpec.inv ≫ U.1.ι).base ⁻¹' (V.1 : Set X)) :=
    (V.1).2.preimage (Scheme.Hom.continuous _)
  have hQV : PrimeSpectrum.comap (starHom U d 𝒰)
      (default : ↥(Spec (CommRingCat.of ((𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)).Germ ℂ))))
      ∈ (U.2.isoSpec.inv ≫ U.1.ι).base ⁻¹' (V.1 : Set X) := by
    have hm := hstarV ⟨default, rfl⟩
    rw [arcPt_chartPoint] at hm
    exact hm
  obtain ⟨f, hf1, hf2⟩ := exists_basicOpen_le_of_mem _ hopen _ hQV
  have hsf : starHom U d 𝒰 f ≠ 0 := (mem_basicOpen_comap_field _ f _).1 hf1
  have hT1 : ∀ᶠ p in (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)), arcHomOf U d p f ≠ 0 :=
    germ_ne_zero_eventually 𝒰 (fun p => arcHomOf U d p f) hsf
  set T : Set (Spec (CommRingCat.of ℂ) ⟶ X) := {p | arcHomOf U d p f ≠ 0} with hTdef
  have hTmem : T ∈ 𝒰 := hT1
  have ht : ∀ t : ↥T,
      Set.range (chartPoint U (CommRingCat.ofHom (arcHomOf U d t.1))).base ⊆ (V.1 : Set X) :=
    fun t => range_chartPoint_subset U V.1 _ f (isUnit_iff_ne_zero.2 t.2) hf2
  have hpi : Set.range (chartPoint U (CommRingCat.ofHom
      (RingHom.pi (fun t : ↥T => arcHomOf U d t.1)))).base ⊆ (V.1 : Set X) :=
    range_chartPoint_subset U V.1 _ f (isUnit_pi_of_ne_zero _ (fun t => t.2)) hf2
  have hgermcond : Set.range (chartPoint U (CommRingCat.ofHom
      ((germRestrict 𝒰 T hTmem).comp (RingHom.pi (fun t : ↥T => arcHomOf U d t.1))))).base
      ⊆ (V.1 : Set X) := by
    rw [← starHom_eq_germRestrict U d 𝒰 T hTmem]
    exact hstarV
  rw [chartHom_congr U V _ _
      (congrArg CommRingCat.ofHom (starHom_eq_germRestrict U d 𝒰 T hTmem)) hstarV hgermcond,
    chartHom_germ_apply U V 𝒰 T hTmem _ hpi hgermcond b, germRestrict_apply, Germ.coe_eq]
  filter_upwards [hTmem] with p hp
  rw [dif_pos hp, chartHom_pi_apply U V _ hpi ⟨p, hp⟩ (ht ⟨p, hp⟩),
    coordFun_eq U V d b p (ht ⟨p, hp⟩)]

/-- ★★超積の座標は、持ち上げの座標の像である。 -/
theorem chartHom_star_eq (U V : X.affineOpens) (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X))
    (l : Spec (CommRingCat.of ↥(finGermSub 𝒰)) ⟶ X)
    (hlV : Set.range l.base ⊆ (V.1 : Set X))
    (hl : Spec.map (CommRingCat.ofHom
        (algebraMap ↥(finGermSub 𝒰)
          (Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ))) ≫ l = starPoint U d 𝒰)
    (hstarV : Set.range (chartPoint U (CommRingCat.ofHom (starHom U d 𝒰))).base
      ⊆ (V.1 : Set X)) (b : Γ(X, V.1)) :
    (chartHom U V (CommRingCat.ofHom (starHom U d 𝒰)) hstarV).hom b
      = algebraMap ↥(finGermSub 𝒰) (Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ)
          ((coordHom V l hlV).hom b) := by
  have hrange : Set.range (Spec.map (CommRingCat.ofHom
      (algebraMap ↥(finGermSub 𝒰)
        (Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ))) ≫ l).base ⊆ (V.1 : Set X) := by
    rintro _ ⟨z, rfl⟩
    exact hlV ⟨_, rfl⟩
  have h1 : chartHom U V (CommRingCat.ofHom (starHom U d 𝒰)) hstarV
      = coordHom V (Spec.map (CommRingCat.ofHom
          (algebraMap ↥(finGermSub 𝒰)
            (Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ))) ≫ l) hrange := by
    rw [chartHom_eq_coordHom]
    exact coordHom_congr V _ _ hl.symm hstarV hrange
  rw [h1, coordHom_natural V l _ hlV hrange]
  rfl

/-- ★★極限点の座標は標準部分である。 -/
theorem coordHom_limit_eq (V : X.affineOpens) {α : Type} (𝒰 : Ultrafilter α)
    (l : Spec (CommRingCat.of ↥(finGermSub 𝒰)) ⟶ X)
    (hlV : Set.range l.base ⊆ (V.1 : Set X))
    (hq : Set.range (Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l).base ⊆ (V.1 : Set X))
    (b : Γ(X, V.1)) :
    (coordHom V (Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l) hq).hom b
      = stHom 𝒰 ((coordHom V l hlV).hom b) := by
  rw [coordHom_natural V l _ hlV hq]
  rfl

/-! ## ★★★★★★★★結論 -/

/-- ★★★★★★★★**固有なら `X^arc` はコンパクトである**——Chow の補題を使わない道。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★超フィルターの極限を**付値判定法**で作る。★★これが (C2) の中身である。 -/
theorem compactSpace_arc {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X)) :
    @CompactSpace (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X) := by
  letI := arcTopology X
  rw [← isCompact_univ_iff, isCompact_iff_ultrafilter_le_nhds]
  intro 𝒰 _
  obtain ⟨U, hU⟩ := exists_affineOpen_mem_ultrafilter 𝒰
  obtain ⟨p₀, hp₀⟩ : {p : Spec (CommRingCat.of ℂ) ⟶ X | arcPt p ∈ U.1}.Nonempty :=
    Filter.nonempty_of_mem hU
  obtain ⟨d, -⟩ := exists_lift_of_mem U.1 p₀ hp₀
  obtain ⟨l, hl⟩ := exists_ultraLift hval U d 𝒰
  refine ⟨Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l, Set.mem_univ _, ?_⟩
  obtain ⟨_, ⟨V0, hVaff, rfl⟩, hqV, -⟩ :=
    (X.isBasis_affineOpens).exists_subset_of_mem_open
      (Set.mem_univ (arcPt (Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l))) isOpen_univ
  set V : X.affineOpens := ⟨V0, hVaff⟩ with hVdef
  have hlV : Set.range l.base ⊆ (V.1 : Set X) := range_ultraLift_subset 𝒰 l V.1 hqV
  have hq : Set.range (Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l).base ⊆ (V.1 : Set X) := by
    rintro _ ⟨z, rfl⟩; exact hlV ⟨_, rfl⟩
  have hstarV : Set.range (chartPoint U (CommRingCat.ofHom (starHom U d 𝒰))).base
      ⊆ (V.1 : Set X) := by
    have hcp : chartPoint U (CommRingCat.ofHom (starHom U d 𝒰)) = starPoint U d 𝒰 := rfl
    rw [hcp, ← hl]
    rintro _ ⟨z, rfl⟩
    exact hlV ⟨_, rfl⟩
  have hUev : ∀ᶠ p in (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)), arcPt p ∈ U.1 := hU
  have hTV : ∀ᶠ p in (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)), arcPt p ∈ V.1 :=
    eventually_arcPt_mem U d 𝒰 V.1 hstarV hUev
  have hrange : Set.range (fun s : Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme => s ≫ V.1.ι) ∈ 𝒰 := by
    filter_upwards [hTV] with p hp
    obtain ⟨s, hs⟩ := exists_lift_of_mem V.1 p hp
    exact ⟨s, hs⟩
  set 𝒱 := Ultrafilter.comap 𝒰 (comp_openImmersion_injective V.1.ι) hrange with h𝒱
  have hmapV : Ultrafilter.map (fun s => s ≫ V.1.ι) 𝒱 = 𝒰 := by
    apply Ultrafilter.coe_injective
    rw [Ultrafilter.coe_map, h𝒱, Ultrafilter.coe_comap, Filter.map_comap_of_mem hrange]
  set qV := liftV V.1 (Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l) hq with hqVdef
  have hqfac : qV ≫ V.1.ι = Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l := liftV_fac _ _ _
  rw [← hqfac]
  refine le_nhds_X V 𝒰 𝒱 hmapV qV (le_nhds_chart V 𝒱 qV ?_)
  intro b
  have hqrange : Set.range (qV ≫ V.1.ι).base ⊆ (V.1 : Set X) := range_comp_ι_subset V.1 qV
  have hrhs : evalAffine Γ(X, V.1) (qV ≫ V.2.isoSpec.hom) b
      = stHom 𝒰 ((coordHom V l hlV).hom b) := by
    rw [evalAffine_eq_coordHom V qV hqrange b,
      show coordHom V (qV ≫ V.1.ι) hqrange
        = coordHom V (Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l) hq from
        coordHom_congr V _ _ hqfac hqrange hq]
    exact coordHom_limit_eq V 𝒰 l hlV hq b
  rw [hrhs]
  have hlhs : ∀ s : Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme,
      evalAffine Γ(X, V.1) (s ≫ V.2.isoSpec.hom) b
        = (coordHom V (s ≫ V.1.ι) (range_comp_ι_subset V.1 s)).hom b :=
    fun s => evalAffine_eq_coordHom V s _ b
  simp only [hlhs]
  have hgerm := germ_coordFun U V d 𝒰 hstarV b
  have hstar_eq := chartHom_star_eq U V d 𝒰 l hlV hl hstarV b
  have hgerm' : ((coordFun U V d b : (Spec (CommRingCat.of ℂ) ⟶ X) → ℂ) :
      Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ)
      = ((coordHom V l hlV).hom b : Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ) := by
    rw [hgerm, hstar_eq]
    rfl
  have htend := tendsto_of_germ_eq 𝒰 (coordFun U V d b) ((coordHom V l hlV).hom b) hgerm'
  have hcoe : (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X))
      = Filter.map (fun s => s ≫ V.1.ι)
        (𝒱 : Filter (Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme)) := by
    rw [← hmapV]
    rfl
  have htend2 : Tendsto ((coordFun U V d b) ∘ (fun s => s ≫ V.1.ι))
      (𝒱 : Filter (Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme))
      (nhds (stHom 𝒰 ((coordHom V l hlV).hom b))) := by
    rw [← tendsto_map'_iff, ← hcoe]
    exact htend
  refine htend2.congr' ?_
  have hUV : ∀ᶠ s in (𝒱 : Filter (Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme)),
      arcPt (s ≫ V.1.ι) ∈ U.1 := by
    rw [hcoe, eventually_map] at hUev
    exact hUev
  filter_upwards [hUV] with s hs
  exact coordFun_eq_coordHom U V d (s ≫ V.1.ι) hs b (range_comp_ι_subset V.1 s)

/-! ## ★出典の紐付け(`.src`) -/

def compactSpace_arc.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——固有性から X^arc のコンパクト性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
