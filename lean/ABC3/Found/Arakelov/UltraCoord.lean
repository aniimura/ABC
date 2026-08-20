import ABC3.Found.Arakelov.UltraBasic

/-!
# Arakelov (C2) 第 87–88 ブロック —— **★★★超積の座標は座標の超積**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★これで収束が言える

`chart` 間の翻訳が**係数環に自然**(第 83)だったので、`ρ` に

| `ρ` | 得られるもの |
|---|---|
| 成分への射影 `(↥T → ℂ) → ℂ` | ★各点の座標 |
| 芽 `(↥T → ℂ) → *ℂ` | ★★超積の座標 |
| `O ↪ *ℂ` | ★★★持ち上げの座標から超積の座標 |
| `st : O → ℂ` | ★★★★持ち上げの座標から極限の座標 |

を入れると、**「各点の座標の超積 = 持ち上げの座標」**が出る。
★★★★★★したがって各点の座標は**標準部分へ収束する**——これが求める収束である。

## ★★積の像を `V` に収める

積 `(↥T → ℂ)` の `Spec` は大きい(超フィルターの分だけ点がある)。
★そのままでは像が `V` に入らない——**基本開集合 `D(f)` で挟む**必要がある。
★★`f` の値が全ての `t ∈ T` で 0 でなければ `φ f` は**単元**であり、
像は `D(f)` に収まる(第 85)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `germRestrict` | ★`T` 上の族から芽を取る環準同型 |
| `chartHom_pi_apply` | ★★積の座標は成分の座標 |
| `chartHom_germ_apply` | ★★★**超積の座標は座標の超積** |
| `starHom_eq_germRestrict` | ★`starHom` は `T` 上の族の芽 |
| `range_chartPoint_subset` | ★★単元条件から像が `V` に入る |
| `coordHom` / `coordHom_natural` | ★★★一般の点の座標とその自然性 |
| `tendsto_of_germ_eq` | ★★★**芽が有限元なら標準部分へ収束** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Filter

/-! ## ★制限からの芽 -/

open scoped Classical in
/-- ★★`T ∈ 𝒰` のとき、`T` 上の族から芽を取る環準同型。 -/
noncomputable def germRestrict {α : Type} (𝒰 : Ultrafilter α) (T : Set α) (hT : T ∈ 𝒰) :
    (↥T → ℂ) →+* Germ (𝒰 : Filter α) ℂ where
  toFun := fun g =>
    (((fun a => if h : a ∈ T then g ⟨a, h⟩ else 0) : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
  map_one' := by
    show ((_ : α → ℂ) : Germ (𝒰 : Filter α) ℂ) = ((1 : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
    rw [Germ.coe_eq]
    filter_upwards [hT] with a ha
    simp [ha]
  map_mul' := by
    intro g h
    show ((_ : α → ℂ) : Germ (𝒰 : Filter α) ℂ) = _
    rw [← Germ.coe_mul, Germ.coe_eq]
    filter_upwards [hT] with a ha
    simp [ha]
  map_zero' := by
    show ((_ : α → ℂ) : Germ (𝒰 : Filter α) ℂ) = ((0 : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
    rw [Germ.coe_eq]
    filter_upwards [hT] with a ha
    simp [ha]
  map_add' := by
    intro g h
    show ((_ : α → ℂ) : Germ (𝒰 : Filter α) ℂ) = _
    rw [← Germ.coe_add, Germ.coe_eq]
    filter_upwards [hT] with a ha
    simp [ha]

/-! ## ★★自然性の適用 -/

variable {X : Scheme.{0}}

/-- ★★積の座標は各成分の座標である。 -/
theorem chartHom_pi_apply (U V : X.affineOpens) {ι : Type} (ψ : ι → (Γ(X, U.1) →+* ℂ))
    (hpi : Set.range (chartPoint U (CommRingCat.ofHom (RingHom.pi ψ))).base ⊆ (V.1 : Set X))
    (t : ι)
    (ht : Set.range (chartPoint U (CommRingCat.ofHom (ψ t))).base ⊆ (V.1 : Set X))
    (x : Γ(X, V.1)) :
    (chartHom U V (CommRingCat.ofHom (RingHom.pi ψ)) hpi).hom x t
      = (chartHom U V (CommRingCat.ofHom (ψ t)) ht).hom x := by
  have hnat := chartHom_natural U V (CommRingCat.ofHom (RingHom.pi ψ))
    (CommRingCat.ofHom (Pi.evalRingHom (fun _ : ι => ℂ) t)) hpi ht
  have hap := congrArg (fun (m : Γ(X, V.1) ⟶ CommRingCat.of ℂ) => m.hom x) hnat
  exact hap.symm

/-- ★★★**超積の座標は座標の超積**——これが要である。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
theorem chartHom_germ_apply (U V : X.affineOpens) {α : Type} (𝒰 : Ultrafilter α)
    (T : Set α) (hT : T ∈ 𝒰) (ψ : ↥T → (Γ(X, U.1) →+* ℂ))
    (hpi : Set.range (chartPoint U (CommRingCat.ofHom (RingHom.pi ψ))).base ⊆ (V.1 : Set X))
    (hgerm : Set.range (chartPoint U
      (CommRingCat.ofHom ((germRestrict 𝒰 T hT).comp (RingHom.pi ψ)))).base ⊆ (V.1 : Set X))
    (x : Γ(X, V.1)) :
    (chartHom U V (CommRingCat.ofHom ((germRestrict 𝒰 T hT).comp (RingHom.pi ψ))) hgerm).hom x
      = germRestrict 𝒰 T hT ((chartHom U V (CommRingCat.ofHom (RingHom.pi ψ)) hpi).hom x) := by
  have hnat := chartHom_natural U V (CommRingCat.ofHom (RingHom.pi ψ))
    (CommRingCat.ofHom (germRestrict 𝒰 T hT)) hpi hgerm
  exact congrArg (fun (m : Γ(X, V.1) ⟶ CommRingCat.of (Germ (𝒰 : Filter α) ℂ)) => m.hom x) hnat

open scoped Classical in
/-- ★★`starHom` は `T` 上の族の芽として書ける(`T ∈ 𝒰` なら)。 -/
theorem starHom_eq_germRestrict (U : X.affineOpens)
    (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X))
    (T : Set (Spec (CommRingCat.of ℂ) ⟶ X)) (hT : T ∈ 𝒰) :
    starHom U d 𝒰
      = (germRestrict 𝒰 T hT).comp (RingHom.pi (fun t : ↥T => arcHomOf U d t.1)) := by
  ext x
  show (((fun p => arcHomOf U d p x) : (Spec (CommRingCat.of ℂ) ⟶ X) → ℂ) :
      Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ) = _
  show _ = (((fun p => if h : p ∈ T then arcHomOf U d p x else 0) :
      (Spec (CommRingCat.of ℂ) ⟶ X) → ℂ) : Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ)
  rw [Germ.coe_eq]
  filter_upwards [hT] with p hp
  simp [hp]

/-- ★★`φ f` が単元で `D(f)` が `V` に入るなら、点は `V` に入る。 -/
theorem range_chartPoint_subset {R : CommRingCat.{0}} (U : X.affineOpens) (V : X.Opens)
    (φ : Γ(X, U.1) ⟶ R) (f : Γ(X, U.1)) (hf : IsUnit (φ.hom f))
    (hsub : (PrimeSpectrum.basicOpen f : Set (PrimeSpectrum Γ(X, U.1)))
      ⊆ (U.2.isoSpec.inv ≫ U.1.ι).base ⁻¹' (V : Set X)) :
    Set.range (chartPoint U φ).base ⊆ (V : Set X) := by
  rintro _ ⟨P, rfl⟩
  rw [arcPt_chartPoint]
  exact hsub (range_comap_subset_basicOpen φ.hom f hf ⟨P, rfl⟩)

/-! ## ★★★一般の点の座標 -/

/-- ★★`V` に入る点の `V` 座標。 -/
noncomputable def coordHom {R : CommRingCat.{0}} (V : X.affineOpens) (g : Spec R ⟶ X)
    (h : Set.range g.base ⊆ (V.1 : Set X)) : Γ(X, V.1) ⟶ R :=
  Spec.preimage (liftV V.1 g h ≫ V.2.isoSpec.hom)

theorem chartHom_eq_coordHom {R : CommRingCat.{0}} (U V : X.affineOpens) (φ : Γ(X, U.1) ⟶ R)
    (h : Set.range (chartPoint U φ).base ⊆ (V.1 : Set X)) :
    chartHom U V φ h = coordHom V (chartPoint U φ) h := rfl

/-- ★★★**座標は係数環に自然である**。 -/
theorem coordHom_natural {R R' : CommRingCat.{0}} (V : X.affineOpens) (g : Spec R ⟶ X)
    (ρ : R ⟶ R') (h : Set.range g.base ⊆ (V.1 : Set X))
    (h' : Set.range (Spec.map ρ ≫ g).base ⊆ (V.1 : Set X)) :
    coordHom V (Spec.map ρ ≫ g) h' = coordHom V g h ≫ ρ := by
  unfold coordHom
  have hlift : liftV V.1 (Spec.map ρ ≫ g) h' = Spec.map ρ ≫ liftV V.1 g h := by
    apply (cancel_mono V.1.ι).1
    rw [liftV_fac, Category.assoc, liftV_fac]
  rw [hlift, Category.assoc, Spec.preimage_comp, Spec.preimage_map]

/-! ## ★★★収束 -/

/-- ★★★**芽が有限元と一致すれば、その標準部分へ収束する**。 -/
theorem tendsto_of_germ_eq {α : Type} (𝒰 : Ultrafilter α) (F : α → ℂ) (c : ↥(finGermSub 𝒰))
    (hc : ((F : α → ℂ) : Germ (𝒰 : Filter α) ℂ) = (c : Germ (𝒰 : Filter α) ℂ)) :
    Tendsto F (𝒰 : Filter α) (nhds (stHom 𝒰 c)) := by
  obtain ⟨C, hC⟩ := c.2
  have hb : ∃ C : ℝ, ∀ᶠ a in (𝒰 : Filter α), ‖F a‖ ≤ C := by
    refine ⟨C, ?_⟩
    rw [← hc] at hC
    rw [Germ.liftPred_coe] at hC
    exact hC
  have htend := tendsto_stOf 𝒰 hb
  have heq : stOf 𝒰 F = stHom 𝒰 c := by
    show stOf 𝒰 F = stG 𝒰 (c : Germ (𝒰 : Filter α) ℂ)
    rw [← hc]
    rfl
  rwa [heq] at htend

/-! ## ★出典の紐付け(`.src`) -/

def chartHom_germ_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——超積の座標が座標の超積であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
