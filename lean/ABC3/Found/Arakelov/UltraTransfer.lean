import ABC3.Found.Arakelov.ChartHom

/-!
# Arakelov (C2) 第 84 ブロック —— **★★★★転送**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★超積の点が `V` に入るなら、ほとんどすべての点も `V` に入る

極限点 `q` の近傍を測る chart `V` は、超積を作った chart `U` と違う。
★そこで「`V` の中で測ってよい」ことを先に言う必要がある。

★★**閉集合は零点集合である**(`PrimeSpectrum.isClosed_iff_zeroLocus`)。
`V` の補集合の引き戻しを `Z = 零点集合 s` と書くと:

    ほとんどすべての `p` で `g(p) = 0` (∀ g ∈ s)
      ⟹ 芽 `g(*) = 0` (∀ g ∈ s)
      ⟹ 超積の点も `Z` に入る

★★★対偶を取ると、超積の点が `V` に入るなら**ほとんどすべての点が `V` に入る**。

★★★★これは超フィルターだから効く——「ほとんど至る所」の否定が「ほとんど至る所」になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `comap_germ_mem_of_eventually` | ★★★閉集合への転送(素スペクトル版) |
| `arcPt_chartPoint` / `chartPoint_arcHomOf` | ★chart と点の往復 |
| `eventually_arcPt_mem` | ★★★★**転送(スキーム版)** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Filter

/-! ## ★★★素スペクトル版 -/

/-- ★★★**転送**——ほとんどすべての点が閉集合に入るなら、超積の点も入る。 -/
theorem comap_germ_mem_of_eventually {A : Type} [CommRing A] {α : Type} (𝒰 : Ultrafilter α)
    (ψ : α → A →+* ℂ) (Z : Set (PrimeSpectrum A)) (hZ : IsClosed Z)
    (P : PrimeSpectrum ℂ) (Q : PrimeSpectrum (Germ (𝒰 : Filter α) ℂ))
    (h : ∀ᶠ a in (𝒰 : Filter α), PrimeSpectrum.comap (ψ a) P ∈ Z) :
    PrimeSpectrum.comap ((Germ.coeRingHom (𝒰 : Filter α)).comp (RingHom.pi ψ)) Q ∈ Z := by
  obtain ⟨s, rfl⟩ := (PrimeSpectrum.isClosed_iff_zeroLocus Z).1 hZ
  rw [PrimeSpectrum.mem_zeroLocus]
  intro g hg
  have hP : P.asIdeal = ⊥ := by
    rcases P.asIdeal.eq_bot_or_top with h1 | h1
    · exact h1
    · exact absurd h1 P.isPrime.ne_top
  have hzero : ∀ᶠ a in (𝒰 : Filter α), ψ a g = 0 := by
    filter_upwards [h] with a ha
    rw [PrimeSpectrum.mem_zeroLocus] at ha
    have hmem := ha hg
    simp only [PrimeSpectrum.comap_asIdeal, hP] at hmem
    exact hmem
  have hg0 : ((Germ.coeRingHom (𝒰 : Filter α)).comp (RingHom.pi ψ)) g = 0 := by
    show ((fun a => ψ a g : α → ℂ) : Germ (𝒰 : Filter α) ℂ) = 0
    have hco : ((fun a => ψ a g : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
        = ((0 : α → ℂ) : Germ (𝒰 : Filter α) ℂ) := by
      rw [Germ.coe_eq]
      filter_upwards [hzero] with a ha
      exact ha
    simpa using hco
  show g ∈ (PrimeSpectrum.comap _ Q).asIdeal
  rw [PrimeSpectrum.comap_asIdeal, Ideal.mem_comap, hg0]
  exact Ideal.zero_mem _

/-! ## ★chart と点の往復 -/

theorem arcPt_chartPoint {X : Scheme.{0}} {R : CommRingCat.{0}} (U : X.affineOpens)
    (φ : Γ(X, U.1) ⟶ R) (P : Spec R) :
    (chartPoint U φ).base P
      = (U.2.isoSpec.inv ≫ U.1.ι).base (PrimeSpectrum.comap φ.hom P) := rfl

/-- ★chart 内の点は、その環準同型から復元できる。 -/
theorem chartPoint_arcHomOf (X : Scheme.{0}) (U : X.affineOpens)
    (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (hp : arcPt p ∈ U.1) :
    chartPoint U (CommRingCat.ofHom (arcHomOf U d p)) = p := by
  obtain ⟨q, hq, hhom⟩ := arcHomOf_of_mem U d p hp
  rw [chartPoint, hhom]
  have hid : CommRingCat.ofHom (Spec.preimage (q ≫ U.2.isoSpec.hom)).hom
      = Spec.preimage (q ≫ U.2.isoSpec.hom) := rfl
  rw [hid, Spec.map_preimage, Category.assoc, Iso.hom_inv_id_assoc, hq]

/-! ## ★★★★スキーム版 -/

/-- ★★★★**転送(スキーム版)**——超積の点が `V` に入るなら、ほとんどすべての点も `V` に入る。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
theorem eventually_arcPt_mem {X : Scheme.{0}} (U : X.affineOpens)
    (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X)) (V : X.Opens)
    (hstar : Set.range (chartPoint U (CommRingCat.ofHom (starHom U d 𝒰))).base ⊆ (V : Set X))
    (hU : ∀ᶠ p in (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)), arcPt p ∈ U.1) :
    ∀ᶠ p in (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)), arcPt p ∈ V := by
  by_contra hcon
  have hne : ∀ᶠ p in (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)), arcPt p ∉ V :=
    Ultrafilter.eventually_not.2 hcon
  have hZ : IsClosed ((U.2.isoSpec.inv ≫ U.1.ι).base ⁻¹' ((V : Set X)ᶜ)) :=
    (V.2.isClosed_compl).preimage (Scheme.Hom.continuous _)
  have hae : ∀ᶠ p in (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)),
      PrimeSpectrum.comap (arcHomOf U d p) (default : ↥(Spec (CommRingCat.of ℂ)))
        ∈ ((U.2.isoSpec.inv ≫ U.1.ι).base ⁻¹' ((V : Set X)ᶜ)) := by
    filter_upwards [hU, hne] with p hpU hpV
    show (U.2.isoSpec.inv ≫ U.1.ι).base _ ∈ ((V : Set X)ᶜ)
    have hp := chartPoint_arcHomOf X U d p hpU
    have hbase : arcPt p
        = (U.2.isoSpec.inv ≫ U.1.ι).base
            (PrimeSpectrum.comap (arcHomOf U d p) (default : ↥(Spec (CommRingCat.of ℂ)))) := by
      conv_lhs => rw [← hp]
      exact arcPt_chartPoint U (CommRingCat.ofHom (arcHomOf U d p))
        (default : ↥(Spec (CommRingCat.of ℂ)))
    rw [← hbase]
    exact hpV
  have hstarZ := comap_germ_mem_of_eventually 𝒰 (fun p => arcHomOf U d p) _ hZ
    (default : ↥(Spec (CommRingCat.of ℂ)))
    (default : ↥(Spec (CommRingCat.of ((𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)).Germ ℂ)))) hae
  have himg : (chartPoint U (CommRingCat.ofHom (starHom U d 𝒰))).base default ∈ (V : Set X) :=
    hstar ⟨default, rfl⟩
  have hout : (chartPoint U (CommRingCat.ofHom (starHom U d 𝒰))).base default
      ∈ ((V : Set X)ᶜ) := by
    rw [arcPt_chartPoint]
    exact Set.mem_preimage.1 hstarZ
  exact hout himg

/-! ## ★出典の紐付け(`.src`) -/

def eventually_arcPt_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——超積の点から近傍への転送)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
