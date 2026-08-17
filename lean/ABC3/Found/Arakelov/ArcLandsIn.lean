import ABC3.Found.Arakelov.ArcBasicOpen

/-!
# Arakelov (C1) の核心 —— **点が基本開集合に落ちる ⟺ 切断が消えない**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★これが `topology_openImmersion` の残り 1 向きの核心である

残るのは「**`(· ≫ f)` が開写像**」。その証明の中心は

    {r : Arc (Spec A) ǀ r の像の点が開集合 O に落ちる}  が開

である。★★これを基本開集合に分解するには

    r の像の点 ∈ D(g)  ⟺  g(r) ≠ 0

が要る。★★★**本ファイルがそれである。**

## ★★機構

`r : Spec ℂ ⟶ Spec A` は `Spec.map (evalHom A r)` に等しい(`ArcEval.lean`)。
★したがって像の点は `PrimeSpectrum.comap (evalHom A r)` を
`Spec ℂ` の**唯一の点**(`⊥`、ℂ が体だから)に当てたもの、すなわち **`ker (evalHom A r)`**。

★★`g ∉ ker φ ⟺ φ g ≠ 0` なので、`mem_basicOpen` と合わせて結論が出る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

/-! ## ★★★像の点と評価 -/

/-- ★★**ℂ-点の像の点は `evalHom` の核である**。 -/
theorem base_default_asIdeal (A : CommRingCat.{0})
    (r : Spec (CommRingCat.of ℂ) ⟶ Spec A) :
    (r.base default).asIdeal = RingHom.ker (evalHom A r).hom := by
  have hdef : (default : ↥(Spec (CommRingCat.of ℂ))).asIdeal = ⊥ := by
    have h : (default : ↥(Spec (CommRingCat.of ℂ)))
        = (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of ℂ)) :=
      Subsingleton.elim _ _
    rw [h]
  conv_lhs => rw [← Spec_map_evalHom A r]
  rw [Spec.map_apply, PrimeSpectrum.comap_asIdeal, hdef, ← RingHom.ker_eq_comap_bot]

/-- ★★★**点が基本開集合に落ちる ⟺ 切断が消えない**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「`X^arc` の位相が代数と噛み合っている」ことの核心である。 -/
theorem mem_basicOpen_base_iff (A : CommRingCat.{0}) (g : A)
    (r : Spec (CommRingCat.of ℂ) ⟶ Spec A) :
    r.base default ∈ PrimeSpectrum.basicOpen g ↔ evalAffine A r g ≠ 0 := by
  constructor
  · intro h hg
    refine ((PrimeSpectrum.mem_basicOpen g (r.base default)).1 h) ?_
    rw [base_default_asIdeal]
    exact RingHom.mem_ker.2 hg
  · intro h
    refine (PrimeSpectrum.mem_basicOpen g (r.base default)).2 ?_
    rw [base_default_asIdeal]
    exact fun hk => h (RingHom.mem_ker.1 hk)

/-- ★★**基本開集合に落ちる点の集合は `arcBasicOpen` と一致する**。 -/
theorem preimage_basicOpen_eq_arcBasicOpen (A : CommRingCat.{0}) (g : A) :
    {r : Spec (CommRingCat.of ℂ) ⟶ Spec A | r.base default ∈ PrimeSpectrum.basicOpen g}
      = arcBasicOpen A g := by
  ext r
  exact mem_basicOpen_base_iff A g r

/-- ★★★**基本開集合に落ちる点の集合は開である**(前の 2 つの帰結)。 -/
theorem isOpen_preimage_basicOpen (A : CommRingCat.{0}) (g : A) :
    @IsOpen _ (arcTopologyAffine A)
      {r : Spec (CommRingCat.of ℂ) ⟶ Spec A | r.base default ∈ PrimeSpectrum.basicOpen g} := by
  rw [preimage_basicOpen_eq_arcBasicOpen]
  exact isOpen_arcBasicOpen A g

/-! ## ★出典の紐付け(`.src`) -/

def mem_basicOpen_base_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——点が基本開集合に落ちることと切断が消えないことの同値)",
    sectionId := "genell-def-1-1-i" }

def isOpen_preimage_basicOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——基本開集合に落ちる点の集合が開であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
