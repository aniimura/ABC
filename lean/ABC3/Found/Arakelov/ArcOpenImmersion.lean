import ABC3.Found.Arakelov.ArcTopologyAffineEq

/-!
# Arakelov (C1) の第七段 —— **開埋め込みに沿った連続性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★C1 に残る最後の条件

`Interface/Arakelov/ArcSpace.lean` の C1 は

    topology X = induced (map f) (topology Y)      (f : X ⟶ Y が開埋め込み)

を要求する。★これは 2 向きある:

| 向き | 意味 | 状態 |
|---|---|---|
| `topology X ≤ induced` | ★**`(· ≫ f)` が連続** | ★★★**本ファイル** |
| `induced ≤ topology X` | 開埋め込みの像で部分空間位相になること | ★未取得 |

## ★★★機構: 開埋め込みはアフィン開をアフィン開に送る

`f` が開埋め込みで `U ⊆ X` がアフィン開なら `f ''ᵁ U ⊆ Y` もアフィン開
(`IsAffineOpen.image_of_isOpenImmersion`、mathlib)。★★さらに

    U.ι ≫ f = e.hom ≫ (f ''ᵁ U).ι

という分解がある(`IsOpenImmersion.isoOfRangeEq_hom_fac`)。
★★★**したがって `X` の chart は `Y` の chart 1 つに収まる**——ここが要点である。
一般の射だと像が複数の chart にまたがるので、この議論は開埋め込みに固有である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory TopologicalSpace

/-! ## ★★アフィン開の像への同型 -/

/-- ★**開埋め込みに沿ったアフィン開の像**。 -/
def imageAffineOpen {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (U : X.affineOpens) : Y.affineOpens :=
  ⟨f ''ᵁ U.1, U.2.image_of_isOpenImmersion f⟩

/-- ★★**`U` とその像は同型**であり、包含と両立する。 -/
noncomputable def imageAffineOpenIso {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (U : X.affineOpens) : U.1.toScheme ≅ (imageAffineOpen f U).1.toScheme :=
  IsOpenImmersion.isoOfRangeEq (U.1.ι ≫ f) (imageAffineOpen f U).1.ι (by
    show Set.range ((U.1.ι ≫ f).base) = _
    rw [Scheme.Opens.range_ι]
    have : ((U.1.ι ≫ f).base : _ → _) = (f.base : _ → _) ∘ (U.1.ι.base : _ → _) := rfl
    rw [this, Set.range_comp, Scheme.Opens.range_ι]
    exact (Scheme.Hom.coe_image (f := f) (U := U.1)).symm)

/-- ★★**分解**: `U.ι ≫ f = iso.hom ≫ (f ''ᵁ U).ι`。 -/
theorem imageAffineOpenIso_fac {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (U : X.affineOpens) :
    (imageAffineOpenIso f U).hom ≫ (imageAffineOpen f U).1.ι = U.1.ι ≫ f :=
  IsOpenImmersion.isoOfRangeEq_hom_fac _ _ _

/-! ## ★★★開埋め込みに沿った合成の連続性 -/

/-- ★★**chart 同士の同型に沿った合成は連続**。

★`iso.hom ≫ isoSpecV.hom = isoSpecU.hom ≫ (isoSpecU.inv ≫ iso.hom ≫ isoSpecV.hom)` と
分解する。★★括弧の中は `Spec Γ(X,U) ⟶ Spec Γ(Y,V)` という**アフィン射**なので
`continuous_comp_affine` が効く。 -/
theorem continuous_comp_imageIso {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (U : X.affineOpens) :
    @Continuous _ _ (arcTopologyOpen U) (arcTopologyOpen (imageAffineOpen f U))
      (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ (imageAffineOpenIso f U).hom) := by
  letI := arcTopologyOpen U
  letI := arcTopologyOpen (imageAffineOpen f U)
  letI := arcTopologyAffine (X.presheaf.obj (Opposite.op U.1))
  letI := arcTopologyAffine (Y.presheaf.obj (Opposite.op (imageAffineOpen f U).1))
  refine continuous_induced_rng.2 ?_
  have h : (fun q : Spec (CommRingCat.of ℂ) ⟶ (imageAffineOpen f U).1.toScheme =>
        q ≫ (imageAffineOpen f U).2.isoSpec.hom)
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme =>
            p ≫ (imageAffineOpenIso f U).hom)
      = (fun q : Spec (CommRingCat.of ℂ) ⟶ Spec (X.presheaf.obj (Opposite.op U.1)) =>
          q ≫ (U.2.isoSpec.inv ≫ (imageAffineOpenIso f U).hom ≫
            (imageAffineOpen f U).2.isoSpec.hom))
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.2.isoSpec.hom) := by
    funext p
    simp only [Function.comp_apply, Category.assoc, Iso.hom_inv_id_assoc]
  rw [h]
  exact Continuous.comp (continuous_comp_affine _) continuous_induced_dom

/-- ★★★**開埋め込みとの合成は連続である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `ArcSpaceData.topology_openImmersion` の**前半**
(`topology X ≤ induced (map f) (topology Y)`)である。

★機構は `imageAffineOpenIso_fac`——**`X` の chart が `Y` の chart 1 つに収まる**。 -/
theorem continuous_comp_openImmersion {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f] :
    @Continuous _ _ (arcTopology X) (arcTopology Y)
      (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ f) := by
  refine continuous_iSup_dom.2 fun U => ?_
  refine continuous_iSup_rng (i := imageAffineOpen f U) ?_
  letI := arcTopologyOpen U
  letI := arcTopologyOpen (imageAffineOpen f U)
  letI : TopologicalSpace (Spec (CommRingCat.of ℂ) ⟶ Y) :=
    TopologicalSpace.coinduced
      (fun q : Spec (CommRingCat.of ℂ) ⟶ (imageAffineOpen f U).1.toScheme =>
        q ≫ (imageAffineOpen f U).1.ι)
      (arcTopologyOpen (imageAffineOpen f U))
  refine continuous_coinduced_dom.2 ?_
  have h : (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ f) ∘
        (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι)
      = (fun q : Spec (CommRingCat.of ℂ) ⟶ (imageAffineOpen f U).1.toScheme =>
          q ≫ (imageAffineOpen f U).1.ι)
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme =>
            p ≫ (imageAffineOpenIso f U).hom) := by
    funext p
    simp only [Function.comp_apply, Category.assoc, imageAffineOpenIso_fac]
  rw [h]
  exact Continuous.comp continuous_coinduced_rng (continuous_comp_imageIso f U)

/-! ## ★★単射性 -/

/-- ★★**開埋め込みとの合成は単射である**(開埋め込みはモノだから)。

★★★`topology_openImmersion` の残り 1 向き
(`induced (· ≫ f) (arcTopology Y) ≤ arcTopology X`)には、
「像が `Arc Y` で開」と合わせてこれが要る——★開集合 `V` に対し
`W := (· ≫ f) '' V` を取り `(· ≫ f) ⁻¹' W = V` を得るのに単射性を使う。 -/
theorem comp_openImmersion_injective {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f] :
    Function.Injective (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ f) := by
  intro p q h
  exact (cancel_mono f).1 h

/-- ★**像の逆像はもとに戻る**(単射性の帰結)。 -/
theorem preimage_image_comp_openImmersion {X Y : Scheme.{0}} (f : X ⟶ Y)
    [IsOpenImmersion f] (V : Set (Spec (CommRingCat.of ℂ) ⟶ X)) :
    (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ f) ⁻¹'
        ((fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ f) '' V) = V :=
  (comp_openImmersion_injective f).preimage_image V

/-! ## ★★★アフィン開の対応(mathlib との接続) -/

/-- ★★★**`imageAffineOpen` は mathlib の `affineOpensEquiv` そのものである**。

`IsOpenImmersion.affineOpensEquiv f : X.affineOpens ≃o {U : Y.affineOpens // U ≤ f.opensRange}`
——★**順序同型**である(mathlib、2026-08-17 実測)。

★★★これが残る 1 向きの鍵である: `X` のアフィン開被覆と
`Y` のアフィン開被覆のうち `f` の像に入るものが**1 対 1 に対応する**。 -/
theorem imageAffineOpen_eq_affineOpensEquiv {X Y : Scheme.{0}} (f : X ⟶ Y)
    [IsOpenImmersion f] (U : X.affineOpens) :
    imageAffineOpen f U = ((IsOpenImmersion.affineOpensEquiv f) U).1 := rfl

/-- ★★**逆向きの対応**——`f` の像に入る `Y` のアフィン開は `X` のアフィン開から来る。 -/
theorem exists_imageAffineOpen {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (V : Y.affineOpens) (hV : V.1 ≤ f.opensRange) :
    ∃ U : X.affineOpens, imageAffineOpen f U = V := by
  refine ⟨(IsOpenImmersion.affineOpensEquiv f).symm ⟨V, hV⟩, ?_⟩
  rw [imageAffineOpen_eq_affineOpensEquiv]
  exact congrArg Subtype.val ((IsOpenImmersion.affineOpensEquiv f).apply_symm_apply ⟨V, hV⟩)

/-! ## ★★★同型の場合の `topology_openImmersion` -/

/-- ★★★**同型に沿っては位相が誘導される**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これは `topology_openImmersion` の**同型の場合**である。
★手順 4b で使う 3 本の同型(`isoSpec` / `Γ(U',D(g)) ≅ Away g` /
`imageAffineOpenIso`)のうち、輸送の骨がこれである。

★機構は簡単: 逆向きの合成 `(· ≫ e.inv)` が連続なので、
像 `(· ≫ e.hom) '' V` は `(· ≫ e.inv) ⁻¹' V` に一致して開である。 -/
theorem arcTopology_iso {X Y : Scheme.{0}} (e : X ≅ Y) :
    arcTopology X
      = TopologicalSpace.induced
          (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ e.hom) (arcTopology Y) := by
  letI := arcTopology X
  letI := arcTopology Y
  refine le_antisymm (continuous_iff_le_induced.1 (continuous_comp_openImmersion e.hom)) ?_
  intro V hV
  refine isOpen_induced_iff.2
    ⟨(fun q : Spec (CommRingCat.of ℂ) ⟶ Y => q ≫ e.inv) ⁻¹' V, ?_, ?_⟩
  · exact (continuous_comp_openImmersion e.inv).isOpen_preimage _ hV
  · ext p
    simp only [Set.mem_preimage, Category.assoc, Iso.hom_inv_id, Category.comp_id]

/-! ## ★★★★chart の位相は内在的である -/

/-- ★★★★**アフィン開の上の位相は、その開部分スキーム自身の `arcTopology` に等しい**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

    arcTopologyOpen U = arcTopology U.1.toScheme

★★★**これは重要な同一視である**——`arcTopologyOpen` は `isoSpec` で
輸送して定義したものだが、それが**開部分スキーム自身の位相**に一致する。
★したがって

    arcTopology X = ⨆ U, coinduced (chart U) (arcTopology U.toScheme)

と書き直せる。★★これが `topology_openImmersion` の一般の場合への還元の骨である。

★機構は 2 つ: `arcTopology_spec`(アフィンでは各点収束)と
`arcTopology_iso`(同型に沿って誘導される)。 -/
theorem arcTopologyOpen_eq_arcTopology {X : Scheme.{0}} (U : X.affineOpens) :
    arcTopologyOpen U = arcTopology U.1.toScheme := by
  rw [arcTopology_iso U.2.isoSpec, arcTopology_spec]

/-! ## ★出典の紐付け(`.src`) -/

def imageAffineOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みがアフィン開を保つこと)",
    sectionId := "genell-def-1-1-i" }

def imageAffineOpenIso_fac.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みの chart の分解)",
    sectionId := "genell-def-1-1-i" }

def continuous_comp_openImmersion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みとの合成の連続性)",
    sectionId := "genell-def-1-1-i" }

def comp_openImmersion_injective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みとの合成が単射であること)",
    sectionId := "genell-def-1-1-i" }

def arcTopologyOpen_eq_arcTopology.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——chart の位相が内在的であること)",
    sectionId := "genell-def-1-1-i" }

def arcTopology_iso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——同型に沿って位相が誘導されること)",
    sectionId := "genell-def-1-1-i" }

def exists_imageAffineOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みでのアフィン開の対応)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
