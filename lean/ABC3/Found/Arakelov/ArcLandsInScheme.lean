import ABC3.Found.Arakelov.ArcLandsIn
import ABC3.Found.Arakelov.ArcTopology
import ABC3.Found.Arakelov.ArcLift

/-!
# Arakelov (C1) の配管 —— **一般スキームでの「点が開集合に落ちる」**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★これが「貼る段」の中核である

`ArcLandsIn.lean` でアフィンの場合(`continuous_base_default`)を取った。
★本ファイルは**一般の `X`** へ上げる。

★★機構は `arcTopology X = ⨆`(アフィン chart)なので、
chart ごとに落とし、各 chart では `isoSpec` で輸送すればよい。
★★★**「像の点を取る」写像が `X^arc` から `X`(Zariski)へ連続**という形に述べると、
`isOpen_landsIn_scheme` はその系として出る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory TopologicalSpace

/-! ## ★★アフィン開の上での連続性(輸送) -/

/-- ★**アフィン開 `U` の上でも「像の点を取る」写像は連続**。

★`arcTopologyOpen U` は `isoSpec` で輸送した位相なので、
アフィンの場合(`continuous_base_default`)から出る。 -/
theorem continuous_base_default_open {X : Scheme.{0}} (U : X.affineOpens) :
    @Continuous _ _ (arcTopologyOpen U) _
      (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p.base default) := by
  letI := arcTopologyOpen U
  letI := arcTopologyAffine (X.presheaf.obj (Opposite.op U.1))
  have h : (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p.base default)
      = (fun x => U.2.isoSpec.inv.base x)
        ∘ (fun r : Spec (CommRingCat.of ℂ) ⟶ Spec (X.presheaf.obj (Opposite.op U.1)) =>
            r.base default)
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.2.isoSpec.hom) := by
    funext p
    show p.base default
      = U.2.isoSpec.inv.base ((p ≫ U.2.isoSpec.hom).base default)
    rw [Scheme.Hom.comp_apply, ← Scheme.Hom.comp_apply, Iso.hom_inv_id]
    rfl
  rw [h]
  exact (U.2.isoSpec.inv.base.hom.continuous_toFun).comp
    ((continuous_base_default _).comp continuous_induced_dom)

/-! ## ★★★一般スキーム版 -/

/-- ★★★**一般の `X` で「像の点を取る」写像は `X^arc` から Zariski 位相へ連続**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★機構は `arcTopology = ⨆` を chart ごとに落とすこと。
★各 chart では `continuous_base_default_open` が効く。 -/
theorem continuous_base_default_scheme {X : Scheme.{0}} :
    @Continuous _ _ (arcTopology X) _
      (fun q : Spec (CommRingCat.of ℂ) ⟶ X => q.base default) := by
  refine continuous_iSup_dom.2 fun U => ?_
  letI := arcTopologyOpen U
  refine continuous_coinduced_dom.2 ?_
  have h : (fun q : Spec (CommRingCat.of ℂ) ⟶ X => q.base default)
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι)
      = (fun x => U.1.ι.base x)
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p.base default) := by
    funext p
    exact Scheme.Hom.comp_apply p U.1.ι default
  rw [h]
  exact (U.1.ι.base.hom.continuous_toFun).comp (continuous_base_default_open U)

/-- ★★★**一般の `X` でも「点が開集合に落ちる」条件は `X^arc` の開集合である**。

★★★これが「貼る段」で要る形そのものである。 -/
theorem isOpen_landsIn_scheme {X : Scheme.{0}} (O : X.Opens) :
    @IsOpen _ (arcTopology X)
      {q : Spec (CommRingCat.of ℂ) ⟶ X | q.base default ∈ O} := by
  letI := arcTopology X
  exact continuous_base_default_scheme.isOpen_preimage _ O.isOpen

/-! ## ★★★アフィン開の上での「落ちる」条件(手順 4 で使う形) -/

/-- ★★★**アフィン開 `U` の上でも「点が開集合に落ちる」条件は開である**。

★`continuous_base_default_open` の直接の系。
★★★これが `topology_openImmersion` の**手順 4** で使う形である
——`U ⊓ U'` に落ちる点の集合が開であることを与える。 -/
theorem isOpen_landsIn_open {X : Scheme.{0}} (U : X.affineOpens)
    (O : U.1.toScheme.Opens) :
    @IsOpen _ (arcTopologyOpen U)
      {p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme | p.base default ∈ O} := by
  letI := arcTopologyOpen U
  exact (continuous_base_default_open U).isOpen_preimage _ O.isOpen

/-- ★★**2 つのアフィン開の共通部分に落ちる点の集合は開である**。

★★★手順 4 では `U ⊓ U'` を扱う——これがその形である。 -/
theorem isOpen_landsIn_inter {X : Scheme.{0}} (U U' : X.affineOpens) :
    @IsOpen _ (arcTopologyOpen U')
      {p : Spec (CommRingCat.of ℂ) ⟶ U'.1.toScheme |
        p.base default ∈ (U'.1.ι ⁻¹ᵁ U.1)} :=
  isOpen_landsIn_open U' _

/-! ## ★★★chart の像は開である -/

/-- ★★★**アフィン chart の像は `X^arc` の開集合である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

    Set.range (· ≫ U.ι) = {q ǀ q の像の点が U に入る}

★★機構は `mem_range_comp_iff_base_default`(像への所属 = 点の条件)と
`isOpen_landsIn_scheme`(点が開集合に落ちる条件は開)。

★★★**これが「chart が開埋め込みである」の半分**である
(残り半分は位相が部分空間位相であること)。 -/
theorem isOpen_range_chart {X : Scheme.{0}} (U : X.affineOpens) :
    @IsOpen _ (arcTopology X)
      (Set.range (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι)) := by
  letI := arcTopology X
  have hset : Set.range (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι)
      = {q : Spec (CommRingCat.of ℂ) ⟶ X | q.base default ∈ U.1} := by
    ext q
    rw [mem_range_comp_iff_base_default U.1.ι q, Set.mem_setOf_eq]
    constructor
    · intro h
      have := Scheme.Opens.range_ι U.1
      rw [this] at h
      exact h
    · intro h
      have := Scheme.Opens.range_ι U.1
      rw [this]
      exact h
  rw [hset]
  exact isOpen_landsIn_scheme U.1

/-! ## ★出典の紐付け(`.src`) -/

def continuous_base_default_scheme.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——像の点を取る写像の連続性、一般のスキーム)",
    sectionId := "genell-def-1-1-i" }

def isOpen_range_chart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィン chart の像が開であること)",
    sectionId := "genell-def-1-1-i" }

def isOpen_landsIn_open.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィン開の上で落ちる条件が開)",
    sectionId := "genell-def-1-1-i" }

def isOpen_landsIn_inter.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——2 つのアフィン開の共通部分に落ちる条件が開)",
    sectionId := "genell-def-1-1-i" }

def isOpen_landsIn_scheme.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——点が開集合に落ちる条件が開、一般のスキーム)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
