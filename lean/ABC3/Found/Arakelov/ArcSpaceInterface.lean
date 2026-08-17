import ABC3.Found.Arakelov.ArcOpenMap
import ABC3.Interface.Arakelov.ArcSpace

/-!
# Arakelov (C1) 段 C —— **`ArcSpaceData` の非空虚 witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`waiting` を外す

`Interface/Arakelov/ArcSpace.lean` の `ArcSpaceData` は 7 要求を課していた。
★★**7 件すべてを `Found/Arakelov/` で構成した**ので、witness を与える。

| 要求 | 実装 |
|---|---|
| `Arc` | `complexPoints X`(= `Spec ℂ ⟶ X`) |
| `topology` | `arcTopology`(アフィン chart の `⨆` coinduced) |
| `conj` | `conjPoint`(`conjSpec ≫ ·`) |
| `conj_involutive` | `conjPoint_conjPoint` |
| `conj_continuous` | `continuous_conjPoint` |
| `equivComplexPoints` | `Equiv.refl`(★台がまさに複素点だから) |
| `evalAffine` | `evalAffine`(`Spec.preimage` の環準同型を切断に当てる) |
| `evalAffine_spec` | ★`rfl`(定義がまさにそれだから) |
| `topology_affine` | `arcTopology_spec`(★`⨆` がアフィンで各点収束位相に潰れる) |
| `map` | `· ≫ f` |
| `topology_openImmersion` | ★★★★★`arcTopology_openImmersion`(段 A + 段 B) |

## ★★退化していないことの根拠

★★★`equivComplexPoints` が `Arc` を**複素点の集合に固定**し、
`evalAffine_spec` が評価を**`Spec` の充満忠実性が与える環準同型に固定**し、
`topology_affine` がアフィンでの位相を**各点収束位相に固定**し、
`topology_openImmersion` が貼り合わせでの位相を**固定する**。
★したがって離散位相・密着位相・1 点集合はいずれもこの構造を満たさない
(負の対照は `Check/Arakelov/ArcSpaceNondegenerate.lean`)。
-/

universe u v

namespace ABC3.Interface.Arakelov

open AlgebraicGeometry CategoryTheory ABC3.Found.Arakelov ABC3.Found.GenEll

/-- ★★★★★**`X^arc` の実装**——複素点の集合に、アフィン chart から貼った位相を載せる。 -/
noncomputable def arcSpaceDataImpl : ArcSpaceData where
  Arc := fun X => complexPoints X
  topology := fun X => arcTopology X
  conj := fun _ p => conjPoint p
  conj_involutive := fun _ p => conjPoint_conjPoint p
  conj_continuous := fun _ => continuous_conjPoint
  equivComplexPoints := fun _ => Equiv.refl _
  evalAffine := fun A p a => evalAffine A p a
  evalAffine_spec := fun _ _ _ => rfl
  topology_affine := fun A => arcTopology_spec A
  map := fun f p => p ≫ f
  topology_openImmersion := fun f hf => @arcTopology_openImmersion _ _ f hf

/-- ★★★★★**`ArcSpaceData` は非空虚である**——C1 達成。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが Arakelov 理論の 9 件のうち **C1** である。 -/
theorem ArcSpaceData.nonvacuous : Nonempty ArcSpaceData :=
  ⟨arcSpaceDataImpl⟩

/-! ## ★出典の紐付け(`.src`) -/

def arcSpaceDataImpl.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——X^arc の構成そのもの)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.Arakelov
