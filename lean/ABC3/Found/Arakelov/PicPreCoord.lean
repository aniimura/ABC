import ABC3.Found.Arakelov.PicHcompatPre

/-!
# Arakelov (B2) 第 230 ブロック —— ★★★**`fromSpec⁻¹` 座標での転送**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-262 で「第 211 の座標変更が要る」と測った所

第 211(`comap_ideal_image`)は `appIso ⊤` の座標で書いてある。
★`hcompat`(第 229)は `fromSpec⁻¹` 座標なので、転送も同じ座標で要る。

## ★★mathlib の補題は**任意のアフィン開**で使える

    ideal_comap_of_isOpenImmersion (I) (f) (U) :
      (I.comap f).ideal U = (I.ideal ⟨f ''ᵁ U, _⟩).comap ((f.appIso U).inv.hom)

★`U` は任意である——第 211 では `U := ⊤` と取ったが、
★★**`U := fromSpec ⁻¹ᵁ A` と取れば `fromSpec⁻¹` 座標になる**。

★★★`fromSpec ⁻¹ᵁ A` がアフィン開であることは
`fromSpec_preimage_self`(`= ⊤`)で転送すれば出る。

## ★★★これで「座標を選べる」ことが確認できた

★§9-262 で「座標に依存しない形で書けるものは書く」と結論したが、
**mathlib の補題自体が座標をパラメータに取っている**ので、
★★**呼び出し側で座標を選ぶだけ**で済んだ——第 211 を書き直す必要はない。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isAffineOpen_pre` | ★`fromSpec ⁻¹ᵁ A` はアフィン開 |
| `comap_fromSpec_pre` | ★★★**`fromSpec⁻¹` 座標での転送** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

/-- ★**`fromSpec ⁻¹ᵁ A` はアフィン開である**。 -/
theorem isAffineOpen_pre {Z : Scheme.{u}} {A : Z.Opens} (hA : IsAffineOpen A) :
    IsAffineOpen (hA.fromSpec ⁻¹ᵁ A) := by
  rw [hA.fromSpec_preimage_self]
  exact isAffineOpen_top _

/-- ★★★**`fromSpec⁻¹` 座標での転送**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★mathlib の `ideal_comap_of_isOpenImmersion` は座標をパラメータに取るので、
**呼び出し側で選ぶだけ**である。 -/
theorem comap_fromSpec_pre {Z : Scheme.{u}} (E : Z.IdealSheafData) {A : Z.Opens}
    (hA : IsAffineOpen A) :
    (E.comap hA.fromSpec).ideal ⟨hA.fromSpec ⁻¹ᵁ A, isAffineOpen_pre hA⟩
      = (E.ideal ⟨hA.fromSpec ''ᵁ (hA.fromSpec ⁻¹ᵁ A),
          (isAffineOpen_pre hA).image_of_isOpenImmersion _⟩).comap
        ((hA.fromSpec.appIso (hA.fromSpec ⁻¹ᵁ A)).inv.hom) :=
  Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion _ _ _

/-! ## ★出典の紐付け(`.src`) -/

def comap_fromSpec_pre.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——fromSpec⁻¹ 座標での転送)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
