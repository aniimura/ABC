import ABC3.Found.Arakelov.PicGammaTensor

/-!
# Arakelov (B1) 第 75 ブロック —— **`tilde` のテンソル比較射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★随伴で転置するだけ

第 74 ブロックの切断の比較射を `tilde ⊣ Γ` の随伴で転置すると

    tilde (M ⊗_R N)  ⟶  tensorModules (tilde M) (tilde N)

が出る。★★**`(Spec R).Modules` はモノイド圏では無い**(mathlib に構造が無く、
我々の `tensorModules` も結合律に局所階数 1 が要る)ので
`Adjunction.leftAdjointOplaxMonoidal` は使えない——★手で転置する。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `tildeTensorCmp` | ★★★★**`tilde (M ⊗ N) ⟶ tilde M ⊗ tilde N`** |
| `tildeUnitIso` | ★★`tilde R ≅ 𝒪_{Spec R}`(mathlib の `tildeSelf` は **`.refl`**) |

## ★★★残り(B1)

| # | 主張 | 状態 |
|---|---|---|
| 1 | 比較射が**同型**であること | ★次(第 29 の器具の 4 度目) |
| 2 | 可逆層が `tilde` の本質像に入ること | ★★`isIso_fromTildeΓ_iff_isLocalizing` から |
| 3 | `equivPicRing` | |

★mathlib は「準連接 ⟹ `fromTildeΓ` 同型」を**まだ持たない**(ファイル内の TODO で実測)。
★★だが可逆層に限れば `IsLocalizing` の判定を通せる。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (R : CommRingCat.{u}) (M N : ModuleCat.{u} (R : Type u))

/-- ★★★★**`tilde` のテンソル比較射**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 74 ブロックの切断の比較射を随伴で転置したものである。 -/
noncomputable def tildeTensorCmp :
    tilde (M ⊗ N) ⟶ tensorModules (tilde M) (tilde N) :=
  ((tilde.adjunction (R := R)).homEquiv _ _).symm
    (((tilde.adjunction (R := R)).unit.app M ⊗ₘ (tilde.adjunction (R := R)).unit.app N)
      ≫ gammaTensorCmp R (tilde M) (tilde N))

/-- ★★**`tilde R ≅ 𝒪_{Spec R}`**——mathlib の `tildeSelf` は **`.refl`** である。 -/
noncomputable abbrev tildeUnitIso :
    tilde (ModuleCat.of (R : Type u) (R : Type u)) ≅ unitModules (Spec R) :=
  tildeSelf

/-! ## ★出典の紐付け(`.src`) -/

def tildeTensorCmp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——tilde のテンソル比較射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
