import ABC3.Found.Arakelov.PicTildeDesc

/-!
# Arakelov (B1) 第 85 ブロック —— **基本開集合での切断は局所化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★比較射が同型であることの土台

第 84 ブロックで層の射

    tildeTensorDesc : tensorModules (tilde M) (tilde N)  ⟶  tilde (M ⊗_R N)

が出た。★これが**同型**であることを示すのが残りである。

## ★★★筋 —— 前層射が基本開集合で同型

★★**第 83 の前層射は、基本開集合 `D(f)` の上で既に同型である**:

    (tilde M)(D f) ⊗_{R_f} (tilde N)(D f)  =  M_f ⊗_{R_f} N_f  ≅  (M ⊗_R N)_f  =  (tilde (M⊗N))(D f)

★★★**可逆性は要らない**——任意の `M`・`N` で成り立つ(第 79 ブロック)。

★基本開集合は基底なので、前層射は**局所全単射**であり、
層化して降ろした `tildeTensorDesc` は同型になる。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `tildeAwaySection` | ★`(tilde M)(D f)` は `M_f`(mathlib の在庫に名を付ける) |
| `tildeAwayEquiv` | ★★**`M_f ≃ₗ[R] (tilde M)(D f)`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace TensorProduct

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (f : (R : Type u))

/-- ★**基本開集合での切断は `M_f` である**——mathlib の在庫に名を付ける。 -/
theorem tildeAwaySection : IsLocalizedModule (Submonoid.powers f)
    (tilde.toOpen M (PrimeSpectrum.basicOpen f)).hom := inferInstance

/-- ★★**`M_f ≃ₗ[R] (tilde M)(D f)`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで前層射の `D(f)` 成分が第 79 ブロックの `localizedTensorEquiv` と
突き合わせられるようになる。 -/
noncomputable def tildeAwayEquiv :
    LocalizedModule (Submonoid.powers f) M
      ≃ₗ[(R : Type u)] ((AlgebraicGeometry.modulesSpecToSheaf.obj (tilde M)).presheaf.obj
        (op (PrimeSpectrum.basicOpen f))) :=
  IsLocalizedModule.iso (Submonoid.powers f) (tilde.toOpen M (PrimeSpectrum.basicOpen f)).hom

/-! ## ★出典の紐付け(`.src`) -/

def tildeAwayEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——基本開集合での切断は局所化)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
