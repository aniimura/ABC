import ABC3.Found.Arakelov.PicAwayScalar

/-!
# Arakelov (B1) 第 125 ブロック —— **`𝒪(D f)` 線型な `M_f ≅ Γ(tilde M, D f)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★2 つの書き方を**使い分ける**のが鍵だった

§9-139 で 5 手詰まった所である。★正体は、前層の切断に

| 書き方 | 持っている構造 |
|---|---|
| `(tilde M).val.obj (op U)` | ★**`Module 𝒪(U)`** |
| `Γ(tilde M, U)`(= `modulesSpecToSheaf` 側) | ★**`Module R`・`IsScalarTower`** |

の**2 通り**があり、台は `rfl` で一致する(第 121)のに
**instance はそれぞれ片方しか持たない**ことだった。

★★★**解決**: `letI` で
**`Module` は `val` 側から、`IsScalarTower` は `Γ` 側から**取る。
★★どちらか一方に揃えようとすると失敗する(実測)。

★★★★これは [[ring-instance-two-paths]] の **4 例目**であり、
**「両方の経路を使う」**という新しい抜け方である。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `tildeAwayEquivScalar` | ★★★★★★**`𝒪(D f)` 線型な `M_f ≅ Γ(tilde M, D f)`** |

## ★★★これで生成元の乗法が `𝒪(D f)` の側で扱える

局所自明性の最後の材料である。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (f : (R : Type u))

/-- ★★★★★★**`𝒪(D f)` 線型な `M_f ≅ Γ(tilde M, D f)`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`Module` は `val` 側から、`IsScalarTower` は `Γ` 側から取るのが要点である
——台は `rfl` で一致するが、instance はそれぞれ片方しか持たない。 -/
noncomputable def tildeAwayEquivScalar :
    LocalizedModule (Submonoid.powers f) M
      ≃ₗ[(Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)]
        (Γ(tilde M, PrimeSpectrum.basicOpen f) : Type u) :=
  letI : Module (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)
      (((modulesSpecToSheaf.obj (tilde M)).presheaf.obj
        (op (PrimeSpectrum.basicOpen f))) : Type u) :=
    inferInstanceAs (Module (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)
      (((tilde M).val.obj (op (PrimeSpectrum.basicOpen f))) : Type u))
  letI : IsScalarTower (R : Type u) (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)
      (((modulesSpecToSheaf.obj (tilde M)).presheaf.obj
        (op (PrimeSpectrum.basicOpen f))) : Type u) :=
    inferInstanceAs (IsScalarTower (R : Type u)
      (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)
      (Γ(tilde M, PrimeSpectrum.basicOpen f) : Type u))
  (tildeAwayEquiv R M f).extendScalarsOfIsLocalization (Submonoid.powers f)
    (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)

/-! ## ★出典の紐付け(`.src`) -/

def tildeAwayEquivScalar.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——𝒪(D f) 線型な M_f ≅ Γ(tilde M, D f))",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
