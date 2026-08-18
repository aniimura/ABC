import ABC3.Found.Arakelov.PicIdealLoc

/-!
# Arakelov (B2) 第 154 ブロック —— **点と素イデアルの対応**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★第 92 ブロックを一般のスキームで使う

`exists_away_linearEquiv`(第 92)は `PrimeSpectrum R` の点について述べてある。
★一般のスキーム `X` のアフィン開 `A` の点 `x` に対しては

    A.2.primeIdealOf x : PrimeSpectrum Γ(X, A)

を通す。★★`fromSpec_preimage_basicOpen` と `fromSpec_primeIdealOf` で
基本開集合の帰属が**翻訳される**:

    x ∈ X.basicOpen g  ↔  primeIdealOf x ∈ PrimeSpectrum.basicOpen g

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mem_basicOpen_iff_primeIdealOf` | ★帰属の翻訳 |
| `exists_basicOpen_free` | ★★**可逆イデアルは各点の近くで基本開集合上自由** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★アフィン開の点と素イデアルの対応で基本開集合の帰属が翻訳される。 -/
theorem mem_basicOpen_iff_primeIdealOf (A : X.affineOpens) (x : A.1)
    (g : (Γ(X, A.1) : Type u)) :
    (x : X) ∈ X.basicOpen g ↔ A.2.primeIdealOf x ∈ PrimeSpectrum.basicOpen g := by
  rw [← A.2.fromSpec_preimage_basicOpen (f := g)]
  show _ ↔ (A.2.fromSpec.base (A.2.primeIdealOf x)) ∈ X.basicOpen g
  rw [A.2.fromSpec_primeIdealOf]

/-- ★★可逆イデアルは各点の近くで基本開集合上自由になる。 -/
theorem exists_basicOpen_free (D : X.IdealSheafData) (A : X.affineOpens)
    [Module.Invertible (Γ(X, A.1) : Type u) (D.ideal A)] (x : A.1) :
    ∃ g : (Γ(X, A.1) : Type u), (x : X) ∈ X.basicOpen g ∧
      Nonempty (LocalizedModule.Away g (D.ideal A)
        ≃ₗ[Localization (Submonoid.powers g)] Localization (Submonoid.powers g)) := by
  obtain ⟨g, hg, he⟩ := exists_away_linearEquiv (Γ(X, A.1) : Type u) (D.ideal A)
    (A.2.primeIdealOf x)
  exact ⟨g, (mem_basicOpen_iff_primeIdealOf A x g).2 hg, he⟩


/-! ## ★出典の紐付け(`.src`) -/

def exists_basicOpen_free.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆イデアルは局所的に自由)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
