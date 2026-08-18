import ABC3.Found.Arakelov.PicComapChain

/-!
# Arakelov (B2) 第 211 ブロック —— **アフィンの対でのイデアルの等式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-230 で足りないと分かった等式(転送つきの形)

    ((D.comap f).ideal A').comap eA
      = ((D.ideal B').comap eB).map (ΓSpecIso.hom ≫ appLE ≫ ΓSpecIso.inv)

ここで `A' = fromSpec ''ᵁ ⊤`、`eA = (fromSpec.appIso ⊤).inv` である。

★★第 210 の 2 本(`comap_ideal_chain` と `comap_fromSpec_gen`)を繋ぐだけで出る。

## ★★★測って分かったこと —— **転送を外した形は `whnf` で落ちる**

    (D.comap f).ideal A' = (D.ideal B').map (eB.hom ≫ … ≫ eA.hom)

という「転送を外した」形を書くと、★**`maxHeartbeats 2000000` でも `whnf` timeout** する
(2026-08-19 実測)。`Γ(X, fromSpec ''ᵁ ⊤)` の綴りが深く、合成の型合わせに
指数的な展開が起きるためである。

★★★**したがって転送つきの形のまま使う**——後段(全射性)でも
`A'` を保ったまま扱い、最後に第 200 と同じ `X.affineOpens` の motive で `cast` すればよい。
★これは [[ring-instance-two-paths]] の「型の綴り」問題が
**`rw` ではなく `whnf` の形で出た**例である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `comap_fromSpec_gen` | ★汎用版の `fromSpec` 転送 |
| `comap_ideal_image` | ★★★★**アフィンの対でのイデアルの等式** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

theorem comap_fromSpec_gen {Z : Scheme.{u}} (E : Z.IdealSheafData) {A : Z.Opens}
    (hA : IsAffineOpen A) :
    (E.comap hA.fromSpec).ideal ⟨⊤, isAffineOpen_top _⟩
      = (E.ideal ⟨hA.fromSpec ''ᵁ (⊤ : (Spec Γ(Z, A)).Opens),
          (isAffineOpen_top _).image_of_isOpenImmersion _⟩).comap
        ((hA.fromSpec.appIso ⊤).inv.hom) :=
  Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion _ _ _

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★★★★アフィンの対でのイデアルの等式。 -/
theorem comap_ideal_image {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B) :
    ((D.comap f).ideal ⟨hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens),
        (isAffineOpen_top _).image_of_isOpenImmersion _⟩).comap
      ((hA.fromSpec.appIso ⊤).inv.hom)
      = ((D.ideal ⟨hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens),
          (isAffineOpen_top _).image_of_isOpenImmersion _⟩).comap
        ((hB.fromSpec.appIso ⊤).inv.hom)).map
        ((Scheme.ΓSpecIso (Γ(Y, B))).hom ≫ f.appLE B A i
          ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv).hom := by
  rw [← comap_fromSpec_gen (D.comap f) hA, ← comap_fromSpec_gen D hB]
  exact comap_ideal_chain f D hA hB i


/-! ## ★出典の紐付け(`.src`) -/

def comap_ideal_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——アフィンの対でのイデアルの等式)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
