import ABC3.Found.Arakelov.PicIdealTransport

/-!
# Arakelov (B2) 第 220 ブロック —— ★★★★★**転送なしのイデアルの等式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-232 で「`whnf` で落ちる」と測った形が出た

    (D.comap f).ideal A' = (D.ideal B').map (appIso ≫ ΓSpecIso ≫ appLE ≫ ΓSpecIso⁻¹ ≫ appIso⁻¹)

★★§9-232 では `maxHeartbeats 2000000` でも落ちた。★★★通った理由は
**第 219(抽象の側で書いた変換)を噛ませた**からである
——スキームの綴りを一度も展開せず、`ideal_of_comap_eq` の**適用 1 回**で出る。

★§9-232 の診断(「綴りが深いから落ちる」)は正しかったが、
**逃げ道(抽象化)を見落としていた**。記録を残す意味はここにある——
測って外したことも、後で読めば逃げ道が見える。

## ★★これで第 218 の仮定が満たせる

第 218 は「イデアルが `map` で書ければ全射」という形だったので、
本ブロックの等式をそのまま渡せばよい。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `comap_ideal_pair` | ★★★★★**転送なしのイデアルの等式** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★★★★★**アフィンの対でのイデアルの等式(転送なし)**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★§9-232 で `whnf` timeout した形が、第 219 を噛ませると**適用 1 回**で出る。 -/
theorem comap_ideal_pair {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B) :
    (D.comap f).ideal ⟨hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens),
        (isAffineOpen_top _).image_of_isOpenImmersion _⟩
      = (D.ideal ⟨hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens),
          (isAffineOpen_top _).image_of_isOpenImmersion _⟩).map
        ((hB.fromSpec.appIso ⊤).hom ≫ ((Scheme.ΓSpecIso (Γ(Y, B))).hom
          ≫ f.appLE B A i ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv)
          ≫ (hA.fromSpec.appIso ⊤).inv).hom :=
  ideal_of_comap_eq (hA.fromSpec.appIso ⊤) (hB.fromSpec.appIso ⊤) _ _ _
    (comap_ideal_image f D hA hB i)


/-! ## ★出典の紐付け(`.src`) -/

def comap_ideal_pair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——転送なしのイデアルの等式)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
