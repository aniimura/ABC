import ABC3.Found.Arakelov.PicComapPair

/-!
# Arakelov (B2) 第 221 ブロック —— ★★★★★★**全射性の接続**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★第 220 を第 218 に渡した

    第 220(転送なしのイデアルの等式)
      → 第 218(イデアルが `map` で書ければ全射)
      → アフィンの対で `pullIdealHom.app` は全射

★★仮定 `hcompat`(`appLE` が `appIso`/`ΓSpecIso` の合成に等しい)は、
第 213–217 で部品が揃っているので**次のブロックで消せる**。
★本ブロックは**仮定として受けて先に接続を通す**——
§9-244 と同じ設計判断(転送に依存しない形にしておく)である。

## ★★これで B2 の残りは「`hcompat` を示す」だけ

    hcompat: appLE (B') (A') = appIso_B ≫ ΓSpecIso_B ≫ appLE B A ≫ ΓSpecIso_A⁻¹ ≫ appIso_A⁻¹

★これは第 214(可換正方形を元に当てる)+ 第 215(`appTop` の分解)
+ 第 217(`appIso` を元で開く)を繋げば出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `surjective_pair` | ★★★★★★**アフィンの対で比較射は全射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★★★★★★**アフィンの対で比較射は全射である**(等式を渡した形)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 220 を第 218 に渡しただけである。 -/
theorem surjective_pair {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (hle : (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens))
      ≤ f ⁻¹ᵁ (hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens)))
    (hcompat : f.appLE (hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens))
        (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens)) hle
      = (hB.fromSpec.appIso ⊤).hom ≫ ((Scheme.ΓSpecIso (Γ(Y, B))).hom
          ≫ f.appLE B A i ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv)
          ≫ (hA.fromSpec.appIso ⊤).inv) :
    Function.Surjective (((pullIdealHom f D).app
      (op (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens)))).hom) := by
  refine surjective_pullIdealHom f D (B' := ⟨hB.fromSpec ''ᵁ ⊤, (isAffineOpen_top _).image_of_isOpenImmersion _⟩)
    (A' := ⟨hA.fromSpec ''ᵁ ⊤, (isAffineOpen_top _).image_of_isOpenImmersion _⟩) hle ?_
  rw [hcompat]
  exact comap_ideal_pair f D hA hB i


/-! ## ★出典の紐付け(`.src`) -/

def surjective_pair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——全射性の接続)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
