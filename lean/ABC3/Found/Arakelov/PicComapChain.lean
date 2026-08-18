import ABC3.Found.Arakelov.PicCartierComap
import ABC3.Found.Arakelov.PicFromSpecTransport

/-!
# Arakelov (B2) 第 210 ブロック —— **`comap` のイデアルの連鎖**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-230 で足りないと分かった等式

第 202 は**可逆性**を連鎖で運ぶが、**イデアルの等式そのもの**は
`fromSpec` の転送に埋もれて取り出せなかった。本ブロックはそれを**等式として**取り出す。

## ★★2 本に分ける

| 定理 | 内容 |
|---|---|
| `comap_ideal_chain` | ★★`Spec` の ⊤ で見た連鎖(`comap_decomp` + `comap_ideal_top` + `appTop_decomp`) |
| `comap_fromSpec_top` | ★`fromSpec` の ⊤ での `comap` は同型に沿った引き戻し |

★★★どちらも既存の補題の**並べ替え**で出る——
第 202 の証明の中で使っていたものを、**可逆性を経由せず**に書いただけである。

★§9-230 の教訓「中間結果を汎用の形で残すかどうかも設計判断」を、ここで**回収**した。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `comap_ideal_chain` | ★★★アフィンの対でのイデアルの連鎖 |
| `comap_fromSpec_top` | ★★`fromSpec` の ⊤ での `comap` |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★★アフィンの対でのイデアルの連鎖(`Spec` の ⊤ で見たもの)。 -/
theorem comap_ideal_chain {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B) :
    ((D.comap f).comap hA.fromSpec).ideal ⟨⊤, isAffineOpen_top _⟩
      = (((D.comap hB.fromSpec).ideal ⟨⊤, isAffineOpen_top _⟩).map
          ((Scheme.ΓSpecIso (Γ(Y, B))).hom ≫ f.appLE B A i
            ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv).hom) := by
  rw [comap_decomp f D hA hB i, comap_ideal_top, appTop_decomp f i]




/-- ★`fromSpec` の ⊤ での comap は同型に沿った引き戻しである。 -/
theorem comap_fromSpec_top {A : X.Opens} (hA : IsAffineOpen A) :
    ((D.comap f).comap hA.fromSpec).ideal ⟨⊤, isAffineOpen_top _⟩
      = ((D.comap f).ideal ⟨hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens),
          (isAffineOpen_top _).image_of_isOpenImmersion _⟩).comap
        ((hA.fromSpec.appIso ⊤).inv.hom) :=
  Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion _ _ _


/-! ## ★出典の紐付け(`.src`) -/

def comap_ideal_chain.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——comap のイデアルの連鎖)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
