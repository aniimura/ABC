import ABC3.Found.Arakelov.PicMapBij

/-!
# Arakelov (B1) 第 113 ブロック —— **生成元を `tilde M` の切断へ運ぶ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★局所自明性の切断

第 76 ブロックで `M_g ≅ R_g` が出た。★その生成元 `e.symm 1`(第 103 の `generatorOf`)を
第 85 ブロックの `tildeAwayEquiv` で **`(tilde M)(D g)` の切断**へ運ぶ。

★★これが第 111 ブロック(層加群版の同型)が要求する切断 `s` である。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `tildeGenSection` | ★★★**生成元に対応する `tilde M` の切断** |

## ★★★残り

    第 111 の仮定「被覆の上で s の倍が全単射」を
    第 100(基底で M_h ≅ R_h)+ 第 112(局所化は全単射を保つ)+ 第 103 で埋める
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (g : (R : Type u))

/-- ★★★**生成元に対応する `tilde M` の切断**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`M_g ≅ R_g` の生成元を、基本開集合での切断へ運んだもの。 -/
noncomputable def tildeGenSection
    (e : LocalizedModule (Submonoid.powers g) M
      ≃ₗ[Localization (Submonoid.powers g)] Localization (Submonoid.powers g)) :
    ((AlgebraicGeometry.modulesSpecToSheaf.obj (tilde M)).presheaf.obj
      (op (PrimeSpectrum.basicOpen g))) :=
  tildeAwayEquiv R M g (generatorOf e)

/-! ## ★出典の紐付け(`.src`) -/

def tildeGenSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元を tilde M の切断へ運ぶ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
