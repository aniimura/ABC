import ABC3.Found.Arakelov.PicIdealRes

/-!
# Arakelov (B2) 第 160 ブロック —— ★★★★★★**制限と局所化の可換図式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★第 120 の `idealSections` 版

    idealResLin ∘ idealAwayEquiv_g = idealAwayEquiv_{t·g} ∘ liftAwayMap

★これが「制限した生成元が生成元である」ことの中身である。

| 段 | 道具 |
|---|---|
| `IsLocalizedModule.ext` の可逆性 | ★第 159 `isUnit_smul_pow` |
| `mk m 1` での値 | ★第 156 `idealAwayEquiv_mk_one` |
| `liftAwayMap` の `mk m 1` | ★第 95 `liftAwayMap_comp` |
| 制限の合成 | ★`Functor.map_comp` + `rfl`(poset なので証明無関係) |

## ★★★逃げ道——**定義と同じファイルに補題を置く**

`idealAwayEquiv_mk_one` を別ファイルで書こうとすると
`show` で `≪≫ₗ` を分解する際に型が決まらず落ちる。
★**定義の直後**(第 156 のファイル)に置くと通る——
そこでは `letI` の文脈と定義本体が見えているからである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData) (A : X.affineOpens)

theorem idealAwayEquiv_res (g t : (Γ(X, A.1) : Type u)) :
    letI := resAlg A g
    letI := resAlg A (t * g)
    (idealResLin D A g t).comp (idealAwayEquiv D A g).toLinearMap
      = (idealAwayEquiv D A (t * g)).toLinearMap.comp
        (liftAwayMap (Γ(X, A.1) : Type u) g t (D.ideal A)) := by
  letI := resAlg A g
  letI := resAlg A (t * g)
  refine IsLocalizedModule.ext (Submonoid.powers g)
    (LocalizedModule.mkLinearMap (Submonoid.powers g) (D.ideal A))
    (isUnit_smul_pow D A g t) ?_
  refine LinearMap.ext fun m => Subtype.ext ?_
  show (X.presheaf.map (homOfLE (boMul_le A g t)).op).hom
      ((idealAwayEquiv D A g (LocalizedModule.mk m 1)) : (Γ(X, X.basicOpen g) : Type u))
    = ((idealAwayEquiv D A (t * g)
        (liftAwayMap (Γ(X, A.1) : Type u) g t (D.ideal A) (LocalizedModule.mk m 1)))
      : (Γ(X, X.basicOpen (t * g)) : Type u))
  rw [show liftAwayMap (Γ(X, A.1) : Type u) g t (D.ideal A) (LocalizedModule.mk m 1)
      = LocalizedModule.mk m 1 from
    DFunLike.congr_fun (liftAwayMap_comp (Γ(X, A.1) : Type u) g t (D.ideal A)) m,
    idealAwayEquiv_mk_one, idealAwayEquiv_mk_one, ← CommRingCat.comp_apply,
    ← Functor.map_comp, ← op_comp]
  rfl


/-! ## ★出典の紐付け(`.src`) -/

def idealAwayEquiv_res.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限と局所化の可換図式)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
