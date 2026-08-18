import ABC3.Found.Arakelov.PicBCMate

/-!
# Arakelov (B1) 第 47 ブロック —— **生成元の制限は生成元である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★Beck–Chevalley を生成元で見るために

第 46 ブロックの mate が**生成元の上で同型**であることを示したい。
★そのためにまず「生成元を制限したものが何か」を知る必要がある:

    (yoneda W)|_V  ≅  yoneda (W ⊓ V)   (`Over V` の中で)

★★理由は 1 行——`Z ≤ V` のとき `Z ≤ W ⟺ Z ≤ W ⊓ V` だからである。

## ★★実装の要点(2026-08-18 実測)

★`Over V` の射は `.left` で決まるので、`Subsingleton` はインスタンスに無い。
★★`Over.OverMorphism.ext` で 1 行作って渡す。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `objOn` | ★`Over V` の中の `W ⊓ V` |
| `yonedaRestrictIso` | ★★★**`(yoneda W)|_V ≅ yoneda (W ⊓ V)`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {Y : Scheme.{u}} (V W : Y.Opens)

/-- ★`Over V` の中の `W ⊓ V`。 -/
noncomputable abbrev objOn : Over V := Over.mk (homOfLE (inf_le_right : W ⊓ V ≤ V))

/-- ★`Over V` の射は `.left` で決まるので subsingleton である。 -/
instance overHomSubsingleton (A B : Over V) : Subsingleton (A ⟶ B) :=
  ⟨fun _ _ => Over.OverMorphism.ext (Subsingleton.elim (α := (A.left ⟶ B.left)) _ _)⟩

/-- ★★★**生成元の制限は生成元である**: `(yoneda W)|_V ≅ yoneda (W ⊓ V)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★理由は 1 行——`Z ≤ V` のとき `Z ≤ W ⟺ Z ≤ W ⊓ V`。

★★★これが Beck–Chevalley を生成元の上で見るための道具である。 -/
noncomputable def yonedaRestrictIso :
    (Over.forget V).op ⋙ yoneda.obj W ≅ yoneda.obj (objOn V W) :=
  NatIso.ofComponents
    (fun Z => Equiv.toIso
      { toFun := fun h => Over.homMk (homOfLE (le_inf (leOfHom h) (leOfHom Z.unop.hom)))
        invFun := fun h => homOfLE (le_trans (leOfHom h.left) inf_le_left)
        left_inv := fun _ => Subsingleton.elim (α := (Z.unop.left ⟶ W)) _ _
        right_inv := fun _ => Subsingleton.elim (α := (Z.unop ⟶ objOn V W)) _ _ })
    (by
      intro Z Z' φ
      ext h
      exact Subsingleton.elim (α := (Z'.unop ⟶ objOn V W)) _ _)

/-! ## ★出典の紐付け(`.src`) -/

def yonedaRestrictIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元の制限が生成元であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
