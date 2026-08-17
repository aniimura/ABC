import ABC3.Found.Arakelov.PicAssoc

/-!
# Arakelov (B1) 第 14 ブロック —— **可逆層の型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★`Pic X` の台

`Pic X` は**可逆層の同型類**である。★本ブロックはその「可逆層」を型にする。

★★★2 条を持たせるのは、**どちらも使う**からである(`Interface/Arakelov/LineBundle.lean`
の `IsInvertibleSheaf` と同じ):

| 条 | 何に効くか |
|---|---|
| 逆の存在 | `Pic X` の**逆元** |
| 局所的に階数 1 自由 | **結合律**(第 13 ブロック `tensorModulesAssoc` の仮定) |

★数学的には同値だが、Lean では片方から他方を出すのに Nakayama 型の議論が要るので、
**両方を持たせるのが正しい**(内容は変わらない)。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

/-- ★★★**可逆層(直線束)**——逆を持ち、局所的に階数 1 自由な層加群。 -/
structure InvertibleSheaf (X : Scheme.{u}) where
  /-- 台となる層加群。 -/
  sheaf : X.Modules
  /-- テンソル積についての逆。 -/
  inv : X.Modules
  /-- ★逆であること。 -/
  isInv : Nonempty (tensorModules sheaf inv ≅ unitModules X)
  /-- ★局所的に階数 1 自由(結合律に要る)。 -/
  rankOne : IsLocallyRankOne X sheaf.val
  /-- ★逆の側も局所的に階数 1 自由。 -/
  invRankOne : IsLocallyRankOne X inv.val

namespace InvertibleSheaf

variable {X : Scheme.{u}}

/-- ★★構造層は可逆層である(`Pic X` の**単位元**)。

★機構は `tensorUnitLeft`(第 3 ブロック)と、
「構造層は各開集合で自分自身を基底に持つ」こと。 -/
noncomputable def one (X : Scheme.{u}) (h : IsLocallyRankOne X (unitModules X).val) :
    InvertibleSheaf X where
  sheaf := unitModules X
  inv := unitModules X
  isInv := ⟨tensorUnitLeft (unitModules X)⟩
  rankOne := h
  invRankOne := h

/-- ★**逆を取る操作**——可逆層の逆もまた可逆層である。

★機構は可換性(第 2 ブロック `tensorModulesComm`)だけ。 -/
noncomputable def symm (L : InvertibleSheaf X) : InvertibleSheaf X where
  sheaf := L.inv
  inv := L.sheaf
  isInv := L.isInv.map fun e => tensorModulesComm L.inv L.sheaf ≪≫ e
  rankOne := L.invRankOne
  invRankOne := L.rankOne

@[simp] theorem symm_sheaf (L : InvertibleSheaf X) : L.symm.sheaf = L.inv := rfl

@[simp] theorem symm_symm (L : InvertibleSheaf X) : L.symm.symm = L := rfl

end InvertibleSheaf

/-! ## ★出典の紐付け(`.src`) -/

def InvertibleSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆層の型)",
    sectionId := "genell-def-1-1-i" }

def InvertibleSheaf.symm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆層の逆もまた可逆層)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
