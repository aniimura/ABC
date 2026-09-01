/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.PointCoordNatural
import ABC3.Meta.Claim

/-!
# 第 1285 ブロック —— **`σ` が固定する曲線の上の点の作用**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——節点 8 の下ごしらえ（2）

☆`σK : K →+* K` が曲線 `W` の係数を固定する（`W.map σK = W`）とき、
点の上の作用 `galAct` を**加法準同型として**定める。

★機構は 2 段:

| 段 | 内容 |
|---|---|
| 1 | `rhPoint` を加法準同型にする（`rhPoint_add`・`rhPoint_zero` は在庫） |
| 2 | 曲線の等式に沿って点を運ぶ（`subst` で `AddEquiv` になる） |

☆座標の言葉では `galAct` はただ `σK` を当てるだけである
——これが第 1283（Tate 一意化の同変性）と突き合わせる形である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Meta

open scoped Classical

variable {F K : Type*} [Field F] [Field K]

/-- ★★★★★★**`rhPoint` を加法準同型として**（第 1285）。 -/
noncomputable def rhPointHom (f : F →+* K) (W : WeierstrassCurve F) :
    W.toAffine.Point →+ (W.map f).toAffine.Point where
  toFun := rhPoint f W
  map_zero' := rhPoint_zero f W
  map_add' := rhPoint_add f W

/-- ★★★★★★**曲線の等式に沿って点を運ぶ**（第 1285）。 -/
noncomputable def pointEquivOfEq {W₁ W₂ : WeierstrassCurve K} (h : W₁ = W₂) :
    W₁.toAffine.Point ≃+ W₂.toAffine.Point := by
  subst h
  exact AddEquiv.refl _

/-- ★★★★**運んでも座標は変わらない**（第 1285）。 -/
theorem pointCoords_pointEquivOfEq {W₁ W₂ : WeierstrassCurve K} (h : W₁ = W₂)
    (P : W₁.toAffine.Point) :
    pointCoords (pointEquivOfEq h P) = pointCoords P := by
  subst h
  rfl

/-- ★★★★★★★★**`σ` が固定する曲線の上の点の作用**（第 1285）。 -/
noncomputable def galAct (σK : K →+* K) (W : WeierstrassCurve K) (hW : W.map σK = W) :
    W.toAffine.Point →+ W.toAffine.Point :=
  (pointEquivOfEq hW).toAddMonoidHom.comp (rhPointHom σK W)

/-- ★★★★★★★★**座標の言葉では `galAct` は `σK` を当てるだけ**——★**無条件**（第 1285）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが第 1283（Tate 一意化の同変性）と突き合わせる形である。 -/
theorem pointCoords_galAct (σK : K →+* K) (W : WeierstrassCurve K) (hW : W.map σK = W)
    (P : W.toAffine.Point) :
    pointCoords (galAct σK W hW P)
      = (σK (pointCoords P).1, σK (pointCoords P).2) := by
  show pointCoords (pointEquivOfEq hW (rhPoint σK W P)) = _
  rw [pointCoords_pointEquivOfEq, pointCoords_rhPoint']

/-! ## ★出典の紐付け(`.src`) -/

def rhPointHom.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(rhPoint を加法準同型として)",
    sectionId := "genell-thm-3-8" }

def pointEquivOfEq.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(曲線の等式に沿って点を運ぶ)",
    sectionId := "genell-thm-3-8" }

def pointCoords_pointEquivOfEq.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(運んでも座標は変わらない。★無条件)",
    sectionId := "genell-thm-3-8" }

def galAct.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ が固定する曲線の上の点の作用)",
    sectionId := "genell-thm-3-8" }

def pointCoords_galAct.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(座標の言葉では galAct は σK を当てるだけ。★無条件)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GenEll
