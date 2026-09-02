/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.LocalInputVc
import ABC3.Meta.Claim

/-!
# 第 1316 ブロック —— **曲線の等式つきの変数変換で移す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——第 1313 を「等しい曲線」で受ける形

第 1315 が与えるのは `C' • W = (Tate モデル)` という**曲線の等式**である。
★第 1313 は `C • W` の上の点をそのまま受けるので、
等式で結ばれた曲線 `V` を受け取れるように一般化しておく（`subst` するだけ）。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Meta

open scoped Classical

variable {F : Type} [Field F]

/-- ★★★★★★**位数 `l` の点は等式つきの変数変換で戻る**——★**無条件**（第 1316）。 -/
theorem exists_point_order_of_vc' (W : WeierstrassCurve F) (C : VariableChange F)
    (V : WeierstrassCurve F) (hV : V = C • W)
    [W.IsElliptic] [(C • W).IsElliptic] {l : ℕ}
    (P : V.toAffine.Point) (hP : addOrderOf P = l) :
    ∃ P₀ : W.toAffine.Point, addOrderOf P₀ = l := by
  subst hV
  exact exists_point_order_of_vc W C P hP

/-- ★★★★★★★★**個数の上界も等式つきの変数変換で戻る**——★**無条件**（第 1316）。 -/
theorem card_le_of_vc' (W : WeierstrassCurve F) (C : VariableChange F)
    (V : WeierstrassCurve F) (hV : V = C • W)
    [W.IsElliptic] [(C • W).IsElliptic] {l : ℕ}
    (h : ∀ T : Finset (V.toAffine.Point), (∀ p ∈ T, l • p = 0) → T.card ≤ l) :
    ∀ T : Finset (W.toAffine.Point), (∀ p ∈ T, l • p = 0) → T.card ≤ l := by
  subst hV
  exact card_le_of_vc W C h

/-! ## ★出典の紐付け(`.src`) -/

def exists_point_order_of_vc'.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(位数 l の点は等式つきの変数変換で戻る。★無条件)",
    sectionId := "genell-thm-3-8" }

def card_le_of_vc'.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(個数の上界も等式つきの変数変換で戻る。★無条件)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GenEll
