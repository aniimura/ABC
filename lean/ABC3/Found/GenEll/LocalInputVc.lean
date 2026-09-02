/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.UnipTransfer
import ABC3.Found.GenEll.PointTransport
import ABC3.Meta.Claim

/-!
# 第 1313 ブロック —— **局所の 2 つの入力は変数変換で移る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——Tate モデルから `E` へ（残り (a) の要）

第 1311・1312 が要る局所の 2 つの入力は

* `P₀ ∈ E(L_v)` が位数 `l`
* `L_v` の `l`-捩れはたかだか `l` 個

であり、どちらも **Tate モデルの側**で第 1297・1304 が与える。
★`E ⊗ L_v` と Tate モデルは変数変換で結ばれる（`exists_tate_model`、在庫）ので、
本ブロックで 2 つとも移す。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {F : Type} [Field F]

/-- ★★★★★★**位数 `l` の点は変数変換で戻る**——★**無条件**（第 1313）。 -/
theorem exists_point_order_of_vc (W : WeierstrassCurve F) (C : VariableChange F)
    [W.IsElliptic] [(C • W).IsElliptic] {l : ℕ}
    (P : (C • W).toAffine.Point) (hP : addOrderOf P = l) :
    ∃ P₀ : W.toAffine.Point, addOrderOf P₀ = l := by
  obtain ⟨P₀, hP₀⟩ := vcPoint_surjective C W P
  refine ⟨P₀, ?_⟩
  have h1 : addOrderOf (vcPoint C W P₀) = addOrderOf P₀ := addOrderOf_vcPoint C W P₀
  have h2 : vcPoint C W P₀ = P := hP₀
  rw [h2, hP] at h1
  exact h1.symm

/-- ★★★★★★★★**`l`-捩れの個数の上界は変数変換で戻る**——★**無条件**（第 1313）。 -/
theorem card_le_of_vc (W : WeierstrassCurve F) (C : VariableChange F)
    [W.IsElliptic] [(C • W).IsElliptic] {l : ℕ}
    (h : ∀ T : Finset ((C • W).toAffine.Point), (∀ p ∈ T, l • p = 0) → T.card ≤ l) :
    ∀ T : Finset (W.toAffine.Point), (∀ p ∈ T, l • p = 0) → T.card ≤ l := by
  intro T hT
  have hinj : Set.InjOn (vcPoint C W) T := by
    intro a _ b _ hab
    exact vcPoint_injective C W hab
  have hcard : (T.image (vcPoint C W)).card = T.card := Finset.card_image_of_injOn hinj
  have htors : ∀ p ∈ T.image (vcPoint C W), l • p = 0 := by
    intro p hp
    obtain ⟨q, hq, rfl⟩ := Finset.mem_image.1 hp
    exact vcPoint_nsmul_eq_zero C W (hT q hq)
  have := h (T.image (vcPoint C W)) htors
  rwa [hcard] at this

/-! ## ★出典の紐付け(`.src`) -/

def exists_point_order_of_vc.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(位数 l の点は変数変換で戻る。★無条件)",
    sectionId := "genell-thm-3-8" }

def card_le_of_vc.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l-捩れの個数の上界は変数変換で戻る。★無条件)",
    sectionId := "genell-thm-3-8" }

def card_le_of_vc.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1313）**——Tate モデルから `E` へ戻す段である。" ++
       "☆`vcPoint` は全単射（在庫）で `l`-捩れを保つので、位数も個数も移る。") 2 ]

end ABC3.Found.GenEll
