/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.AlphaMemImage
import ABC3.Found.GaloisRep.GalRepWitness
import ABC3.Meta.Claim

/-!
# 第 1321 ブロック —— **`h2`・`h1` から `ImageContainsSL2J` へ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——II 側の**出口**

第 1249（`imageContainsSL2J_of_galTate′`）は基底 `e` を 1 つ受け取る。
★基底は `nonempty_tate_basis`（在庫、`T_l E ≅ ℤ_l²`）が与えるので、
**`h2`・`h1` を満たす `σ` があれば `ImageContainsSL2J` が出る**。

☆`h2`・`h1` は基底に依らない言明なので、`σ` だけあればよい。
-/

namespace ABC3.Found.GenEll

open ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★★★★★
**`h2`・`h1` を満たす `σ` があれば `ImageContainsSL2J`**——★**無条件**（第 1321）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆基底は `nonempty_tate_basis`（在庫）が与える。 -/
theorem imageContainsSL2J_of_h2_h1 (E : SSCurve) (l : ℕ) [Fact l.Prime] (hl5 : 5 ≤ l)
    (hno : ¬ HasLCyclicJ E l)
    (h : ∃ σ : E.alg ≃ₐ[E.fld] E.alg,
      (∀ x : E.tate l, ∃ u : E.tate l,
          galTate E.W l σ (galTate E.W l σ x) + x
            = galTate E.W l σ x + galTate E.W l σ x + l • u) ∧
        (∃ x : E.tate l, ∀ u : E.tate l, galTate E.W l σ x ≠ x + l • u)) :
    ImageContainsSL2J E l := by
  obtain ⟨σ, h2, h1⟩ := h
  obtain ⟨e⟩ := nonempty_tate_basis (L := E.alg) E.W E.isEll l
  exact imageContainsSL2J_of_galTate' E l hl5 e σ h2 h1 hno

/-! ## ★出典の紐付け(`.src`) -/

def imageContainsSL2J_of_h2_h1.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(h2・h1 を満たす σ があれば ImageContainsSL2J。★無条件)",
    sectionId := "genell-thm-3-8" }

def imageContainsSL2J_of_h2_h1.needs : List ProofObligation :=
  [ .citation "[ABC3]" "imageContainsSL2J_of_galTate′(第 1249、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.imageContainsSL2J_of_galTate'") 1,
    .citation "[ABC3]" "nonempty_tate_basis(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.nonempty_tate_basis") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1321）**——II 側の**出口**である。" ++
       "☆`h2`・`h1` は基底に依らない言明なので、`σ` だけあればよい。") 2 ]

end ABC3.Found.GenEll
