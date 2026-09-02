/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38PiFromPhi
import ABC3.Meta.Claim

/-!
# 第 1304 ブロック —— **`l ∤ v(Q)` なら基礎体の `l`-捩れは `l` 個以下**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`h1` の入力（`E[l] ⊄ E(K₀)`）の本体

第 1279 は「`l`-捩れが `l` 個より多ければ `Q` は `l` 乗になる」と言う。
★その**対偶**が本ブロックである:

> `l ∤ v(Q)` なら、基礎体の上の `l`-捩れはたかだか `l` 個。

☆`E[l]` は代数閉体の上では `l²` 個あるので、`l ≥ 2` なら
**`E[l]` のどれかは基礎体に無い**——これが `h1`（第 1301）の入力である。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Meta

variable {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★★★★★★★★★★★
**`l ∤ v(Q)` なら基礎体の `l`-捩れは `l` 個以下**——★**無条件**（第 1304）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1279 の対偶である——`l` 個より多ければ `π^l = Q` なる `π` が取れ、
`v(Q) = l·v(π)` になってしまう。 -/
theorem card_torsion_le_of_not_dvd (S : TateSetup R I K) {l : ℕ} (hl : l.Prime)
    {ζ : Kˣ} (hζ : IsPrimitiveRoot ((ζ : K)) l)
    {P : Type*} [AddCommGroup P] (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (hnd : ¬ ((l : ℤ) ∣ vAdd S.v S.Q))
    (T : Finset P) (hT : ∀ p ∈ T, l • p = 0) : T.card ≤ l := by
  by_contra hcon
  push_neg at hcon
  obtain ⟨π, hπ⟩ := exists_pi_of_phi S hl hζ Φ T hT hcon
  refine hnd ⟨vAdd S.v π, ?_⟩
  rw [← hπ, ← zpow_natCast π l, vAdd_zpow]

/-! ## ★出典の紐付け(`.src`) -/

def card_torsion_le_of_not_dvd.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l ∤ v(Q) なら基礎体の l-捩れは l 個以下。★無条件)",
    sectionId := "genell-thm-3-8" }

def card_torsion_le_of_not_dvd.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_pi_of_phi(第 1279、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pi_of_phi") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1304）**——`h1` の入力（`E[l] ⊄ E(K₀)`）の本体である。" ++
       "☆`E[l]` は代数閉体の上では `l²` 個あるので、`l ≥ 2` なら" ++
       "**どれかは基礎体に無い**。") 2 ]

end ABC3.Found.GenEll
