/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.TateMuPoint
import ABC3.Found.GenEll.TorsionCardLocal
import ABC3.Meta.Claim

/-!
# 第 1314 ブロック —— **Tate 一意化が与える局所の 2 つの入力**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——II 側の**局所の入力を 1 本に**

第 1311・1312 が要る局所の 2 つの入力を、Tate 一意化から**同時に**出す:

| 入力 | 出どころ |
|---|---|
| 位数 `l` の点 | `Φ[ζ]`（第 1297） |
| `l`-捩れはたかだか `l` 個 | `l ∤ v(Q)`（第 1304） |

☆受け取るのは `TateSetup`・`Φ`・`ζ`（原始 `l` 乗根）・`l ∤ v(Q)` の 4 つだけ。
★★★どれも `mkTateSetup`・`dvrTatePhiAddEquiv`（在庫、**無条件**）と
`PrimeToLocalHeights` から出る。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Meta

variable {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★★★★★★★★★★★
**Tate 一意化が局所の 2 つの入力を与える**——★**無条件**（第 1314）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1297（`Φ[ζ]` は位数 `l`）と第 1304（`l ∤ v(Q)` なら `l` 個以下）を並べただけである。 -/
theorem tate_local_inputs (S : TateSetup R I K) {l : ℕ} (hl : l.Prime)
    {ζ : Kˣ} (hζ : IsPrimitiveRoot ((ζ : K)) l)
    {P : Type*} [AddCommGroup P] (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (hnd : ¬ ((l : ℤ) ∣ vAdd S.v S.Q)) :
    (∃ p : P, addOrderOf p = l) ∧
      (∀ T : Finset P, (∀ p ∈ T, l • p = 0) → T.card ≤ l) :=
  ⟨⟨Φ (Additive.ofMul (QuotientGroup.mk ζ)), addOrderOf_tatePhi_zeta S hl Φ ζ hζ⟩,
    card_torsion_le_of_not_dvd S hl hζ Φ hnd⟩

/-! ## ★出典の紐付け(`.src`) -/

def tate_local_inputs.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 一意化が局所の 2 つの入力を与える。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_local_inputs.needs : List ProofObligation :=
  [ .citation "[ABC3]" "addOrderOf_tatePhi_zeta(第 1297、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.addOrderOf_tatePhi_zeta") 1,
    .citation "[ABC3]" "card_torsion_le_of_not_dvd(第 1304、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.card_torsion_le_of_not_dvd") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1314）**——II 側の**局所の入力を 1 本に**まとめた形である。" ++
       "☆受け取るのは `TateSetup`・`Φ`・`ζ`・`l ∤ v(Q)` の 4 つだけで、" ++
       "どれも `mkTateSetup`・`dvrTatePhiAddEquiv`（在庫、無条件）と" ++
       "`PrimeToLocalHeights` から出る。") 2 ]

end ABC3.Found.GenEll
