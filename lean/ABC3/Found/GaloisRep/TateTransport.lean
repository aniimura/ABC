/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.QTorsion
import ABC3.Meta.Claim

/-!
# 第 1227 ブロック —— **Tate 一意化で `T_l E ≅ ℤ_l²` を移す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——α 橋の II 側の配管

`tateModule_qQuot`（在庫）は `Kˣ/qℤ` の側で `T_l ≅ ℤ_l²` を与える。
☆Tate 一意化 `Φ : Additive (Kˣ/qℤ) ≃+ E(K)` があれば、
`E` の側でも同じことが言える。

★`addEquiv_limTors` が要るのは「`#A[m] = m²`」だけなので、
`Φ` で捩れ部分群を移せば済む。

★★★これが `α` を `galRep` の行列として読む段の**基底の側**の配管である。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Meta

universe u v

/-- ☆加法同型は `m`-捩れを移す。 -/
def kerNsmulEquiv {A : Type u} {B : Type v} [AddCommGroup A] [AddCommGroup B]
    (Φ : A ≃+ B) (m : ℕ) : (nsmulHom A m).ker ≃ (nsmulHom B m).ker :=
  Equiv.subtypeEquiv Φ.toEquiv (fun a => by
    show m • a = 0 ↔ m • (Φ a) = 0
    constructor
    · intro h
      rw [← map_nsmul, h, map_zero]
    · intro h
      have h2 : Φ (m • a) = Φ 0 := by rw [map_nsmul, h, map_zero]
      exact Φ.injective h2)

/-- ★★★★★★★★★★★★
**Tate 一意化で `T_l ≅ ℤ_l²` を移す**——★**無条件**（第 1227）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`addEquiv_limTors` が要るのは「`#A[m] = m²`」だけなので、
加法同型 `Φ` で捩れ部分群を移せば済む。

★★★これが `α` を `galRep` の行列として読む段の**基底の側**の配管である。 -/
theorem addEquiv_limTors_of_addEquiv {A : Type u} {B : Type u}
    [AddCommGroup A] [AddCommGroup B] (Φ : A ≃+ B)
    (hfin : ∀ m : ℕ, 1 ≤ m → Finite (nsmulHom A m).ker)
    (hcard : ∀ m : ℕ, 1 ≤ m → Nat.card (nsmulHom A m).ker = m ^ 2)
    (l : ℕ) [Fact l.Prime] :
    Nonempty (limTors B l ≃+ (ℤ_[l] × ℤ_[l])) := by
  refine addEquiv_limTors (A := B) ?_ ?_ l
  · intro m hm
    haveI := hfin m hm
    exact Finite.of_equiv _ (kerNsmulEquiv Φ m)
  · intro m hm
    rw [← Nat.card_congr (kerNsmulEquiv Φ m)]
    exact hcard m hm

/-! ## ★出典の紐付け(`.src`) -/

def kerNsmulEquiv.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(加法同型は m-捩れを移す。★無条件)",
    sectionId := "genell-thm-3-8" }

def addEquiv_limTors_of_addEquiv.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 一意化で T_l ≅ Z_l² を移す。★無条件)",
    sectionId := "genell-thm-3-8" }

def addEquiv_limTors_of_addEquiv.needs : List ProofObligation :=
  [ .citation "[ABC3]" "addEquiv_limTors(証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.addEquiv_limTors") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1227）**——`α` を `galRep` の行列として読む段の" ++
       "**基底の側**の配管である。☆`addEquiv_limTors` が要るのは" ++
       "「`#A[m] = m²`」だけなので、加法同型で捩れ部分群を移せば済む。") 2 ]

end ABC3.Found.GaloisRep
