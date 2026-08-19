import ABC3.Found.GaloisRep.EdsTarget

/-!
# Galois (G1) 第 23 ブロック —— **偶数の帰納段**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★§9-359 の証明証明書をそのまま渡す

`T(j−1)`, `T(j)`, `T(j+1)` から `T(2j)` が出る。
★係数は Python で計算し(`ResearchPaper/eds-certificates.json`)、
**`linear_combination` にそのまま渡す**——`ring` に探させない。
-/

namespace ABC3.Found.GaloisRep

open scoped Classical

universe u

variable {R : Type u} [CommRing R] [CharP R 2] (b c d A : R)

/-- ★★★★**偶数の帰納段**。 -/
theorem edsTarget_even_step (j : ℤ)
    (h1 : EdsTarget b c d A j) (h2 : EdsTarget b c d A (j + 1))
    (h3 : EdsTarget b c d A (j - 1)) :
    EdsTarget b c d A (2 * j) := by
  have hz : (2 : R) = 0 := CharP.cast_eq_zero R 2
  simp only [EdsTarget, edsEps, edsP] at h1 h2 h3 ⊢
  rw [show 2 * j - 1 = 2 * (j - 1) + 1 by ring, show 2 * j + 2 = 2 * (j + 1) by ring,
    show 2 * j - 2 = 2 * (j - 1) by ring]
  rw [preNormEDS_odd, preNormEDS_even, preNormEDS_even, preNormEDS_even]
  rw [show j - 1 - 1 = j - 2 by ring, show j - 1 + 2 = j + 1 by ring,
    show j - 1 - 2 = j - 3 by ring, show j - 1 + 1 = j by ring,
    show j + 1 - 1 = j by ring, show j + 1 + 2 = j + 3 by ring,
    show j + 1 - 2 = j - 1 by ring, show j + 1 + 1 = j + 2 by ring]
  rw [if_pos (even_two_mul j)]
  by_cases hj : Even j
  · rw [if_pos hj] at h1
    rw [if_neg (by simpa [Int.even_add_one] using hj)] at h2
    rw [if_neg (by simpa [Int.even_sub] using hj)] at h3
    rw [if_pos hj, if_neg (by simpa [Int.even_sub] using hj)]
  linear_combination
    (b^4*edsP b c d (j - 2)^2*edsP b c d j^3*edsP b c d (j + 1)^4 + b^4*edsP b c d (j - 1)^4*edsP b c d j^3*edsP b c d (j + 2)^2) * h1
      + (b^8*edsP b c d (j - 2)^2*edsP b c d j^6*edsP b c d (j + 1) + edsP b c d (j - 1)^6*edsP b c d (j + 1)^3) * h2
      + (b^8*edsP b c d (j - 1)*edsP b c d j^6*edsP b c d (j + 2)^2 + edsP b c d (j - 1)^3*edsP b c d (j + 1)^6) * h3
      + (A*b^9*edsP b c d (j - 2)*edsP b c d (j - 1)^2*edsP b c d j^7*edsP b c d (j + 2)^2 + A*b^5*edsP b c d (j - 2)^2*edsP b c d (j - 1)*edsP b c d j^4*edsP b c d (j + 1)^5 + A*b*edsP b c d (j - 1)^6*edsP b c d j*edsP b c d (j + 1)^4*edsP b c d (j + 2) - b^8*edsP b c d (j - 2)^2*edsP b c d (j - 1)*edsP b c d j^6*edsP b c d (j + 1)*edsP b c d (j + 2)^2 + b^8*edsP b c d (j - 2)^2*edsP b c d j^6*edsP b c d (j + 1)^4 + b^8*edsP b c d (j - 1)^4*edsP b c d j^6*edsP b c d (j + 2)^2 - b^8*edsP b c d (j - 1)*edsP b c d j^8*edsP b c d (j + 2)^2*edsP b c d (j - 3) - 3*b^4*edsP b c d (j - 2)^2*edsP b c d (j - 1)^2*edsP b c d j^3*edsP b c d (j + 1)^4*edsP b c d (j + 2) + 2*b^4*edsP b c d (j - 2)*edsP b c d (j - 1)^4*edsP b c d j^3*edsP b c d (j + 1)^2*edsP b c d (j + 2)^2 - b^4*edsP b c d (j - 2)*edsP b c d (j - 1)^3*edsP b c d j^5*edsP b c d (j + 1)^2*edsP b c d (j + 3) - b^4*edsP b c d (j - 1)^6*edsP b c d j^3*edsP b c d (j + 2)^3 + b^4*edsP b c d (j - 1)^2*edsP b c d j^5*edsP b c d (j + 1)^3*edsP b c d (j + 2)*edsP b c d (j - 3) - edsP b c d (j - 1)^7*edsP b c d (j + 1)^3*edsP b c d (j + 2)^2 + edsP b c d (j - 1)^6*edsP b c d (j + 1)^6 - edsP b c d (j - 1)^3*edsP b c d j^2*edsP b c d (j + 1)^6*edsP b c d (j - 3)) * hz
  · rw [if_neg hj] at h1
    rw [if_pos (by simpa [Int.even_add_one] using hj)] at h2
    rw [if_pos (by simpa [Int.even_sub] using hj)] at h3
    rw [if_neg hj, if_pos (by simpa [Int.even_sub] using hj)]
  linear_combination
    (A*b*edsP b c d (j - 2)^2*edsP b c d j^4*edsP b c d (j + 1)^2*edsP b c d (j + 2) + A*b*edsP b c d (j - 2)*edsP b c d (j - 1)^2*edsP b c d j^4*edsP b c d (j + 2)^2 + edsP b c d (j - 2)^2*edsP b c d j^5*edsP b c d (j + 1)*edsP b c d (j + 3) + edsP b c d (j - 1)*edsP b c d j^5*edsP b c d (j + 2)^2*edsP b c d (j - 3)) * h1
      + (A^2*b^2*edsP b c d (j - 2)*edsP b c d (j - 1)^4*edsP b c d j^2*edsP b c d (j + 1)*edsP b c d (j + 2) + A*b^5*edsP b c d (j - 2)*edsP b c d (j - 1)^4*edsP b c d j*edsP b c d (j + 1)^3 + A*b^5*edsP b c d (j - 1)^6*edsP b c d j*edsP b c d (j + 1)*edsP b c d (j + 2) + A*b*edsP b c d (j - 2)^2*edsP b c d (j - 1)^3*edsP b c d j*edsP b c d (j + 1)^2*edsP b c d (j + 2) + A*b*edsP b c d (j - 2)^2*edsP b c d (j - 1)*edsP b c d j^4*edsP b c d (j + 1)^2 + A*b*edsP b c d (j - 2)*edsP b c d (j - 1)^5*edsP b c d j*edsP b c d (j + 2)^2 + A*b*edsP b c d (j - 2)*edsP b c d (j - 1)^4*edsP b c d j^3*edsP b c d (j + 3) + A*b*edsP b c d (j - 1)^3*edsP b c d j^3*edsP b c d (j + 1)*edsP b c d (j + 2)*edsP b c d (j - 3) + b^4*edsP b c d (j - 2)^2*edsP b c d (j - 1)^3*edsP b c d (j + 1)^4 + b^4*edsP b c d (j - 1)^7*edsP b c d (j + 2)^2 + b^4*edsP b c d (j - 1)^6*edsP b c d j^2*edsP b c d (j + 3) + b^4*edsP b c d (j - 1)^3*edsP b c d j^2*edsP b c d (j + 1)^3*edsP b c d (j - 3) + edsP b c d (j - 2)^3*edsP b c d j^3*edsP b c d (j + 1)^3 + edsP b c d (j - 2)^2*edsP b c d (j - 1)^4*edsP b c d (j + 1)*edsP b c d (j + 2)^2 + edsP b c d (j - 2)^2*edsP b c d (j - 1)^3*edsP b c d j^2*edsP b c d (j + 1)*edsP b c d (j + 3) + edsP b c d (j - 2)^2*edsP b c d (j - 1)^2*edsP b c d j^3*edsP b c d (j + 1)*edsP b c d (j + 2) + edsP b c d (j - 1)^4*edsP b c d j^2*edsP b c d (j + 2)^2*edsP b c d (j - 3) + edsP b c d (j - 1)^3*edsP b c d j^4*edsP b c d (j + 3)*edsP b c d (j - 3)) * h2
      + (A^2*b^2*edsP b c d (j - 1)^3*edsP b c d j^2*edsP b c d (j + 1)^2*edsP b c d (j + 2)^2 + A*b*edsP b c d (j - 1)^2*edsP b c d j^4*edsP b c d (j + 1)*edsP b c d (j + 2)^2 + edsP b c d (j - 2)*edsP b c d (j - 1)*edsP b c d j^3*edsP b c d (j + 1)^2*edsP b c d (j + 2)^2 + edsP b c d (j - 1)^5*edsP b c d (j + 2)^4 + edsP b c d (j - 1)^3*edsP b c d j^4*edsP b c d (j + 3)^2 + edsP b c d (j - 1)^3*edsP b c d j^3*edsP b c d (j + 2)^3) * h3
      + (A^3*b^3*edsP b c d (j - 2)*edsP b c d (j - 1)^4*edsP b c d j^3*edsP b c d (j + 1)^2*edsP b c d (j + 2)^2 + A^2*b^6*edsP b c d (j - 2)*edsP b c d (j - 1)^4*edsP b c d j^2*edsP b c d (j + 1)^4*edsP b c d (j + 2) + A^2*b^6*edsP b c d (j - 1)^6*edsP b c d j^2*edsP b c d (j + 1)^2*edsP b c d (j + 2)^2 + A^2*b^2*edsP b c d (j - 2)^2*edsP b c d (j - 1)*edsP b c d j^5*edsP b c d (j + 1)^3*edsP b c d (j + 2) + A^2*b^2*edsP b c d (j - 2)*edsP b c d (j - 1)^3*edsP b c d j^5*edsP b c d (j + 1)*edsP b c d (j + 2)^2 + A*b^9*edsP b c d (j - 1)^6*edsP b c d j*edsP b c d (j + 1)^4*edsP b c d (j + 2) + A*b^5*edsP b c d (j - 2)^2*edsP b c d (j - 1)^3*edsP b c d j*edsP b c d (j + 1)^5*edsP b c d (j + 2) + A*b^5*edsP b c d (j - 2)^2*edsP b c d (j - 1)*edsP b c d j^4*edsP b c d (j + 1)^5 + A*b^5*edsP b c d (j - 1)^3*edsP b c d j^3*edsP b c d (j + 1)^4*edsP b c d (j + 2)*edsP b c d (j - 3) - A*b*edsP b c d (j - 2)^2*edsP b c d (j - 1)^2*edsP b c d j^4*edsP b c d (j + 1)^2*edsP b c d (j + 2)^2 - A*b*edsP b c d (j - 2)*edsP b c d (j - 1)^5*edsP b c d j^3*edsP b c d (j + 2)^2*edsP b c d (j + 3) + A*b*edsP b c d (j - 2)*edsP b c d (j - 1)^2*edsP b c d j^7*edsP b c d (j + 2)^2 + b^8*edsP b c d (j - 2)^2*edsP b c d (j - 1)^3*edsP b c d (j + 1)^7 + b^8*edsP b c d (j - 1)^6*edsP b c d j^2*edsP b c d (j + 1)^3*edsP b c d (j + 3) + b^4*edsP b c d (j - 2)^3*edsP b c d j^3*edsP b c d (j + 1)^6 - 2*b^4*edsP b c d (j - 2)^2*edsP b c d (j - 1)^2*edsP b c d j^3*edsP b c d (j + 1)^4*edsP b c d (j + 2) + 3*b^4*edsP b c d (j - 2)*edsP b c d (j - 1)^4*edsP b c d j^3*edsP b c d (j + 1)^2*edsP b c d (j + 2)^2 - b^4*edsP b c d (j - 2)*edsP b c d (j - 1)^3*edsP b c d j^5*edsP b c d (j + 1)^2*edsP b c d (j + 3) - b^4*edsP b c d (j - 1)^7*edsP b c d j^2*edsP b c d (j + 2)^2*edsP b c d (j + 3) + b^4*edsP b c d (j - 1)^2*edsP b c d j^5*edsP b c d (j + 1)^3*edsP b c d (j + 2)*edsP b c d (j - 3) - edsP b c d (j - 2)^3*edsP b c d (j - 1)*edsP b c d j^3*edsP b c d (j + 1)^3*edsP b c d (j + 2)^2 - edsP b c d (j - 2)^3*edsP b c d j^5*edsP b c d (j + 1)^3*edsP b c d (j + 3) - edsP b c d (j - 2)^2*edsP b c d (j - 1)^5*edsP b c d (j + 1)*edsP b c d (j + 2)^4 - edsP b c d (j - 2)^2*edsP b c d (j - 1)^4*edsP b c d j^2*edsP b c d (j + 1)*edsP b c d (j + 2)^2*edsP b c d (j + 3) - edsP b c d (j - 2)^2*edsP b c d (j - 1)^3*edsP b c d j^4*edsP b c d (j + 1)*edsP b c d (j + 3)^2 - edsP b c d (j - 2)^2*edsP b c d (j - 1)^3*edsP b c d j^3*edsP b c d (j + 1)*edsP b c d (j + 2)^3 - edsP b c d (j - 2)^2*edsP b c d (j - 1)^2*edsP b c d j^5*edsP b c d (j + 1)*edsP b c d (j + 2)*edsP b c d (j + 3) + edsP b c d (j - 2)^2*edsP b c d j^8*edsP b c d (j + 1)*edsP b c d (j + 3) - edsP b c d (j - 2)*edsP b c d (j - 1)*edsP b c d j^5*edsP b c d (j + 1)^2*edsP b c d (j + 2)^2*edsP b c d (j - 3) - edsP b c d (j - 1)^5*edsP b c d j^2*edsP b c d (j + 2)^4*edsP b c d (j - 3) - edsP b c d (j - 1)^4*edsP b c d j^4*edsP b c d (j + 2)^2*edsP b c d (j + 3)*edsP b c d (j - 3) - edsP b c d (j - 1)^3*edsP b c d j^6*edsP b c d (j + 3)^2*edsP b c d (j - 3) - edsP b c d (j - 1)^3*edsP b c d j^5*edsP b c d (j + 2)^3*edsP b c d (j - 3)) * hz

/-! ## ★出典の紐付け(`.src`) -/

def edsTarget_even_step.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ω_n の分子——偶数の帰納段)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
