/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateDeepPoints
import ABC3.Found.GaloisRep.TatePhiTwo
import ABC3.Found.GaloisRep.TateQUnique
import ABC3.Found.GaloisRep.TateVeluMu
import ABC3.Found.GaloisRep.TateLinearize
import ABC3.Found.GaloisRep.TateMultRed
import ABC3.Meta.Claim

/-!
# TateVeluDeep —— `[GenEll] Definition 3.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GenEll Finset
open scoped Classical
variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-! ## ★★深い類の代表は `𝔪` に入る -/

/-- ★★**`v(Q) ∤ i·v(y)` なら `[yⁱ]` の代表は `I` に入る**。

☆`normRep_vAdd_pos_of_not_dvd`（第 1411）で付値が正になり、
`tateAOf_mem_of_pos`（在庫）で `I` の元になる。 -/
theorem tateAOf_pow_mem (S : TateSetup R I K) (y : Kˣ) (i : ℕ)
    (hi : ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y)) :
    tateAOf S (QuotientGroup.mk (y ^ i)) ∈ I := by
  refine tateAOf_mem_of_pos S _ ?_
  refine normRep_vAdd_pos_of_not_dvd S.v S.Q S.hQ (y ^ i) ?_
  rwa [show vAdd S.v (y ^ i) = (i : ℤ) * vAdd S.v y by rw [← zpow_natCast y i, vAdd_zpow]]

/-! ## ★出典の紐付け(`.src`) -/

def tateAOf_pow_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(深い類の代表は 𝔪 に入る)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
