/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluGoodPrime
import ABC3.Found.GaloisRep.VeluIntegralFromCoords
import ABC3.Meta.Claim

/-!
# 第 1432 ブロック —— **良い素点：核の座標が整なら `p ∣ l` でも半安定**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1387 から `hlu` を外す

`semistableAt_veluQuot_good`（第 1387＋1402）は
`hlu : IsUnit (l : primeSubring p)`（すなわち `p ∤ l`）を仮定していたが、
★使いどころは**核の座標の整性 2 か所だけ**である:

1. `isIntegral_veluQuotientFull_of_addOrderOf_prime`（第 1074）——商の整性
2. `veluKernel_coords_mem`（第 1386）——`N = ∏(2y+a₁x+a₃)` の整性

☆どちらも「核の座標が `primeSubring p` に入る」で置き換えられる
（(1) は第 1420、(2) は `veluKernelNorm_mem_primeSubring` にそのまま渡すだけ）。

★★★したがって **`p ∣ l` でも、核に深い点が無ければ良い素点は閉じる**。

☆残るのは「良い素点・`p ∣ l`・核に**深い点がある**」場合だけである
——それは形式群の `l`-捩れであり、局所体の分岐指数が `e ≥ l−1` のときにしか起きない。
★第 1431 でその場合も **`jExp p E′ ≥ 0` 1 本**に帰着させてある。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GaloisRep ABC3.Skeleton.GenEll ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★
**良い素点では Vélu の商は半安定**——★**核の座標が整でありさえすればよい**（第 1432）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1387 の `hlu`（`p ∤ l`）を仮説 `hmem` に置き換えた形である。
★`p ∣ l` でも、核に深い点（形式群の点）が無ければこれで閉じる。 -/
theorem semistableAt_veluQuot_good_of_mem [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    [WeierstrassCurve.IsIntegral (primeSubring p) E]
    (hΔ0 : valAdd p (Units.mk0 E.Δ E.isUnit_Δ.ne_zero) = 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hmem : ∀ k ∈ (range l).erase 0,
      (pointCoords (k • Q)).1 ∈ primeSubring p ∧ (pointCoords (k • Q)).2 ∈ primeSubring p)
    (hΔ' : (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ ≠ 0) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  set S := ((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)) with hS
  have hid := disc_pow_eq_veluQuot_mul E hl hodd Q hQ
  haveI hintE' : WeierstrassCurve.IsIntegral (primeSubring p) (veluQuotientFull E S) :=
    isIntegral_veluQuotientFull_of_pointCoords_mem p E hl Q hQ hmem
  have hSmem : ∀ q ∈ S, q.1 ∈ primeSubring p ∧ q.2 ∈ primeSubring p := by
    intro q hq
    rw [hS] at hq
    obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 hq
    exact hmem k hk
  have hNint : veluKernelNorm E S ∈ primeSubring p :=
    veluKernelNorm_mem_primeSubring p E S hSmem
  have hNne : veluKernelNorm E S ≠ 0 := by
    intro h0
    rw [h0] at hid
    simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero] at hid
    exact pow_ne_zero l E.isUnit_Δ.ne_zero hid
  exact semistableAt_of_disc_pow_eq p E (veluQuotientFull E S)
    E.isUnit_Δ.ne_zero hΔ' hintE' hΔ0 hNne hNint hid

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuot_good_of_mem.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い素点では核の座標が整なら Vélu の商は半安定。★p ∤ l 不要)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_good_of_mem.needs : List ProofObligation :=
  [ .citation "[ABC3]" "disc_pow_eq_veluQuot_mul(第 1402、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.disc_pow_eq_veluQuot_mul") 1,
    .citation "[ABC3]" "isIntegral_veluQuotientFull_of_pointCoords_mem(第 1420、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isIntegral_veluQuotientFull_of_pointCoords_mem") 1,
    .citation "[ABC3]" "semistableAt_of_disc_pow_eq(第 1385、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_of_disc_pow_eq") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1432）**——第 1387 が `hlu`（`p ∤ l`）を使うのは" ++
       "**核の座標の整性 2 か所だけ**だと測った。☆それを仮説にした。" ++
       "★★★これで `p ∣ l` でも、核に深い点（形式群の点）が無ければ良い素点は閉じる。" ++
       "☆残るのは「良い素点・`p ∣ l`・核に深い点がある」場合だけであり、" ++
       "そこは第 1431 で `jExp p E′ ≥ 0` 1 本に帰着させてある。") 17 ]

end ABC3.Found.GenEll
