/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.VeluDiscIdentity
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GenEll.VeluQuotElliptic
import ABC3.Meta.Claim

/-!
# 第 1403 ブロック —— **良い素点（`p ∤ l`）では Vélu の商は半安定**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——**恒等式が閉じたので繋がった**

第 1402 で判別式の恒等式が閉じ、`semistableAt_veluQuot_good`（良い素点、極小モデルの上）が
無条件になった。★本ブロックは**極小モデルへの正規化**を足して、
`SemistableAt p E` と `0 ≤ jExp p E`（＝良い素点）と `p ∤ l` だけから

    SemistableAt p (veluQuotientFull E ⟨Q⟩)

を出す。

☆道は 3 行である:

1. 半安定なので `minDeltaExp p E = max(0, −jExp p E) = 0`
   （`minDeltaExp_eq_maxJ_of_semistable`、在庫）
2. 極小モデル `C • E` を取ると `v_p(Δ(C•E)) = minDeltaExp p E = 0`（`minDeltaExp_eq`、在庫）
   ——`IsMinimal` から `IsIntegral` はインスタンスで出る
3. `semistableAt_veluQuot_good` を `C • E` に当て、
   `veluQuotientFull_vcPoint_eq`（在庫）で `C` を外へ出し、
   `semistableAt_variableChange`（在庫）で `C⁻¹` を掛けて戻す

★★★これで `semistableAt_veluQuotientFull` の**良い素点（`p ∤ l`）側が閉じた**。
☆残るのは**悪い素点**（Tate モデルの配管）と **`p ∣ l`**（形式群）である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GaloisRep ABC3.Skeleton.GenEll ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 良い素点（`p ∤ l`）では Vélu の商は半安定**——★（第 1403）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆仮定は `SemistableAt p E`・`0 ≤ jExp p E`（良い素点）・`l` は奇素数・`p ∤ l` だけである。

★★★第 1402 で判別式の恒等式が閉じたので、あとは極小モデルへの正規化だけだった。 -/
theorem semistableAt_veluQuot_goodPrime [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E) (hj : 0 ≤ jExp p E)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : primeSubring p)))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  have hΔ : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) E
  haveI := hC
  have hmd : minDeltaExp p E = 0 := by
    rw [minDeltaExp_eq_maxJ_of_semistable p E hss]
    exact max_eq_left (by omega)
  haveI hCell : (C • E).IsElliptic :=
    ⟨isUnit_iff_ne_zero.2 (variableChange_Delta_ne_zero E hΔ C)⟩
  have hΔC : valAdd p (Units.mk0 ((C • E).Δ) (C • E).isUnit_Δ.ne_zero) = 0 := by
    have h := minDeltaExp_eq p E hΔ C hC
    rw [hmd] at h
    rw [← h]
  have hQ' : addOrderOf (vcPoint C E Q) = l := by rw [addOrderOf_vcPoint, hQ]
  haveI hell' := isElliptic_veluQuotientFull_nsmul_nf' L (C • E) hQ'
  have hΔ' : (veluQuotientFull (C • E)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C E Q)))).Δ ≠ 0 :=
    (veluQuotientFull (C • E)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C E Q)))).isUnit_Δ.ne_zero
  have hgood := semistableAt_veluQuot_good p (C • E) hΔC hl hodd hlu (vcPoint C E Q) hQ' hΔ'
  have heq := veluQuotientFull_vcPoint_eq C E _ hQ two_ne_zero rfl
  rw [← heq] at hgood
  have hres := semistableAt_variableChange p _ C⁻¹ hgood
  rwa [inv_smul_smul] at hres

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuot_goodPrime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い素点 p ∤ l では Vélu の商は半安定。★半安定＋jExp ≥ 0＋l 奇素数＋p ∤ l)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_goodPrime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuot_good(第 1387＋第 1402、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.semistableAt_veluQuot_good") 1,
    .citation "[ABC3]" "minDeltaExp_eq_maxJ_of_semistable(第 738、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_maxJ_of_semistable") 1,
    .citation "[ABC3]" "veluQuotientFull_vcPoint_eq・semistableAt_variableChange(在庫)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1403）**——第 1402 で判別式の恒等式が閉じたので繋がった。" ++
       "☆あとは極小モデルへの正規化（3 行）だけだった。" ++
       "★★★これで `semistableAt_veluQuotientFull` の**良い素点（`p ∤ l`）側が閉じた**。" ++
       "残るのは**悪い素点**（Tate モデルの配管、第 1388 の分岐版）と" ++
       "**`p ∣ l`**（形式群）である。") 17 ]

end ABC3.Found.GenEll
