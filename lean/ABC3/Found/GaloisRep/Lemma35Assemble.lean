/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TowerInstances
import ABC3.Found.GenEll.DegreeLine
import ABC3.Meta.Claim

/-!
# 第 1223 ブロック —— **`Lemma 3.5` を安定直線の側で述べる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——`L(H)` の道の完成

| 段 | 定理 | 第 |
|---|---|---|
| 安定直線 → 位数 `l` の点 | `exists_stablePoint_of_hasLCyclicJ` | 1205 |
| 点を次数 `≤ l−1` の `M` へ降ろす | `exists_point_descent_of_stable` | 1219 |
| インスタンス | `numberField_and_towers` | 1222 |
| 仮説を `L` の側だけに | `lemma_3_5_height_ineq_descend` | 1221 |
| 有限拡大で確かめれば足りる | `lemma_3_5_height_ineq_over_extension` | 1199 |

★★★本ブロックはこれらを 1 本に繋ぐ。

## ★★逸脱の記録

☆Vélu の商が**楕円かつ半安定**であることは仮説として受ける
（原文は「同種なので自動」と括弧で述べる。
`Found/GaloisRep/Lemma35Unconditional.lean` の逸脱表と同じ扱い）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Interface.GaloisRep ABC3.Meta
open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 安定直線の側**（第 1223）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`E[l]` の中に `Gal`-安定な直線があれば、その生成元は次数 `≤ l−1` の
拡大 `M` で有理になり（第 1219）、`Lemma 3.5` の不等式が `L` の上で出る
（第 1199・第 1221）。

★★逸脱: Vélu の商が楕円かつ半安定であることは仮説として受ける
（原文は「同種なので自動」と括弧で述べる）。 -/
theorem lemma_3_5_height_ineq_stableLine (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic]
      [(E.baseChange (AlgebraicClosure L)).IsElliptic] (l : ℕ), l.Prime →
      (∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E) →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      ∀ Q : (E.baseChange (AlgebraicClosure L)).toAffine.Point, addOrderOf Q = l →
      (∀ σ : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L,
        ∃ k : ℕ, galPoint E σ Q = k • Q) →
      (∀ (M : IntermediateField L (AlgebraicClosure L)) [FiniteDimensional L M]
         (Q' : (E.baseChange M).toAffine.Point),
         letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
         letI : NumberField (M : Type) := NumberField.of_module_finite L M
         letI : IsScalarTower (𝓞 L) L M := isScalarTower_ringOfIntegers_base L M
         letI : IsScalarTower (𝓞 L) (𝓞 (M : Type)) M := isScalarTower_ringOfIntegers_top L M
         addOrderOf Q' = l →
         (veluQuotientFull (E.baseChange M)
            (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q')))).IsElliptic ∧
         ∀ P : HeightOneSpectrum (𝓞 (M : Type)),
           SemistableAt P (veluQuotientFull (E.baseChange M)
             (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'))))) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_height_ineq_descend eps heps
  refine ⟨C, ?_⟩
  intro L _ _ E _ _ l hl hss hcop Q hQ hst hquot
  obtain ⟨M, hfin, hdeg, Q', hord, _himg⟩ := exists_point_descent_of_stable E hl Q hQ hst
  haveI := hfin
  letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
  letI : NumberField (M : Type) := NumberField.of_module_finite L M
  letI : IsScalarTower (𝓞 L) L M := isScalarTower_ringOfIntegers_base L M
  letI : IsScalarTower (𝓞 L) (𝓞 (M : Type)) M := isScalarTower_ringOfIntegers_top L M
  haveI hE0 : (E.baseChange (M : Type)).IsElliptic := by
    show (E.map (algebraMap L M)).IsElliptic
    infer_instance
  obtain ⟨hell, hssq⟩ := hquot M Q' hord
  have hl2 : 2 ≤ l := hl.two_le
  have hdlt : Module.finrank L (M : Type) < l := by omega
  exact hC L E l hl hss hcop (M : Type) hdlt _ hE0 hell Q' hord rfl hssq

/-! ## ★出典の紐付け(`.src`) -/

def lemma_3_5_height_ineq_stableLine.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(安定直線の側——商の楕円性と半安定性は仮説として受ける)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_height_ineq_stableLine.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_point_descent_of_stable(第 1219、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_point_descent_of_stable") 1,
    .citation "[ABC3]" "lemma_3_5_height_ineq_descend(第 1221、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_height_ineq_descend") 1,
    .citation "[ABC3]" "numberField_and_towers(第 1222、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.numberField_and_towers") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1223）**——`L(H)` の道の**完成**である。" ++
       "☆逸脱: Vélu の商が楕円かつ半安定であることは仮説として受ける" ++
       "（原文は「同種なので自動」と括弧で述べる。" ++
       "`Lemma35Unconditional` の逸脱表と同じ扱い）。") 3 ]

end ABC3.Found.GaloisRep
