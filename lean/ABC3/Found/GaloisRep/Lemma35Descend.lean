/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma35StableLine
import ABC3.Found.GaloisRep.JExpCoprime
import ABC3.Meta.Claim

/-!
# 第 1221 ブロック —— **`Lemma 3.5` の仮説を `L` の側だけにする**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か

第 1199（`lemma_3_5_height_ineq_over_extension`）は `L''` の側で

* `∀ P, SemistableAt P (E ⊗ L'')`
* `∀ P, jExp P (E ⊗ L'') < 0 → ¬ (l ∣ jExp P (E ⊗ L''))`

を要求していた。☆第 1220 でこの 2 つが `L` の側の仮定

* `∀ p, SemistableAt p E`
* `∀ p, ¬ (l ∣ jExp p E)`（＝ `PrimeToLocalHeights`）
* `[L'':L] < l`

から出ることが分かったので、それを取り込む。

★★★これで `Lemma 3.5` の仮説は
**`L` の側のもの ＋ `L''` の上の点と商だけ**になった。
☆残るのは Vélu の商が楕円かつ半安定であること（原文が「同種なので自動」と
括弧で述べる段。`Lemma35Unconditional` の逸脱表と同じ扱い）である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta
open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 仮説は `L` の側だけ**（第 1221）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`L''` の側の半安定性と素性は、`L` の側の同じ仮定と `[L'':L] < l` から出る
（第 1220）。

★★★これが `Lemma 3.5` を安定直線の側で述べるための**受け口**である
——第 1219 が `[L'':L] ≤ l−1` の `L''` と点 `Q'` を与える。 -/
theorem lemma_3_5_height_ineq_descend (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic]
      (l : ℕ), l.Prime →
      (∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E) →
      (∀ p : HeightOneSpectrum (𝓞 L), ¬ ((l : ℤ) ∣ jExp p E)) →
      ∀ (L'' : Type) [Field L''] [NumberField L''] [Algebra L L''] [IsScalarTower ℚ L L'']
        [Algebra (𝓞 L) L''] [IsScalarTower (𝓞 L) L L'']
        [IsScalarTower (𝓞 L) (𝓞 L'') L''] [Module.Finite (𝓞 L) (𝓞 L'')]
        [Algebra.IsIntegral (𝓞 L) (𝓞 L'')],
      Module.finrank L L'' < l →
      ∀ (E'' : WeierstrassCurve L''),
      ∀ _hE0 : (E.baseChange L'').IsElliptic, ∀ _hE1 : E''.IsElliptic,
      ∀ Q : (E.baseChange L'').toAffine.Point, addOrderOf Q = l →
      E'' = veluQuotientFull (E.baseChange L'')
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) →
      (∀ P : HeightOneSpectrum (𝓞 L''), SemistableAt P E'') →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_height_ineq_over_extension eps heps
  refine ⟨C, ?_⟩
  intro L _ _ E _ l hl hss hcop L'' _ _ _ _ _ _ _ _ _ hdeg E'' hE0 hE1 Q hQ hE' hssE'
  haveI := hE0
  haveI := hE1
  exact hC L E l hl hss L'' E'' hE0 hE1 Q hQ hE'
    (semistableAt_baseChange_all L L'' E hss) hssE'
    (fun P _ => not_dvd_jExp_baseChange_all L L'' E hl hdeg hcop P)

/-! ## ★出典の紐付け(`.src`) -/

def lemma_3_5_height_ineq_descend.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(仮説は L の側だけ——L'' の側の半安定性と素性は降ろせる)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_height_ineq_descend.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_height_ineq_over_extension(第 1199、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_height_ineq_over_extension") 1,
    .citation "[ABC3]" "semistableAt_baseChange_all(第 1220、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_baseChange_all") 1,
    .citation "[ABC3]" "not_dvd_jExp_baseChange_all(第 1220、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.not_dvd_jExp_baseChange_all") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1221）**——`Lemma 3.5` の仮説が" ++
       "**`L` の側のもの ＋ `L''` の上の点と商だけ**になった。" ++
       "☆残るのは Vélu の商が楕円かつ半安定であること" ++
       "（原文が「同種なので自動」と括弧で述べる段。" ++
       "`Lemma35Unconditional` の逸脱表と同じ扱い）である。") 2 ]

end ABC3.Found.GaloisRep
