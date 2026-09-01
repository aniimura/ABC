/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma37CondB
import ABC3.Found.GaloisRep.Lemma35Unconditional
import ABC3.Meta.Claim

/-!
# 第 1150 ブロック —— **`Lemma 3.7` の条件 (b) からも `(†)` が消えた**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★これは何か

第 1088 は**条件 (a)** の側で `(†)`（＝ `Lemma 3.5` の結論）を埋めた。
☆そのとき条件 (b) の側は埋められなかった——理由は 1 つだけで、

    条件 (b) が `l` について言うのは **`2 ≤ l`** だけであり、
    `Lemma 3.5` が置いていた `l ≠ 2` を吸収できない

からである（`Lemma35Unconditional.lean` の「残る `l ≠ 2` がどこまで吸収されるか」）。

★★★★**2026-09-01（第 1149）で `Lemma 3.5` から `l ≠ 2` が外れた**。
☆したがって条件 (b) でも `(†)` はそのまま埋まる。本ファイルがそれである。

## ★★★★条件 (a) との対比

| | `l` について仮定すること | `ht^Falt` の上界の出どころ |
|---|---|---|
| (a) | `l ≥ 100d·(ht^Falt + C·d^ϵ)`（★`l` が**大きい**） | `l` の大きさ |
| (b) | `2 ≤ l`（★それだけ） | `[E_L] ∈ K_V`、すなわち `h_∞(j) ≤ M` |

★どちらの側でも `ht^Falt` が `l` に依らない定数で抑えられ、
`Northcott` で類が有限個になる。☆これが `Lemma 3.7` の第 3 の主張の中身である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open ABC3.Found.GenEll ABC3.Meta
open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★
**条件 (b) だけで `ht^Falt ≤ C″`**（第 1150）——`(†)` は要らない。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★★★★**2026-09-01（第 1150）**——`lemma_3_5_height_ineq`（第 1149 で条なしに
なったもの）を `htFalt_le_of_condB` に代入した形である。
☆外から受けるのは `h_∞(j) ≤ M`（＝ `[E_L] ∈ K_V` の内容）と
半安定性・`l`-cyclic の情報だけで、**原典の `Lemma 3.7` が置くものだけ**である。 -/
theorem htFalt_le_of_condB_lcyclic :
    ∃ C₂ C' : ℝ, 0 ≤ C' ∧
      ∀ (M : ℝ) (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
        [E.IsElliptic] [E'.IsElliptic] (l : ℕ),
        l.Prime →
        (∀ r : HeightOneSpectrum (𝓞 L), SemistableAt r E) →
        (∀ r : HeightOneSpectrum (𝓞 L), SemistableAt r E') →
        htArchJ L E ≤ M →
        ∀ Q : E.toAffine.Point, addOrderOf Q = l →
        E' = veluQuotientFull E (((Finset.range l).erase 0).image
            (fun k : ℕ => pointCoords (k • Q))) →
        (∀ r : HeightOneSpectrum (𝓞 L), jExp r E < 0 → ¬ ((l : ℤ) ∣ jExp r E)) →
        htFaltOf L E ≤ |M + C₂| / 5 + 28 / 5 + 1.4 * C' := by
  obtain ⟨C₂, hC₂⟩ := htFalt_le_of_condB
  obtain ⟨C₅, h₅⟩ := lemma_3_5_height_ineq (1/6) (by norm_num)
  refine ⟨C₂, max C₅ 0, le_max_right _ _, ?_⟩
  intro M L _ _ E E' _ _ l hl hssE hssE' harch Q hQ hE' hcop
  -- ☆`Lemma 3.5` で `(†)` を埋める（★`l ≠ 2` は第 1149 で外れた）
  have hdag : ((l : ℝ) / 14) * degInfOf L E
      ≤ htFaltOf L E + 2 * Real.log l + max C₅ 0 := by
    have h := h₅ L E E' l hl Q hQ hE' hssE hssE' hcop
    have hcoef : (1 : ℝ) / (12 * (1 + 1/6)) = 1 / 14 := by norm_num
    rw [hcoef] at h
    have hrw : (1 : ℝ) / 14 * (l : ℝ) * degInfOf L E
        = ((l : ℝ) / 14) * degInfOf L E := by ring
    rw [hrw] at h
    exact le_trans h (by linarith [le_max_left C₅ (0:ℝ)])
  exact hC₂ M (max C₅ 0) (le_max_right _ _) L E hssE harch l hl.two_le hdag

/-! ## ★出典の紐付け(`.src`) -/

def htFalt_le_of_condB_lcyclic.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 3 の主張・条件 (b)——(†) を埋めた形。★外部入力なし)",
    sectionId := "genell-lemma-3-7" }

def htFalt_le_of_condB_lcyclic.needs : List ProofObligation :=
  [ .citation "[ABC3]" "htFalt_le_of_condB(数値の核、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_le_of_condB") 1,
    .citation "[ABC3]" "lemma_3_5_height_ineq((†)、第 1149 で条なし、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_height_ineq") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1150）の記録**——本定理は第 1088 の時点では " ++
       "**書けなかった**。☆条件 (b) が `l` について言うのは `2 ≤ l` だけで、" ++
       "`Lemma 3.5` が置いていた `l ≠ 2` を吸収できなかったからである。" ++
       "★第 1149 でその仮説が外れたので、条件 (a) と同じ 1 行で埋まった。") 1 ]

end ABC3.Found.GaloisRep
