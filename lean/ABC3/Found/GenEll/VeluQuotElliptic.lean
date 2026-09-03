/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Lemma37Hdag
import ABC3.Found.GenEll.VeluEllipticNF
import ABC3.Meta.Claim

/-!
# 第 1336 ブロック —— **`VeluQuotOK` の楕円性の穴が閉じる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か

`VeluQuotOK`（第 1225）は **2 つ**を要求していた:

| 穴 | 内容 | 状態 |
|---|---|---|
| 楕円性 | `veluQuotientFull (E ⊗ M) ⟨Q'⟩` が楕円 | ★**本ブロックで閉じる**（第 1335） |
| 半安定性 | 各素点で半安定 | ☆悪い素点は第 1327 で閉じた。良い素点が残る |

★★★したがって `Lemma 3.7`（安定直線の側）が外部に受けている既知数学は
**「Vélu の商の半安定性」ただ 1 つ**になった。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve
open IsDedekindDomain NumberField Finset ABC3.Meta
open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`VeluQuotOK` は半安定性だけに落ちる**——★**無条件**（第 1336）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆楕円性は第 1335（`isElliptic_veluQuotientFull_nsmul_nf'`）が無条件で与える。 -/
theorem veluQuotOK_of_semistable (E : SSCurve) (l : ℕ)
    (hss : ∀ (M : IntermediateField E.fld E.alg) [FiniteDimensional E.fld M]
      (Q' : (E.W.baseChange M).toAffine.Point),
      letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
      letI : NumberField (M : Type) := NumberField.of_module_finite E.fld M
      letI : IsScalarTower (𝓞 E.fld) E.fld M := isScalarTower_ringOfIntegers_base E.fld M
      letI : IsScalarTower (𝓞 E.fld) (𝓞 (M : Type)) M := isScalarTower_ringOfIntegers_top E.fld M
      addOrderOf Q' = l →
      ∀ P : HeightOneSpectrum (𝓞 (M : Type)),
        SemistableAt P (veluQuotientFull (E.W.baseChange M)
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'))))) :
    VeluQuotOK E l := by
  intro M _ Q' hQ'
  letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
  letI : NumberField (M : Type) := NumberField.of_module_finite E.fld M
  haveI hEell : (E.W.baseChange (M : Type)).IsElliptic := by
    show (E.W.map (algebraMap E.fld M)).IsElliptic
    infer_instance
  exact ⟨isElliptic_veluQuotientFull_nsmul_nf' (M : Type) _ hQ', hss M Q' hQ'⟩

/-! ## ★出典の紐付け(`.src`) -/

def veluQuotOK_of_semistable.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(VeluQuotOK は半安定性だけに落ちる——楕円性は無条件。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotOK_of_semistable.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isElliptic_veluQuotientFull_nsmul_nf'(第 1335、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isElliptic_veluQuotientFull_nsmul_nf'") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1336）**——`Lemma 3.7`（安定直線の側）が外部に受けている" ++
       "既知数学は**「Vélu の商の半安定性」ただ 1 つ**になった。" ++
       "☆その悪い素点側は第 1327 で閉じている。") 3 ]

end ABC3.Found.GenEll
