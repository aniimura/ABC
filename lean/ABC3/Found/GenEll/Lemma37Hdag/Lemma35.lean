/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma35Assemble
import ABC3.Found.GenEll.LCyclicPoint
import ABC3.Found.GenEll.Lemma37StableLineCop
import ABC3.Meta.Claim

/-!
# Lemma37Hdag —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve
open IsDedekindDomain NumberField Finset ABC3.Meta
open scoped Classical

/-- ☆Vélu の商が楕円かつ半安定であること（原文が「同種なので自動」と括弧で述べる段）。

★★★★★★★★☆**2026-09-02（第 1322-1323）——悪い素点側の見通し**

☆半安定性は 2 つに分かれる:

* **悪い素点**（`jExp < 0`）——★**核は取れた**。
  `isUnit_c4_velu_tate`（第 1323）が「Tate 曲線の Vélu の商は `c₄` が単元」を、
  `semistableAt_of_c4_valAdd_zero`（第 1322）が「`c₄` が単元の整モデルは半安定」を与える。
  ☆残るのは商のモデルを Tate モデルに移す変数変換
  （`veluQuotientFull_vcPoint_eq` 第 969・`vAdd_tateModel_u_eq_zero` 第 1056、在庫）と
  整性（`veluIntegralClosed`、在庫）の配管である。
* **良い素点**（`jExp ≥ 0`）——☆**同種で良還元が保たれる**（Néron–Ogg–Shafarevich）が要る。

★これと Vélu の商の楕円性（Vélu の定理）が `VeluQuotOK` に残る 2 つの既知数学である。

★★☆**2026-09-02（第 1327-1328）——悪い素点側は閉じ、残りは 1 本に絞れた**

☆第 1327（`semistableAt_velu_of_veluCurve_eq`）で**悪い素点の半安定性は閉じた**。

★良い素点側は `SemistableAt` の第 1 の選択肢（`minDeltaExp = 0`＝商も良還元）が要るが、
その道筋は
(1) `E` の整モデルは `v(Δ) = 0`、
(2) `l`-捩れの点は整で還元は単射（`l ≠ char k`）、
(3) Vélu の `v`・`w` は整なので商の還元は剰余体の上の Vélu の商、
(4) 剰余体の上で商が楕円なら `Δ` は単元、
であり、**根は「Vélu の定理（商の `Δ ≠ 0`）」1 本**である。

★★★したがって `VeluQuotOK` の 2 つの穴（楕円性・良い素点の半安定性）は
**どちらも Vélu の定理に落ちる**

★★★★★★★★★★★★★★★★**2026-09-02（第 1330-1336）——楕円性の穴は閉じた**

☆解析側の在庫（`exists_velu_rep`・`latticeCurve_eq_veluQuotientFull`）が**Vélu の定理そのもの**を持っていた。
★第 1333（`ℂ` の格子曲線）→ 第 1334（`ℂ` への埋め込みがある体）→ 第 1335（数体）と降り、第 1336 の `veluQuotOK_of_semistable` で
**`VeluQuotOK` は半安定性だけに落ちた**。

★★★残る既知数学は**良い素点の半安定性 1 本**である。 -/
def VeluQuotOK (E : SSCurve) (l : ℕ) : Prop :=
  ∀ (M : IntermediateField E.fld E.alg) [FiniteDimensional E.fld M]
     (Q' : (E.W.baseChange M).toAffine.Point),
     letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
     letI : NumberField (M : Type) := NumberField.of_module_finite E.fld M
     letI : IsScalarTower (𝓞 E.fld) E.fld M := isScalarTower_ringOfIntegers_base E.fld M
     letI : IsScalarTower (𝓞 E.fld) (𝓞 (M : Type)) M := isScalarTower_ringOfIntegers_top E.fld M
     addOrderOf Q' = l →
     (veluQuotientFull (E.W.baseChange M)
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q')))).IsElliptic ∧
     ∀ P : HeightOneSpectrum (𝓞 (M : Type)),
       SemistableAt P (veluQuotientFull (E.W.baseChange M)
         (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'))))

def VeluQuotOK.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商が楕円かつ半安定であること——原文が括弧で述べる段)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
