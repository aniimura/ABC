/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VAddC4Velu
import ABC3.Found.GaloisRep.SemistableFromC4
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Meta.Claim

/-!
# 第 1327 ブロック —— **Vélu の商は悪い素点で半安定**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——`VeluQuotOK` の**悪い素点側**

☆5 段を繋ぐ:

| 段 | 内容 | 第 |
|---|---|---|
| 1 | Tate 曲線の Vélu の商は `c₄` が単元 | 1323 |
| 2 | 曲線の等式から `v(c₄(X)) = 0` | 1326 |
| 3 | `vAdd`（局所）→ `valAdd`（大域） | `vAdd_algebraMap_eq_valAdd`（在庫） |
| 4 | `c₄` が単元の整モデルは半安定 | 1322 |

★★★これで**「原文が『同種なので自動』と括弧で述べた段」の悪い素点側が閉じる**。
☆残るのは良い素点側（同種で良還元が保たれる＝Néron–Ogg–Shafarevich）である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll ABC3.Meta

variable {L : Type} [Field L] [NumberField L]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**Vélu の商は悪い素点で半安定**——★**無条件**（第 1327）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`E′ ⊗ L_v` が Tate 曲線の Vélu の商と変数変換（`u` の付値 `0`）で結ばれていれば、
`c₄` の付値が `0` になり、整性と合わせて半安定である。 -/
theorem semistableAt_velu_of_veluCurve_eq {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E' : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) E']
    (hΔ : E'.Δ ≠ 0) (hc4 : E'.c₄ ≠ 0)
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (v w : R)
    (hunit : IsUnit ((tateCurveAt q hq).c₄ + 240 * v))
    (C₀ : WeierstrassCurve.VariableChange R)
    (hEq : (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv)
      = (veluCurve (tateCurveAt q hq) v w).map (algebraMap R Lv))
    (hu : vAdd (tateDvrVal R Lv) ((C₀.map (algebraMap R Lv)).u) = 0) :
    SemistableAt p E' := by
  have hinj : Function.Injective (algebraMap L Lv) := (algebraMap L Lv).injective
  have hc4' : (E'.baseChange Lv).c₄ ≠ 0 := by
    rw [WeierstrassCurve.baseChange, WeierstrassCurve.map_c₄]
    exact (map_ne_zero_iff _ hinj).2 hc4
  have hloc := vAdd_c4_of_veluCurve_eq (R := R) (K := Lv) q hq v w hunit C₀
    (E'.baseChange Lv) hc4' hEq hu
  -- 局所の付値を大域の `valAdd` に直す
  have hne2 : algebraMap L Lv E'.c₄ ≠ 0 := (map_ne_zero_iff _ hinj).2 hc4
  have hEqU : (Units.mk0 ((E'.baseChange Lv).c₄) hc4')
      = Units.mk0 (algebraMap L Lv E'.c₄) hne2 := by
    refine Units.ext ?_
    show (E'.baseChange Lv).c₄ = algebraMap L Lv E'.c₄
    rw [WeierstrassCurve.baseChange, WeierstrassCurve.map_c₄]
  rw [hEqU, vAdd_algebraMap_eq_valAdd (R := R) p hp E'.c₄ hc4 hne2] at hloc
  exact semistableAt_of_c4_valAdd_zero p E' hΔ hc4 hloc

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_velu_of_veluCurve_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商は悪い素点で半安定。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_velu_of_veluCurve_eq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "vAdd_c4_of_veluCurve_eq(第 1326、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_c4_of_veluCurve_eq") 1,
    .citation "[ABC3]" "semistableAt_of_c4_valAdd_zero(第 1322、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_of_c4_valAdd_zero") 1,
    .citation "[ABC3]" "vAdd_algebraMap_eq_valAdd(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_algebraMap_eq_valAdd") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1327）**——`VeluQuotOK` の**悪い素点側**である。" ++
       "☆原文が「同種なので自動」と括弧で述べた段のうち、悪い素点の半安定性がこれで閉じる。" ++
       "★残るのは良い素点側（同種で良還元が保たれる＝Néron–Ogg–Shafarevich）と" ++
       "Vélu の商の楕円性である。") 3 ]

end ABC3.Found.GaloisRep
