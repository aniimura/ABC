/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateCurveWitness
import ABC3.Meta.Claim

/-!
# 第 1315 ブロック —— **`E` と Tate モデルを結ぶ変数変換（`K` の上で）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か——残り (a) の第 1 段

`exists_tate_model`（在庫）は `R` の上の変数変換 `C` で
`C • integralModel R W = tateCurveAt q` を与える。
★これを `K` の上に降ろすと

> `(C.map (algebraMap R K)) • W = (tateCurveAt q).map (algebraMap R K)`

になる（`baseChange_integralModel_eq`）。
☆この形が第 1313（局所の 2 入力を変数変換で戻す段）が受け取るものである。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Found.GaloisRep ABC3.Meta

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★★★★★★★★★
**`E` と Tate モデルを結ぶ変数変換（`K` の上で）**——★**無条件**（第 1315）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

☆在庫の `tateParamR_spec` を `K` へ降ろしただけである。 -/
theorem exists_vc_tateModel (W : WeierstrassCurve K) [hell : W.IsElliptic]
    [WeierstrassCurve.IsMinimal R W] (h : W.HasSplitMultiplicativeReduction R) :
    ∃ C' : VariableChange K,
      C' • W = (tateCurveAt (tateParamR W h) (tateParamR_mem W h)).map (algebraMap R K) := by
  obtain ⟨hq, C, hne, hCE⟩ := tateParamR_spec W h
  refine ⟨C.map (algebraMap R K), ?_⟩
  have hWmap : (WeierstrassCurve.integralModel R W).map (algebraMap R K) = W :=
    WeierstrassCurve.baseChange_integralModel_eq R W
  conv_lhs => rw [← hWmap]
  rw [WeierstrassCurve.map_variableChange, hCE]

/-! ## ★出典の紐付け(`.src`) -/

def exists_vc_tateModel.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(E と Tate モデルを結ぶ変数変換——K の上で。★無条件)",
    sectionId := "genell-def-3-3" }

def exists_vc_tateModel.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateParamR_spec(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateParamR_spec") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1315）**——残り (a) の第 1 段である。" ++
       "☆この形が第 1313（局所の 2 入力を変数変換で戻す段）が受け取るものである。") 2 ]

end ABC3.Found.GenEll
