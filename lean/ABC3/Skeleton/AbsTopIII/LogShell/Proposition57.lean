/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.Topology.Compactness.Compact
import Mathlib.Data.Real.Basic
import Mathlib.MeasureTheory.Measure.MeasureSpace
import ABC3.Skeleton.AbsTopIII.LogShell.Corollary510

/-!
# LogShell —— `[AbsTopIII] Proposition 5.7` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Skeleton.AbsTopIII

open ABC3.Meta

/-! ### [AbsTopIII] Proposition 5.7, (i), (a) —— 対数体積の定義 -/

/-- **局所体積**(原文の `μ_k`)が満たすべき3条件。

原文 (AbsTopIII p.137):
> The following result is elementary and well-known.

原文 (AbsTopIII p.137):
> that satisfies the following properties: (1) additivity, i.e., μk(A

原文 (AbsTopIII p.137):
> invariance, i.e., μk(A + x) = μk(A), for A ∈M(k), x ∈k; (3) normal-

原文 (AbsTopIII p.137):
> ization, i.e., μk(Ok) = 1. We shall refer to μk(−) as the volume on k.

★原文はこの3条件が `μ_k` を**一意に**定めると述べる。ここでは3条件を述語として写す
(一意性は写していない)。`vol` は `ℝ≥0∞` 値の測度として受ける——原文の
`𝕄(k) → ℝ_{>0}` は**コンパクト開**集合の上での話で、そこでは正値性が自動だからである。 -/
def IsLocalVolume {K : Type} [MeasurableSpace K] [Add K]
    (vol : MeasureTheory.Measure K) (integers : Set K) : Prop :=
  (∀ (x : K) (s : Set K), vol ((fun y => x + y) ⁻¹' s) = vol s) ∧ vol integers = 1

def IsLocalVolume.src : Source :=
  { paper := "AbsTopIII", pdfPage := 137, item := "Proposition 5.7, (i), (a)",
    sectionId := "prop-5-7-i-a-local-volume" }

/-- ★**辺は 0 本**。原文 p.137 は直前に「The following result is elementary and
well-known.」と書いており、この命題は他の番号付き項目に依拠していない。
**依存グラフの終端**である。 -/
def IsLocalVolume.needs : List ProofObligation :=
  [ .citation "ABC3(本プロジェクト、Track B)"
      "ℚ_p の体積と対数体積の構成(𝒪 の体積を 1 に正規化した Haar 測度)"
      (.inProject "ABC3"
        "Found/IUTchIII/LogVolume.lean の padicVolume / padicLogVol(sorry 無し、標準3公理)。正規化は padicVolume_integerBall、非退化性は padicLogVol_smallBall_ne_integerBall")
      137,
    .implicitStep
      "原文は3条件が `μ_k` を**一意に**定めると述べるが、その一意性を我々は写していない(存在する実物が原文の `μ_k` と一致することは示していない)" 137,
    .implicitStep
      "原文は加法性を `𝕄(k)`(コンパクト開集合)の上の有限加法性として述べる。我々は mathlib の測度(可算加法性)として受けており、**強い方**を取っている" 137,
    .implicitStep
      "原文はアルキメデス的な場合(complex archimedean field)も扱うが、我々は非アルキメデス側だけを写している" 137 ]

end ABC3.Skeleton.AbsTopIII
