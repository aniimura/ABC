/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.Topology.Compactness.Compact
import Mathlib.Data.Real.Basic
import Mathlib.MeasureTheory.Measure.MeasureSpace

/-!
# LogShell —— `[AbsTopIII] Corollary 5.10` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Skeleton.AbsTopIII

open ABC3.Meta

/-- **(L1)** の転写。`shell` がコンパクトであり、かつその log-volume が有限である。

原文 (AbsTopIII p.5):
> (L1) a log-shell is compact and hence of finite “log-volume” [cf. Corollary

★`shell` と `logVol` を引数として受けるのは、log-shell という対象を我々が
一般には持っていないためである(上の docstring)。 -/
def LogShellCompactAndFinite {K : Type} (τ : TopologicalSpace K)
    (shell : Set K) (logVol : Set K → WithTop ℝ) : Prop :=
  @IsCompact K τ shell ∧ logVol shell ≠ ⊤

def LogShellCompactAndFinite.src : Source :=
  { paper := "AbsTopIII", pdfPage := 145, item := "Corollary 5.10, (i)",
    sectionId := "cor-5-10-log-shell-properties" }

/-- (L1) が要求するもの。

★1件目は**この計画で初めての着地**である——外部の補題ではなく、
我々が `Found/` に作った実物を指している。`inProject` の本来の用途は
「mathlib 外の公開プロジェクト」だが、ここでは**本プロジェクト内**を指している。
移植は不要なので判断は軽いが、用途を広げていることは明示しておく。 -/
def LogShellCompactAndFinite.needs : List ProofObligation :=
  [ .otherPaper "[AbsTopIII]" "Proposition 5.7, (i), (ii)(局所体積・対数体積)" 137,
    .citation "ABC3(本プロジェクト、Track B)"
      "ℚ_p の log-shell の構成とそのコンパクト性"
      (.inProject "ABC3"
        "Found/IUTchIII/PadicLog.lean の logShell / isCompact_logShell'(sorry 無し、標準3公理)。非退化性は logOneAdd_ne_zero と logShell_ne_singleton_zero、加法で閉じることは Found/IUTchIII/PadicLogMul.lean の logShell_add_mem")
      5,
    .implicitStep
      "★原文 (L1) は一般の p 進局所体 k の log_{k̄}(𝒪^×_k) について述べる。我々が構成したのは k = ℚ_p の場合だけであり、しかも 𝒪^× ではなく 1 + 𝔪 の像である(𝒪^× ≅ μ × (1 + 𝔪) は mathlib に無いと 2026-08-15 に実測)"
      5,
    .implicitStep
      "「hence of finite log-volume」の hence は、対数体積が 𝔐(−)(空でないコンパクト開集合の集まり)の上で ℝ 値に定まること([IUTchIII] Proposition 3.9, (i))を使う。我々の Interface ではこれを logVolCompact の型で受けている"
      5 ]

end ABC3.Skeleton.AbsTopIII
