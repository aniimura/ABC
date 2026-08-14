import ABC3.Meta.Claim
import Mathlib.Topology.Compactness.Compact
import Mathlib.Data.Real.Basic

/-!
# [AbsTopIII] (L1) — log-shell はコンパクトで、ゆえに有限な log-volume を持つ

原典: S. Mochizuki, *Topics in Absolute Anabelian Geometry III: Global Reconstruction
Algorithms* [AbsTopIII]、物理 p.5(全 164 ページ。**400 dpi 目視確認済み**)。

原文 (AbsTopIII p.5):
> (L1) a log-shell is compact and hence of finite “log-volume” [cf. Corollary

## なぜここを張るか

`tools/check.mjs` の依存グラフで `[AbsTopIII] (L1)` が**未展開の葉**として現れた
(`Skeleton/IUTchIII/Cor312Claims.lean` の `LogShellPacketCompact` から出た辺)。

★そして**これは現時点で我々が着地させられる唯一の葉**である——
`Found/IUTchIII/PadicLog.lean` に ℚ_p の log-shell を `sorry` 無しで構成し、
コンパクト性(`isCompact_logShell'`)も非退化性(`logShell_ne_singleton_zero`)も
証明済みだから。

## ★何を写し、何を写していないか

写したのは主張の**形**である——「ある集合がコンパクトであり、かつその log-volume が
有限である」。log-shell という**対象そのもの**は引数として受ける。理由は2つ:

1. 一般の p 進局所体 `k` の `log_{k̄}(𝒪^×_k)` を我々は持っていない
   (作れているのは `k = ℚ_p` の場合だけ)。
2. `Skeleton` は `Found` を import しない(2トラック構成)。実物との接続は
   `.needs` の `.inProject` として記録し、`Check/` 側で確かめる。
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
  { paper := "AbsTopIII", pdfPage := 5, item := "(L1)",
    sectionId := "l1-log-shell-compact" }

/-- (L1) が要求するもの。

★1件目は**この計画で初めての着地**である——外部の補題ではなく、
我々が `Found/` に作った実物を指している。`inProject` の本来の用途は
「mathlib 外の公開プロジェクト」だが、ここでは**本プロジェクト内**を指している。
移植は不要なので判断は軽いが、用途を広げていることは明示しておく。 -/
def LogShellCompactAndFinite.needs : List ProofObligation :=
  [ .citation "ABC3(本プロジェクト、Track B)"
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
