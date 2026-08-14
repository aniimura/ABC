import ABC3.Skeleton.AbsTopIII.LogShell
import ABC3.Found.IUTchIII.PadicLog
import ABC3.Found.IUTchIII.PadicLogMul

/-!
# ★着地 — `Found/` の ℚ_p の log-shell が [AbsTopIII] (L1) を満たすこと

依存グラフの葉 `[AbsTopIII] (L1)` は、この計画で**初めて実物に着地する葉**である。
`Skeleton/AbsTopIII/LogShell.lean` が原文の主張の形を持ち、
`Found/IUTchIII/PadicLog.lean` がその実物(ℚ_p の場合)を持つ。
ここは**その2つを突き合わせる場所**である。

★`Skeleton` は `Found` を import しない(2トラック構成)ので、突き合わせは
`Check/` でしか書けない。`Check/` は「我々のモデルについての検査」の置き場である。

## ★★何が着地して、何が着地していないか

- **着地した**: コンパクト性。`isCompact_logShell'` は `sorry` 無し・標準3公理。
- **着地していない**: (a) 有限な log-volume——**対数体積を持っていない**ので、
  `logVol` とその有限性を仮定として受け取る形でしか書けない。
  (b) 一般の p 進局所体 `k`——作れているのは `k = ℚ_p` だけ。
  (c) `𝒪^×` ではなく `1 + 𝔪` の像である(`𝒪^× ≅ μ × (1 + 𝔪)` は
  mathlib に無いと 2026-08-15 に実測)。
-/

namespace ABC3.Check.AbsTopIII

open ABC3.Skeleton.AbsTopIII ABC3.Found.IUTchIII

variable {p : ℕ} [Fact p.Prime]

/-- ★**(L1) の前半が着地した** —— `Found/` の ℚ_p の log-shell はコンパクトである。

原文 (AbsTopIII p.5):
> (L1) a log-shell is compact and hence of finite “log-volume” [cf. Corollary
-/
theorem padicLogShell_isCompact : IsCompact (logShell (p := p)) :=
  isCompact_logShell'

/-- ★**(L1) は非自明に着地している** —— log-shell は `{0}` ではない。

コンパクト性だけなら `logShell := ∅` でも満たせるので、**非退化性が対で要る**。 -/
theorem padicLogShell_nondegenerate : logShell (p := p) ≠ {0} :=
  logShell_ne_singleton_zero

/-- ★log-shell は**加法で閉じる** —— [AbsTopIII] p.3 が log-shell を
「canonical rigid **integral structure**」と呼ぶときの structure の一部。 -/
theorem padicLogShell_add_closed {a b : ℚ_[p]}
    (ha : a ∈ logShell (p := p)) (hb : b ∈ logShell (p := p)) :
    a + b ∈ logShell (p := p) :=
  logShell_add_mem ha hb

/-- ★**(L1) 全体**。ただし後半は仮定として受け取っている——
我々は対数体積を持っていないので、これ以上は書けない。 -/
theorem padicLogShell_satisfies_L1 (logVol : Set ℚ_[p] → WithTop ℝ)
    (hfin : ∀ U : Set ℚ_[p], IsCompact U → logVol U ≠ ⊤) :
    LogShellCompactAndFinite (inferInstance : TopologicalSpace ℚ_[p])
      (logShell (p := p)) logVol :=
  ⟨isCompact_logShell', hfin _ isCompact_logShell'⟩

end ABC3.Check.AbsTopIII
