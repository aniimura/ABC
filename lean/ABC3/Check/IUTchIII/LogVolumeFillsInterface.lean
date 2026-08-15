import ABC3.Interface.IUTchIII.PilotObjects
import ABC3.Found.IUTchIII.LogVolume
import ABC3.Skeleton.AbsTopIII.LogShell

/-!
# ★`Interface` の posit が初めて `Found/` の実物に置き換わる

`Interface/IUTchIII/PilotObjects.lean` の

```
logVolCompact : ∀ U : Set Amb, @IsCompact Amb topology U → ℝ
```

は、`Amb := ℚ_[p]` について **`Found/IUTchIII/LogVolume.lean` の実物で埋まる**。

★`PilotObjectData` **全体**を埋めてはいない。他のフィールド(Θ/q パイロット対象・
正則包・可能な像)には中身が無いことは 2026-08-15 に確認済みで、
埋めれば `Check/IUTchIII/Cor312Degenerate.lean` の退化 witness と同じものになる。
ここで示すのは「**このフィールドは実物で埋まる**」ことだけである。

## ★★見つかった不整合(この作業の副産物)

原文 (AbsTopIII p.137, Proposition 5.7, (i)):
> for the maximal ideal of Ok and M(k) for the set of compact open subsets

原文の `μ_k : 𝕄(k) → ℝ_{>0}` は **コンパクト開**集合の上でしか定義されておらず、
値域は **正の実数**である。ところが `Interface` のフィールドは
**すべてのコンパクト集合**に実数値を要求している。

コンパクトだが開でない集合(例: `∅`)は体積 0 で、対数体積は本来 `−∞` である。
Lean の `Real.log 0` は **0**(ゴミ値)なので型は付くが、それは原文の値ではない。
これを `padicLogVolCompact_empty_is_junk` で**実際に示す**。

★すなわち `Interface` のこのフィールドは**原文より広い定義域を要求している**。
`tools/check.mjs` 冒頭 A6(仮説の強化)の鏡像で、**実装側への要求の強化**にあたる。
-/

namespace ABC3.Check.IUTchIII

open ABC3.Interface.IUTchIII ABC3.Found.IUTchIII ABC3.Skeleton.AbsTopIII MeasureTheory Metric

variable {p : ℕ} [Fact p.Prime]

/-- ★**`PilotObjectData.logVolCompact` の型が ℚ_p の実物で埋まる**。 -/
noncomputable def padicLogVolCompact (p : ℕ) [Fact p.Prime] :
    ∀ U : Set ℚ_[p], @IsCompact ℚ_[p] inferInstance U → ℝ :=
  fun U _ => padicLogVol U

/-- 埋めたものは**原文の正規化**を満たす —— `𝒪_k` の対数体積は 0
(体積が 1 だから)。

原文 (AbsTopIII p.137):
> ization, i.e., μk(Ok) = 1. We shall refer to μk(−) as the volume on k.
-/
theorem padicLogVolCompact_integerBall :
    padicLogVolCompact p (closedBall (0 : ℚ_[p]) 1) (isCompact_closedBall 0 1) = 0 :=
  padicLogVol_integerBall

/-- ★**埋めたものは非退化である** —— 単位球と `1/p` 球で値が異なる。

これが無ければ「実物で埋めた」とは言えない(`fun _ _ => 0` でも型は付く)。 -/
theorem padicLogVolCompact_nondegenerate :
    padicLogVolCompact p (closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹))
        (isCompact_closedBall 0 _)
      ≠ padicLogVolCompact p (closedBall (0 : ℚ_[p]) 1) (isCompact_closedBall 0 1) :=
  padicLogVol_smallBall_ne_integerBall

/-- ★★**`Interface` のフィールドは原文より広い定義域を要求している**。

`∅` はコンパクトだが開ではなく、体積 0、対数体積は本来 `−∞`。
Lean の `Real.log 0` は 0 なので型は付くが、それは原文の値ではない
——原文の `μ_k` は `ℝ_{>0}` に値を取る。

**これは「型が付くが何も言っていない」の実例であり、我々の `Interface` の側の欠陥である。** -/
theorem padicLogVolCompact_empty_is_junk :
    padicLogVolCompact p (∅ : Set ℚ_[p]) isCompact_empty = 0 := by
  simp [padicLogVolCompact, padicLogVol]

/-- 原文の定義域(コンパクト**開**かつ空でない)では、体積は正で有限であり、
対数体積は本物の実数である。 -/
theorem padicLogVolCompact_genuine_on_compactOpen {U : Set ℚ_[p]}
    (hc : IsCompact U) (ho : IsOpen U) (hne : U.Nonempty) :
    0 < padicVolume p U ∧ padicVolume p U < ⊤ :=
  padicVolume_pos_lt_top_of_isCompact_isOpen hc ho hne

/-- ★[AbsTopIII] (L1) の未着地部分が1つ埋まった —— `Found/` の log-shell は
**有限体積**を持つ。

原文 (AbsTopIII p.5):
> (L1) a log-shell is compact and hence of finite “log-volume” [cf. Corollary

★まだ埋まっていないのは「対数体積が実数である」の側で、それには体積が**正**である
ことが要る(= `logShell` が開であること)。我々はそれを示していない。 -/
theorem padicLogShell_finite_volume : padicVolume p (logShell (p := p)) < ⊤ :=
  padicVolume_logShell_lt_top

/-- ★★**2件目の着地** —— `Found/` の `padicVolume` は
[AbsTopIII] Proposition 5.7 が要求する条件(平行移動不変性・正規化)を満たす。

原文 (AbsTopIII p.137):
> invariance, i.e., μk(A + x) = μk(A), for A ∈M(k), x ∈k; (3) normal-

★これで依存グラフに **Corollary 3.12 から実物まで届く鎖**が1本できた:
`cor_3_12` → `[IUTchIII] Theorem 3.11` → `[AbsTopIII] Corollary 5.10` →
`[AbsTopIII] Proposition 5.7` → `Found/IUTchIII/LogVolume.lean`(実物)。 -/
theorem padicVolume_isLocalVolume :
    IsLocalVolume (padicVolume p) (closedBall (0 : ℚ_[p]) 1) :=
  ⟨fun x s => measure_preimage_add _ _ _, padicVolume_integerBall⟩

end ABC3.Check.IUTchIII
