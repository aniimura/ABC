import Mathlib.MeasureTheory.Measure.Haar.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.NumberTheory.Padics.ProperSpace
import ABC3.Found.IUTchIII.PadicLog

/-!
# ℚ_p の**対数体積**(log-volume)——[AbsTopIII] Proposition 5.7 の転写

原典: S. Mochizuki, *Topics in Absolute Anabelian Geometry III* [AbsTopIII]、
物理 p.137(全 164 ページ。**400 dpi 目視確認 2026-08-15**)。

原文 (AbsTopIII p.137):
> The following result is elementary and well-known.

原文 (AbsTopIII p.137, Proposition 5.7):
> (Local Volumes) Let k be either a mixed-characteristic

原文 (AbsTopIII p.137, Proposition 5.7, (i)):
> for the maximal ideal of Ok and M(k) for the set of compact open subsets

★**定義は完全に具体的である**。原文が要求するのは次の3つだけ:

原文 (AbsTopIII p.137):
> that satisfies the following properties: (1) additivity, i.e., μk(A

原文 (AbsTopIII p.137):
> invariance, i.e., μk(A + x) = μk(A), for A ∈M(k), x ∈k; (3) normal-

原文 (AbsTopIII p.137):
> ization, i.e., μk(Ok) = 1. We shall refer to μk(−) as the volume on k.

そして対数体積は音量の自然対数である:

原文 (AbsTopIII p.137):
> logarithm R>0 →R] and refer to μlog

## ★この3条件は mathlib の **Haar 測度** そのものである

(1) 加法性 = 測度、(2) 平行移動不変 = 左不変、(3) 正規化 `μ(𝒪_k) = 1` =
`addHaarMeasure K₀` の `K₀ := 𝒪_k` における正規化(`addHaarMeasure_self`)。
すなわち**原文の μ_k は「𝒪_k の体積を 1 に正規化した Haar 測度」**である。

## ★★測定(2026-08-15、S1-S4)

- **(S1) 原文側の呼称**: “volume” / “log-volume” / “compact open subsets” /
  “additivity” / “translation invariance” / “normalization μ_k(𝒪_k) = 1”。
- **(S2) 全出現の列挙**: `ls Mathlib/MeasureTheory/Measure/Haar/` は 10 ファイル
  (Basic, Disintegration, DistribChar, Extension, InnerProductSpace,
  MulEquivHaarChar, NormedSpace, OfBasis, Quotient, Unique)。
  `Measure.addHaarMeasure` と `Measure.addHaarMeasure_self` がある。
- **(S3) ℚ_p 側**: `grep -rn "Padic|ℤ_[|ℚ_[" Mathlib/MeasureTheory/` → **0 件**。
  逆向きに `grep -rniE "measure|haar|volume" Mathlib/NumberTheory/Padics/*.lean` →
  実質 0 件(`AddChar.lean` の散文に "measures" が1回出るだけ)。
  すなわち **ℚ_p に測度空間のインスタンスは無い**。
- **(S4) 判定**: 一般論(局所コンパクト群の Haar)は**ある**。ℚ_p への**適用が無い**だけ。
  Lean で実測したところ `LocallyCompactSpace ℚ_[p]` / `T2Space ℚ_[p]` /
  `IsTopologicalAddGroup ℚ_[p]` / `SecondCountableTopology ℚ_[p]` はすべて
  `infer_instance` で通る。足りないのは `MeasurableSpace` / `BorelSpace` の
  インスタンスだけで、これは `borel` で入る。

## ★`MeasurableSpace ℚ_[p]` を大域インスタンスとして置くことについて

mathlib に ℚ_p の測度空間インスタンスは無い(上の実測)ので、衝突しない。
`borel` を入れて `BorelSpace` を宣言する。
-/

namespace ABC3.Found.IUTchIII

open MeasureTheory Metric

variable {p : ℕ} [Fact p.Prime]

noncomputable instance : MeasurableSpace ℚ_[p] := borel _

instance : BorelSpace ℚ_[p] := ⟨rfl⟩

/-! ### 𝒪_k —— 正規化の基準 -/

/-- `𝒪_k`(ℚ_p の整数環)を `ℚ_p` の中で単位閉球として見たもの。

原文 (AbsTopIII p.137):
> ization, i.e., μk(Ok) = 1. We shall refer to μk(−) as the volume on k.

★これが正規化の基準である。 -/
noncomputable def integerBall (p : ℕ) [Fact p.Prime] : TopologicalSpace.PositiveCompacts ℚ_[p] where
  carrier := closedBall 0 1
  isCompact' := isCompact_closedBall 0 1
  interior_nonempty' := by
    refine ⟨0, ?_⟩
    exact mem_interior_iff_mem_nhds.mpr (closedBall_mem_nhds 0 one_pos)

/-- **体積** —— 原文の `μ_k`。`𝒪_k` の体積を 1 に正規化した Haar 測度。 -/
noncomputable def padicVolume (p : ℕ) [Fact p.Prime] : Measure ℚ_[p] :=
  Measure.addHaarMeasure (integerBall p)

/-- **正規化**(原文の条件 (3))。 -/
@[simp]
theorem padicVolume_integerBall : padicVolume p (closedBall (0 : ℚ_[p]) 1) = 1 :=
  Measure.addHaarMeasure_self

instance : (padicVolume p).IsAddLeftInvariant :=
  inferInstanceAs (Measure.addHaarMeasure (integerBall p)).IsAddLeftInvariant

instance : IsFiniteMeasureOnCompacts (padicVolume p) :=
  inferInstanceAs (IsFiniteMeasureOnCompacts (Measure.addHaarMeasure (integerBall p)))

instance : (padicVolume p).IsOpenPosMeasure :=
  inferInstanceAs (Measure.addHaarMeasure (integerBall p)).IsOpenPosMeasure

/-! ### 対数体積 -/

/-- **対数体積** —— 原文の `μ_k^log = log ∘ μ_k`。

原文 (AbsTopIII p.137):
> logarithm R>0 →R] and refer to μlog

★原文の定義域は**コンパクト開集合** `𝕄(k)` である。ここでは全体で定義した関数として
置き、コンパクト開集合の上で原文と一致することを別に示す
(`padicLogVol_wellDefined_of_isCompact_isOpen`)。 -/
noncomputable def padicLogVol (U : Set ℚ_[p]) : ℝ :=
  Real.log (padicVolume p U).toReal

@[simp]
theorem padicLogVol_integerBall : padicLogVol (closedBall (0 : ℚ_[p]) 1) = 0 := by
  simp [padicLogVol]

/-! ### コンパクト集合の体積は有限 —— [AbsTopIII] (L1) の「hence of finite」 -/

/-- コンパクト集合の体積は有限。 -/
theorem padicVolume_lt_top_of_isCompact {U : Set ℚ_[p]} (hU : IsCompact U) :
    padicVolume p U < ⊤ :=
  hU.measure_lt_top

/-- ★**コンパクト開集合の上で対数体積は本物の実数である**——
体積が `0 < μ < ∞` を満たすので `log` が意味を持つ。

原文 (AbsTopIII p.137) が `μ_k : 𝕄(k) → ℝ_{>0}` と書き、`𝕄(k)` を
**コンパクト開集合**の集まりとしているのはこのためである。 -/
theorem padicVolume_pos_lt_top_of_isCompact_isOpen {U : Set ℚ_[p]}
    (hc : IsCompact U) (ho : IsOpen U) (hne : U.Nonempty) :
    0 < padicVolume p U ∧ padicVolume p U < ⊤ :=
  ⟨ho.measure_pos _ hne, hc.measure_lt_top⟩

/-- 上の系: コンパクト開かつ空でなければ `toReal` は正。 -/
theorem padicLogVol_wellDefined_of_isCompact_isOpen {U : Set ℚ_[p]}
    (hc : IsCompact U) (ho : IsOpen U) (hne : U.Nonempty) :
    0 < (padicVolume p U).toReal := by
  obtain ⟨hpos, hlt⟩ := padicVolume_pos_lt_top_of_isCompact_isOpen hc ho hne
  exact ENNReal.toReal_pos hpos.ne' hlt.ne

/-! ### ★非退化性 -/

/-- `‖x‖ ≤ 1/p` と `‖x - 1‖ ≤ 1/p` は両立しない(超距離)。 -/
private theorem disjoint_small_balls :
    Disjoint (closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹))
      ((fun x : ℚ_[p] => (-1) + x) ⁻¹' closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹)) := by
  have hp1 : (1 : ℝ) < p := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hinv : ((p : ℝ)⁻¹) < 1 := by
    rw [inv_lt_one_iff₀]; right; exact hp1
  rw [Set.disjoint_left]
  intro x hx hx2
  simp only [Set.mem_preimage, mem_closedBall, dist_zero_right] at hx hx2
  have h1 : ‖(1 : ℚ_[p])‖ ≤ ((p : ℝ)⁻¹) := by
    have hrw : (1 : ℚ_[p]) = x + (-((-1) + x)) := by ring
    rw [hrw]
    refine le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le hx ?_)
    rw [norm_neg]; exact hx2
  rw [norm_one] at h1
  linarith

/-- ★★**非退化性**: 単位球と `1/p` 球で対数体積が**異なる**。

すなわち `padicLogVol` は定数関数ではない。**これが無ければ
「対数体積を作った」とは言えない**——`logVol ≡ 0` でも型は付くからである。

証明は原文の3条件だけを使う: `1/p` 球とその `1` 平行移動は**互いに素**で、
どちらも単位球に入るので、加法性と平行移動不変性から
`2 · μ(1/p 球) ≤ μ(単位球) = 1`。 -/
theorem padicVolume_smallBall_le_half :
    2 * padicVolume p (closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹)) ≤ 1 := by
  set A : Set ℚ_[p] := closedBall 0 ((p : ℝ)⁻¹) with hA
  set B : Set ℚ_[p] := (fun x : ℚ_[p] => (-1) + x) ⁻¹' A with hB
  have hp1 : (1 : ℝ) < p := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hinv : ((p : ℝ)⁻¹) ≤ 1 := by
    rw [inv_le_one_iff₀]; right; linarith
  have hmeas : padicVolume p B = padicVolume p A := measure_preimage_add _ _ _
  have hAsub : A ⊆ closedBall (0 : ℚ_[p]) 1 := closedBall_subset_closedBall hinv
  have hBsub : B ⊆ closedBall (0 : ℚ_[p]) 1 := by
    intro x hx
    simp only [hB, hA, Set.mem_preimage, mem_closedBall, dist_zero_right] at hx
    simp only [mem_closedBall, dist_zero_right]
    have hrw : x = 1 + ((-1) + x) := by ring
    rw [hrw]
    refine le_trans (IsUltrametricDist.norm_add_le_max _ _) ?_
    rw [norm_one]
    exact max_le le_rfl (le_trans hx hinv)
  have hBmeas : MeasurableSet B :=
    (measurable_const_add (-1 : ℚ_[p])) measurableSet_closedBall
  have hunion : padicVolume p (A ∪ B) = padicVolume p A + padicVolume p B :=
    measure_union disjoint_small_balls hBmeas
  have hle : padicVolume p (A ∪ B) ≤ 1 := by
    rw [← padicVolume_integerBall (p := p)]
    exact measure_mono (Set.union_subset hAsub hBsub)
  rw [hunion, hmeas] at hle
  rw [two_mul]
  exact hle

/-- `1/p` 球の体積は**正**(空でない開球を含むから)。 -/
theorem padicVolume_smallBall_pos :
    0 < padicVolume p (closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹)) := by
  have hp0 : (0 : ℝ) < (p : ℝ)⁻¹ := by
    have : (0 : ℝ) < p := by exact_mod_cast (Fact.out : p.Prime).pos
    positivity
  exact lt_of_lt_of_le (measure_ball_pos _ 0 hp0) (measure_mono ball_subset_closedBall)

/-- ★★**非退化性(本体)**: 単位球と `1/p` 球で**対数体積が異なる**。

`padicLogVol` は定数関数ではない。★これが無ければ「対数体積を作った」とは言えない
——`logVol ≡ 0` でも型は付くからである。 -/
theorem padicLogVol_smallBall_ne_integerBall :
    padicLogVol (closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹))
      ≠ padicLogVol (closedBall (0 : ℚ_[p]) 1) := by
  rw [padicLogVol_integerBall, padicLogVol]
  have hfin : padicVolume p (closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹)) ≠ ⊤ :=
    (isCompact_closedBall 0 _).measure_lt_top.ne
  have hpos : 0 < (padicVolume p (closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹))).toReal :=
    ENNReal.toReal_pos padicVolume_smallBall_pos.ne' hfin
  have h2 : 2 * (padicVolume p (closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹))).toReal ≤ 1 := by
    have h := padicVolume_smallBall_le_half (p := p)
    have := ENNReal.toReal_mono (by simp) h
    rwa [ENNReal.toReal_mul, ENNReal.toReal_ofNat, ENNReal.toReal_one] at this
  intro hlog
  have h1 : (padicVolume p (closedBall (0 : ℚ_[p]) ((p : ℝ)⁻¹))).toReal = 1 := by
    rw [← Real.exp_log hpos, hlog, Real.exp_zero]
  rw [h1] at h2
  linarith

/-! ### [AbsTopIII] (L1) の「hence of finite」—— `Found/` の log-shell に当てる -/

/-- ★`Found/IUTchIII/PadicLog.lean` の log-shell は**有限体積**を持つ。

原文 (AbsTopIII p.5):
> (L1) a log-shell is compact and hence of finite “log-volume” [cf. Corollary

★原文の `hence` はコンパクト性から来る。ここではまさにそれを使っている。
**ただし対数体積が実数になるには体積が正であることも要る**(原文の定義域が
**コンパクト開**集合であるのはそのため)。`logShell` が開であることは示していない。 -/
theorem padicVolume_logShell_lt_top : padicVolume p (logShell (p := p)) < ⊤ :=
  isCompact_logShell'.measure_lt_top

end ABC3.Found.IUTchIII
