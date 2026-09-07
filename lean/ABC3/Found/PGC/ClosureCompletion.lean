import ABC3.Found.PGC.UnramifiedCompletion
import ABC3.Found.PGC.DworkFixedRing
import ABC3.Found.PGC.LubinTateDistinguishedSeparable
import Mathlib.Analysis.Normed.Algebra.Ultra
import Mathlib.Analysis.Normed.Module.Completion
import Mathlib.Topology.Algebra.Valued.NormedValued
import Mathlib.Topology.Algebra.Valued.ValuedField
import Mathlib.RingTheory.PowerSeries.PiTopology
import Mathlib.RingTheory.MvPowerSeries.PiTopology

/-!
# `ℂ_K := Completion K^al` —— `θ(α)` の住処(`sorry` 無し)

経路 Λ の節点 **M3**。Λ6(Dwork の補題)が作る形式冪級数 `θ` を、
`K^al` の元 `α`(`‖α‖ < 1`)で**実際に評価する**ための器を建てる:

```
ℂ_K := (K^al のスペクトルノルムによる完備化)   𝒪_{ℂ_K}   ι_ur : K̂^{ur} ↪ ℂ_K
```

## 典拠(Milne, Class Field Theory)

原典 Milne CFT の地の文(物理ページ 49):

> We write `K̂^un` for its completion and `B` for the valuation ring of `K̂^un`.

同(物理ページ 51):

> for any `α ∈ K^al` … `g(θ(α)) = 0`

★**原典は `θ(α)` の住処に名前を付けていない**。`θ` の係数は `B = 𝒪_{K̂^un}` に
あり、`α` は `K^al` にある。両者を同時に受け止める環が要る——本ファイルは
その器を `𝒪_{ℂ_K}`(`ℂ_K = K^al の完備化` の整数環)として明示する。
省略の合図語(`immediately` / `formally` 等)も無いので `hedge-index` にも
出ない、「原典が名前を付け忘れた節点」である。

★`.src` は**書いていない**。MilneCFT は `ResearchPaper/1_Structured/` に
無く、`sectionId` を正直に埋められないため。Λ6a′
(`LubinTateFieldFIndependent.lean`)と同じ扱いで、逐語引用と物理ページを
この docstring に置くに留める。**嘘の `sectionId` は書かない。**

## ★作らなかったもの(「偽」と確定しているので作ってはならない)

1. **`IsAdicComplete (maximalIdeal 𝒪_{ℂ_K}) 𝒪_{ℂ_K}` は偽**。
   `K^al` の値群 `ℚ` は稠密なので `𝔪² = 𝔪`、したがって
   `⋂ₙ 𝔪ⁿ = 𝔪 ≠ 0` で **`IsHausdorff` が成り立たない**。
   ★本ファイルは代わりに **`CompleteSpace` + `IsLinearTopology`** を置く。
   `PowerSeries.aeval` が要求するのはこの 2 つであって `IsAdicComplete` では
   ない。「adically complete」という言い回しをそのまま Lean に写すと
   **偽の statement になる**ので、写していない。
2. **`IsDiscreteValuationRing 𝒪_{ℂ_K}` も偽**(同じ理由——値群が稠密)。
   `UnramifiedCompletionDVR.lean` には対応物があるが、あれは値群が離散な
   `K̂^{ur}` の話であって、`ℂ_K` に転写してはならない。
3. **`IsAlgClosed ℂ_K`(Ax–Sen–Tate)は入れていない**。要らないから。
   本ファイルの用途(冪級数の評価)に必要なのは
   「完備な超距離ノルム体」だけである。

## 設計の決定(記録)

* **`𝒪_{ℂ_K}` は `ValuationSubring` で作った**(素朴な `Subring` ではない)。
  理由: (a) `unramifiedCompletionInt`(`K̂^{ur}` 側、`UnramifiedCompletion.lean`)
  が `ValuationSubring` なので、`ι_ur` の整数環への制限が同じ形で書ける;
  (b) `IsLocalRing` / `ValuationRing` / `maximalIdeal` が無料で付く
  (`AdjoinIntegers.lean` の素朴な `Subring` 版は `ValuationRing` を手で
  証明していた)。代償は `AdjoinIntegers.lean` で `rfl` だった所属判定が
  `mem_closureCompletionInt` の書き換え 1 行になること——実測で問題なし。
* **`Valued` を `K^al` の instance にしていない**。instance にすると
  `Valued.valuedCompletion` が `ℂ_K` に第二の `Valued` を作り、
  `NormedField.toValued` と菱形になる(`UnramifiedCompletion.lean:69–74`、
  `lean-idioms.md` の「完備化に `Valued` が二重に付く」)。
  `CompletableTopField` は `letI` の中で `Valued.completable` を借りる。
  この選び方のおかげで `Valued.v z = ‖z‖₊` が `ℂ_K` 上で `rfl` になる
  (`valuation_closureCompletion`)。
* **結論に自由なパラメータを出していない**——`ℂ_K` も `𝒪_{ℂ_K}` も `K` に
  しか依存しない。

## 抽象核(§1)

`§1` は**原典の設定に一切依らない**。可換環でも群作用でもなく
「超距離ノルム体 `A` の、ノルム `≤ 1` の元からなる部分環」だけの話であり、
分岐・付値・Galois の語彙が 1 つも出てこない:

* `normBallIdeal` —— 半径 `ε` の閉球が部分環のイデアルになる。
* `isLinearTopology_of_norm_le_one` —— したがって線形位相環。
* `completeSpace_of_norm_le_one_iff` —— 閉単位球は完備。

`§4` はこれに `A := ℂ_K`, `S := 𝒪_{ℂ_K}` を代入するだけである。
同じ構成が `AdjoinIntegers.lean`(`K⟮x⟯` の整数環)にも重複して
書かれているが、あちらは既に着地しているので**書き換えていない**
(本ファイルの `§1` は将来そちらを差し替えられる形にしてある)。

## 退化の自己検査

* `ι_ur` の等長性(`norm_unramifiedToClosureCompletion`)を落とすと
  `‖θ の係数‖ ≤ 1` が `ℂ_K` 側で保証されず、`𝒪_{K̂^{ur}} → 𝒪_{ℂ_K}` が
  そもそも定義できない(`unramifiedIntToClosureCompletionInt` の
  codRestrict の条件が破れる)。すると `HasEval` を持っていても
  `PowerSeries.aeval` の始域が作れず、主張が空虚化する。
* `IsLinearTopology`(`isLinearTopology_closureCompletionInt`)を落とすと
  `PowerSeries.aeval` が**型検査を通らない**。すなわち
  `closureCompletionEval` の存在そのものがこの条件の必要性の証拠である。
* `hasEval_of_norm_lt_one` の `‖z‖ < 1` を `≤ 1` に弱めると偽——
  `z = 1` で `zⁿ = 1 ↛ 0`。
* `𝒪_{ℂ_K}` を `ℂ_K` 全体に置き換えると `IsLinearTopology` が偽になる
  (体の `0` の近傍にイデアルは `0` しか無い)。「整数環に降りる」ことは
  省略できない。

## 逸脱(記録)

* 古典的には `ℂ_p` は「`ℚ̄_p` の完備化」として作り、そのあと Ax–Sen–Tate で
  代数閉性を示す。本ファイルは**完備化までで止めている**(代数閉性は
  作らない)。`Λ_n` はモニック多項式の根集合なので、下流で必要になるのは
  `∏(θλ − ι r) = 0` から所属を出す**整域性だけ**であり、代数閉性は要らない。
* mathlib には `NumberTheory/Padics/Complex.lean` に `ℂ_[p]`
  (`PadicComplex`、`IsAlgClosed` つき)がある。**代案として記録するに留める**:
  本ファイルは `K` を基点にしており(`ℚ_[p]` ではない)、また
  `UnramifiedCompletion.lean` の構成の逐語転写になるぶん予測可能なので、
  `Completion K.closure` を自前で作る道を選んだ。
* `K^ur` の完備化 `K̂^{ur}` と `K^al` の完備化 `ℂ_K` の関係は、本ファイルでは
  **等長埋め込み `ι_ur`(単射)まで**しか主張しない。`ℂ_K` が `K̂^{ur}` 上
  どういう拡大かは射程外。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal Valued

/-! ## 1. 抽象核 —— 超距離ノルム体の「ノルム ≤ 1 の部分環」

★ここには `PAdicLocalField` も分岐も付値も出てこない。純粋に
「超距離ノルム体 `A` と、その中の部分環 `S` で全元のノルムが `≤ 1` のもの」
だけの話である。§4 はこれに代入するだけ。 -/

section UnitBallCore

variable {A : Type*} [NormedField A] [IsUltrametricDist A]
  {σ : Type*} [SetLike σ A] [SubringClass σ A]

/-- **抽象核**——超距離ノルム体 `A` の部分環 `S`(全元のノルムが `≤ 1`)の中で、
半径 `ε` の閉球は `S` の**イデアル**になる。

非アルキメデス三角不等式(`add_mem'`)と「ノルム `≤ 1` の元を掛けても
ノルムは増えない」(`smul_mem'`)だけを使う。`Valued`・
`IsDiscreteValuationRing`・adic 位相を一切経由しない。

`ε ≤ 0` のときは `max ε 0 = 0` に退化して `{y | ‖y‖ ≤ 0}` になるが、
使うのは `ε > 0` の場合だけなので問題ない。 -/
def normBallIdeal (S : σ) (hS : ∀ x : A, x ∈ S → ‖x‖ ≤ 1) (ε : ℝ) : Ideal ↥S where
  carrier := {y | ‖(y : A)‖ ≤ max ε 0}
  zero_mem' := by simp
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    calc ‖(a : A) + (b : A)‖ ≤ max ‖(a : A)‖ ‖(b : A)‖ := IsUltrametricDist.norm_add_le_max _ _
      _ ≤ max ε 0 := max_le ha hb
  smul_mem' := by
    intro c a ha
    simp only [Set.mem_setOf_eq] at *
    show ‖(c : A) * (a : A)‖ ≤ max ε 0
    rw [norm_mul]
    calc ‖(c : A)‖ * ‖(a : A)‖ ≤ 1 * max ε 0 :=
          mul_le_mul (hS _ c.2) ha (norm_nonneg _) zero_le_one
      _ = max ε 0 := one_mul _

/-- `ε > 0` のとき `normBallIdeal S hS ε` は(集合として)ちょうど半径 `ε` の
閉球である。 -/
theorem coe_normBallIdeal (S : σ) (hS : ∀ x : A, x ∈ S → ‖x‖ ≤ 1) {ε : ℝ} (hε : 0 < ε) :
    (↑(normBallIdeal S hS ε) : Set ↥S) = Metric.closedBall 0 ε := by
  ext y
  simp only [SetLike.mem_coe, normBallIdeal, Submodule.mem_mk, AddSubmonoid.mem_mk,
    AddSubsemigroup.mem_mk, Set.mem_setOf_eq, max_eq_left hε.le, Metric.mem_closedBall,
    dist_eq_norm, sub_zero]
  show ‖(y : A)‖ ≤ ε ↔ ‖y‖ ≤ ε
  rfl

/-- `S` における `0` の近傍フィルターは `normBallIdeal`(`ε > 0`)を基底に持つ。 -/
theorem hasBasis_nhds_zero_normBallIdeal (S : σ) (hS : ∀ x : A, x ∈ S → ‖x‖ ≤ 1) :
    (nhds (0 : ↥S)).HasBasis (fun ε : ℝ => 0 < ε)
      (fun ε => (↑(normBallIdeal S hS ε) : Set ↥S)) :=
  (Metric.nhds_basis_closedBall (x := (0 : ↥S))).congr (fun _ => Iff.rfl)
    (fun _ hε => (coe_normBallIdeal S hS hε).symm)

/-- **抽象核**——超距離ノルム体の「ノルム `≤ 1` の元からなる部分環」は
**線形位相環**。`PowerSeries.aeval` が要求する条件のうち、いちばん
配線が面倒な一つがこれだけで出る。 -/
theorem isLinearTopology_of_norm_le_one (S : σ) (hS : ∀ x : A, x ∈ S → ‖x‖ ≤ 1) :
    IsLinearTopology ↥S ↥S :=
  IsLinearTopology.mk_of_hasBasis _ (hasBasis_nhds_zero_normBallIdeal S hS)

omit [IsUltrametricDist A] [SubringClass σ A] in
/-- **抽象核**——「ノルム `≤ 1`」でちょうど特徴づけられる部分集合は閉。 -/
theorem isClosed_of_norm_le_one_iff (S : σ) (hS : ∀ x : A, x ∈ S ↔ ‖x‖ ≤ 1) :
    IsClosed ((S : Set A)) := by
  have heq : (S : Set A) = Metric.closedBall 0 1 := by
    ext y
    simp only [SetLike.mem_coe, hS, Metric.mem_closedBall, dist_eq_norm, sub_zero]
  rw [heq]
  exact Metric.isClosed_closedBall

omit [IsUltrametricDist A] [SubringClass σ A] in
/-- **抽象核**——完備なノルム体の閉単位球は完備。 -/
theorem completeSpace_of_norm_le_one_iff [CompleteSpace A] (S : σ)
    (hS : ∀ x : A, x ∈ S ↔ ‖x‖ ≤ 1) : CompleteSpace ↥S :=
  haveI := isClosed_of_norm_le_one_iff S hS
  IsClosed.completeSpace_coe

end UnitBallCore

variable {p : ℕ} [Fact p.Prime]

/-! ## 2. `ℂ_K` —— `K^al` のスペクトルノルムによる完備化 -/

/-- `K^al` は completable(付値体だから)。

★`Valued` を `K.closure` の instance にせず `letI` で借りているのは、
`Valued.valuedCompletion` が `ℂ_K` に第二の `Valued` を作るのを防ぐため
(モジュール docstring 参照)。 -/
noncomputable scoped instance completableTopField_closure (K : PAdicLocalField p) :
    CompletableTopField K.closure :=
  letI : Valued K.closure NNReal := NormedField.toValued
  Valued.completable

/-- **`ℂ_K`** —— 代数閉包 `K^al` のスペクトルノルムによる完備化。

原典 Milne CFT p.51 の `θ(α)`(`α ∈ K^al`)が住む体。 -/
abbrev closureCompletion (K : PAdicLocalField p) : Type :=
  UniformSpace.Completion K.closure

/-- `K^al` のノルムは `K` 上のスペクトルノルムそのもの(`closureNormedField` の定義)。 -/
theorem norm_eq_spectralNorm_closure (K : PAdicLocalField p) (x : K.closure) :
    ‖x‖ = spectralNorm K.carrier K.closure x := rfl

theorem norm_algebraMap_closure (K : PAdicLocalField p) (x : K.carrier) :
    ‖algebraMap K.carrier K.closure x‖ = ‖x‖ := spectralNorm_extends x

noncomputable scoped instance isUltrametricDist_closureCompletion (K : PAdicLocalField p) :
    IsUltrametricDist (closureCompletion K) := IsUltrametricDist.of_normedAlgebra K.carrier

/-- `ℂ_K` の付値。**`Valued.v z = ‖z‖₊` が `rfl`** になるように、
完備化の `Valued.valuedCompletion` ではなくノルムから作る。 -/
noncomputable scoped instance valued_closureCompletion (K : PAdicLocalField p) :
    Valued (closureCompletion K) NNReal := NormedField.toValued

theorem norm_coe_closureCompletion (K : PAdicLocalField p) (x : K.closure) :
    ‖(x : closureCompletion K)‖ = ‖x‖ := UniformSpace.Completion.norm_coe x

theorem norm_algebraMap_closureCompletion (K : PAdicLocalField p) (x : K.carrier) :
    ‖algebraMap K.carrier (closureCompletion K) x‖ = ‖x‖ := by
  rw [UniformSpace.Completion.algebraMap_def, norm_coe_closureCompletion, norm_algebraMap_closure]

theorem valuation_closureCompletion (K : PAdicLocalField p) (z : closureCompletion K) :
    Valued.v z = ‖z‖₊ := rfl

theorem valuation_coe_closureCompletion (K : PAdicLocalField p) (x : K.closure) :
    Valued.v (x : closureCompletion K) = ‖x‖₊ := by
  rw [valuation_closureCompletion]
  ext
  exact norm_coe_closureCompletion K x

noncomputable scoped instance nontriviallyNormedField_closureCompletion
    (K : PAdicLocalField p) : NontriviallyNormedField (closureCompletion K) where
  non_trivial := by
    obtain ⟨x, hx⟩ := NontriviallyNormedField.non_trivial (α := K.carrier)
    exact ⟨algebraMap K.carrier (closureCompletion K) x,
      by rwa [norm_algebraMap_closureCompletion]⟩

noncomputable scoped instance charZero_closureCompletion (K : PAdicLocalField p) :
    CharZero (closureCompletion K) := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  exact (RingHom.charZero_iff
    (algebraMap K.carrier (closureCompletion K)).injective).mp inferInstance

/-- `K^al` は `ℂ_K` で稠密。 -/
theorem denseRange_coe_closureCompletion (K : PAdicLocalField p) :
    DenseRange (fun x : K.closure => (x : closureCompletion K)) :=
  UniformSpace.Completion.denseRange_coe

/-- **任意精度の近似**。 -/
theorem exists_norm_sub_lt_closureCompletion (K : PAdicLocalField p) (z : closureCompletion K)
    {ε : ℝ} (hε : 0 < ε) : ∃ x : K.closure, ‖z - (x : closureCompletion K)‖ < ε := by
  obtain ⟨-, ⟨x, rfl⟩, hx⟩ :=
    Metric.mem_closure_iff.mp
      ((denseRange_coe_closureCompletion K).closure_eq ▸ Set.mem_univ z) ε hε
  exact ⟨x, by rwa [← dist_eq_norm]⟩

/-- **整数の近似**——`‖z‖ ≤ 1` なら近似元も `𝒪_{K^al}` の中に取れる。

退化の自己検査:`hz` を落とすと `‖x‖ ≤ 1` は**偽**——`‖z‖ > 1` なら
`ε` を小さく取ったとき `‖x‖ = ‖z‖ > 1` になる。 -/
theorem exists_norm_le_one_norm_sub_lt_closureCompletion (K : PAdicLocalField p)
    (z : closureCompletion K) (hz : ‖z‖ ≤ 1) {ε : ℝ} (hε : 0 < ε) :
    ∃ x : K.closure, ‖x‖ ≤ 1 ∧ ‖z - (x : closureCompletion K)‖ < ε := by
  obtain ⟨x, hx⟩ := exists_norm_sub_lt_closureCompletion K z (lt_min hε one_pos)
  refine ⟨x, ?_, lt_of_lt_of_le hx (min_le_left _ _)⟩
  have h1 : ‖(x : closureCompletion K)‖ ≤ 1 := by
    have hmax := IsUltrametricDist.norm_add_le_max z (-(z - (x : closureCompletion K)))
    simp only [neg_sub, add_sub_cancel] at hmax
    refine hmax.trans (max_le hz ?_)
    rw [norm_sub_rev]
    exact le_of_lt (lt_of_lt_of_le hx (min_le_right _ _))
  rwa [norm_coe_closureCompletion] at h1

/-! ## 3. `𝒪_{ℂ_K}` —— 原典 p.49 の `B` にあたるものの `K^al` 版 -/

/-- **`𝒪_{ℂ_K}`** —— `ℂ_K` の整数環(付値 `≤ 1` の元)。

★`ValuationSubring` で作る(モジュール docstring の「設計の決定」参照)。
これで `IsLocalRing` / `ValuationRing` / `maximalIdeal` が無料で付き、
`unramifiedCompletionInt`(`K̂^{ur}` 側)と同じ形で扱える。 -/
noncomputable def closureCompletionInt (K : PAdicLocalField p) :
    ValuationSubring (closureCompletion K) :=
  (Valued.v : Valuation (closureCompletion K) NNReal).valuationSubring

@[simp] theorem mem_closureCompletionInt (K : PAdicLocalField p) (z : closureCompletion K) :
    z ∈ closureCompletionInt K ↔ ‖z‖ ≤ 1 := by
  rw [closureCompletionInt, Valuation.mem_valuationSubring_iff, valuation_closureCompletion]
  exact ⟨fun h => by exact_mod_cast h, fun h => by exact_mod_cast h⟩

theorem closureCompletionInt_integers (K : PAdicLocalField p) :
    (Valued.v : Valuation (closureCompletion K) NNReal).Integers (closureCompletionInt K) :=
  Valuation.integer.integers _

theorem norm_coe_closureCompletionInt (K : PAdicLocalField p) (w : ↥(closureCompletionInt K)) :
    ‖w‖ = ‖(w : closureCompletion K)‖ := rfl

theorem norm_le_one_closureCompletionInt (K : PAdicLocalField p)
    (w : ↥(closureCompletionInt K)) : ‖(w : closureCompletion K)‖ ≤ 1 :=
  (mem_closureCompletionInt K _).mp w.2

theorem not_isUnit_closureCompletionInt (K : PAdicLocalField p) (w : ↥(closureCompletionInt K)) :
    ¬ IsUnit w ↔ ‖(w : closureCompletion K)‖ < 1 := by
  have h := Valuation.Integer.not_isUnit_iff_valuation_lt_one
    (v := (Valued.v : Valuation (closureCompletion K) NNReal)) (x := w)
  rw [valuation_closureCompletion] at h
  exact h.trans ⟨fun hh => by exact_mod_cast hh, fun hh => by exact_mod_cast hh⟩

/-- `𝒪_{ℂ_K}` の極大イデアルはちょうど `‖·‖ < 1`。

★ここが `IsAdicComplete` が**偽**になる理由の所在でもある: `K^al` の
値群は稠密なので `𝔪² = 𝔪` となり、`⋂ₙ 𝔪ⁿ = 𝔪 ≠ 0`。本ファイルは
`IsAdicComplete` を主張しない(モジュール docstring 参照)。 -/
theorem mem_maximalIdeal_closureCompletionInt (K : PAdicLocalField p)
    (w : ↥(closureCompletionInt K)) :
    w ∈ IsLocalRing.maximalIdeal ↥(closureCompletionInt K)
      ↔ ‖(w : closureCompletion K)‖ < 1 :=
  (IsLocalRing.mem_maximalIdeal w).trans (not_isUnit_closureCompletionInt K w)

/-! ## 4. `𝒪_{ℂ_K}` の位相 —— §1 の抽象核に代入するだけ -/

/-- `𝒪_{ℂ_K}` は `ℂ_K` の中で閉。 -/
theorem isClosed_closureCompletionInt (K : PAdicLocalField p) :
    IsClosed ((↑(closureCompletionInt K) : Set (closureCompletion K))) :=
  isClosed_of_norm_le_one_iff (closureCompletionInt K) (mem_closureCompletionInt K)

/-- **`𝒪_{ℂ_K}` は完備**——`ℂ_K` が完備で `𝒪_{ℂ_K}` がその中で閉だから。
★`PowerSeries.aeval` が要求する条件その 1。 -/
noncomputable scoped instance completeSpace_closureCompletionInt (K : PAdicLocalField p) :
    CompleteSpace ↥(closureCompletionInt K) :=
  completeSpace_of_norm_le_one_iff (closureCompletionInt K) (mem_closureCompletionInt K)

/-- 半径 `ε` の閉球(`𝒪_{ℂ_K}` のイデアル)。§1 の `normBallIdeal` の specialisation。 -/
noncomputable def closureCompletionBall (K : PAdicLocalField p) (ε : ℝ) :
    Ideal ↥(closureCompletionInt K) :=
  normBallIdeal (closureCompletionInt K) (fun _ hz => (mem_closureCompletionInt K _).mp hz) ε

theorem coe_closureCompletionBall (K : PAdicLocalField p) {ε : ℝ} (hε : 0 < ε) :
    (↑(closureCompletionBall K ε) : Set ↥(closureCompletionInt K)) = Metric.closedBall 0 ε :=
  coe_normBallIdeal _ _ hε

theorem hasBasis_nhds_closureCompletionInt (K : PAdicLocalField p) :
    (nhds (0 : ↥(closureCompletionInt K))).HasBasis (fun ε : ℝ => 0 < ε)
      (fun ε => (↑(closureCompletionBall K ε) : Set ↥(closureCompletionInt K))) :=
  hasBasis_nhds_zero_normBallIdeal _ _

/-- **`𝒪_{ℂ_K}` は線形位相環**。★`PowerSeries.aeval` が要求する条件その 2。

★`IsDiscreteValuationRing` も `Ideal.isLinearTopology`(adic 位相)も
経由していない——`ℂ_K` の値群は稠密なのでそれらは**使えない**。 -/
theorem isLinearTopology_closureCompletionInt (K : PAdicLocalField p) :
    IsLinearTopology ↥(closureCompletionInt K) ↥(closureCompletionInt K) :=
  isLinearTopology_of_norm_le_one (closureCompletionInt K)
    (fun _ hz => (mem_closureCompletionInt K _).mp hz)

/-! ## 5. `ι_al : K^al ↪ ℂ_K` -/

/-- **`ι_al : K^al →+* ℂ_K`**。 -/
noncomputable def closureCompletionCoe (K : PAdicLocalField p) :
    K.closure →+* closureCompletion K :=
  UniformSpace.Completion.coeRingHom

@[simp] theorem closureCompletionCoe_apply (K : PAdicLocalField p) (x : K.closure) :
    closureCompletionCoe K x = (x : closureCompletion K) := rfl

theorem norm_closureCompletionCoe (K : PAdicLocalField p) (x : K.closure) :
    ‖closureCompletionCoe K x‖ = ‖x‖ := norm_coe_closureCompletion K x

theorem isometry_closureCompletionCoe (K : PAdicLocalField p) :
    Isometry (closureCompletionCoe K) :=
  Isometry.of_dist_eq fun x y => by
    rw [dist_eq_norm, dist_eq_norm, ← map_sub, norm_closureCompletionCoe]

theorem closureCompletionCoe_injective (K : PAdicLocalField p) :
    Function.Injective (closureCompletionCoe K) :=
  UniformSpace.Completion.coe_injective K.closure

/-! ## 6. `ι_ur : K̂^{ur} ↪ ℂ_K`

`unramifiedClosure K` は `IntermediateField K.carrier K.closure` なので、
包含 `K^ur → K^al` は**定義からノルムを保つ**(どちらも `K` 上の
スペクトルノルムの制限)。そこへ mathlib の `Isometry.extensionHom` を当てる。 -/

/-- `K^ur ⊆ K^al` はノルムを保つ(どちらも `K` 上のスペクトルノルム)。 -/
theorem norm_coe_unramifiedClosure_closure (K : PAdicLocalField p) (x : ↥(unramifiedClosure K)) :
    ‖(x : K.closure)‖ = ‖x‖ := rfl

/-- `K^ur ↪ ℂ_K`(体の包含と `ι_al` の合成)。 -/
noncomputable def unramifiedClosureToClosureCompletion (K : PAdicLocalField p) :
    ↥(unramifiedClosure K) →+* closureCompletion K :=
  (closureCompletionCoe K).comp (SubringClass.subtype (unramifiedClosure K))

@[simp] theorem unramifiedClosureToClosureCompletion_apply (K : PAdicLocalField p)
    (x : ↥(unramifiedClosure K)) :
    unramifiedClosureToClosureCompletion K x = ((x : K.closure) : closureCompletion K) := rfl

theorem norm_unramifiedClosureToClosureCompletion (K : PAdicLocalField p)
    (x : ↥(unramifiedClosure K)) : ‖unramifiedClosureToClosureCompletion K x‖ = ‖x‖ :=
  norm_coe_closureCompletion K (x : K.closure)

theorem isometry_unramifiedClosureToClosureCompletion (K : PAdicLocalField p) :
    Isometry (unramifiedClosureToClosureCompletion K) :=
  Isometry.of_dist_eq fun x y => by
    rw [dist_eq_norm, dist_eq_norm, ← map_sub, norm_unramifiedClosureToClosureCompletion]

/-- **`ι_ur : K̂^{ur} →+* ℂ_K`** —— 等長写像 `K^ur → ℂ_K` の完備化への一意な
連続延長(`Isometry.extensionHom`)。

原典 Milne CFT p.49 の `K̂^un` と、p.51 の `θ(α)` の住む体を繋ぐ写像。 -/
noncomputable def unramifiedToClosureCompletion (K : PAdicLocalField p) :
    unramifiedCompletion K →+* closureCompletion K :=
  (isometry_unramifiedClosureToClosureCompletion K).extensionHom

@[simp] theorem unramifiedToClosureCompletion_coe (K : PAdicLocalField p)
    (x : ↥(unramifiedClosure K)) :
    unramifiedToClosureCompletion K (x : unramifiedCompletion K)
      = ((x : K.closure) : closureCompletion K) :=
  UniformSpace.Completion.extension_coe
    (isometry_unramifiedClosureToClosureCompletion K).uniformContinuous x

theorem continuous_unramifiedToClosureCompletion (K : PAdicLocalField p) :
    Continuous (unramifiedToClosureCompletion K) :=
  UniformSpace.Completion.continuous_extension

/-- **`ι_ur` は等長**——`K^ur` 上で等長で、両辺が連続だから稠密性で全体へ。

退化の自己検査: これを落とすと `𝒪_{K̂^{ur}} → 𝒪_{ℂ_K}` が定義できず、
`θ` の係数を `ℂ_K` 側に持って行けない(モジュール docstring 参照)。 -/
theorem norm_unramifiedToClosureCompletion (K : PAdicLocalField p) (z : unramifiedCompletion K) :
    ‖unramifiedToClosureCompletion K z‖ = ‖z‖ := by
  have h : (fun w : unramifiedCompletion K => ‖unramifiedToClosureCompletion K w‖)
      = (fun w : unramifiedCompletion K => ‖w‖) :=
    UniformSpace.Completion.ext
      (continuous_norm.comp (continuous_unramifiedToClosureCompletion K)) continuous_norm
      (fun x => by
        rw [unramifiedToClosureCompletion_coe, norm_coe_closureCompletion,
          norm_coe_unramifiedCompletion, norm_coe_unramifiedClosure_closure])
  exact congrFun h z

theorem isometry_unramifiedToClosureCompletion (K : PAdicLocalField p) :
    Isometry (unramifiedToClosureCompletion K) :=
  Isometry.of_dist_eq fun x y => by
    rw [dist_eq_norm, dist_eq_norm, ← map_sub, norm_unramifiedToClosureCompletion]

/-- **`ι_ur` は単射**(等長だから)。 -/
theorem unramifiedToClosureCompletion_injective (K : PAdicLocalField p) :
    Function.Injective (unramifiedToClosureCompletion K) :=
  (isometry_unramifiedToClosureCompletion K).injective

theorem coe_algebraMap_unramifiedClosure (K : PAdicLocalField p) (a : K.carrier) :
    ((algebraMap K.carrier ↥(unramifiedClosure K) a : ↥(unramifiedClosure K)) : K.closure)
      = algebraMap K.carrier K.closure a := rfl

theorem unramifiedToClosureCompletion_algebraMap (K : PAdicLocalField p) (a : K.carrier) :
    unramifiedToClosureCompletion K (algebraMap K.carrier (unramifiedCompletion K) a)
      = algebraMap K.carrier (closureCompletion K) a := by
  rw [UniformSpace.Completion.algebraMap_def, unramifiedToClosureCompletion_coe,
    coe_algebraMap_unramifiedClosure, ← UniformSpace.Completion.algebraMap_def]

/-- **`ι_ur` は `K`-代数準同型**。 -/
noncomputable def unramifiedToClosureCompletionAlg (K : PAdicLocalField p) :
    unramifiedCompletion K →ₐ[K.carrier] closureCompletion K :=
  { unramifiedToClosureCompletion K with
    commutes' := unramifiedToClosureCompletion_algebraMap K }

/-! ## 7. 整数環のレベル —— `𝒪_K → 𝒪_{K̂^{ur}} → 𝒪_{ℂ_K}` -/

/-- `ι_ur` の整数環への制限 **`𝒪_{K̂^{ur}} → 𝒪_{ℂ_K}`**。
`ι_ur` が等長であること(`norm_unramifiedToClosureCompletion`)が
codRestrict の条件を保証する。 -/
noncomputable def unramifiedIntToClosureCompletionInt (K : PAdicLocalField p) :
    ↥(unramifiedCompletionInt K) →+* ↥(closureCompletionInt K) :=
  RingHom.codRestrict
    ((unramifiedToClosureCompletion K).comp (SubringClass.subtype (unramifiedCompletionInt K)))
    (closureCompletionInt K)
    (fun w => by
      rw [mem_closureCompletionInt]
      show ‖unramifiedToClosureCompletion K (w : unramifiedCompletion K)‖ ≤ 1
      rw [norm_unramifiedToClosureCompletion]
      exact (mem_unramifiedCompletionInt K _).mp w.2)

@[simp] theorem unramifiedIntToClosureCompletionInt_coe (K : PAdicLocalField p)
    (w : ↥(unramifiedCompletionInt K)) :
    ((unramifiedIntToClosureCompletionInt K w : ↥(closureCompletionInt K)) : closureCompletion K)
      = unramifiedToClosureCompletion K (w : unramifiedCompletion K) := rfl

theorem unramifiedIntToClosureCompletionInt_injective (K : PAdicLocalField p) :
    Function.Injective (unramifiedIntToClosureCompletionInt K) := by
  intro a b hab
  apply Subtype.ext
  apply unramifiedToClosureCompletion_injective K
  rw [← unramifiedIntToClosureCompletionInt_coe, ← unramifiedIntToClosureCompletionInt_coe, hab]

/-- **`𝒪_{ℂ_K}` は `𝒪_{K̂^{ur}}`-代数**。★`PowerSeries.aeval` が要求する条件その 3
(`θ` の係数環をスカラーにする)。 -/
noncomputable scoped instance algebra_unramifiedInt_closureCompletionInt
    (K : PAdicLocalField p) :
    Algebra ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K) :=
  (unramifiedIntToClosureCompletionInt K).toAlgebra

theorem algebraMap_unramifiedInt_closureCompletionInt_coe (K : PAdicLocalField p)
    (w : ↥(unramifiedCompletionInt K)) :
    ((algebraMap ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K) w
        : ↥(closureCompletionInt K)) : closureCompletion K)
      = unramifiedToClosureCompletion K (w : unramifiedCompletion K) := rfl

theorem continuous_algebraMap_closureCompletionInt (K : PAdicLocalField p) :
    Continuous (algebraMap ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K)) := by
  apply Continuous.subtype_mk
  exact (continuous_unramifiedToClosureCompletion K).comp continuous_subtype_val

/-- ★`PowerSeries.aeval` が要求する条件その 4。 -/
theorem continuousSMul_closureCompletionInt (K : PAdicLocalField p) :
    ContinuousSMul ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K) :=
  continuousSMul_of_algebraMap _ _ (continuous_algebraMap_closureCompletionInt K)

/-- `𝒪_{K̂^{ur}}` は `𝒪_K`-代数(`baseIntHom`、`DworkFixedRing.lean`)。 -/
noncomputable scoped instance algebra_base_unramifiedCompletionInt (K : PAdicLocalField p) :
    Algebra (𝒪[K.carrier]) ↥(unramifiedCompletionInt K) :=
  (baseIntHom K).toAlgebra

/-- `𝒪_{ℂ_K}` は `𝒪_K`-代数。★**合成として定義する**ので
`IsScalarTower` が `rfl` で出る。 -/
noncomputable scoped instance algebra_base_closureCompletionInt (K : PAdicLocalField p) :
    Algebra (𝒪[K.carrier]) ↥(closureCompletionInt K) :=
  (((unramifiedIntToClosureCompletionInt K).comp (baseIntHom K))).toAlgebra

/-- **`𝒪_K → 𝒪_{K̂^{ur}} → 𝒪_{ℂ_K}` はスカラー塔**。 -/
scoped instance isScalarTower_base_closureCompletionInt (K : PAdicLocalField p) :
    IsScalarTower (𝒪[K.carrier]) ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K) :=
  IsScalarTower.of_algebraMap_eq' rfl

/-- `𝒪_K → 𝒪_{ℂ_K}` は `algebraMap K.carrier ℂ_K` の制限にほかならない。 -/
theorem algebraMap_base_closureCompletionInt_coe (K : PAdicLocalField p) (a : 𝒪[K.carrier]) :
    ((algebraMap (𝒪[K.carrier]) ↥(closureCompletionInt K) a : ↥(closureCompletionInt K))
        : closureCompletion K)
      = algebraMap K.carrier (closureCompletion K) (a : K.carrier) := by
  show unramifiedToClosureCompletion K ((baseIntHom K a : ↥(unramifiedCompletionInt K))
      : unramifiedCompletion K) = _
  rw [baseIntHom_coe, unramifiedToClosureCompletion_algebraMap]

/-! ## 8. 冪級数の評価 —— `θ(α)` -/

/-- `‖z‖ < 1` なら `z` は `𝒪_{ℂ_K}` の位相的冪零元。

退化の自己検査: `<` を `≤` に弱めると偽(`z = 1` で `zⁿ = 1 ↛ 0`)。 -/
theorem hasEval_of_norm_lt_one (K : PAdicLocalField p) {z : closureCompletion K} (hz : ‖z‖ < 1) :
    PowerSeries.HasEval
      (⟨z, (mem_closureCompletionInt K z).mpr hz.le⟩ : ↥(closureCompletionInt K)) := by
  apply tendsto_pow_atTop_nhds_zero_of_norm_lt_one
  exact hz

/-- `α ∈ K^al` が `‖α‖ < 1` なら、その `ℂ_K` での像は `𝒪_{ℂ_K}` の
位相的冪零元。 -/
theorem hasEval_coe_of_norm_lt_one (K : PAdicLocalField p) {lam : K.closure} (hlam : ‖lam‖ < 1) :
    PowerSeries.HasEval
      (⟨(lam : closureCompletion K),
        (mem_closureCompletionInt K _).mpr
          (by rw [norm_coe_closureCompletion]; exact hlam.le)⟩ : ↥(closureCompletionInt K)) :=
  hasEval_of_norm_lt_one K (by rw [norm_coe_closureCompletion]; exact hlam)

/-- ★★★★★★★★★★★★★★**(M3)`θ ↦ θ(α)`** —— 係数が `𝒪_{K̂^{ur}}` にある
形式冪級数を、`‖α‖ < 1` なる `α ∈ K^al` で評価する `𝒪_{K̂^{ur}}`-代数準同型。

原典 Milne CFT p.51 の `θ(α)` の**器**。`PowerSeries.aeval` が要求する
4 条件(`CompleteSpace`・`IsLinearTopology`・`Algebra`・`ContinuousSMul`)は
`§4` と `§7` で揃えてある。

退化の自己検査: `IsLinearTopology` を落とすとこの `def` は**型検査を
通らない**。すなわち本宣言の存在が `§4` の必要性の証拠である。 -/
noncomputable def closureCompletionEval (K : PAdicLocalField p) {lam : K.closure}
    (hlam : ‖lam‖ < 1) :
    PowerSeries ↥(unramifiedCompletionInt K)
      →ₐ[↥(unramifiedCompletionInt K)] ↥(closureCompletionInt K) :=
  haveI := isLinearTopology_closureCompletionInt K
  haveI := continuousSMul_closureCompletionInt K
  PowerSeries.aeval (hasEval_coe_of_norm_lt_one K hlam)

/-- **`θ(α) ∈ 𝒪_{ℂ_K}`** —— 完了判定そのもの。 -/
noncomputable def evalAt (K : PAdicLocalField p)
    (theta : PowerSeries ↥(unramifiedCompletionInt K))
    {lam : K.closure} (hlam : ‖lam‖ < 1) : ↥(closureCompletionInt K) :=
  closureCompletionEval K hlam theta

/-! ### `Λ_n` の元での評価(Λ6 の消費側が実際に使う形) -/

/-- `Λ_n` の元はノルムが `1` 未満(`LubinTateDistinguishedSeparable.lean` の
`spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints` を `‖·‖` の
言葉に直しただけ——両者は `rfl` で繋がる)。 -/
theorem norm_lt_one_of_mem_iteratedLubinTateTorsionPoints (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    ‖x‖ < 1 :=
  spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n x hx

/-- ★★★★★★★★★★★★★★★★★★**(M3)`θ(λ)` を `Λ_n` の元 `λ` で実際に作る**。

`θ : PowerSeries 𝒪_{K̂^{ur}}` は Λ6(Dwork)が出す形式冪級数、
`λ ∈ Λ_n ⊆ K^al` は Lubin–Tate の捩れ点。値は `𝒪_{ℂ_K}` に落ちる。

★これが M3 の完了判定である: 「`θ` の係数環(`K̂^{ur}` 側)」と
「`λ` の住む体(`K^al` 側)」という**別々の完備化**を、
`ι_ur` の等長性を通して 1 つの環 `𝒪_{ℂ_K}` の上で会わせている。 -/
noncomputable def evalAtTorsionPoint (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (theta : PowerSeries ↥(unramifiedCompletionInt K)) (lam : K.closure)
    (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    ↥(closureCompletionInt K) :=
  evalAt K theta
    (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n lam hlam)

end ABC3.Found.PGC
