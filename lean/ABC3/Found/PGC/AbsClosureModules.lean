import ABC3.Found.PGC.ClosureCompletion
import ABC3.Found.PGC.RamificationFiltrationZero
import ABC3.Skeleton.PGC.Section2

/-!
# `𝒪_{K̄}` と `K̄^∧` を `Γ_K`-加群として構成する —— [pGC] Proposition 2.2 の主語

`Skeleton/PGC/Section2.lean` の `prop_2_2` は

```lean
theorem prop_2_2 (_RF : RamificationFiltration p)
    (IntKbar CompKbar : PAdicLocalField p → Type*)
    [∀ K, AddCommGroup (IntKbar K)] [∀ K, DistribMulAction K.absGal (IntKbar K)]
    [∀ K, AddCommGroup (CompKbar K)] [∀ K, DistribMulAction K.absGal (CompKbar K)] :
    RecoverableAsAddModule IntKbar ∧ RecoverableAsAddModule CompKbar
```

と、原文の `𝒪_{K̄}`・`K̄^∧` を**未構築の対象**として抽象化していた
(「具体的な構成は別途 `Found/` の課題として残す」)。本ファイルがその課題を埋める:

```
IntKbar K := absClosureInt K    = 𝒪_{K̄}   (K^al のスペクトルノルムの付値環)
CompKbar K := closureCompletion K = K̄^∧ = ℂ_K
```

の 2 つを**実際に構成し**、`AddCommGroup` と `DistribMulAction K.absGal`
(実際にはより強い `MulSemiringAction K.absGal`)を付ける。

## 典拠

原文 (pGC 物理 p.5, Proposition 2.2):

> Suppose that we are given the following group-theoretic data: the topological group Γ_K,
> together with the indexed filtration Γ_K^v for all v > 0. Then the Γ_K-modules O[scr]_K[bar], and
> K[bar]∧ can be recovered group-theoretically from this group-theoretic data.

同, 直前の setup (`#setup-2-completion`):

> Now let us denote by K^{∧} the p-adic completion of the field K. Thus, K^{∧} is the quotient
> field of the p-adic completion of the ring of integers O_K.

同, Proof(`0_Source` の `.txt` 204–219 行を `brief.mjs` が切り出したもの。
◆抽出器の字形の潰れは直していない):

> To obtain K
> ∧, we note that K
> ∧is just the p-adic completion of OK tensored over Zp
> with Qp.

★**本ファイルは Proposition 2.2 そのものを証明していない**(群論的回復の主張には
分岐フィルトレーションが要る)。埋めたのは**主語の構成**、すなわち
「`𝒪_{K̄}` と `K̄^∧` が実在し、`Γ_K` が加法群の自己同型として(実際には環の自己同型
として)作用する」ことである。

## ★docstring の前提の訂正(実測)

`prop_2_2` の docstring は「`Found/PGC/LocalFieldNorm.lean` のスペクトルノルム機構
(有限次拡大にのみ適用)をそのままでは使えない」と書いていたが、これは**古い**。
`LocalFieldNorm.lean:132–141` は既に

> `spectralNorm.normedField`/`normedAlgebra` は基点が完備でありさえすれば
> 有限次拡大に限らず任意の代数拡大に対して働くので、`K.closure` が
> `K.carrier` 上有限次でなくても(実際、代数閉包は一般に無限次)問題ない。

と述べ、`closureNormedField : NormedField K.closure` を無条件に与えている。
したがって「無限次拡大 `K̄` 上の付値の構成」は既に済んでいた。

## ★原文の「p 進完備化」= スペクトルノルムによる完備化 であること

原文は `K̄^∧` を「`𝒪_{K̄}` の p 進完備化の商体」と定義する。一方
`ClosureCompletion.lean` の `closureCompletion K` は `K^al` の
**スペクトルノルムによる完備化**である。両者が一致する理由は

```
p^n · 𝒪_{K̄}  =  {x ∈ K̄ : ‖x‖ ≤ ‖p‖^n}
```

——`𝒪_{K̄}` の p 進位相はノルム位相にほかならない、ということである。
この等式は本ファイルの `mem_pow_p_mul_absClosureInt` で **`sorry` 無しに**示した
(`‖x/p^n‖ = ‖x‖/‖p‖^n` だけで出る。値群の可除性は要らない)。
★これは逸脱ではなく**同定**であり、原文の定義を我々の構成へ橋渡しする。

## 抽象核(§1)—— 原典の設定に一切依らない

`§1` には分岐・付値・Galois・`PAdicLocalField` の語彙が **1 つも出てこない**。
中身は「位相環の完備化」と「モノイド作用」だけである:

| 宣言 | 主張 |
|---|---|
| `isometry_smul_of_norm_smul` | ノルムを保つ加法的作用は等長 |
| `uniformContinuousConstSMul_of_norm_smul` | したがって一様連続 |
| `completionMulSemiringAction` | ★一様連続な**環**作用は完備化へ延びる(mathlib 不在) |
| `norm_smul_completion` | 延ばした作用もノルムを保つ |
| `faithfulSMul_completion` | 忠実性は完備化へ上がる |
| `mulSemiringActionOfSubringClass` | 不変部分環への作用の制限 |
| `eq_of_smul_eq_smul_of_scaling` | 忠実性が部分集合へ降りるための十分条件 |

`§2`–`§4` はこれに `M := K.absGal`, `A := K.closure` を**代入するだけ**である。

★mathlib には `DistribMulAction M (Completion α)`
(`Topology/Algebra/GroupCompletion.lean:165`)は在るが、
`MulSemiringAction M (Completion α)` は**無い**(実測: `.cache/mathlib-index.txt` を
`UniformSpace.Completion` で引いた 45 個の instance に該当なし)。
`completionMulSemiringAction` はそれを補うもので、mathlib の
`DistribMulAction` インスタンスを `with` で受け継いでいる(菱形を作らない)。

## ★退化の自己検査

* **`absClosureInt`(`𝒪_{K̄}`)と `closureCompletionInt`(`𝒪_{ℂ_K}`)は別物**である。
  前者は `K^al` の(完備でない)付値環、後者はその完備化の付値環。
  本ファイルは両者を結ぶ**単射・等長・`Γ_K` 同変**の環準同型
  `absClosureIntToClosureCompletionInt` を作るが、**全射だとは主張しない**
  (`K̄ ⊊ ℂ_K` は `K̄` が完備でないことによる。それ自体は本ファイルの射程外)。
  対応して、`ClosureCompletion.lean` の `completeSpace_of_norm_le_one_iff` は
  `CompleteSpace A` を要求するので `𝒪_{K̄}` には**適用できない**——
  すなわち `CompleteSpace ↥(absClosureInt K)` は本ファイルにはない(意図的)。
* **作用がノルムを保つことを落とすと完備化へ延びない**。
  `uniformContinuousConstSMul_closure` を落とすと mathlib の
  `DistribMulAction M (Completion α)` の instance 条件が破れ、
  `closureCompletionMulSemiringAction` は**型検査を通らない**。
* ★★**自明な作用で型クラスを満たしていない**ことの証拠は 2 段構えである:
  1. `smul_coe` 系(`coe_smul_absClosureInt`・`smul_coe_closureCompletion`)が
     作用を `σ • x = σ x`(共役そのもの)と**同定**する。`rfl` で閉じる。
  2. `FaithfulSMul K.absGal ↥(absClosureInt K)` と
     `FaithfulSMul K.absGal (closureCompletion K)`。自明な作用
     `σ • x := x` は `Γ_K` が非自明なとき **`FaithfulSMul` を満たさない**
     (`Found/PGC/QpNonAbelian.lean` により `Γ_{ℚ_p}` は非可換、特に非自明)。
     ★`𝒪_{K̄}` 側の忠実性は自明ではない——`K̄` 上の忠実性から
     「係数 `p^n` を掛けて単位球に落とす」ことで初めて降りる
     (`eq_of_smul_eq_smul_of_scaling`)。
* `DistribMulAction`(原文の「Γ_K-加群」)より強い `MulSemiringAction` まで
  出している。下流(`prop_2_2` を消費する側)が環構造を使えるようにするため。

## 逸脱の記録

* 原文は `K̄^∧` を「`𝒪_{K̄}` の p 進完備化を `ℤ_p` 上 `ℚ_p` とテンソルしたもの」と
  書くが、本ファイルは `K^al` のスペクトルノルムによる完備化として構成した。
  一致の根拠は上記(`mem_pow_p_mul_absClosureInt`)。★テンソル積による記述は
  **写していない**(同じ体を作る 2 通りの言い方であり、後続が使うのは体としての
  `ℂ_K` だけである)。
* `prop_2_2` 本体(群論的回復)は**証明していない**。本ファイルは主語の構成のみ。
* `Skeleton/PGC/Section2.lean` は**書き換えていない**(import のみ)。
  `IntKbarRecoverable` / `CompKbarRecoverable` は、`prop_2_2` が要求する
  4 つの型クラスがすべて我々の構成で満たされることを**型で**確認するために
  `RecoverableAsAddModule` を具体化した `Prop` である(証明はしていない)。
* `Valued K.closure NNReal` を instance に**していない**——`ClosureCompletion.lean`
  の設計判断(完備化に `Valued` が二重に付くのを防ぐ)をそのまま踏襲し、
  `closureValuation` という項として持つ。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal

/-! ## §1 抽象核 —— 位相環の完備化と、その上のモノイド作用

★ここには `PAdicLocalField` も代数閉包も分岐も付値も Galois も出てこない。 -/

section AbstractCore

open UniformSpace

variable {M X : Type*}

/-- **抽象核 1(a)**。ノルムを保つ加法的作用は等長写像で与えられる。 -/
theorem isometry_smul_of_norm_smul [Monoid M] [SeminormedAddCommGroup X] [DistribMulAction M X]
    (h : ∀ (c : M) (x : X), ‖c • x‖ = ‖x‖) (c : M) : Isometry (fun x : X => c • x) :=
  Isometry.of_dist_eq fun x y => by simp only [dist_eq_norm, ← smul_sub, h]

/-- **抽象核 1(b)**。したがってその作用は一様連続。 -/
theorem uniformContinuousConstSMul_of_norm_smul [Monoid M] [SeminormedAddCommGroup X]
    [DistribMulAction M X] (h : ∀ (c : M) (x : X), ‖c • x‖ = ‖x‖) :
    UniformContinuousConstSMul M X :=
  ⟨fun c => (isometry_smul_of_norm_smul h c).uniformContinuous⟩

/-- **抽象核 2**——★**一様連続な環作用は完備化へ一意に延びる**。

mathlib には `DistribMulAction M (Completion α)`
(`Topology/Algebra/GroupCompletion.lean:165`)しかなく、
`MulSemiringAction`(乗法まで込み)は無い。ここではその `DistribMulAction`
インスタンスを `with` でそのまま受け継ぐので、**菱形は生じない**。

`smul_mul` は稠密性(`Completion.induction_on₂`)と、両辺の連続性
(`ContinuousConstSMul` + `IsTopologicalRing`)だけで出る。 -/
@[implicit_reducible] noncomputable def completionMulSemiringAction {A : Type*} [Monoid M]
    [Ring A] [UniformSpace A] [IsTopologicalRing A] [IsUniformAddGroup A]
    [MulSemiringAction M A] [UniformContinuousConstSMul M A] :
    MulSemiringAction M (Completion A) :=
  { (inferInstance : DistribMulAction M (Completion A)) with
    smul_one := fun c => by
      rw [← Completion.coe_one, ← Completion.coe_smul, smul_one]
    smul_mul := fun c x y => by
      have hc : Continuous (fun z : Completion A => c • z) := continuous_const_smul c
      refine Completion.induction_on₂ x y (isClosed_eq ?_ ?_) ?_
      · exact hc.comp (continuous_fst.mul continuous_snd)
      · exact (hc.comp continuous_fst).mul (hc.comp continuous_snd)
      · intro a b
        rw [← Completion.coe_mul, ← Completion.coe_smul, ← Completion.coe_smul,
            ← Completion.coe_smul, ← Completion.coe_mul, smul_mul'] }

/-- **抽象核 3**——完備化へ延ばした作用もノルムを保つ。 -/
theorem norm_smul_completion [Monoid M] [SeminormedAddCommGroup X] [DistribMulAction M X]
    [UniformContinuousConstSMul M X] (h : ∀ (c : M) (x : X), ‖c • x‖ = ‖x‖) (c : M)
    (z : Completion X) : ‖c • z‖ = ‖z‖ := by
  refine Completion.induction_on z (isClosed_eq ?_ ?_) ?_
  · exact continuous_norm.comp (continuous_const_smul c)
  · exact continuous_norm
  · intro a
    rw [← Completion.coe_smul, Completion.norm_coe, Completion.norm_coe, h]

/-- **抽象核 4**——忠実性は完備化へ上がる(`T0Space` なら埋め込みが単射だから)。 -/
theorem faithfulSMul_completion [Monoid M] [UniformSpace X] [T0Space X] [MulAction M X]
    [UniformContinuousConstSMul M X] [FaithfulSMul M X] : FaithfulSMul M (Completion X) where
  eq_of_smul_eq_smul {c d} h := by
    refine eq_of_smul_eq_smul (α := X) fun x => ?_
    have hx := h (x : Completion X)
    rwa [← Completion.coe_smul, ← Completion.coe_smul, Completion.coe_inj] at hx

/-- **抽象核 5**——作用で不変な部分環には環作用が制限される。

`SubringClass` で書いてあるので `Subring` にも `ValuationSubring` にも当たる。 -/
@[implicit_reducible] def mulSemiringActionOfSubringClass {A S : Type*} [Monoid M] [Ring A]
    [SetLike S A] [SubringClass S A] [MulSemiringAction M A] (s : S)
    (hs : ∀ (c : M) (x : A), x ∈ s → c • x ∈ s) : MulSemiringAction M ↥s where
  smul c x := ⟨c • (x : A), hs c x x.2⟩
  one_smul x := Subtype.ext (one_smul M (x : A))
  mul_smul c d x := Subtype.ext (mul_smul c d (x : A))
  smul_zero c := Subtype.ext (smul_zero c)
  smul_add c x y := Subtype.ext (smul_add c (x : A) (y : A))
  smul_one c := Subtype.ext (smul_one c)
  smul_mul c x y := Subtype.ext (smul_mul' c (x : A) (y : A))

/-- **抽象核 6**——忠実性が**部分集合へ降りる**ための十分条件。

体 `A` 上の忠実な環作用に対し、「どの `x` も、作用で固定される `0` でない係数
`a` を掛ければ `s` に入る」なら、`s` の上だけで一致する 2 元は等しい。

★これが無いと「`𝒪_{K̄}` 上の作用が自明でない」ことが言えない
——`K̄` 上の忠実性だけでは単位球の中身について何も言えないからである。 -/
theorem eq_of_smul_eq_smul_of_scaling {A : Type*} [Monoid M] [Field A] [MulSemiringAction M A]
    [FaithfulSMul M A] (s : Set A)
    (hscale : ∀ x : A, ∃ a : A, a ≠ 0 ∧ (∀ c : M, c • a = a) ∧ a * x ∈ s)
    {c d : M} (h : ∀ x ∈ s, c • x = d • x) : c = d := by
  refine eq_of_smul_eq_smul (α := A) fun x => ?_
  obtain ⟨a, ha0, hafix, hax⟩ := hscale x
  have hcd := h _ hax
  rw [smul_mul', smul_mul', hafix c, hafix d] at hcd
  exact mul_left_cancel₀ ha0 hcd

end AbstractCore

/-! ## §2 `Γ_K` は `K^al` に等長な環自己同型として作用する -/

variable {p : ℕ} [Fact p.Prime]

/-- `Γ_K` の `K^al` への作用は共役そのもの(`AlgEquiv.applyMulSemiringAction`)。

★これが「自明な作用ではない」ことの第 1 の証拠。 -/
theorem smul_closure_def (K : PAdicLocalField p) (σ : K.absGal) (x : K.closure) :
    σ • x = σ x := rfl

/-- ★**Galois 作用はスペクトルノルムを保つ**。

mathlib の `spectralNorm_eq_of_equiv`
(`Analysis/Normed/Unbundled/SpectralNorm.lean:449`)がそのまま使える。
★実測: この補題は `FiniteNormal` セクションの**前**にあり、
`[NormedField K] [Field L] [Algebra K L]` しか要求しない
(有限次性も正規性も不要)。証明は `minpoly.algEquiv_eq` 一行である。 -/
theorem norm_smul_closure (K : PAdicLocalField p) (σ : K.absGal) (x : K.closure) :
    ‖σ • x‖ = ‖x‖ := by
  rw [smul_closure_def, norm_eq_spectralNorm_closure, norm_eq_spectralNorm_closure,
    ← spectralNorm_eq_of_equiv σ x]

/-- `Γ_K` は `K^al` の等長変換として作用する。 -/
theorem isometry_smul_closure (K : PAdicLocalField p) (σ : K.absGal) :
    Isometry (fun x : K.closure => σ • x) :=
  isometry_smul_of_norm_smul (norm_smul_closure K) σ

/-- したがって作用は一様連続——★これが完備化へ延ばすための唯一の入力である。 -/
scoped instance uniformContinuousConstSMul_closure (K : PAdicLocalField p) :
    UniformContinuousConstSMul K.absGal K.closure :=
  uniformContinuousConstSMul_of_norm_smul (norm_smul_closure K)

/-! ### `p` のノルム(`𝒪_{K̄}` へ落とすときの係数) -/

/-- `‖p‖ = 1/p`(`K` 上)。 -/
theorem norm_natCast_p_carrier (K : PAdicLocalField p) :
    ‖((p : ℕ) : K.carrier)‖ = ((p : ℕ) : ℝ)⁻¹ := by
  rw [show ((p : ℕ) : K.carrier) = algebraMap ℚ_[p] K.carrier ((p : ℕ) : ℚ_[p]) from
        (map_natCast _ p).symm, norm_algebraMap, Padic.norm_p]

/-- `‖p‖ = 1/p`(`K^al` 上)。スペクトルノルムは `K` のノルムを延長するから。 -/
theorem norm_natCast_p_closure (K : PAdicLocalField p) :
    ‖((p : ℕ) : K.closure)‖ = ((p : ℕ) : ℝ)⁻¹ := by
  rw [show ((p : ℕ) : K.closure) = algebraMap K.carrier K.closure ((p : ℕ) : K.carrier) from
        (map_natCast _ p).symm, norm_algebraMap_closure, norm_natCast_p_carrier]

theorem norm_natCast_p_closure_pos (K : PAdicLocalField p) :
    0 < ‖((p : ℕ) : K.closure)‖ := by
  rw [norm_natCast_p_closure]
  exact inv_pos.mpr (by exact_mod_cast (Fact.out : p.Prime).pos)

theorem norm_natCast_p_closure_lt_one (K : PAdicLocalField p) :
    ‖((p : ℕ) : K.closure)‖ < 1 := by
  rw [norm_natCast_p_closure]
  exact inv_lt_one_of_one_lt₀ (by exact_mod_cast (Fact.out : p.Prime).one_lt)

theorem natCast_p_closure_ne_zero (K : PAdicLocalField p) : ((p : ℕ) : K.closure) ≠ 0 :=
  norm_pos_iff.mp (norm_natCast_p_closure_pos K)

/-! ## §3 `𝒪_{K̄}` —— `K^al` の付値環 -/

/-- `K^al` 上のスペクトルノルムが定める付値。

★`Valued K.closure NNReal` を **instance にしていない**——`ClosureCompletion.lean`
の設計判断(完備化に第二の `Valued` が付くのを防ぐ)をそのまま踏襲する。 -/
noncomputable def closureValuation (K : PAdicLocalField p) : Valuation K.closure NNReal :=
  letI : Valued K.closure NNReal := NormedField.toValued
  Valued.v

theorem closureValuation_apply (K : PAdicLocalField p) (x : K.closure) :
    closureValuation K x = ‖x‖₊ := rfl

/-- **`𝒪_{K̄}`** —— 原文の `O[scr]_{K[bar]}`。`K^al` のノルム `≤ 1` の元のなす環。

★`ClosureCompletion.lean` の `closureCompletionInt`(`𝒪_{ℂ_K}`)とは**別物**。
こちらは完備化する**前**の `K^al` の付値環である。 -/
noncomputable def absClosureInt (K : PAdicLocalField p) : ValuationSubring K.closure :=
  (closureValuation K).valuationSubring

@[simp] theorem mem_absClosureInt (K : PAdicLocalField p) (x : K.closure) :
    x ∈ absClosureInt K ↔ ‖x‖ ≤ 1 := by
  rw [absClosureInt, Valuation.mem_valuationSubring_iff, closureValuation_apply]
  exact ⟨fun h => by exact_mod_cast h, fun h => by exact_mod_cast h⟩

theorem norm_le_one_absClosureInt (K : PAdicLocalField p) (w : ↥(absClosureInt K)) :
    ‖(w : K.closure)‖ ≤ 1 := (mem_absClosureInt K _).mp w.2

/-- **`𝒪_{K̄}` は `K^al` の中で閉**。§1(`ClosureCompletion.lean`)の抽象核の代入。 -/
theorem isClosed_absClosureInt (K : PAdicLocalField p) :
    IsClosed ((↑(absClosureInt K) : Set K.closure)) :=
  isClosed_of_norm_le_one_iff (absClosureInt K) (mem_absClosureInt K)

/-- **`𝒪_{K̄}` は線形位相環**——原文の言う「p 進位相」がここに入っている。

★`ClosureCompletion.lean` の `completeSpace_of_norm_le_one_iff` は `CompleteSpace A`
を要求するので `𝒪_{K̄}` には**当たらない**(`K^al` は完備でない)。
`CompleteSpace ↥(absClosureInt K)` は**意図的に作っていない**。 -/
theorem isLinearTopology_absClosureInt (K : PAdicLocalField p) :
    IsLinearTopology ↥(absClosureInt K) ↥(absClosureInt K) :=
  isLinearTopology_of_norm_le_one (absClosureInt K)
    (fun _ hx => (mem_absClosureInt K _).mp hx)

/-- ★★**原文の「p 進完備化」= ノルムによる完備化 であることの根拠**。

`p^n · 𝒪_{K̄}` はちょうどノルム球 `{x : ‖x‖ ≤ ‖p‖^n}` である。したがって
`𝒪_{K̄}` の p 進位相はスペクトルノルムの位相と一致し、
「`𝒪_{K̄}` の p 進完備化」は `ClosureCompletion.lean` の `closureCompletionInt`、
その商体は `closureCompletion K = ℂ_K` にほかならない。

★値群の可除性は使っていない——`‖x/p^n‖ = ‖x‖/‖p‖^n` だけである。 -/
theorem mem_pow_p_mul_absClosureInt (K : PAdicLocalField p) (n : ℕ) (x : K.closure) :
    (∃ y : K.closure, y ∈ absClosureInt K ∧ x = ((p : ℕ) : K.closure) ^ n * y)
      ↔ ‖x‖ ≤ ‖((p : ℕ) : K.closure)‖ ^ n := by
  have hp0 : ((p : ℕ) : K.closure) ≠ 0 := natCast_p_closure_ne_zero K
  have hpn0 : ((p : ℕ) : K.closure) ^ n ≠ 0 := pow_ne_zero _ hp0
  constructor
  · rintro ⟨y, hy, rfl⟩
    rw [norm_mul, norm_pow]
    calc ‖((p : ℕ) : K.closure)‖ ^ n * ‖y‖
        ≤ ‖((p : ℕ) : K.closure)‖ ^ n * 1 :=
          mul_le_mul_of_nonneg_left ((mem_absClosureInt K _).mp hy) (by positivity)
      _ = ‖((p : ℕ) : K.closure)‖ ^ n := mul_one _
  · intro hx
    refine ⟨x / ((p : ℕ) : K.closure) ^ n, ?_, by field_simp⟩
    rw [mem_absClosureInt, norm_div, norm_pow, div_le_one (by positivity)]
    exact hx

/-! ### `Γ_K` の `𝒪_{K̄}` への作用 -/

theorem smul_mem_absClosureInt (K : PAdicLocalField p) (σ : K.absGal) {x : K.closure}
    (hx : x ∈ absClosureInt K) : σ • x ∈ absClosureInt K := by
  rw [mem_absClosureInt] at hx ⊢
  rwa [norm_smul_closure]

/-- ★★**`𝒪_{K̄}` は `Γ_K`-環**。§1 の抽象核 5 の代入。

`DistribMulAction`(原文の「Γ_K-加群」)より強い `MulSemiringAction` を出している。 -/
noncomputable scoped instance absClosureIntMulSemiringAction (K : PAdicLocalField p) :
    MulSemiringAction K.absGal ↥(absClosureInt K) :=
  mulSemiringActionOfSubringClass (absClosureInt K)
    (fun σ _ hx => smul_mem_absClosureInt K σ hx)

/-- ★作用は共役そのもの(`rfl`)——自明な作用ではないことの第 1 の証拠。 -/
@[simp] theorem coe_smul_absClosureInt (K : PAdicLocalField p) (σ : K.absGal)
    (w : ↥(absClosureInt K)) : ((σ • w : ↥(absClosureInt K)) : K.closure) = σ w := rfl

theorem norm_smul_absClosureInt (K : PAdicLocalField p) (σ : K.absGal)
    (w : ↥(absClosureInt K)) : ‖((σ • w : ↥(absClosureInt K)) : K.closure)‖ = ‖(w : K.closure)‖ :=
  norm_smul_closure K σ (w : K.closure)

/-- 任意の `x ∈ K^al` は、`Γ_K` が固定する `0` でない係数 `p^n` を掛ければ
`𝒪_{K̄}` に入る。§1 抽象核 6 の仮説を満たすための具体層。 -/
theorem exists_smul_invariant_mul_mem_absClosureInt (K : PAdicLocalField p) (x : K.closure) :
    ∃ a : K.closure, a ≠ 0 ∧ (∀ σ : K.absGal, σ • a = a) ∧
      a * x ∈ (absClosureInt K : Set K.closure) := by
  by_cases hx : x = 0
  · refine ⟨1, one_ne_zero, fun σ => by rw [smul_closure_def, map_one], ?_⟩
    simp [hx]
  · have hx0 : 0 < ‖x‖ := norm_pos_iff.mpr hx
    obtain ⟨n, hn⟩ := exists_pow_lt_of_lt_one (inv_pos.mpr hx0) (norm_natCast_p_closure_lt_one K)
    refine ⟨((p : ℕ) : K.closure) ^ n, pow_ne_zero _ (natCast_p_closure_ne_zero K), ?_, ?_⟩
    · intro σ
      rw [smul_closure_def, map_pow, map_natCast]
    · have : ‖((p : ℕ) : K.closure) ^ n * x‖ < 1 := by
        rw [norm_mul, norm_pow]
        calc ‖((p : ℕ) : K.closure)‖ ^ n * ‖x‖ < ‖x‖⁻¹ * ‖x‖ := by
              exact mul_lt_mul_of_pos_right hn hx0
          _ = 1 := inv_mul_cancel₀ (ne_of_gt hx0)
      simpa using this.le

/-- ★★**`𝒪_{K̄}` への作用は忠実**——自明な作用ではないことの第 2 の証拠。

自明な作用 `σ • w := w` は `Γ_K` が非自明なら `FaithfulSMul` を満たさない。
`Γ_{ℚ_p}` が非自明(実際は非可換)であることは
`Found/PGC/QpNonAbelian.lean` に在る。 -/
scoped instance absClosureIntFaithfulSMul (K : PAdicLocalField p) :
    FaithfulSMul K.absGal ↥(absClosureInt K) where
  eq_of_smul_eq_smul {σ τ} h := by
    refine eq_of_smul_eq_smul_of_scaling (A := K.closure) ((absClosureInt K : Set K.closure))
      (exists_smul_invariant_mul_mem_absClosureInt K) ?_
    intro x hx
    exact congrArg Subtype.val (h ⟨x, hx⟩)

/-! ## §4 `K̄^∧ = ℂ_K` —— 完備化への作用の延長 -/

/-- ★★**`ℂ_K` は `Γ_K`-環**。§1 の抽象核 2(完備化への延長)の代入。

`DistribMulAction K.absGal (closureCompletion K)` は mathlib の
`Topology/Algebra/GroupCompletion.lean:165` から自動で付くが、
乗法まで込みの `MulSemiringAction` は mathlib に無いので抽象核 2 で作る。 -/
noncomputable scoped instance closureCompletionMulSemiringAction (K : PAdicLocalField p) :
    MulSemiringAction K.absGal (closureCompletion K) :=
  completionMulSemiringAction

/-- 埋め込み `K^al ↪ ℂ_K` は `Γ_K` 同変。 -/
@[simp] theorem smul_coe_closureCompletion (K : PAdicLocalField p) (σ : K.absGal)
    (x : K.closure) :
    σ • ((x : closureCompletion K)) = ((σ • x : K.closure) : closureCompletion K) :=
  (UniformSpace.Completion.coe_smul σ x).symm

/-- ★延長した作用もノルムを保つ(抽象核 3 の代入)。 -/
theorem norm_smul_closureCompletion (K : PAdicLocalField p) (σ : K.absGal)
    (z : closureCompletion K) : ‖σ • z‖ = ‖z‖ :=
  norm_smul_completion (norm_smul_closure K) σ z

theorem isometry_smul_closureCompletion (K : PAdicLocalField p) (σ : K.absGal) :
    Isometry (fun z : closureCompletion K => σ • z) :=
  isometry_smul_of_norm_smul (norm_smul_closureCompletion K) σ

theorem continuous_smul_closureCompletion (K : PAdicLocalField p) (σ : K.absGal) :
    Continuous (fun z : closureCompletion K => σ • z) :=
  (isometry_smul_closureCompletion K σ).continuous

/-- ★★**`ℂ_K` への作用も忠実**(抽象核 4 の代入)——自明な作用ではない。 -/
scoped instance closureCompletionFaithfulSMul (K : PAdicLocalField p) :
    FaithfulSMul K.absGal (closureCompletion K) :=
  faithfulSMul_completion

/-- `𝒪_{ℂ_K}` も `Γ_K` で不変。 -/
theorem smul_mem_closureCompletionInt (K : PAdicLocalField p) (σ : K.absGal)
    {z : closureCompletion K} (hz : z ∈ closureCompletionInt K) :
    σ • z ∈ closureCompletionInt K := by
  rw [mem_closureCompletionInt] at hz ⊢
  rwa [norm_smul_closureCompletion]

/-- `𝒪_{ℂ_K}` も `Γ_K`-環。 -/
noncomputable scoped instance closureCompletionIntMulSemiringAction (K : PAdicLocalField p) :
    MulSemiringAction K.absGal ↥(closureCompletionInt K) :=
  mulSemiringActionOfSubringClass (closureCompletionInt K)
    (fun σ _ hz => smul_mem_closureCompletionInt K σ hz)

@[simp] theorem coe_smul_closureCompletionInt (K : PAdicLocalField p) (σ : K.absGal)
    (w : ↥(closureCompletionInt K)) :
    ((σ • w : ↥(closureCompletionInt K)) : closureCompletion K) = σ • (w : closureCompletion K) :=
  rfl

/-- ★**`𝒪_{K̄} ↪ 𝒪_{ℂ_K}`** —— 単射・等長・`Γ_K` 同変の環準同型。

★★**全射だとは主張しない**(`K̄ ⊊ ℂ_K`)。2 つの整数環が別物であることの
明示的な記録である。 -/
noncomputable def absClosureIntToClosureCompletionInt (K : PAdicLocalField p) :
    ↥(absClosureInt K) →+* ↥(closureCompletionInt K) :=
  RingHom.codRestrict ((closureCompletionCoe K).comp (absClosureInt K).subtype)
    (closureCompletionInt K)
    (fun w => by
      show ((w : K.closure) : closureCompletion K) ∈ closureCompletionInt K
      rw [mem_closureCompletionInt, norm_coe_closureCompletion]
      exact norm_le_one_absClosureInt K w)

@[simp] theorem coe_absClosureIntToClosureCompletionInt (K : PAdicLocalField p)
    (w : ↥(absClosureInt K)) :
    ((absClosureIntToClosureCompletionInt K w : ↥(closureCompletionInt K)) :
        closureCompletion K) = ((w : K.closure) : closureCompletion K) := rfl

theorem absClosureIntToClosureCompletionInt_injective (K : PAdicLocalField p) :
    Function.Injective (absClosureIntToClosureCompletionInt K) := by
  intro w₁ w₂ hw
  have := congrArg (fun v : ↥(closureCompletionInt K) => (v : closureCompletion K)) hw
  simp only [coe_absClosureIntToClosureCompletionInt] at this
  exact Subtype.ext (closureCompletionCoe_injective K this)

theorem absClosureIntToClosureCompletionInt_smul (K : PAdicLocalField p) (σ : K.absGal)
    (w : ↥(absClosureInt K)) :
    absClosureIntToClosureCompletionInt K (σ • w) = σ • absClosureIntToClosureCompletionInt K w :=
  Subtype.ext (by
    rw [coe_absClosureIntToClosureCompletionInt, coe_smul_closureCompletionInt,
      coe_absClosureIntToClosureCompletionInt, coe_smul_absClosureInt, smul_coe_closureCompletion,
      smul_closure_def])

/-! ## §5 `prop_2_2` が要求する型クラスをすべて満たすことの確認

★`Skeleton/PGC/Section2.lean` は**書き換えていない**。ここでは
`RecoverableAsAddModule` を我々の構成で具体化した `Prop` を**書き下す**ことで、
`prop_2_2` の 4 つの instance 引数

```
[∀ K, AddCommGroup (IntKbar K)] [∀ K, DistribMulAction K.absGal (IntKbar K)]
[∀ K, AddCommGroup (CompKbar K)] [∀ K, DistribMulAction K.absGal (CompKbar K)]
```

がすべて合成できることを**型で**検査する(命題自体は証明していない)。 -/

/-- 原文の `O[scr]_{K[bar]}` の実体。 -/
abbrev IntKbar (K : PAdicLocalField p) : Type := ↥(absClosureInt K)

/-- 原文の `K[bar]^∧` の実体。 -/
abbrev CompKbar (K : PAdicLocalField p) : Type := closureCompletion K

/-- `prop_2_2` の `IntKbar` 側を我々の構成で具体化した `Prop`。

★これが**型検査を通ること自体**が、`AddCommGroup (𝒪_{K̄})` と
`DistribMulAction Γ_K (𝒪_{K̄})` が揃っていることの証明である。 -/
def IntKbarRecoverable : Prop :=
  RecoverableAsAddModule (p := p) (fun K => IntKbar K)

def IntKbarRecoverable.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- `prop_2_2` の `CompKbar` 側を我々の構成で具体化した `Prop`。 -/
def CompKbarRecoverable : Prop :=
  RecoverableAsAddModule (p := p) (fun K => CompKbar K)

def CompKbarRecoverable.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-! ### 型クラスが揃っていることの明示的な証拠 -/

/-- `𝒪_{K̄}` は加法群。 -/
@[implicit_reducible] noncomputable def intKbarAddCommGroup (K : PAdicLocalField p) :
    AddCommGroup (IntKbar K) := inferInstance

/-- `𝒪_{K̄}` は `Γ_K`-加群。 -/
@[implicit_reducible] noncomputable def intKbarDistribMulAction (K : PAdicLocalField p) :
    DistribMulAction K.absGal (IntKbar K) := inferInstance

/-- `ℂ_K` は加法群。 -/
@[implicit_reducible] noncomputable def compKbarAddCommGroup (K : PAdicLocalField p) :
    AddCommGroup (CompKbar K) := inferInstance

/-- `ℂ_K` は `Γ_K`-加群。 -/
@[implicit_reducible] noncomputable def compKbarDistribMulAction (K : PAdicLocalField p) :
    DistribMulAction K.absGal (CompKbar K) := inferInstance

/-! ### 退化していないことの型による記録 -/

/-- ★★`𝒪_{K̄}` への作用は**自明ではない**: `σ` が `𝒪_{K̄}` 上恒等なら `σ = 1`。 -/
theorem eq_one_of_smul_eq_self_absClosureInt (K : PAdicLocalField p) {σ : K.absGal}
    (h : ∀ w : ↥(absClosureInt K), σ • w = w) : σ = 1 :=
  eq_of_smul_eq_smul (α := ↥(absClosureInt K)) fun w => by rw [h w, one_smul]

/-- ★★`ℂ_K` への作用は**自明ではない**: `σ` が `ℂ_K` 上恒等なら `σ = 1`。 -/
theorem eq_one_of_smul_eq_self_closureCompletion (K : PAdicLocalField p) {σ : K.absGal}
    (h : ∀ z : closureCompletion K, σ • z = z) : σ = 1 :=
  eq_of_smul_eq_smul (α := closureCompletion K) fun z => by rw [h z, one_smul]

end ABC3.Found.PGC
