import ABC3.Found.GenEll.Sl2Level
import Mathlib.NumberTheory.Padics.RingHoms
import Mathlib.Topology.Algebra.Group.Matrix
import Mathlib.Topology.MetricSpace.Basic
import Mathlib.Analysis.SpecificLimits.Basic

/-!
# [GenEll] Lemma 3.1, (iv) の**位相段**(段 B)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.14。

原文 (GenEll p.14):
> (iv) Let J ⊆ GL2(Zl) be a closed subgroup whose image HJ in GL2(Fl) contains

## ★★これが「位相が要る最後の 1 歩」である

`Sl2Level.lean` が段 (A)(有限段、位相不要)を終えている:

> `l ≥ 5` 素数、`H ≤ SL₂(ℤ/l^{k+1})` が `SL₂(𝔽_l)` へ全射なら `H = ⊤`。

★本ファイルはそれを `ℤ_l` へ持ち上げる:

> **閉部分群 `J ≤ SL₂(ℤ_l)` が `SL₂(𝔽_l)` へ全射なら `J = ⊤`。**

これが原文の引く **[Serre] Chapter IV, §3.4, Lemma 3** そのものである
(その [Serre] は `0_Source` に無い)。

## ★機構は「点列を作って収束させる」だけ

各 `n` について、段 (A) から `mod l^n` の像は全体である。
ゆえに `g ∈ SL₂(ℤ_l)` に対し `j_n ∈ J` を `j_n ≡ g (mod l^n)` に取れる。
すると各成分で `‖j_n − g‖ ≤ l^{-n} → 0` なので **`j_n → g`**。
`J` は閉だから `g ∈ J`。

★`PadicInt.ker_toZModPow`(核が `span {l^n}`)と
`PadicInt.norm_le_pow_iff_mem_span_pow`(そこから距離評価)がそのまま効く。
-/

namespace ABC3.Found.GenEll

open Matrix
open scoped MatrixGroups

section Bridge

variable (l : ℕ) [Fact (Nat.Prime l)]

/-- `SL₂(ℤ_l) → SL₂(ℤ/l^n)` の還元。 -/
noncomputable def redPadicPow (n : ℕ) : SL(2, ℤ_[l]) →* SL(2, ZMod (l ^ n)) :=
  Matrix.SpecialLinearGroup.map (PadicInt.toZModPow n)

/-- ★`ZMod (l^n) → ZMod l` と `toZModPow n` の合成は `toZMod` である。

★**証明は「`x = appr x 1 + l·y` と書く」だけ**——両辺とも
`(appr x 1 : ZMod l)` になる(`l` は `ZMod l` で消えるから)。 -/
theorem castHom_toZModPow_eq_toZMod {n : ℕ} (hn : n ≠ 0) (x : ℤ_[l]) :
    ZMod.castHom (dvd_pow_self l hn) (ZMod l) (PadicInt.toZModPow n x)
      = PadicInt.toZMod x := by
  obtain ⟨y, hy⟩ := Ideal.mem_span_singleton.1 (PadicInt.appr_spec 1 x)
  have hx : x = ((x.appr 1 : ℕ) : ℤ_[l]) + (l : ℤ_[l]) * y := by
    have : (l : ℤ_[l]) ^ 1 * y = (l : ℤ_[l]) * y := by ring
    rw [← this, ← hy]
    ring
  rw [hx]
  simp only [map_add, map_mul, map_natCast]
  rw [ZMod.natCast_self, zero_mul, add_zero, zero_mul, add_zero]

/-- ★`redPow ∘ redPadicPow = redPadic`。 -/
theorem redPow_comp_redPadicPow {n : ℕ} (hn : n ≠ 0) (j : SL(2, ℤ_[l])) :
    redPow l n hn (redPadicPow l n j) = redPadic l j := by
  apply Subtype.ext
  ext i j'
  exact castHom_toZModPow_eq_toZMod l hn _

end Bridge

section Density

variable (l : ℕ) [Fact (Nat.Prime l)]

/-- ★各段の像は全体である(段 (A) を `ℤ_l` の像に適用する)。 -/
theorem redPadicPow_map_eq_top (h5 : 5 ≤ l) (J : Subgroup SL(2, ℤ_[l]))
    (hJ : ∀ g : SL(2, ZMod l), ∃ j ∈ J, redPadic l j = g) (k : ℕ) :
    J.map (redPadicPow l (k + 1)) = ⊤ := by
  refine subgroup_eq_top_of_redPow_surj l h5 k _ ?_
  intro g
  obtain ⟨j, hj, hjg⟩ := hJ g
  refine ⟨redPadicPow l (k + 1) j, ⟨j, hj, rfl⟩, ?_⟩
  rw [redPow_comp_redPadicPow l (Nat.succ_ne_zero k) j, hjg]

/-- ★`mod l^n` が一致すれば、成分の差のノルムは `l^{-n}` 以下。 -/
theorem norm_sub_le_of_redPadicPow_eq {n : ℕ} (x y : SL(2, ℤ_[l]))
    (h : redPadicPow l n x = redPadicPow l n y) (i j : Fin 2) :
    ‖(x : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j - (y : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j‖
      ≤ (l : ℝ) ^ (-(n : ℤ)) := by
  have hker : PadicInt.toZModPow n
      ((x : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j - (y : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j) = 0 := by
    have := congrFun (congrFun (congrArg
      (fun z : SL(2, ZMod (l ^ n)) => (z : Matrix (Fin 2) (Fin 2) (ZMod (l ^ n)))) h) i) j
    rw [map_sub]
    simpa [redPadicPow, Matrix.SpecialLinearGroup.map] using sub_eq_zero.2 this
  have hmem : (x : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j - (y : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j
      ∈ Ideal.span {(l : ℤ_[l]) ^ n} := by
    rw [← PadicInt.ker_toZModPow]
    exact hker
  exact (PadicInt.norm_le_pow_iff_mem_span_pow _ n).2 hmem

end Density

/-! ## ★段 B の本体 —— 点列を作って収束させる -/

section SerreLemma

variable (l : ℕ) [Fact (Nat.Prime l)]

theorem tendsto_inv_pow_atTop_zero (hl : 1 < l) :
    Filter.Tendsto (fun k : ℕ => ((l : ℝ)⁻¹) ^ (k + 1)) Filter.atTop (nhds 0) := by
  have h0 : (0 : ℝ) ≤ (l : ℝ)⁻¹ := by positivity
  have h1 : (l : ℝ)⁻¹ < 1 := by
    rw [inv_lt_one_iff₀]
    right
    exact_mod_cast hl
  exact (tendsto_pow_atTop_nhds_zero_of_lt_one h0 h1).comp (Filter.tendsto_add_atTop_nat 1)

/-- ★★**[Serre] Chapter IV, §3.4, Lemma 3**(`SL₂` 版)。

`l ≥ 5` 素数、`J ≤ SL₂(ℤ_l)` が**閉**で `SL₂(𝔽_l)` へ全射なら **`J = ⊤`**。

★原文はこれを [Serre] から引くが、`0_Source` に無いので**自分で証明した**。
段 (A)(`Sl2Level.lean`、有限群論)と本ファイル(位相)の合成である。 -/
theorem sl2_padic_eq_top_of_isClosed (h5 : 5 ≤ l) (J : Subgroup SL(2, ℤ_[l]))
    (hclosed : IsClosed (J : Set SL(2, ℤ_[l])))
    (hJ : ∀ g : SL(2, ZMod l), ∃ j ∈ J, redPadic l j = g) :
    J = ⊤ := by
  have hl : Nat.Prime l := Fact.out
  refine (Subgroup.eq_top_iff' J).2 fun g => ?_
  -- 各段で `g` に一致する `J` の元を選ぶ
  have hpick : ∀ k : ℕ, ∃ x, x ∈ J ∧ redPadicPow l (k + 1) x = redPadicPow l (k + 1) g := by
    intro k
    have htop := redPadicPow_map_eq_top l h5 J hJ k
    have hmem : redPadicPow l (k + 1) g ∈ J.map (redPadicPow l (k + 1)) := by
      rw [htop]; trivial
    obtain ⟨x, hx, hxg⟩ := hmem
    exact ⟨x, hx, hxg⟩
  choose j hjmem hjeq using hpick
  -- 各成分で `‖j k − g‖ ≤ l^{-(k+1)} → 0`
  have hentry : ∀ i j' : Fin 2, Filter.Tendsto
      (fun k => (j k : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j') Filter.atTop
      (nhds ((g : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j')) := by
    intro i j'
    rw [tendsto_iff_norm_sub_tendsto_zero]
    refine squeeze_zero (fun k => norm_nonneg _) (fun k => ?_)
      (tendsto_inv_pow_atTop_zero l hl.one_lt)
    have hb := norm_sub_le_of_redPadicPow_eq l (j k) g (hjeq k) i j'
    have hconv : ((l : ℝ)) ^ (-((k + 1 : ℕ) : ℤ)) = ((l : ℝ)⁻¹) ^ (k + 1) := by
      rw [zpow_neg, zpow_natCast, inv_pow]
    rwa [hconv] at hb
  -- 部分型・行列の位相に持ち上げる
  have hind : Topology.IsInducing
      (fun x : SL(2, ℤ_[l]) => (x : Matrix (Fin 2) (Fin 2) ℤ_[l])) := ⟨rfl⟩
  have htend : Filter.Tendsto j Filter.atTop (nhds g) := by
    rw [hind.tendsto_nhds_iff]
    exact tendsto_pi_nhds.2 fun i => tendsto_pi_nhds.2 fun j' => hentry i j'
  exact hclosed.mem_of_tendsto htend (Filter.Eventually.of_forall hjmem)

end SerreLemma

/-! ## ★おまけ —— `SL₂(𝔽_l)` の完全性

★`layer0-triage.md` の実測では「群の完全性は mathlib に**無い**」
(`perfect` の hit は**すべて perfect matching**(グラフ理論))。
`lemma_3_1_ii` からただちに出るので、ここで取っておく。 -/

section Perfect

variable (l : ℕ) [Fact (Nat.Prime l)]

/-- ★**`SL₂(𝔽_l)` は完全群**(`l ≥ 5`)。

`Lemma 3.1, (ii)`(「正規部分群で商が可換なら全体」)を交換子群に適用するだけである。
★これが `Lemma 3.1, (iv)` の `GL → SL` の還元(交換子群を取る段)で効く。 -/
theorem sl2_commutator_eq_top (h5 : 5 ≤ l) :
    commutator SL(2, ZMod l) = ⊤ := by
  refine lemma_3_1_ii l h5 (commutator SL(2, ZMod l)) ?_
  intro x y
  exact mul_comm (G := Abelianization SL(2, ZMod l)) x y

end Perfect

end ABC3.Found.GenEll
