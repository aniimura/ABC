import ABC3.Found.GenEll.Sl2Level
import Mathlib.NumberTheory.Padics.RingHoms
import Mathlib.Topology.Algebra.Group.Matrix
import Mathlib.Topology.MetricSpace.Basic
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.GroupTheory.Commutator.Basic
import Mathlib.GroupTheory.Abelianization.Defs

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

/-! ## ★★`Lemma 3.1, (iv)` —— `GL` から `SL` へ降ろす

★**位相的閉包は要らなかった。** `J.comap toGL` が既に閉である
(`toGL` は連続、`mathlib` の `continuous_toGL`)。

`GL₂(𝔽_l)` での像が `SL₂(𝔽_l)` を含むことから、**交換子を取って行列式を 1 に落とす**。
`⁅SL₂(𝔽_l), SL₂(𝔽_l)⁆ = SL₂(𝔽_l)`(上の `sl2_commutator_eq_top`)がここで効く。 -/

section GLtoSL

variable (l : ℕ) [Fact (Nat.Prime l)]

open Matrix.SpecialLinearGroup

/-- `GL₂(ℤ_l) → GL₂(𝔽_l)` の還元。★原文 (iv) の「image `H_J` in `GL₂(𝔽_l)`」。 -/
noncomputable def glRedPadic : GL (Fin 2) ℤ_[l] →* GL (Fin 2) (ZMod l) :=
  Units.map (PadicInt.toZMod.mapMatrix (m := Fin 2)).toMonoidHom

/-- ★還元と包含は可換。 -/
theorem glRedPadic_toGL (x : SL(2, ℤ_[l])) :
    glRedPadic l (toGL x) = toGL (redPadic l x) := by
  apply Units.ext
  rfl

/-- ★★**[GenEll] Lemma 3.1, (iv)**。

原文 (GenEll p.14):
> (iv) Let J ⊆ GL2(Zl) be a closed subgroup whose image HJ in GL2(Fl) contains

`GL₂(ℤ_l)` の**閉**部分群 `J` の `GL₂(𝔽_l)` における像が
`α = (1 1 / 0 1)` と**非上三角**な行列を含むなら、**`SL₂(ℤ_l) ⊆ J`**。

★原文は [Serre] Chapter IV, §3.4, Lemma 3 を引くが `0_Source` に無いので、
段 (A)(`Sl2Level.lean`)＋段 (B)(本ファイル)で**自分で証明した**。 -/
theorem lemma_3_1_iv (h5 : 5 ≤ l) (J : Subgroup (GL (Fin 2) ℤ_[l]))
    (hclosed : IsClosed (J : Set (GL (Fin 2) ℤ_[l])))
    (hα : (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l)) ∈ J.map (glRedPadic l))
    (hMex : ∃ M ∈ J.map (glRedPadic l), (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 1 0 ≠ 0)
    (g : SL(2, ℤ_[l])) : (toGL g : GL (Fin 2) ℤ_[l]) ∈ J := by
  classical
  -- (1) `J₁ := J.comap toGL` は閉部分群
  have hclosed₁ : IsClosed ((J.comap (toGL : SL(2, ℤ_[l]) →* GL (Fin 2) ℤ_[l]))
      : Set SL(2, ℤ_[l])) :=
    hclosed.preimage continuous_toGL
  -- (2) `H_J` は `toGL SL₂(𝔽_l)` を含む —— これは (iii)
  obtain ⟨M, hM1, hM2⟩ := hMex
  have hHJ : ∀ s : SL(2, ZMod l),
      (toGL s : GL (Fin 2) (ZMod l)) ∈ J.map (glRedPadic l) :=
    fun s => lemma_3_1_iii l (J.map (glRedPadic l)) hα M hM1 hM2 s
  -- (3) `J₁` の `mod l` 像は `SL₂(𝔽_l)` 全体
  have hcomm : commutator SL(2, ZMod l)
      ≤ (J.comap (toGL : SL(2, ℤ_[l]) →* GL (Fin 2) ℤ_[l])).map (redPadic l) := by
    rw [commutator_def]
    refine Subgroup.commutator_le.2 fun a _ b _ => ?_
    obtain ⟨α, hαJ, hαa⟩ := hHJ a
    obtain ⟨β, hβJ, hβb⟩ := hHJ b
    have hγJ : α * β * α⁻¹ * β⁻¹ ∈ J :=
      J.mul_mem (J.mul_mem (J.mul_mem hαJ hβJ) (J.inv_mem hαJ)) (J.inv_mem hβJ)
    have hdet : Matrix.GeneralLinearGroup.det (α * β * α⁻¹ * β⁻¹) = 1 := by
      simp only [map_mul, map_inv]
      rw [mul_comm (Matrix.GeneralLinearGroup.det α) (Matrix.GeneralLinearGroup.det β)]
      group
    obtain ⟨x, hx⟩ : α * β * α⁻¹ * β⁻¹
        ∈ Set.range (toGL : SL(2, ℤ_[l]) → GL (Fin 2) ℤ_[l]) := by
      rw [range_toGL]; exact hdet
    refine ⟨x, ?_, ?_⟩
    · show toGL x ∈ J
      rw [hx]; exact hγJ
    · apply toGL_injective
      show (toGL (redPadic l x) : GL (Fin 2) (ZMod l)) = toGL (a * b * a⁻¹ * b⁻¹)
      rw [← glRedPadic_toGL, hx]
      simp only [map_mul, map_inv, hαa, hβb]
  -- (4) Serre の補題を適用する
  have hsurj : ∀ s : SL(2, ZMod l),
      ∃ x ∈ J.comap (toGL : SL(2, ℤ_[l]) →* GL (Fin 2) ℤ_[l]), redPadic l x = s := by
    intro s
    have hmem : s ∈ (J.comap (toGL : SL(2, ℤ_[l]) →* GL (Fin 2) ℤ_[l])).map (redPadic l) := by
      refine hcomm ?_
      rw [sl2_commutator_eq_top l h5]
      trivial
    obtain ⟨x, hx, hxs⟩ := hmem
    exact ⟨x, hx, hxs⟩
  have htop := sl2_padic_eq_top_of_isClosed l h5 _ hclosed₁ hsurj
  have hgJ : g ∈ J.comap (toGL : SL(2, ℤ_[l]) →* GL (Fin 2) ℤ_[l]) := by
    rw [htop]; trivial
  exact hgJ

end GLtoSL

/-! ## ★★`Lemma 3.1` 全体(i)–(iv)

★(i)(ii)(iii) は `Lemma31.lean`、(iv) は本ファイル。**4 条すべてが `sorry` 無しで揃った。**
ここで初めて **条なしの `.src`** を付けられる(2 値規則)。 -/

section Whole

variable (l : ℕ) [Fact (Nat.Prime l)]

open Matrix.SpecialLinearGroup

/-- **[GenEll] Lemma 3.1**(The Structure of SL2)—— **4 条すべて**。

原文 (GenEll p.13):
> Lemma 3.1. (The Structure of SL2) Let l ≥ 5 be a prime number. Then:

(i) `SL₂(𝔽_l)` は `α = (1 1 / 0 1)` と `β = (1 0 / 1 1)` で生成される。
(ii) `SL₂(𝔽_l)` の正規部分群で商が可換なものは全体。
(iii) `GL₂(𝔽_l)` の部分群が `α` と**非上三角**な行列を含めば `SL₂(𝔽_l)` を含む。
(iv) `GL₂(ℤ_l)` の**閉**部分群 `J` の `GL₂(𝔽_l)` での像が同じ条件を満たせば `SL₂(ℤ_l) ⊆ J`。

★★**原文が (iv) で引く [Serre] Chapter IV, §3.4, Lemma 3 は `0_Source` に無い。**
段 (A)(`Sl2Level.lean`、有限群論 813 行)と段 (B)(本ファイル、位相)で
**自分で証明した**——これが本項目が「well-known facts の復習」でありながら
最も手数を要した理由である。 -/
theorem lemma_3_1 (hl : 5 ≤ l) :
    (Subgroup.closure ({upper (1 : ZMod l), lower (1 : ZMod l)} : Set SL(2, ZMod l)) = ⊤)
  ∧ (∀ (N : Subgroup SL(2, ZMod l)) [N.Normal],
        (∀ x y : SL(2, ZMod l) ⧸ N, x * y = y * x) → N = ⊤)
  ∧ (∀ H : Subgroup (GL (Fin 2) (ZMod l)),
        (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l)) ∈ H →
        (∃ M ∈ H, (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 1 0 ≠ 0) →
        ∀ g : SL(2, ZMod l), (toGL g : GL (Fin 2) (ZMod l)) ∈ H)
  ∧ (∀ J : Subgroup (GL (Fin 2) ℤ_[l]), IsClosed (J : Set (GL (Fin 2) ℤ_[l])) →
        (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l)) ∈ J.map (glRedPadic l) →
        (∃ M ∈ J.map (glRedPadic l), (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 1 0 ≠ 0) →
        ∀ g : SL(2, ℤ_[l]), (toGL g : GL (Fin 2) ℤ_[l]) ∈ J) := by
  refine ⟨lemma_3_1_i l, fun N _ habel => lemma_3_1_ii l hl N habel, ?_, ?_⟩
  · rintro H hα ⟨M, hM, hMnt⟩ g
    exact lemma_3_1_iii l H hα M hM hMnt g
  · intro J hclosed hα hMex g
    exact lemma_3_1_iv l hl J hclosed hα hMex g

/-! ### ★出典の紐付け -/

/-- ★**条なし** —— `Lemma 3.1` は 4 条すべてが揃った。 -/
def lemma_3_1.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 13, item := "Lemma 3.1",
    sectionId := "genell-lemma-3-1" }

def lemma_3_1_iv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 14, item := "Lemma 3.1, (iv)",
    sectionId := "genell-lemma-3-1-iv" }

end Whole

end ABC3.Found.GenEll
