import Mathlib.LinearAlgebra.Matrix.SpecialLinearGroup
import Mathlib.LinearAlgebra.Matrix.GeneralLinearGroup.Defs
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Field.ZMod
import Mathlib.RingTheory.LocalRing.Basic
import Mathlib.NumberTheory.Padics.PadicIntegers
import Mathlib.NumberTheory.Padics.RingHoms
import Mathlib.Tactic.LinearCombination
import Mathlib.GroupTheory.QuotientGroup.Basic
import ABC3.Meta.Claim

/-!
# [GenEll] Lemma 3.1 —— The Structure of SL₂

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.13(宣言)– p.14(4 条と証明)。**200 dpi 目視確認 2026-08-16**。

原文 (GenEll p.13):
> Lemma 3.1. (The Structure of SL_2) Let l ≥ 5 be a prime number. Then:

## ★なぜこの項目から始めるか

長期ゴールは `ResearchPaper/genell-goal.md`。要点:

- [GenEll] の同一論文内グラフで **`Lemma 3.1` は到達 0 件の真の葉**である(実測)。
- **`Theorem 3.8` の証明が使う群論は `Lemma 3.1, (iv)` ただ 1 つ**である(原文 p.20 を通読して確認)。
- ★**(iii) は [IUTchIV] が名指しで使う**。[IUTchIV] 物理 p.46 は
  「(P6) follows formally from (P2), (P4), and [GenEll], Lemma 3.1, (iii)」と書き、
  その (P6)(= l 捩れ点への Galois 表現の像が `SL₂(𝔽_l)` を含む)が揃うことで
  **[IUTchI] Definition 3.1 の初期 Θ-データが構成される**。

## ★この項目の記法上の危険(実測)

`pdftotext` の**既定モード**は 2×2 行列の組版を潰すとき**順序を入れ替える**。
原文 p.14 の (i) は `α ≝ (1 1 / 0 1)`(上三角)・`β ≝ (1 0 / 1 1)`(下三角)の順だが、
既定モードは `1 0 / 1 1` を先に出し、しかも α・β のギリシャ文字自体が空白に脱落する。
★**実害**: (iii) の仮定は目視では `α`(上三角)だが、既定モードの出力だけを読むと
**β(下三角)についての主張に読める**——命題が別物になる。
`-layout` モードでは順序が保たれる。記録は `papers.json` の `GenEll.notationNotes` と
`1_Structured/Arithmetic Elliptic Curves in General Position/lemma-3-1.html`。

## ★規模と、原文との道の違い

| 条 | 主張 | 原文の道 | 本ファイルの道 |
|---|---|---|---|
| (i) | `⟨α, β⟩ = SL₂(𝔽_l)` | ベクトル `v = (0,1)` の軌道 | ★**分解を明示**(`matrix_factor`) |
| (ii) | 可換な商は自明のみ | `ε·α·ε⁻¹·α⁻¹ = α^{λ²−1}` | **同じ** |
| (iii) | `α` と非上三角 1 つ ⟹ `SL₂ ⊆ H` | l-Sylow の個数 `n_H = l+1` | ★**共役で `β` を作る**(`gl_conj_key`) |

★**(i) で道を変えた理由**: 原文は `γ ≝ (0 1 / 1 0)` を経由するが、
**det γ = −1 なので `γ ∉ SL₂(𝔽_l)`** である(原文は注記していない)。
軌道の議論は「群の外の行列」を要し、Lean では型を跨ぐ。
一方 `c` が可逆なときの分解 `g = upper((a−1)e) · lower(c) · upper((d−1)e)`(`ce = 1`)は
**座標だけの計算**で閉じ、`c` が可逆でない場合は `lower 1 · g` に帰着する。

★★**その結果、(i) は体でなく局所環で成り立つ形になった**(`mem_of_upper_lower_mem`)。
使うのは「`c` が可逆」か「`a` が可逆」の二択だけで、それは**局所環の定義そのもの**である。
★これは **`Lemma 3.1, (iv)` の 2 本柱の 1 本**を先に取ったことになる——
`ℤ_[l]` は局所環なので `SL(2, ℤ_[l])` は基本行列で**代数的に**生成される
(`closure_elementary_padic`。**位相を一切使っていない**)。

★**(iii) で道を変えた理由**: 原文の道は
**`|GL₂(𝔽_l)|` の計算**と**正規化群 = 上三角全体の同定**の 2 つを前提にしており、
原文自身は「the general theory of Sylow subgroups」の一言で済ませている。
mathlib には Sylow 論(272 件)はあるが、この 2 つの具体形は無い
(探索範囲: `Mathlib` 全体を `Sylow` / `GL2` で grep、および
`Mathlib/LinearAlgebra/Matrix/SpecialLinearGroup.lean`・`GeneralLinearGroup/Defs.lean` の全宣言)。
代わりに `W ≝ upper(−a/c) · M` が左上成分 0 を持つことを使い、
`W · upper(u) · W⁻¹ = lower(1)` を**行列の等式 1 本**にした。

## ★`l ≥ 5` はどこで効くか(実測)

**(ii) だけである。** 原文の証明も `λ² ≠ 1` なる `λ` の存在にのみ使う。
(i) の原文の証明は `l ≥ 5` を一度も使わない。
負の対照は `lemma_3_1_ii_neg_control`(l = 2, 3 では `λ² = 1` が全ての `λ` で成り立つ)。

## ★(iv) を含めない理由

(iv) は **[Serre], Chapter IV, §3.4, Lemma 3** に依存し、
その [Serre](*Abelian l-adic Representations and Elliptic Curves*, Benjamin 1968)は
`0_Source` に**無い**。したがって (iv) は引用ではなく**自分で証明する**ことになる。

★**(iv) に残っている仕事の切り分け**(2026-08-16 に測った):

| 柱 | 内容 | 状態 |
|---|---|---|
| 1 | `SL(2, ℤ_[l])` が `upper t` / `lower t` で生成される | ★**済**(`closure_elementary_padic`。位相不要) |
| 2 | 閉部分群 `J` が `upper t` / `lower t` を**全部**含むこと | **未**。[Serre] の補題に当たる |

柱 2 には pro-l 群の位相が要る。mathlib の実測(2026-08-16、探索範囲つき):
`ProfiniteGrp` は **135 件**あるが、**pro-p 群は 2 件**で実質不在。
`Matrix.SpecialLinearGroup.map (f : R →+* S)` はあるので還元 `SL₂(ℤ_l) → SL₂(𝔽_l)` は書ける。
合同部分群の族は `SL(2, ℤ)` について `NumberTheory/ModularForms/CongruenceSubgroups.lean` に
あるだけで、`ℤ_[l]` については無い。

★**したがって本ファイルに `Lemma 3.1` 全体の `.src` は付けない**——
条つき `.src` は locator の記録であり、命題全体の完全実装の主張ではない。

## ★器具についての発見(2026-08-16、本ファイルの作業中に実測)

`first | ring | …` は**期待どおりに動かない**。mathlib の `ring` は失敗すると
`ring_nf` にフォールバックして**成功扱いになる**ため、`first` が次の候補へ進まない。
`ring1` を使えばよい。実際 `matrix_factor` の 1 ゴールがこれで 2 回失敗した。
-/

namespace ABC3.Found.GenEll

open Matrix
open scoped MatrixGroups

section Elementary

variable {R : Type*} [CommRing R]

/-- 原文 `α` の族: `!![1, t; 0, 1]`。原文の `α` は `upper 1`。 -/
def upper (t : R) : SL(2, R) := ⟨!![1, t; 0, 1], by simp [Matrix.det_fin_two_of]⟩

/-- 原文 `β` の族: `!![1, 0; t, 1]`。原文の `β` は `lower 1`。 -/
def lower (t : R) : SL(2, R) := ⟨!![1, 0; t, 1], by simp [Matrix.det_fin_two_of]⟩

@[simp] theorem coe_upper (t : R) :
    ((upper t : SL(2, R)) : Matrix (Fin 2) (Fin 2) R) = !![1, t; 0, 1] := rfl

@[simp] theorem coe_lower (t : R) :
    ((lower t : SL(2, R)) : Matrix (Fin 2) (Fin 2) R) = !![1, 0; t, 1] := rfl

theorem upper_mul_upper (s t : R) : upper s * upper t = upper (s + t) := by
  apply Subtype.ext
  simp [Matrix.SpecialLinearGroup.coe_mul, add_comm]

theorem lower_mul_lower (s t : R) : lower s * lower t = lower (s + t) := by
  apply Subtype.ext
  simp [Matrix.SpecialLinearGroup.coe_mul]

@[simp] theorem upper_zero : upper (0 : R) = 1 := by
  apply Subtype.ext
  simp [Matrix.one_fin_two]

@[simp] theorem lower_zero : lower (0 : R) = 1 := by
  apply Subtype.ext
  simp [Matrix.one_fin_two]

theorem upper_inv (t : R) : (upper t)⁻¹ = upper (-t) :=
  inv_eq_of_mul_eq_one_right (by rw [upper_mul_upper]; simp)

theorem lower_inv (t : R) : (lower t)⁻¹ = lower (-t) :=
  inv_eq_of_mul_eq_one_right (by rw [lower_mul_lower]; simp)

/-- ★原文が `α^λ`(λ は `𝔽_l` の元)と書けるのは `α^l = 1` だからである。
Lean では指数を `ℕ` に取り、`upper 1 ^ n = upper n` として同じことを言う。 -/
theorem upper_pow (t : R) (n : ℕ) : upper t ^ n = upper (n * t) := by
  induction n with
  | zero => simp
  | succ k ih => rw [pow_succ, ih, upper_mul_upper]; push_cast; ring_nf

theorem lower_pow (t : R) (n : ℕ) : lower t ^ n = lower (n * t) := by
  induction n with
  | zero => simp
  | succ k ih => rw [pow_succ, ih, lower_mul_lower]; push_cast; ring_nf

end Elementary

section GeneratedByElementary

variable {R : Type*} [CommRing R]

/-- `SL(2, R)` の元の行列式の展開。 -/
theorem det_entries (g : SL(2, R)) :
    (g : Matrix (Fin 2) (Fin 2) R) 0 0 * (g : Matrix (Fin 2) (Fin 2) R) 1 1
      - (g : Matrix (Fin 2) (Fin 2) R) 0 1 * (g : Matrix (Fin 2) (Fin 2) R) 1 0 = 1 := by
  have h := g.2
  rwa [Matrix.det_fin_two] at h

/-- ★分解を座標だけで書いたもの(`SL` の型は関係しない)。

**左下成分 `c` が可逆**(`c · e = 1`)なら
`!![a,b;c,d] = upper((a−1)·e) · lower(c) · upper((d−1)·e)`。
★体でなく **`c` の可逆性だけ**を使うので、`ℤ_[l]` のような局所環でもそのまま効く。
`(0,1)` 成分でだけ**行列式の条件 `ad − bc = 1` を使う**——他の 3 成分は使わない。 -/
theorem matrix_factor (a b c d e : R) (hdet : a * d - b * c = 1) (he : c * e = 1) :
    !![a, b; c, d]
      = (!![1, (a - 1) * e; 0, 1] * !![1, 0; c, 1]) * !![1, (d - 1) * e; 0, 1] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp <;>
    first
      | linear_combination (1 - a) * he
      | linear_combination (-(b + (d - 1) * (a - 1) * e)) * he + (-e) * hdet
      | linear_combination (1 - d) * he

/-- ★左下成分が可逆なときの明示的な分解。 -/
theorem eq_upper_mul_lower_mul_upper (g : SL(2, R)) {e : R}
    (he : (g : Matrix (Fin 2) (Fin 2) R) 1 0 * e = 1) :
    g = upper (((g : Matrix (Fin 2) (Fin 2) R) 0 0 - 1) * e)
        * lower ((g : Matrix (Fin 2) (Fin 2) R) 1 0)
        * upper (((g : Matrix (Fin 2) (Fin 2) R) 1 1 - 1) * e) := by
  apply Subtype.ext
  rw [Matrix.SpecialLinearGroup.coe_mul, Matrix.SpecialLinearGroup.coe_mul,
    coe_upper, coe_lower, coe_upper]
  conv_lhs => rw [Matrix.eta_fin_two ((g : Matrix (Fin 2) (Fin 2) R))]
  exact matrix_factor _ _ _ _ _ (det_entries g) he

variable {S : Subgroup SL(2, R)} (hu : ∀ t : R, upper t ∈ S) (hl : ∀ t : R, lower t ∈ S)

include hu hl in
theorem mem_of_isUnit_lowerLeft (g : SL(2, R))
    (hc : IsUnit ((g : Matrix (Fin 2) (Fin 2) R) 1 0)) : g ∈ S := by
  obtain ⟨e, he⟩ := hc.exists_right_inv
  rw [eq_upper_mul_lower_mul_upper g he]
  exact mul_mem (mul_mem (hu _) (hl _)) (hu _)

end GeneratedByElementary

section LocalRing

variable {R : Type*} [CommRing R] [IsLocalRing R]

/-- 局所環では「単元 + 非単元 = 単元」。 -/
theorem isUnit_add_of_not_isUnit {a c : R} (ha : IsUnit a) (hc : ¬ IsUnit c) :
    IsUnit (a + c) := by
  have h : IsUnit ((a + c) + (-c)) := by simpa using ha
  rcases IsLocalRing.isUnit_or_isUnit_of_isUnit_add h with h1 | h2
  · exact h1
  · exact absurd (by simpa using h2.neg) hc

variable {S : Subgroup SL(2, R)} (hu : ∀ t : R, upper t ∈ S) (hl : ∀ t : R, lower t ∈ S)

include hu hl in
/-- ★★**局所環 `R` 上では、基本行列を全部含む部分群は `SL(2, R)` 全体である。**

`c` が可逆なら分解がそのまま効く。可逆でないときは `ad − bc = 1` から
**`a` が可逆**で(局所環だから `ad` か `bc` の一方は可逆、`bc` が可逆なら `c` も可逆)、
`lower 1 · g` の左下成分が `a + c` = 単元 + 非単元 = 単元になる。

★**これは `Lemma 3.1, (iv)` の 2 本柱の 1 本**である——`ℤ_[l]` は局所環なので、
`SL(2, ℤ_[l])` は基本行列で(位相を使わず)**代数的に**生成される。
残る 1 本(閉部分群が基本行列を全部含むこと)は [Serre] Ch IV §3.4 Lemma 3 に当たる。 -/
theorem mem_of_upper_lower_mem (g : SL(2, R)) : g ∈ S := by
  by_cases hc : IsUnit ((g : Matrix (Fin 2) (Fin 2) R) 1 0)
  · exact mem_of_isUnit_lowerLeft hu hl g hc
  · -- `a` が可逆であること
    have hdet := det_entries g
    have hsum : IsUnit ((g : Matrix (Fin 2) (Fin 2) R) 0 0 * (g : Matrix (Fin 2) (Fin 2) R) 1 1
        + -((g : Matrix (Fin 2) (Fin 2) R) 0 1 * (g : Matrix (Fin 2) (Fin 2) R) 1 0)) := by
      rw [← sub_eq_add_neg, hdet]; exact isUnit_one
    have ha : IsUnit ((g : Matrix (Fin 2) (Fin 2) R) 0 0) := by
      rcases IsLocalRing.isUnit_or_isUnit_of_isUnit_add hsum with h1 | h2
      · exact isUnit_of_mul_isUnit_left h1
      · exact absurd (isUnit_of_mul_isUnit_right (by simpa using h2.neg)) hc
    -- `lower 1 * g` の左下成分は `a + c`
    have hentry : ((lower 1 * g : SL(2, R)) : Matrix (Fin 2) (Fin 2) R) 1 0
        = (g : Matrix (Fin 2) (Fin 2) R) 0 0 + (g : Matrix (Fin 2) (Fin 2) R) 1 0 := by
      rw [Matrix.SpecialLinearGroup.coe_mul]
      simp [Matrix.mul_apply, Fin.sum_univ_succ]
    have hstep : lower 1 * g ∈ S :=
      mem_of_isUnit_lowerLeft hu hl _ (by rw [hentry]; exact isUnit_add_of_not_isUnit ha hc)
    have hg : g = lower (-1) * (lower 1 * g) := by
      rw [← mul_assoc, lower_mul_lower]; simp
    rw [hg]
    exact mul_mem (hl _) hstep

/-- 基本行列の生成する部分群は `SL(2, R)` 全体(局所環 `R`)。 -/
theorem closure_elementary_eq_top :
    Subgroup.closure ((Set.range (upper : R → SL(2, R))) ∪ (Set.range (lower : R → SL(2, R))))
      = ⊤ :=
  eq_top_iff.2 fun g _ =>
    mem_of_upper_lower_mem (fun t => Subgroup.subset_closure (Or.inl ⟨t, rfl⟩))
      (fun t => Subgroup.subset_closure (Or.inr ⟨t, rfl⟩)) g

end LocalRing

section PadicPillar

/-- ★★**`SL(2, ℤ_[l])` は基本行列で代数的に生成される。**

`Lemma 3.1, (iv)` の 2 本柱の 1 本。★**位相を一切使っていない**——
`ℤ_[l]` が局所環であることだけで出る。
残る 1 本は「閉部分群 `J` が `upper t`・`lower t` を**全部**含むこと」であり、
そこが [Serre] Ch IV §3.4 Lemma 3(`0_Source` に無い)に当たる。 -/
theorem closure_elementary_padic (l : ℕ) [Fact (Nat.Prime l)] :
    Subgroup.closure
        ((Set.range (upper : ℤ_[l] → SL(2, ℤ_[l]))) ∪ (Set.range (lower : ℤ_[l] → SL(2, ℤ_[l]))))
      = ⊤ :=
  closure_elementary_eq_top

end PadicPillar

section Reduction

variable {R : Type*} [CommRing R] {S : Type*} [CommRing S]

@[simp] theorem map_upper (f : R →+* S) (t : R) :
    Matrix.SpecialLinearGroup.map f (upper t) = upper (f t) := by
  apply Subtype.ext
  ext i j
  fin_cases i <;> fin_cases j <;> simp [Matrix.SpecialLinearGroup.map, upper]

@[simp] theorem map_lower (f : R →+* S) (t : R) :
    Matrix.SpecialLinearGroup.map f (lower t) = lower (f t) := by
  apply Subtype.ext
  ext i j
  fin_cases i <;> fin_cases j <;> simp [Matrix.SpecialLinearGroup.map, lower]

end Reduction


section Diagonal

variable {K : Type*} [Field K]

/-- 原文 `ε ≝ diag(λ, λ⁻¹)`。 -/
def diagUnit (u : Kˣ) : SL(2, K) :=
  ⟨!![(u : K), 0; 0, ((u : K))⁻¹], by
    rw [Matrix.det_fin_two_of, mul_inv_cancel₀ u.ne_zero, mul_zero, sub_zero]⟩

@[simp] theorem coe_diagUnit (u : Kˣ) :
    ((diagUnit u : SL(2, K)) : Matrix (Fin 2) (Fin 2) K)
      = !![(u : K), 0; 0, ((u : K))⁻¹] := rfl

theorem diag_mul_upper (u : Kˣ) (t : K) :
    diagUnit u * upper t = upper ((u : K) ^ 2 * t) * diagUnit u := by
  have hu : (u : K) ≠ 0 := u.ne_zero
  apply Subtype.ext
  rw [Matrix.SpecialLinearGroup.coe_mul, Matrix.SpecialLinearGroup.coe_mul,
    coe_upper, coe_upper, coe_diagUnit]
  ext i j
  fin_cases i <;> fin_cases j <;> (simp; try field_simp)

theorem diag_mul_lower (u : Kˣ) (t : K) :
    diagUnit u * lower ((u : K) ^ 2 * t) = lower t * diagUnit u := by
  have hu : (u : K) ≠ 0 := u.ne_zero
  apply Subtype.ext
  rw [Matrix.SpecialLinearGroup.coe_mul, Matrix.SpecialLinearGroup.coe_mul,
    coe_lower, coe_lower, coe_diagUnit]
  ext i j
  fin_cases i <;> fin_cases j <;> (simp; try field_simp)

/-- ★原文の `ε · α · ε⁻¹ · α⁻¹ = α^{λ²−1}` そのもの。 -/
theorem commutator_diag_upper (u : Kˣ) (t : K) :
    diagUnit u * upper t * (diagUnit u)⁻¹ * (upper t)⁻¹ = upper (((u : K) ^ 2 - 1) * t) := by
  have h : diagUnit u * upper t * (diagUnit u)⁻¹ = upper ((u : K) ^ 2 * t) := by
    rw [diag_mul_upper, mul_assoc, mul_inv_cancel, mul_one]
  rw [h, upper_inv, upper_mul_upper]
  congr 1
  ring

/-- ★下三角側。原文は「(and, similarly, β)」の一言で済ませている段。 -/
theorem commutator_diag_lower (u : Kˣ) (t : K) :
    diagUnit u * lower ((u : K) ^ 2 * t) * (diagUnit u)⁻¹ * (lower ((u : K) ^ 2 * t))⁻¹
      = lower ((1 - (u : K) ^ 2) * t) := by
  have h : diagUnit u * lower ((u : K) ^ 2 * t) * (diagUnit u)⁻¹ = lower t := by
    rw [diag_mul_lower, mul_assoc, mul_inv_cancel, mul_one]
  rw [h, lower_inv, lower_mul_lower]
  congr 1
  ring

end Diagonal

section PrimeField

variable (l : ℕ) [Fact (Nat.Prime l)]

theorem upper_mem_closure (t : ZMod l) :
    upper t ∈ Subgroup.closure ({upper 1, lower 1} : Set SL(2, ZMod l)) := by
  haveI : NeZero l := ⟨(Fact.out : Nat.Prime l).ne_zero⟩
  have ht : upper t = upper (1 : ZMod l) ^ t.val := by
    rw [upper_pow, mul_one, ZMod.natCast_val, ZMod.cast_id]
  rw [ht]
  exact pow_mem (Subgroup.subset_closure (by simp)) _

theorem lower_mem_closure (t : ZMod l) :
    lower t ∈ Subgroup.closure ({upper 1, lower 1} : Set SL(2, ZMod l)) := by
  haveI : NeZero l := ⟨(Fact.out : Nat.Prime l).ne_zero⟩
  have ht : lower t = lower (1 : ZMod l) ^ t.val := by
    rw [lower_pow, mul_one, ZMod.natCast_val, ZMod.cast_id]
  rw [ht]
  exact pow_mem (Subgroup.subset_closure (by simp)) _

/-! ### ★(i) -/

/-- **[GenEll] Lemma 3.1, (i)**。

原文 (GenEll p.14):
> (i) Let G ⊆ SL_2(F[bb]_l) be the subgroup generated by the matrices α ≝

原文 (GenEll p.14):
> Then G = SL_2(F[bb]_l).

★原文の `l ≥ 5` は (i) では使わない(原文の証明も使っていない)。 -/
theorem lemma_3_1_i :
    Subgroup.closure ({upper 1, lower 1} : Set SL(2, ZMod l)) = ⊤ :=
  eq_top_iff.2 fun g _ => mem_of_upper_lower_mem (upper_mem_closure l) (lower_mem_closure l) g

/-! ### ★(ii) -/

/-- ★`l ≥ 5` が効く唯一の箇所: `λ² ≠ 1` なる `λ ∈ 𝔽_l^×` の存在。

原文「Note that such a λ exists since l ≥ 5」の中身。
`{0, 1, −1}` は高々 3 元、`|𝔽_l| = l ≥ 5` なので、その外に元がある。 -/
theorem exists_unit_sq_ne_one (hl : 5 ≤ l) : ∃ u : (ZMod l)ˣ, ((u : ZMod l)) ^ 2 ≠ 1 := by
  haveI : NeZero l := ⟨(Fact.out : Nat.Prime l).ne_zero⟩
  by_contra hcon
  simp only [not_exists, not_not] at hcon
  have hsub : (Finset.univ : Finset (ZMod l)) ⊆ ({0, 1, -1} : Finset (ZMod l)) := by
    intro a _
    by_cases ha : a = 0
    · simp [ha]
    · have hsq : a ^ 2 = 1 := by simpa using hcon (Units.mk0 a ha)
      have hfac : (a - 1) * (a + 1) = 0 := by linear_combination hsq
      rcases mul_eq_zero.1 hfac with h1 | h2
      · have : a = 1 := by linear_combination h1
        simp [this]
      · have : a = -1 := by linear_combination h2
        simp [this]
  have hcard := Finset.card_le_card hsub
  have h3 : ({0, 1, -1} : Finset (ZMod l)).card ≤ 3 := by
    refine le_trans (Finset.card_insert_le _ _) ?_
    have h2 := Finset.card_insert_le (1 : ZMod l) ({-1} : Finset (ZMod l))
    simp only [Finset.card_singleton] at h2
    omega
  rw [Finset.card_univ, ZMod.card] at hcard
  omega

/-- **[GenEll] Lemma 3.1, (ii)**。

原文 (GenEll p.14):
> (ii) The finite group SL_2(F[bb]_l) has no nontrivial abelian quotients.

★「非自明な可換商を持たない」を**そのままの形**で書いた——
「正規部分群 `N` の商が可換なら `N = ⊤`」。
交換子群が全体であることと同値だが、原文の主語は商の側なので商の側で述べる。 -/
theorem lemma_3_1_ii (hl : 5 ≤ l) (N : Subgroup SL(2, ZMod l)) [N.Normal]
    (habel : ∀ x y : SL(2, ZMod l) ⧸ N, x * y = y * x) : N = ⊤ := by
  obtain ⟨u, hu⟩ := exists_unit_sq_ne_one l hl
  have hne : ((u : ZMod l)) ^ 2 - 1 ≠ 0 := sub_ne_zero.2 hu
  have hne' : (1 : ZMod l) - ((u : ZMod l)) ^ 2 ≠ 0 := by
    intro h
    exact hne (by linear_combination -h)
  have hcomm : ∀ x y : SL(2, ZMod l), x * y * x⁻¹ * y⁻¹ ∈ N := by
    intro x y
    have h1 : (QuotientGroup.mk' N) (x * y * x⁻¹ * y⁻¹) = 1 := by
      simp only [map_mul, map_inv]
      rw [habel ((QuotientGroup.mk' N) x) ((QuotientGroup.mk' N) y)]
      group
    simpa using (QuotientGroup.eq_one_iff _).1 h1
  have hupper : ∀ s : ZMod l, upper s ∈ N := by
    intro s
    have h := hcomm (diagUnit u) (upper (s / (((u : ZMod l)) ^ 2 - 1)))
    rw [commutator_diag_upper] at h
    have hs : (((u : ZMod l)) ^ 2 - 1) * (s / (((u : ZMod l)) ^ 2 - 1)) = s := by field_simp
    rwa [hs] at h
  have hlower : ∀ s : ZMod l, lower s ∈ N := by
    intro s
    have h := hcomm (diagUnit u)
      (lower (((u : ZMod l)) ^ 2 * (s / (1 - ((u : ZMod l)) ^ 2))))
    rw [commutator_diag_lower] at h
    have hs : (1 - ((u : ZMod l)) ^ 2) * (s / (1 - ((u : ZMod l)) ^ 2)) = s := by field_simp
    rwa [hs] at h
  have hle : Subgroup.closure ({upper 1, lower 1} : Set SL(2, ZMod l)) ≤ N := by
    refine (Subgroup.closure_le _).2 ?_
    rintro x (rfl | rfl)
    · exact hupper 1
    · exact hlower 1
  rw [lemma_3_1_i] at hle
  exact top_le_iff.1 hle

end PrimeField

/-! ### ★(ii) の負の対照 —— `l ≥ 5` は飾りではない -/

/-- ★`l = 2, 3` では **すべての** `λ ∈ 𝔽_l^×` が `λ² = 1` を満たす。

すなわち `exists_unit_sq_ne_one` は `l ≥ 5` 無しでは**偽**であり、
(ii) の証明の要になる `ε` は作れない。
(★これは「(ii) 自身が偽」の証明ではない。(ii) が実際に `l = 3` で偽であること
——`SL₂(𝔽_3)` が位数 3 の可換商を持つこと——は本ファイルでは示していない。) -/
theorem lemma_3_1_ii_neg_control :
    (∀ u : (ZMod 2)ˣ, ((u : ZMod 2)) ^ 2 = 1) ∧ (∀ u : (ZMod 3)ˣ, ((u : ZMod 3)) ^ 2 = 1) :=
  ⟨by decide, by decide⟩

section GeneralLinear

open Matrix.SpecialLinearGroup

variable {K : Type*} [Field K]

@[simp] theorem coe_toGL (g : SL(2, K)) :
    ((toGL g : GL (Fin 2) K) : Matrix (Fin 2) (Fin 2) K)
      = (g : Matrix (Fin 2) (Fin 2) K) := rfl

/-- ★(iii) の核となる行列の等式。

`W ≝ upper(−a/c) · M` は**左上成分が 0** になる(`a + (−a/c)·c = 0`)。
その `W` で `upper` を共役すると**下三角の基本行列になる**——
これが「上三角でない行列が 1 つあれば `β` が作れる」の中身である。
`c ≠ 0` は `−a/c` を書くためにも、結論の `lower 1` を得るためにも要る。 -/
theorem gl_conj_key (a b c d : K) (hc : c ≠ 0) :
    (!![1, -a / c; 0, 1] * !![a, b; c, d]) * !![1, -(a * d - b * c) / c ^ 2; 0, 1]
      = !![1, 0; 1, 1] * (!![1, -a / c; 0, 1] * !![a, b; c, d]) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp <;> field_simp <;>
    first
      | ring1
      | (refine Or.inl ?_; ring1)

end GeneralLinear

section PrimeFieldGL

open Matrix.SpecialLinearGroup

variable (l : ℕ) [Fact (Nat.Prime l)]

/-! ### ★(iii) -/

/-- **[GenEll] Lemma 3.1, (iii)**。

原文 (GenEll p.14):
> (iii) Let H ⊆ GL_2(F[bb]_l) be a subgroup that contains the matrix α ≝

原文 (GenEll p.14):
> well as at least one matrix that is not upper triangular. Then SL_2(F[bb]_l) ⊆ H.

★「上三角でない」は 2×2 では **`(1,0)` 成分 ≠ 0** のことである。
★結論の `SL₂(𝔽_l) ⊆ H` は、`toGL` による像がすべて `H` に入ることとして書いた。 -/
theorem lemma_3_1_iii (H : Subgroup (GL (Fin 2) (ZMod l)))
    (hα : (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l)) ∈ H)
    (M : GL (Fin 2) (ZMod l)) (hM : M ∈ H)
    (hMnt : (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 1 0 ≠ 0) (g : SL(2, ZMod l)) :
    (toGL g : GL (Fin 2) (ZMod l)) ∈ H := by
  haveI : NeZero l := ⟨(Fact.out : Nat.Prime l).ne_zero⟩
  -- `α` の冪で `upper` 全体が入る
  have hupper : ∀ t : ZMod l, (toGL (upper t) : GL (Fin 2) (ZMod l)) ∈ H := by
    intro t
    have ht : upper t = upper (1 : ZMod l) ^ t.val := by
      rw [upper_pow, mul_one, ZMod.natCast_val, ZMod.cast_id]
    rw [ht, map_pow]
    exact pow_mem hα _
  set a := (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 0 0 with ha
  set b := (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 0 1 with hb
  set c := (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 1 0 with hcdef
  set d := (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 1 1 with hd
  set W : GL (Fin 2) (ZMod l) := (toGL (upper (-a / c))) * M with hW
  have hWmem : W ∈ H := mul_mem (hupper _) hM
  have hUmem : (toGL (upper (-(a * d - b * c) / c ^ 2)) : GL (Fin 2) (ZMod l)) ∈ H := hupper _
  have hkey : W * (toGL (upper (-(a * d - b * c) / c ^ 2)))
      = (toGL (lower (1 : ZMod l))) * W := by
    apply Units.ext
    simp only [hW, Units.val_mul, coe_toGL, coe_upper, coe_lower]
    conv_lhs => rw [Matrix.eta_fin_two ((M : Matrix (Fin 2) (Fin 2) (ZMod l)))]
    conv_rhs => rw [Matrix.eta_fin_two ((M : Matrix (Fin 2) (Fin 2) (ZMod l)))]
    exact gl_conj_key a b c d hMnt
  have hlower1 : (toGL (lower (1 : ZMod l)) : GL (Fin 2) (ZMod l)) ∈ H := by
    have hrw : (toGL (lower (1 : ZMod l)) : GL (Fin 2) (ZMod l))
        = W * (toGL (upper (-(a * d - b * c) / c ^ 2))) * W⁻¹ :=
      eq_mul_inv_of_mul_eq hkey.symm
    rw [hrw]
    exact mul_mem (mul_mem hWmem hUmem) (inv_mem hWmem)
  -- (i) から `SL₂` 全体が入る
  have hcomap : (⊤ : Subgroup SL(2, ZMod l)) ≤ Subgroup.comap (toGL : SL(2, ZMod l) →* _) H := by
    rw [← lemma_3_1_i l]
    refine (Subgroup.closure_le _).2 ?_
    rintro x (rfl | rfl)
    · exact hupper 1
    · exact hlower1
  exact hcomap (Subgroup.mem_top g)

end PrimeFieldGL

section PadicReduction

variable (l : ℕ) [Fact (Nat.Prime l)]

/-- 還元 `SL(2, ℤ_l) → SL(2, 𝔽_l)`。(iv) はこの写像の**切断**についての主張である。 -/
noncomputable def redPadic : SL(2, ℤ_[l]) →* SL(2, ZMod l) :=
  Matrix.SpecialLinearGroup.map (PadicInt.toZMod)

/-- ★還元 `SL(2, ℤ_l) → SL(2, 𝔽_l)` は**全射**である。

★証明は (i) をそのまま使う——`SL(2, 𝔽_l)` は `upper 1` と `lower 1` で生成され、
その 2 つは明らかに持ち上がるから。**位相も Hensel も要らない。** -/
theorem redPadic_surjective : Function.Surjective (redPadic l) := by
  rw [← MonoidHom.range_eq_top, eq_top_iff, ← lemma_3_1_i l]
  refine (Subgroup.closure_le _).2 ?_
  rintro x (rfl | rfl)
  · exact ⟨upper 1, by simp [redPadic]⟩
  · exact ⟨lower 1, by simp [redPadic]⟩

end PadicReduction

section FiniteReduction

open Matrix.SpecialLinearGroup

variable (l n : ℕ) [Fact (Nat.Prime l)]

/-- 有限段の還元 `SL(2, ℤ/l^n) → SL(2, 𝔽_l)`。

★`Lemma 3.1, (iv)` の証明は「閉部分群」の主張だが、
**位相が要るのは逆極限の 1 歩だけ**で、帰納法の本体はこの有限段で回る
(`ResearchPaper/genell-goal.md` の (A) / (B) の切り分け)。 -/
def redPow (hn : n ≠ 0) : SL(2, ZMod (l ^ n)) →* SL(2, ZMod l) :=
  Matrix.SpecialLinearGroup.map (ZMod.castHom (dvd_pow_self l hn) (ZMod l))

/-- ★有限段の還元も**全射**である。証明は `ℤ_l` 版と同じ——(i) をそのまま使う。 -/
theorem redPow_surjective (hn : n ≠ 0) : Function.Surjective (redPow l n hn) := by
  rw [← MonoidHom.range_eq_top, eq_top_iff, ← lemma_3_1_i l]
  refine (Subgroup.closure_le _).2 ?_
  rintro x (rfl | rfl)
  · exact ⟨upper 1, by rw [redPow, map_upper, map_one]⟩
  · exact ⟨lower 1, by rw [redPow, map_lower, map_one]⟩

end FiniteReduction

/-! ## ★出典の紐付け(`.src`)

★条つきの `.src` は **locator の記録**であって、命題全体を完全に実装したという主張ではない。
`Lemma 3.1` は (iv) が未実装なので、**条なしの `.src` は付けない**
(`tools/frdi-progress.mjs` と同じ規則: 条つきは数に入らない)。 -/

def lemma_3_1_i.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 14, item := "Lemma 3.1, (i)",
    sectionId := "genell-lemma-3-1-i" }

def lemma_3_1_ii.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 14, item := "Lemma 3.1, (ii)",
    sectionId := "genell-lemma-3-1-ii" }

def lemma_3_1_iii.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 14, item := "Lemma 3.1, (iii)",
    sectionId := "genell-lemma-3-1-iii" }

end ABC3.Found.GenEll
