/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Data.NNReal.Basic
import Mathlib.Algebra.Module.Rat
import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Tactic.FieldSimp
import ABC3.Meta.Claim

/-!
# `ℝ≥0` の加法自己同型は正の実数倍 —— `Theorem 6.4, (ii)` の核(`Found`、`sorry` 無し)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.114):
> 5.3] of C1, C2. Then there exists an element deg(Ψrlf) ∈ R>0 such that for all

## ★★なにを閉じたか

原文 `Theorem 6.4, (ii)` は「`Φ^rlf₁(A₁) ≅ Φ^rlf₂(A₂)` は単系の同型だから
順序を保つ」と 1 行で書き、そこから `deg(Ψ^rlf) ∈ ℝ>0` を出す。

★単系の同型が順序を保つのは自明である(`ℝ≥0` では `a ≤ b ↔ ∃ c, a + c = b`)。
★★中身は**その先**、すなわち

```
加法自己同型 f : ℝ≥0 ≃+ ℝ≥0  ⟹  ∃ c ≠ 0, ∀ x, f x = c · x
```

である。★これが `deg(Ψ^rlf) ∈ ℝ>0` の実体で、**mathlib に無かった**
(2026-08-25 実測。`AddMonoidHom.toRealLinearMap` は**連続性**を要求し、
単調性から連続性を出す補題は無い)。

## ★★★中身は 3 段

| 段 | 中身 | 在庫 |
|---|---|---|
| 1 | `f` は単調 | `le_iff_exists_add`(自明) |
| 2 | `f (r) = f 1 · r`(`r : ℚ≥0`) | `map_nnrat_smul`(mathlib) |
| 3 | 単調＋稠密性で全部の `x` に広げる | ★本ファイル |

★段 2 が **`ℚ` でなく `ℚ≥0`** で済むのが `ℝ≥0` の得なところである
(`ℝ≥0` は群ではないので `ℚ`-加群にはならないが、`ℚ≥0`-加群にはなる)。

★段 3 は「`f x < c·x` なら `f x / c < r < x` なる有理数を取って矛盾」の 2 回。
そのために `ℚ≥0` が `ℝ≥0` で稠密であること(`nnreal_exists_nnrat_btwn`)を
`exists_rat_btwn` から作る —— ★これも mathlib に `ℝ≥0` 版が無かった。
-/

namespace ABC3.Found

open scoped NNReal

/-! ## ★1. `ℚ≥0` は `ℝ≥0` で稠密 -/

/-- ★★**`ℚ≥0` は `ℝ≥0` で稠密** —— `ℝ` の `exists_rat_btwn` から作る。

★mathlib に `ℝ≥0` 版は無い(2026-08-25 実測)。 -/
theorem nnreal_exists_nnrat_btwn (a b : ℝ≥0) (h : a < b) :
    ∃ r : ℚ≥0, a < (r : ℝ≥0) ∧ (r : ℝ≥0) < b := by
  obtain ⟨q, hq1, hq2⟩ := exists_rat_btwn (NNReal.coe_lt_coe.mpr h)
  have hq0 : (0 : ℝ) ≤ (q : ℝ) := le_trans a.coe_nonneg hq1.le
  have hq0' : (0 : ℚ) ≤ q := by exact_mod_cast hq0
  have hcast : (((q.toNNRat : ℚ≥0) : ℝ≥0) : ℝ) = (q : ℝ) := by
    rw [show (((q.toNNRat : ℚ≥0) : ℝ≥0) : ℝ) = ((q.toNNRat : ℚ) : ℝ) by norm_cast,
      Rat.coe_toNNRat _ hq0']
  exact ⟨q.toNNRat, by rw [← NNReal.coe_lt_coe, hcast]; exact hq1,
    by rw [← NNReal.coe_lt_coe, hcast]; exact hq2⟩

/-! ## ★2. 有理数点での値 -/

/-- ★★**加法同型は `ℚ≥0` 倍を保つ** —— `f r = f 1 · r`。

★`ℝ≥0` は `ℚ≥0`-加群なので `map_nnrat_smul` がそのまま当たる。 -/
theorem addEquiv_nnrat_cast (f : ℝ≥0 ≃+ ℝ≥0) (r : ℚ≥0) :
    f (r : ℝ≥0) = f 1 * (r : ℝ≥0) := by
  have h1 : ((r : ℚ≥0) • (1 : ℝ≥0)) = (r : ℝ≥0) := by simp
  have h2 := map_nnrat_smul f r (1 : ℝ≥0)
  rw [h1, NNRat.smul_def] at h2
  rw [h2, mul_comm]

/-! ## ★3. 本体 -/

/-- ★★**加法同型は順序を保つ**(原文が 1 行で書いている段)。

★`ℝ≥0` では `a ≤ b ↔ ∃ c, a + c = b` なので、加法性だけで単調になる。 -/
theorem monotone_of_addEquiv (f : ℝ≥0 ≃+ ℝ≥0) : Monotone f := by
  intro a b hab
  obtain ⟨d, rfl⟩ := le_iff_exists_add.mp hab
  simp [map_add]

/-- ★★★★★★**`ℝ≥0` の加法自己同型は正の実数倍**。

★★これが `Theorem 6.4, (ii)` の `deg(Ψ^rlf) ∈ ℝ>0` の実体である。 -/
theorem exists_smul_eq_of_addEquiv (f : ℝ≥0 ≃+ ℝ≥0) :
    ∃ c : ℝ≥0, c ≠ 0 ∧ ∀ x : ℝ≥0, f x = c * x := by
  have hmono : Monotone f := monotone_of_addEquiv f
  have hc : f 1 ≠ 0 := by
    intro h
    have h0 : f (1 : ℝ≥0) = f 0 := by rw [h, map_zero]
    exact one_ne_zero (f.injective h0)
  have hcpos : 0 < f 1 := pos_iff_ne_zero.mpr hc
  refine ⟨f 1, hc, fun x => le_antisymm ?_ ?_⟩
  · -- ★`f 1 · x < f x` なら `x < f x / f 1` なる有理数を挟んで矛盾
    by_contra hlt
    have hlt' : f 1 * x < f x := not_le.mp hlt
    have hx : x < f x / f 1 := by rw [lt_div_iff₀ hcpos, mul_comm]; exact hlt'
    obtain ⟨r, hr1, hr2⟩ := nnreal_exists_nnrat_btwn x (f x / f 1) hx
    have hle : f x ≤ f (r : ℝ≥0) := hmono hr1.le
    rw [addEquiv_nnrat_cast f r] at hle
    have hlt2 : f 1 * (r : ℝ≥0) < f x :=
      calc f 1 * (r : ℝ≥0) < f 1 * (f x / f 1) := mul_lt_mul_of_pos_left hr2 hcpos
        _ = f x := by field_simp
    exact absurd hle (not_le.mpr hlt2)
  · -- ★`f x < f 1 · x` なら `f x / f 1 < r < x` なる有理数を挟んで矛盾
    by_contra hlt
    have hlt' : f x < f 1 * x := not_le.mp hlt
    have hx : f x / f 1 < x := by rw [div_lt_iff₀ hcpos, mul_comm]; exact hlt'
    obtain ⟨r, hr1, hr2⟩ := nnreal_exists_nnrat_btwn (f x / f 1) x hx
    have hle : f (r : ℝ≥0) ≤ f x := hmono hr2.le
    rw [addEquiv_nnrat_cast f r] at hle
    have hlt2 : f x < f 1 * (r : ℝ≥0) :=
      calc f x = f 1 * (f x / f 1) := by field_simp
        _ < f 1 * (r : ℝ≥0) := mul_lt_mul_of_pos_left hr1 hcpos
    exact absurd hle (not_le.mpr hlt2)

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

open ABC3.Meta in
def exists_smul_eq_of_addEquiv.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (ii) — deg(Ψ^rlf) ∈ ℝ>0(加法自己同型は正数倍)",
    sectionId := "frdi-thm-6-4" }

open ABC3.Meta in
def exists_smul_eq_of_addEquiv.needs : List ProofObligation :=
  [ .citation "[mathlib]" "map_nnrat_smul(加法準同型は ℚ≥0 倍を保つ)"
      (.inMathlib "map_nnrat_smul") 114,
    .citation "[mathlib]" "exists_rat_btwn(ℝ における ℚ の稠密性)"
      (.inMathlib "exists_rat_btwn") 114,
    .citation "[mathlib]" "単調な加法準同型が線型であること"
      (.absent "AddMonoidHom.toRealLinearMap は連続性を要求する。単調性から線型性を出す宣言は無い(2026-08-25 実測)") 114,
    .derivation "単調＋有理点での値＋稠密性で挟み撃ちにする(両向き)" 114,
    .implicitStep
      "★原文は「単系の同型だから順序を保つ」の 1 行で deg(Ψ^rlf) ∈ ℝ>0 の存在まで畳んでいる" 114 ]

open ABC3.Meta in
def nnreal_exists_nnrat_btwn.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (ii) — ℚ≥0 は ℝ≥0 で稠密",
    sectionId := "frdi-thm-6-4" }

open ABC3.Meta in
def nnreal_exists_nnrat_btwn.needs : List ProofObligation :=
  [ .citation "[mathlib]" "exists_rat_btwn"
      (.inMathlib "exists_rat_btwn") 114,
    .derivation "ℝ の側で取って Rat.toNNRat で落とす" 114 ]

open ABC3.Meta in
def addEquiv_nnrat_cast.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (ii) — 加法同型は ℚ≥0 倍を保つ",
    sectionId := "frdi-thm-6-4" }

open ABC3.Meta in
def addEquiv_nnrat_cast.needs : List ProofObligation :=
  [ .citation "[mathlib]" "map_nnrat_smul"
      (.inMathlib "map_nnrat_smul") 114 ]

open ABC3.Meta in
def monotone_of_addEquiv.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (ii) — 単系の同型は順序を保つ",
    sectionId := "frdi-thm-6-4" }

open ABC3.Meta in
def monotone_of_addEquiv.needs : List ProofObligation :=
  [ .citation "[mathlib]" "le_iff_exists_add"
      (.inMathlib "le_iff_exists_add") 114,
    .implicitStep "★原文は「単系の同型だから順序を保つ」と 1 行で置いている" 114 ]

end ABC3.Found
