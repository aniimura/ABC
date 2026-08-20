/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithDivisor
import Mathlib.NumberTheory.NumberField.Completion.InfinitePlace

/-!
# `ord(F_v)` と `Φ(F) = ⊕_v ord(O_v^▷)`(鎖 `arith` の `ord-mon` / `arith-phi`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.112–113。

原文 (FrdI p.112):
> for the multiplicative monoid of elements of valuation ≤1; μ(F) ⊆O×

## ★★原文の言葉

原文は `ord(F_v) := F_v^× / O_v^×`、`ord(O_v^▷) := O_v^▷ / O_v^×` と置き、

* `v` が非アルキメデス: `ord(F_v) ≅ ℤ`、`ord(O_v^▷) ≅ ℤ≥0`
* `v` が**アルキメデス**: `ord(F_v) ≅ ℝ`、`ord(O_v^▷) ≅ ℝ≥0`

と述べる。★**分母も分子も完備化 `F_v` の側**であることに注意する ——
アルキメデスで `ℝ` 全体が出るのはそのためである(`L` 自身の絶対値の像は
`ℝ_{>0}` の稠密な真部分群にすぎない)。

## ★本ファイルの 2 段

| 段 | 中身 |
|---|---|
| §1 | アルキメデス側 `ord(F_v) ≅ ℝ`(`ordArch`、全射性は完備化が `ℝ` か `ℂ` であることから) |
| §2 | `Φ(F)^gp ≅ (V_fin →₀ ℤ) × (V_∞ →₀ ℝ)` と、その有効部分が `Φ(F)`(`arithOfParts`) |

★非アルキメデス側(`ordFin`、`ord(F_v) ≅ ℤ`)は `ArithDivisor.lean` にある。
-/

namespace ABC3.Found.Divisor

open NumberField

universe u

variable {L : Type u} [Field L] [NumberField L]

/-! ## ★1. アルキメデス素点 —— `ord(F_v) ≅ ℝ` -/

section Arch

open NumberField.InfinitePlace NumberField.InfinitePlace.Completion

/-- ★完備化への埋め込みはノルムを保つ(等長かつ `0 ↦ 0`)。 -/
theorem norm_extensionEmbedding (v : InfinitePlace L) (x : v.Completion) :
    ‖extensionEmbedding v x‖ = ‖x‖ :=
  (isometry_extensionEmbedding v).norm_map_of_map_zero (map_zero _) x

/-- ★★**完備化のノルムは `[0, ∞)` を全部取る** —— `F_v` は `ℝ` か `ℂ` だからである。

★これが原文の「`ord(F_v) ≅ ℝ`(`v` がアルキメデス)」の中身である。 -/
theorem exists_norm_eq (v : InfinitePlace L) {r : ℝ} (hr : 0 ≤ r) :
    ∃ x : v.Completion, ‖x‖ = r := by
  obtain ⟨y, hy⟩ : ∃ y : ℂ, y = (r : ℂ) := ⟨_, rfl⟩
  rcases v.isReal_or_isComplex with hv | hv
  · obtain ⟨x, hx⟩ := surjective_extensionEmbeddingOfIsReal hv r
    refine ⟨x, ?_⟩
    have h := (isometry_extensionEmbeddingOfIsReal hv).norm_map_of_map_zero (map_zero _) x
    rw [← h, hx, Real.norm_eq_abs, abs_of_nonneg hr]
  · obtain ⟨x, hx⟩ := surjective_extensionEmbedding_of_isComplex hv (r : ℂ)
    refine ⟨x, ?_⟩
    rw [← norm_extensionEmbedding v x, hx, Complex.norm_real, Real.norm_eq_abs,
      abs_of_nonneg hr]

/-- ★★**アルキメデス素点の `ord`** —— `ord_v(x) = -log ‖x‖`。

★`ord(O_v^▷) = {x | ‖x‖ ≤ 1}` の像はちょうど `ℝ≥0` である(`ordArch_surjective`)。 -/
noncomputable def ordArch (v : InfinitePlace L) (x : v.Completion) : ℝ := -Real.log ‖x‖

/-- ★`ord` は乗法を加法に移す。 -/
theorem ordArch_mul (v : InfinitePlace L) {x y : v.Completion} (hx : x ≠ 0) (hy : y ≠ 0) :
    ordArch v (x * y) = ordArch v x + ordArch v y := by
  have hnx : ‖x‖ ≠ 0 := norm_ne_zero_iff.mpr hx
  have hny : ‖y‖ ≠ 0 := norm_ne_zero_iff.mpr hy
  simp only [ordArch, norm_mul, Real.log_mul hnx hny]
  ring

/-- ★★**`ord : F_v^× → ℝ` は全射** —— すなわち `ord(F_v) ≅ ℝ`。 -/
theorem ordArch_surjective (v : InfinitePlace L) (r : ℝ) :
    ∃ x : v.Completion, x ≠ 0 ∧ ordArch v x = r := by
  obtain ⟨x, hx⟩ := exists_norm_eq v (r := Real.exp (-r)) (Real.exp_nonneg _)
  refine ⟨x, ?_, ?_⟩
  · intro h
    rw [h, norm_zero] at hx
    exact absurd hx.symm (Real.exp_ne_zero _)
  · rw [ordArch, hx, Real.log_exp]
    ring

/-- ★★**`ord(O_v^▷) = ℝ≥0`** —— `‖x‖ ≤ 1` はちょうど `0 ≤ ord_v(x)` である。 -/
theorem ordArch_nonneg_iff (v : InfinitePlace L) {x : v.Completion} (hx : x ≠ 0) :
    0 ≤ ordArch v x ↔ ‖x‖ ≤ 1 := by
  have hpos : 0 < ‖x‖ := norm_pos_iff.mpr hx
  rw [ordArch, neg_nonneg]
  exact Real.log_nonpos_iff hpos.le

/-- ★**`ord_v(x) = 0` はちょうど `x ∈ O_v^×`**。 -/
theorem ordArch_eq_zero_iff (v : InfinitePlace L) {x : v.Completion} (hx : x ≠ 0) :
    ordArch v x = 0 ↔ ‖x‖ = 1 := by
  have hpos : 0 < ‖x‖ := norm_pos_iff.mpr hx
  constructor
  · intro h
    have h0 : Real.log ‖x‖ = 0 := by simpa [ordArch, neg_eq_zero] using h
    have := Real.exp_log hpos
    rw [h0, Real.exp_zero] at this
    exact this.symm
  · intro h
    simp [ordArch, h]

end Arch

/-! ## ★2. `Φ(F)^gp = ⊕_v ord(F_v)` —— 直和分解 -/

/-- ★非アルキメデス素点の絶対ノルムの対数(正)。 -/
noncomputable def logAbsNorm (w : FinitePlace L) : ℝ :=
  Real.log (Ideal.absNorm (FinitePlace.maximalIdeal w).asIdeal : ℝ)

theorem logAbsNorm_pos (w : FinitePlace L) : 0 < logAbsNorm w :=
  log_absNorm_pos (FinitePlace.maximalIdeal w)

theorem logAbsNorm_ne_zero (w : FinitePlace L) : logAbsNorm w ≠ 0 :=
  (logAbsNorm_pos w).ne'

/-- ★`ℤ` 係数を `ord(F_v) = ℤ·log(N v)` に埋め込む。 -/
noncomputable def scaleFin (a : FinitePlace L →₀ ℤ) : FinitePlace L →₀ ℝ where
  support := a.support
  toFun w := (a w : ℝ) * logAbsNorm w
  mem_support_toFun w := by
    constructor
    · intro hw
      have ha : a w ≠ 0 := Finsupp.mem_support_iff.mp hw
      exact mul_ne_zero (Int.cast_ne_zero.mpr ha) (logAbsNorm_ne_zero w)
    · intro hne
      refine Finsupp.mem_support_iff.mpr (fun h => hne ?_)
      rw [h]
      simp

@[simp] theorem scaleFin_apply (a : FinitePlace L →₀ ℤ) (w : FinitePlace L) :
    scaleFin a w = (a w : ℝ) * logAbsNorm w := rfl

theorem scaleFin_add (a b : FinitePlace L →₀ ℤ) :
    scaleFin (a + b) = scaleFin a + scaleFin b := by
  ext w
  simp only [scaleFin_apply, Finsupp.add_apply, Int.cast_add]
  ring

/-- ★★**`(V_fin →₀ ℤ) × (V_∞ →₀ ℝ) → Φ(F)^gp`** —— 原文の `⊕_v ord(F_v)` からの写像。 -/
noncomputable def arithOfParts :
    ((FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)) →+ (ArithPlace L →₀ ℝ) where
  toFun p := (Finsupp.sumFinsuppAddEquivProdFinsupp (M := ℝ)).symm (scaleFin p.1, p.2)
  map_zero' := by
    have h : scaleFin (0 : FinitePlace L →₀ ℤ) = 0 := by ext w; simp
    rw [show ((0 : (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)).1) = 0 from rfl, h]
    exact map_zero _
  map_add' p q := by
    show (Finsupp.sumFinsuppAddEquivProdFinsupp (M := ℝ)).symm (scaleFin (p.1 + q.1), p.2 + q.2)
      = _
    rw [scaleFin_add]
    exact map_add _ (scaleFin p.1, p.2) (scaleFin q.1, q.2)

@[simp] theorem arithOfParts_inl (p : (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ))
    (w : FinitePlace L) : arithOfParts p (Sum.inl w) = (p.1 w : ℝ) * logAbsNorm w := by
  have h := Finsupp.fst_sumFinsuppAddEquivProdFinsupp
    ((Finsupp.sumFinsuppAddEquivProdFinsupp (M := ℝ)).symm (scaleFin p.1, p.2)) w
  rw [AddEquiv.apply_symm_apply] at h
  exact h.symm

@[simp] theorem arithOfParts_inr (p : (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ))
    (w : InfinitePlace L) : arithOfParts p (Sum.inr w) = p.2 w := by
  have h := Finsupp.snd_sumFinsuppAddEquivProdFinsupp
    ((Finsupp.sumFinsuppAddEquivProdFinsupp (M := ℝ)).symm (scaleFin p.1, p.2)) w
  rw [AddEquiv.apply_symm_apply] at h
  exact h.symm

/-- ★`arithOfParts` は単射(`log(N v) ≠ 0` だから成分ごとに戻せる)。 -/
theorem arithOfParts_injective :
    Function.Injective (arithOfParts (L := L)) := by
  intro p q h
  have hfin : ∀ w : FinitePlace L, p.1 w = q.1 w := by
    intro w
    have h1 := congrArg (fun d : ArithPlace L →₀ ℝ => d (Sum.inl w)) h
    simp only [arithOfParts_inl] at h1
    have := mul_right_cancel₀ (logAbsNorm_ne_zero w) h1
    exact_mod_cast this
  have hinf : ∀ w : InfinitePlace L, p.2 w = q.2 w := by
    intro w
    have h1 := congrArg (fun d : ArithPlace L →₀ ℝ => d (Sum.inr w)) h
    simpa only [arithOfParts_inr] using h1
  refine Prod.ext ?_ ?_
  · ext w; exact hfin w
  · ext w; exact hinf w

/-- ★`arithOfParts` の像は `Φ(F)^gp` に入る。 -/
theorem arithOfParts_mem (p : (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)) :
    arithOfParts p ∈ arithDivGroup L := by
  intro w
  exact ⟨p.1 w, arithOfParts_inl p w⟩

/-- ★★**`arithOfParts` の像はちょうど `Φ(F)^gp`**。

★`d ∈ Φ(F)^gp` は「各非アルキメデス成分が `log(N v)` の整数倍」であり、
その整数は `log(N v) ≠ 0` から一意に定まる。 -/
theorem arithOfParts_surjective {d : ArithPlace L →₀ ℝ} (hd : d ∈ arithDivGroup L) :
    ∃ p, arithOfParts p = d := by
  classical
  -- ★非アルキメデス成分の整数を選ぶ
  have hchoice : ∀ w : FinitePlace L, ∃ n : ℤ, d (Sum.inl w) = (n : ℝ) * logAbsNorm w := hd
  choose n hn using hchoice
  -- ★`n` は有限台
  have hfin : (Function.support n).Finite := by
    refine Set.Finite.subset (d.support.finite_toSet.preimage
      (Set.injOn_of_injective Sum.inl_injective)) ?_
    intro w hw
    refine Finsupp.mem_support_iff.mpr ?_
    rw [hn w]
    exact mul_ne_zero (Int.cast_ne_zero.mpr hw) (logAbsNorm_ne_zero w)
  let a : FinitePlace L →₀ ℤ :=
    { support := hfin.toFinset, toFun := n,
      mem_support_toFun := fun w => by simp [Set.Finite.mem_toFinset, Function.mem_support] }
  let b : InfinitePlace L →₀ ℝ :=
    { support := d.support.preimage Sum.inr (Set.injOn_of_injective Sum.inr_injective),
      toFun := fun w => d (Sum.inr w),
      mem_support_toFun := fun w => by
        simp only [Finset.mem_preimage, Finsupp.mem_support_iff] }
  have hav : ∀ w, a w = n w := fun _ => rfl
  have hbv : ∀ w, b w = d (Sum.inr w) := fun _ => rfl
  refine ⟨(a, b), ?_⟩
  ext x
  cases x with
  | inl w =>
      rw [arithOfParts_inl]
      show ((a w : ℤ) : ℝ) * logAbsNorm w = d (Sum.inl w)
      rw [hav w, ← hn w]
  | inr w =>
      rw [arithOfParts_inr]
      show b w = d (Sum.inr w)
      rw [hbv w]

/-- ★★★★★**[FrdI] Example 6.3 —— `Φ(F)^gp = ⊕_v ord(F_v)`**。

★非アルキメデス成分は `ord(F_v) ≅ ℤ`(`ArithDivisor.lean` の `ordFin`)、
アルキメデス成分は `ord(F_v) ≅ ℝ`(本ファイル §1 の `ordArch`)。

原文 (FrdI p.113):
> as an arithmetic divisor on F. Thus, there is a natural homomorphism of groups -/
noncomputable def arithDivGroupEquiv :
    ((FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)) ≃+ arithDivGroup L :=
  AddEquiv.ofBijective
    ((arithOfParts (L := L)).codRestrict (arithDivGroup L).toAddSubmonoid arithOfParts_mem)
    ⟨fun p q h => arithOfParts_injective (congrArg Subtype.val h),
     fun d => by
       obtain ⟨p, hp⟩ := arithOfParts_surjective d.2
       exact ⟨p, Subtype.ext hp⟩⟩

/-- ★★★★★**[FrdI] Example 6.3 —— `Φ(F) = ⊕_v ord(O_v^▷)`**。

★有効算術因子はちょうど「非アルキメデス成分が `ℤ≥0`、アルキメデス成分が `ℝ≥0`」である。
`ord(O_v^▷)` は非アルキメデスで `ℤ≥0`(`ordFin ≥ 0`)、
アルキメデスで `ℝ≥0`(`ordArch_nonneg_iff`)だから、これが原文の直和分解である。

原文 (FrdI p.113):
> as an effective arithmetic divisor on F, and to an element of the group -/
theorem arithOfParts_mem_arithEff_iff (p : (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)) :
    arithOfParts p ∈ arithEff L ↔ ((∀ w, 0 ≤ p.1 w) ∧ ∀ w, 0 ≤ p.2 w) := by
  constructor
  · intro h
    refine ⟨fun w => ?_, fun w => ?_⟩
    · have h1 := h.2 (Sum.inl w)
      rw [arithOfParts_inl] at h1
      by_contra hneg
      push_neg at hneg
      have hlt : ((p.1 w : ℤ) : ℝ) < 0 := by exact_mod_cast hneg
      nlinarith [logAbsNorm_pos w]
    · have h1 := h.2 (Sum.inr w)
      rwa [arithOfParts_inr] at h1
  · rintro ⟨h1, h2⟩
    refine ⟨arithOfParts_mem p, fun x => ?_⟩
    cases x with
    | inl w =>
        rw [arithOfParts_inl]
        exact mul_nonneg (by exact_mod_cast h1 w) (logAbsNorm_pos w).le
    | inr w =>
        rw [arithOfParts_inr]
        exact h2 w

/-! ### ★出典の紐付け -/

/-- ★★locator —— `Example 6.3` のアルキメデス側 `ord(F_v) ≅ ℝ`、`ord(O_v^▷) ≅ ℝ≥0`。 -/
def ordArch.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 112,
    item := "Example 6.3 — ord(F_v) ≅ ℝ(アルキメデス)",
    sectionId := "frdi-example-6-3" }

/-- ★★★locator —— `Example 6.3` の `Φ(F)^gp = ⊕_v ord(F_v)`。 -/
def arithDivGroupEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(F)^gp = ⊕_v ord(F_v)(直和分解)",
    sectionId := "frdi-example-6-3" }

/-- ★★★locator —— `Example 6.3` の `Φ(F) = ⊕_v ord(O_v^▷)`。 -/
def arithOfParts_mem_arithEff_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(F) = ⊕_v ord(O_v^▷)(有効部分)",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
