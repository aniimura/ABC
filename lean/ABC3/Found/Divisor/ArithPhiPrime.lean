/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithMonoprime
import ABC3.Found.Divisor.ArithOrd

/-!
# `Φ(L)` が実係数の枠に収まること —— `Prime(Φ(L)) ≃ V(L)`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★何をするか

`ArithPhiMonoid.lean` / `ArithMonoprime.lean` は
「`Γ ⊆ (S →₀ ℝ)` が各素点で正の元を持ち、局所群が離散か全体」という
**抽象的な枠**で議論した。ここではその枠に

  `S = V(L) = ArithPlace L`、`Γ = arithDivGroup L`

を当てはめる。★確かめるのは 2 つだけである:

* `IsGenSubgroupR` —— 非アルキメデス `w` で `log(N w) > 0`、アルキメデス `w` で `1`。
* `IsLocallyMonoprimeR` —— 非アルキメデス `w` で `Γ_w = ℤ·log(N w)`(**離散**)、
  アルキメデス `w` で `Γ_w = ℝ`(**全体**)。

★★これで [FrdI] `Example 6.3` の

* 「there is a natural bijection `Prime(Φ(L)) → V(L)`」
* 「the supports of elements of `Φ(L)` are precisely the finite subsets of `V(L)`」
* `Definition 2.4, (i)` の条件 (a)(divisorial)・(b)(`M_p` は monoprime)

が閉じる。★残るのは (c)(d)(因子分解写像)である。
-/

namespace ABC3.Found.Divisor

open NumberField Finsupp ABC3.Found.FrdI

universe u

variable {L : Type u} [Field L] [NumberField L]

/-! ## ★1. `Φ(L) = effR (arithDivGroup L)` -/

/-- ★`arithEff L` は抽象版 `effR` そのものである。 -/
theorem arithEff_eq_effR : arithEff L = effR (arithDivGroup L) := by
  ext d
  constructor
  · rintro ⟨h1, h2⟩
    exact mem_effR.mpr ⟨h1, h2⟩
  · intro h
    obtain ⟨h1, h2⟩ := mem_effR.mp h
    exact ⟨h1, h2⟩

/-! ## ★2. 枠の仮定を確かめる -/

/-- ★台が `{Sum.inl w}` の元の非アルキメデス条件。 -/
theorem single_inl_mem_arithDivGroup (w : FinitePlace L) {c : ℝ} {n : ℤ}
    (hc : c = (n : ℝ) * logAbsNorm w) :
    (single (Sum.inl w : ArithPlace L) c) ∈ arithDivGroup L := by
  classical
  intro w'
  rcases eq_or_ne w w' with rfl | hne
  · exact ⟨n, by simpa [logAbsNorm] using hc⟩
  · refine ⟨0, ?_⟩
    have : (single (Sum.inl w : ArithPlace L) c) (Sum.inl w') = 0 := by
      rw [Finsupp.single_apply, if_neg]
      exact fun h => hne (Sum.inl_injective h)
    rw [this]
    simp

/-- ★アルキメデス素点では条件は空虚に成り立つ。 -/
theorem single_inr_mem_arithDivGroup (w : InfinitePlace L) (c : ℝ) :
    (single (Sum.inr w : ArithPlace L) c) ∈ arithDivGroup L := by
  classical
  intro w'
  refine ⟨0, ?_⟩
  have : (single (Sum.inr w : ArithPlace L) c) (Sum.inl w') = 0 := by
    rw [Finsupp.single_apply, if_neg]
    intro h
    exact absurd h (by simp)
  rw [this]
  simp

/-- ★★**`arithDivGroup L` は各素点で正の元を持つ**。 -/
theorem isGenSubgroupR_arithDivGroup :
    IsGenSubgroupR (arithDivGroup L) := by
  intro s
  cases s with
  | inl w => exact ⟨logAbsNorm w, logAbsNorm_pos w, single_inl_mem_arithDivGroup w (n := 1) (by ring)⟩
  | inr w => exact ⟨1, one_pos, single_inr_mem_arithDivGroup w 1⟩

/-- ★★★**局所群は非アルキメデスで離散、アルキメデスで全体**。

★これが原文の `ord(O_v^▷) ≅ ℤ≥0`(非アルキメデス)/ `≅ ℝ≥0`(アルキメデス)である。 -/
theorem isLocallyMonoprimeR_arithDivGroup :
    IsLocallyMonoprimeR (arithDivGroup L) := by
  intro s
  cases s with
  | inl w =>
      refine Or.inl ⟨logAbsNorm w, logAbsNorm_pos w, fun c => ?_⟩
      constructor
      · intro hmem
        obtain ⟨n, hn⟩ := hmem w
        refine ⟨n, ?_⟩
        rw [Finsupp.single_eq_same] at hn
        simpa [logAbsNorm] using hn
      · rintro ⟨n, hn⟩
        exact single_inl_mem_arithDivGroup w hn
  | inr w => exact Or.inr (single_inr_mem_arithDivGroup w)

/-! ## ★3. 帰結 -/

/-- ★★★**[FrdI] Example 6.3 —— `Φ(L)` は divisorial**(`Definition 2.4, (i)` の (a))。 -/
theorem isDivisorial_arithEff : IsDivisorial (effR (arithDivGroup L)) :=
  isDivisorial_effR _

/-- ★★★★★**[FrdI] Example 6.3 —— `Prime(Φ(L)) ≃ V(L)`**。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model] -/
noncomputable def arithPrimeEquiv :
    Prime (effR (arithDivGroup L)) ≃ ArithPlace L :=
  effRPrimeEquiv isGenSubgroupR_arithDivGroup

/-- ★★★★**[FrdI] Example 6.3 —— 台は `V(L)` の有限部分集合ちょうど**。 -/
theorem exists_arithEff_support_eq (T : Finset (ArithPlace L)) :
    ∃ a : effR (arithDivGroup L), (a : ArithPlace L →₀ ℝ).support = T :=
  exists_effR_support_eq isGenSubgroupR_arithDivGroup T

/-- ★★★★★**[FrdI] Example 6.3 —— `M_p` は monoprime**
(`Definition 2.4, (i)` の (b))。

★非アルキメデス素点で `ℤ`-monoprime(`≃+ ℕ`)、
アルキメデス素点で **`ℝ`-monoprime**(`≃+ ℝ≥0`)。 -/
theorem isMonoprime_Mp_arithEff (p : Prime (effR (arithDivGroup L))) :
    IsMonoprime (Mp (effR (arithDivGroup L)) p) :=
  isMonoprime_Mp_effR isLocallyMonoprimeR_arithDivGroup p

/-- ★★**`Φ(L) ≠ 0`** —— 原文の「`Φ(L) ̸= 0`」。

★素点が 1 つでもあれば(数体には必ずアルキメデス素点がある)非零元がある。 -/
theorem exists_ne_zero_arithEff (w : InfinitePlace L) :
    ∃ a : effR (arithDivGroup L), a ≠ 0 := by
  refine ⟨genR (isGenSubgroupR_arithDivGroup (L := L)) (Sum.inr w), ?_⟩
  intro hz
  have hsup := genR_support (isGenSubgroupR_arithDivGroup (L := L)) (Sum.inr w)
  rw [(effR_eq_zero_iff _).mp hz] at hsup
  simp at hsup

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.3` の `Prime(Φ(L)) ≃ V(L)`。 -/
def arithPrimeEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Prime(Φ(L)) ≃ V(L)",
    sectionId := "frdi-example-6-3" }

/-- ★★locator —— `Example 6.3` の「台は `V(L)` の有限部分集合ちょうど」。 -/
def exists_arithEff_support_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(L) の元の台は V(L) の有限部分集合ちょうど",
    sectionId := "frdi-example-6-3" }

/-- ★★★locator —— `Example 6.3` の `M_p` は monoprime。 -/
def isMonoprime_Mp_arithEff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(L) の M_p は monoprime (Definition 2.4, (i) の (b))",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
