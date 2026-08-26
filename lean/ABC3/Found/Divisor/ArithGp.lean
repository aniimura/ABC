/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithPhiMonoid
import ABC3.Found.Divisor.CartierMonoid

/-!
# `Φ(L)^gp ≃ Γ`(実係数版)—— `Example 6.3` の `divB` の受け皿

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★幾何版との対応

`CartierMonoid.lean` は `ℤ` 係数で
「`Q`-Cartier なら `Φ^gp ≃ Γ`」(`effSubGpEquiv`)を示した。
★`Example 6.3` は係数が `ℝ` なので、同じことを `effR` について書く必要がある。

## ★★★`ℤ` 係数との違いは 1 か所だけ

`γ ∈ Γ` を有効な 2 元の差に書くとき、`ℤ` 係数では
「各素点の生成元 `n_s·s` を `max 0 (−γ s)` 倍する」で足りた。
★★**実係数では `Γ` が `ℝ`-スカラー倍で閉じていない**(部分**群**なので `ℤ`-倍だけ)ので、
**天井関数**を挟む:

  `b := Σ_{s ∈ supp γ} ⌈max 0 (−γ s) / c_s⌉ · (c_s · s)`

★`⌈x⌉ ≥ x` だから `b s ≥ max 0 (−γ s)`、したがって `γ + b ≥ 0` になる。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `exists_effR_sub` | ★`Γ` の元は有効な 2 元の差 |
| `effRGpHom` / `effRGpHom_injective` | `Φ^gp → Γ` は単射 |
| `effRGpEquiv` | ★★★**`Φ(L)^gp ≃ Γ`**(実係数版) |
-/

namespace ABC3.Found.Divisor

open ABC3.Found.FrdI Finsupp

variable {S : Type*}

/-! ## ★1. `Γ` の元は有効な 2 元の差 -/

/-- ★★★**`Γ` が各素点で正の元を持てば、`Γ` の任意の元は有効な 2 元の差**(実係数版)。

★★`ℤ` 係数版(`exists_effSub_sub`)との違いは**天井関数**だけである ——
`Γ` は `ℝ`-スカラー倍で閉じていないので、生成元 `c_s·s` の**整数**倍を取る。 -/
theorem exists_effR_sub (Γ : AddSubgroup (S →₀ ℝ)) (hG : IsGenSubgroupR Γ)
    (γ : S →₀ ℝ) (hγ : γ ∈ Γ) :
    ∃ a b : S →₀ ℝ, a ∈ effR Γ ∧ b ∈ effR Γ ∧ γ = a - b := by
  classical
  set c : S → ℝ := fun s => (hG s).choose with hc
  have hcpos : ∀ s, 0 < c s := fun s => (hG s).choose_spec.1
  have hcmem : ∀ s, (single s (c s)) ∈ Γ := fun s => (hG s).choose_spec.2
  set n : S → ℕ := fun s => ⌈max 0 (- γ s) / c s⌉₊ with hn
  set b : S →₀ ℝ := γ.support.sum (fun s => single s ((n s : ℝ) * c s)) with hb
  have hbapp : ∀ t, b t = if t ∈ γ.support then (n t : ℝ) * c t else 0 := by
    intro t
    rw [hb, Finsupp.finsetSum_apply]
    simp only [Finsupp.single_apply]
    exact Finset.sum_ite_eq' γ.support t (fun s => (n s : ℝ) * c s)
  have hbmem : b ∈ Γ := by
    rw [hb]
    refine AddSubgroup.sum_mem _ (fun s _ => ?_)
    have hsm : single s ((n s : ℝ) * c s) = (n s) • (single s (c s)) := by
      rw [Finsupp.smul_single, nsmul_eq_mul]
    rw [hsm]
    exact AddSubgroup.nsmul_mem _ (hcmem s) _
  have hkey : ∀ t, max 0 (- γ t) ≤ (n t : ℝ) * c t := by
    intro t
    have h1 : max 0 (- γ t) / c t ≤ (n t : ℝ) := Nat.le_ceil _
    exact (div_le_iff₀ (hcpos t)).mp h1
  have hbnn : ∀ t, 0 ≤ b t := by
    intro t
    rw [hbapp t]
    split
    · exact le_trans (le_max_left _ _) (hkey t)
    · exact le_refl 0
  have hgbnn : ∀ t, 0 ≤ γ t + b t := by
    intro t
    by_cases ht : t ∈ γ.support
    · have hbt : b t = (n t : ℝ) * c t := by simp [hbapp t, ht]
      have h2 : - γ t ≤ (n t : ℝ) * c t := le_trans (le_max_right _ _) (hkey t)
      rw [hbt]
      linarith
    · have h0 : γ t = 0 := Finsupp.notMem_support_iff.mp ht
      have h2 := hbnn t
      rw [h0]
      linarith
  refine ⟨γ + b, b, mem_effR.mpr ⟨Γ.add_mem hγ hbmem, ?_⟩, mem_effR.mpr ⟨hbmem, hbnn⟩,
    (add_sub_cancel_right γ b).symm⟩
  intro t
  simpa using hgbnn t

/-! ## ★2. `Φ^gp ≃ Γ` -/

/-- ★`Φ = Γ ∩ ℝ≥0[S] ↪ Γ`。 -/
def effRIncl (Γ : AddSubgroup (S →₀ ℝ)) : effR Γ →+ Γ where
  toFun x := ⟨(x : S →₀ ℝ), (mem_effR.mp x.2).1⟩
  map_zero' := rfl
  map_add' _ _ := rfl

/-- ★★`Φ^gp →+ Γ`。 -/
noncomputable def effRGpHom (Γ : AddSubgroup (S →₀ ℝ)) : Gp (effR Γ) →+ Γ :=
  gpLift (effRIncl Γ)

theorem effRGpHom_injective (Γ : AddSubgroup (S →₀ ℝ)) :
    Function.Injective (effRGpHom Γ) := by
  rw [injective_iff_map_eq_zero]
  intro z hz
  obtain ⟨a, b, rfl⟩ := gp_eq_sub z
  rw [effRGpHom, map_sub, gpLift_toGp, gpLift_toGp, sub_eq_zero] at hz
  have h1 : (a : S →₀ ℝ) = (b : S →₀ ℝ) := congrArg (fun t : Γ => (t : S →₀ ℝ)) hz
  rw [Subtype.ext h1, sub_self]

/-- ★★★★**`Γ` が各素点で正の元を持てば `Φ^gp ≃ Γ`**(実係数版)。

★★これが `Example 6.3` の `divB : B(L) → Φ(L)^gp` の**受け皿**である ——
`arithDiv : L^× → arithDivGroup L` をこの同型の逆で `Φ^gp` へ運ぶ。 -/
noncomputable def effRGpEquiv (Γ : AddSubgroup (S →₀ ℝ)) (hG : IsGenSubgroupR Γ) :
    Gp (effR Γ) ≃+ Γ :=
  AddEquiv.ofBijective (effRGpHom Γ) ⟨effRGpHom_injective Γ, by
    rintro ⟨γ, hγ⟩
    obtain ⟨a, b, ha, hb, hab⟩ := exists_effR_sub Γ hG γ hγ
    refine ⟨toGp _ ⟨a, ha⟩ - toGp _ ⟨b, hb⟩, ?_⟩
    rw [effRGpHom, map_sub, gpLift_toGp, gpLift_toGp]
    exact Subtype.ext hab.symm⟩

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.3` の `Φ(L)^gp` と算術因子の群の同一視。 -/
def effRGpEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(L)^gp は算術因子の群と同一視される",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
