/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.CartierPrime
import Mathlib.Data.Rat.Cast.Order

/-!
# `Φ(L)^pf = ℚ≥0[D_L]`

★★★[FrdI] `Example 6.1` は

> Φ(L)pf = Q≥0[DL] ⊆ Q[DL] = (Φ(L)pf)gp
> [since DK is K-Q-Cartier]

と書く。★**これをここで閉じる。**

`Φ := Γ ∩ ℤ≥0[S]` の完全化 `Φ^pf`(形式的な分数 `m/a`)から `ℚ[S]` への写像

* `pfCoeffHom : Φ^pf →+ ℚ[S]` —— `m/a ↦ (1/a)·m`

を作り、**単射**であること(`ℤ[S]` に捻れが無い)と、**像がちょうど `ℚ≥0[S]`** で
あること(`Q`-Cartier 性)を示す。

## ★像が `ℚ≥0[S]` 全体になる理由

`Q`-Cartier だから各素因子 `s` について `n_s·s ∈ Γ`。したがって

* `pfCoeff (qcGen s / n_s) = 1·s` —— **素因子そのものが `Φ^pf` の中にある**
* `q = p/r ≥ 0` なら `pfCoeff ((p·qcGen s) / (n_s·r)) = q·s`
* 一般の `x ≥ 0` は `x = Σ_{s ∈ supp x} (x s)·s` と有限和に分かれ、
  像は部分単系だから閉じている

★逆向き(像が `ℚ≥0[S]` に入る)は `m ≥ 0` と `a > 0` から成分ごとに出る。
-/

namespace ABC3.Found.FrdI

open Finsupp

variable {S : Type*}

/-! ## ★1. `ℤ[S] → ℚ[S]` -/

/-- ★`ℤ[S] → ℚ[S]`。 -/
noncomputable def intToRatFs : (S →₀ ℤ) →+ (S →₀ ℚ) :=
  Finsupp.mapRange.addMonoidHom (Int.castAddHom ℚ)

@[simp] theorem intToRatFs_apply (x : S →₀ ℤ) (s : S) :
    (intToRatFs x) s = ((x s : ℤ) : ℚ) := rfl

theorem intToRatFs_injective : Function.Injective (intToRatFs (S := S)) := by
  intro x y h
  ext s
  have h1 := congrArg (fun t : S →₀ ℚ => t s) h
  simp only [intToRatFs_apply] at h1
  exact_mod_cast h1

variable {Γ : AddSubgroup (S →₀ ℤ)}

/-! ## ★2. `Φ^pf → ℚ[S]` -/

/-- ★★**`Φ^pf → ℚ[S]`** —— `m/a ↦ (1/a)·m`。

★well-defined 性: `k·a'·m = k·a·m'` の両辺を `ℚ` に送り、`k` を消して
`m/a = m'/a'` にする。 -/
noncomputable def pfCoeff : Pf (effSub Γ) → (S →₀ ℚ) :=
  Quotient.lift (fun p : effSub Γ × ℕ+ =>
      (((p.2 : ℕ) : ℚ))⁻¹ • intToRatFs ((p.1 : effSub Γ) : S →₀ ℤ)) (by
    rintro ⟨m, a⟩ ⟨m', a'⟩ ⟨k, hk⟩
    ext s
    simp only [Finsupp.smul_apply, smul_eq_mul, intToRatFs_apply]
    have hkZ : ((k : ℕ) * (a' : ℕ)) • ((m : effSub Γ) : S →₀ ℤ)
        = ((k : ℕ) * (a : ℕ)) • ((m' : effSub Γ) : S →₀ ℤ) :=
      congrArg (fun t : effSub Γ => (t : S →₀ ℤ)) hk
    have h1 := congrArg (fun t : S →₀ ℤ => t s) hkZ
    simp only [Finsupp.smul_apply, nsmul_eq_mul] at h1
    have ha : ((a : ℕ) : ℚ) ≠ 0 := by positivity
    have ha' : ((a' : ℕ) : ℚ) ≠ 0 := by positivity
    have hkne : ((k : ℕ) : ℚ) ≠ 0 := by positivity
    have h1Q : ((k : ℕ) : ℚ) * (((a' : ℕ) : ℚ) * ((((m : effSub Γ) : S →₀ ℤ) s : ℤ) : ℚ))
        = ((k : ℕ) : ℚ) * (((a : ℕ) : ℚ) * ((((m' : effSub Γ) : S →₀ ℤ) s : ℤ) : ℚ)) := by
      have h2 := congrArg (fun z : ℤ => (z : ℚ)) h1
      push_cast at h2 ⊢
      rw [← mul_assoc, ← mul_assoc]
      exact h2
    have h2 := mul_left_cancel₀ hkne h1Q
    rw [inv_mul_eq_div, inv_mul_eq_div, div_eq_div_iff ha ha']
    rw [mul_comm] at h2
    rw [h2, mul_comm])

@[simp] theorem pfCoeff_mk (m : effSub Γ) (a : ℕ+) :
    pfCoeff (Pf.mk m a) = (((a : ℕ) : ℚ))⁻¹ • intToRatFs ((m : effSub Γ) : S →₀ ℤ) := rfl

/-- ★★`Φ^pf →+ ℚ[S]`。 -/
noncomputable def pfCoeffHom : Pf (effSub Γ) →+ (S →₀ ℚ) where
  toFun := pfCoeff
  map_zero' := by
    show pfCoeff (Pf.mk (0 : effSub Γ) 1) = 0
    rw [pfCoeff_mk]
    simp
  map_add' x y := by
    induction x using Pf.inductionOn with | _ m a =>
    induction y using Pf.inductionOn with | _ m' a' =>
    rw [Pf.mk_add_mk, pfCoeff_mk, pfCoeff_mk, pfCoeff_mk]
    ext s
    have ha : ((a : ℕ) : ℚ) ≠ 0 := by positivity
    have ha' : ((a' : ℕ) : ℚ) ≠ 0 := by positivity
    have hcoe : (((((a' : ℕ) • m + (a : ℕ) • m' : effSub Γ)) : S →₀ ℤ))
        = (a' : ℕ) • ((m : effSub Γ) : S →₀ ℤ) + (a : ℕ) • ((m' : effSub Γ) : S →₀ ℤ) := rfl
    simp only [Finsupp.add_apply, Finsupp.smul_apply, smul_eq_mul, intToRatFs_apply, hcoe,
      PNat.mul_coe, nsmul_eq_mul]
    push_cast
    field_simp

@[simp] theorem pfCoeffHom_mk (m : effSub Γ) (a : ℕ+) :
    (pfCoeffHom (Pf.mk m a) : S →₀ ℚ)
      = (((a : ℕ) : ℚ))⁻¹ • intToRatFs ((m : effSub Γ) : S →₀ ℤ) := rfl

/-- ★★`Φ^pf → ℚ[S]` は**単射** —— `ℤ[S]` に捻れが無いから。 -/
theorem pfCoeffHom_injective : Function.Injective (pfCoeffHom (Γ := Γ)) := by
  intro x y h
  induction x using Pf.inductionOn with | _ m a =>
  induction y using Pf.inductionOn with | _ m' a' =>
  rw [pfCoeffHom_mk, pfCoeffHom_mk] at h
  have ha : ((a : ℕ) : ℚ) ≠ 0 := by positivity
  have ha' : ((a' : ℕ) : ℚ) ≠ 0 := by positivity
  refine Pf.sound 1 (Subtype.ext ?_)
  have hcoeL : ((((1 : ℕ+) : ℕ) * (a' : ℕ)) • m : effSub Γ).1
      = (((1 : ℕ+) : ℕ) * (a' : ℕ)) • ((m : effSub Γ) : S →₀ ℤ) := rfl
  have hcoeR : ((((1 : ℕ+) : ℕ) * (a : ℕ)) • m' : effSub Γ).1
      = (((1 : ℕ+) : ℕ) * (a : ℕ)) • ((m' : effSub Γ) : S →₀ ℤ) := rfl
  rw [hcoeL, hcoeR]
  ext s
  have h1 := congrArg (fun t : S →₀ ℚ => t s) h
  simp only [Finsupp.smul_apply, smul_eq_mul, intToRatFs_apply] at h1
  rw [inv_mul_eq_div, inv_mul_eq_div, div_eq_div_iff ha ha'] at h1
  simp only [Finsupp.smul_apply, nsmul_eq_mul, PNat.one_coe, one_mul]
  have h2 : ((((a' : ℕ) : ℤ) * ((m : effSub Γ) : S →₀ ℤ) s : ℤ) : ℚ)
      = ((((a : ℕ) : ℤ) * ((m' : effSub Γ) : S →₀ ℤ) s : ℤ) : ℚ) := by
    push_cast
    rw [mul_comm, mul_comm (((a : ℕ) : ℚ)) _]
    exact h1
  exact_mod_cast h2

/-! ## ★3. 像はちょうど `ℚ≥0[S]` -/

/-- ★`qcGen` の係数(正の自然数)。 -/
noncomputable def qcDen (hQ : IsQCartierSubgroup Γ) (s : S) : ℕ+ :=
  ⟨(hQ s).choose, (hQ s).choose_spec.1⟩

theorem qcGen_coe (hQ : IsQCartierSubgroup Γ) (s : S) :
    ((qcGen hQ s : effSub Γ) : S →₀ ℤ) = single s ((qcDen hQ s : ℕ) : ℤ) := rfl

/-- ★★**素因子 `s` は `Φ^pf` の中で「1 倍の `s`」として実現される**。 -/
theorem pfCoeff_qcGen (hQ : IsQCartierSubgroup Γ) (s : S) :
    pfCoeffHom (Pf.mk (qcGen hQ s) (qcDen hQ s)) = single s (1 : ℚ) := by
  rw [pfCoeffHom_mk, qcGen_coe]
  have hne : (((qcDen hQ s : ℕ)) : ℚ) ≠ 0 := by positivity
  ext t
  simp only [Finsupp.smul_apply, smul_eq_mul, intToRatFs_apply]
  rcases eq_or_ne s t with rfl | hst
  · simp [hne]
  · simp [hst]

/-- ★★`q ≥ 0` なら `q·s` は `Φ^pf` の像に入る。 -/
theorem single_mem_range_pfCoeff (hQ : IsQCartierSubgroup Γ) (s : S) {q : ℚ} (hq : 0 ≤ q) :
    single s q ∈ Set.range (pfCoeffHom (Γ := Γ)) := by
  set p : ℕ := q.num.toNat with hp
  set r : ℕ+ := ⟨q.den, q.pos⟩ with hr
  refine ⟨Pf.mk (p • qcGen hQ s) (qcDen hQ s * r), ?_⟩
  rw [pfCoeffHom_mk]
  have hcoe : ((p • qcGen hQ s : effSub Γ) : S →₀ ℤ)
      = single s ((p : ℤ) * ((qcDen hQ s : ℕ) : ℤ)) := by
    show p • ((qcGen hQ s : effSub Γ) : S →₀ ℤ) = _
    rw [qcGen_coe, Finsupp.smul_single]
    simp
  rw [hcoe]
  have hN : (((qcDen hQ s : ℕ)) : ℚ) ≠ 0 := by positivity
  have hrne : ((q.den : ℚ)) ≠ 0 := by
    have := q.pos
    positivity
  have hpq : ((p : ℚ)) = (q.num : ℚ) := by
    rw [hp]
    have hnn : (0 : ℤ) ≤ q.num := Rat.num_nonneg.mpr hq
    exact_mod_cast Int.toNat_of_nonneg hnn
  have hnd : (q.num : ℚ) = q * (q.den : ℚ) := (div_eq_iff hrne).mp (Rat.num_div_den q)
  ext t
  simp only [Finsupp.smul_apply, smul_eq_mul, intToRatFs_apply]
  rcases eq_or_ne s t with rfl | hst
  · simp only [Finsupp.single_eq_same]
    push_cast
    rw [hpq]
    show (((qcDen hQ s : ℕ) : ℚ) * (q.den : ℚ))⁻¹ * ((q.num : ℚ) * ((qcDen hQ s : ℕ) : ℚ)) = q
    field_simp
    rw [hnd, mul_comm]
  · simp [hst]

theorem pfCoeff_nonneg (x : Pf (effSub Γ)) : 0 ≤ pfCoeffHom x := by
  induction x using Pf.inductionOn with | _ m a =>
  rw [pfCoeffHom_mk]
  refine Finsupp.le_def.mpr (fun s => ?_)
  have h1 : (0 : ℤ) ≤ ((m : effSub Γ) : S →₀ ℤ) s := effSub_nonneg m s
  have h2 : (0 : ℚ) ≤ ((((m : effSub Γ) : S →₀ ℤ) s : ℤ) : ℚ) := by exact_mod_cast h1
  have h3 : (0 : ℚ) ≤ (((a : ℕ)) : ℚ)⁻¹ := by positivity
  simp only [Finsupp.coe_zero, Pi.zero_apply, Finsupp.smul_apply, smul_eq_mul, intToRatFs_apply]
  exact mul_nonneg h3 h2

/-- ★★★★**`Φ^pf` の像はちょうど `ℚ≥0[D_L]`** —— [FrdI] `Example 6.1` の
`Φ(L)^pf = ℚ≥0[D_L] ⊆ ℚ[D_L] = (Φ(L)^pf)^gp`。 -/
theorem pfCoeffHom_range (hQ : IsQCartierSubgroup Γ) :
    Set.range (pfCoeffHom (Γ := Γ)) = {x : S →₀ ℚ | 0 ≤ x} := by
  classical
  refine Set.Subset.antisymm ?_ ?_
  · rintro _ ⟨x, rfl⟩
    exact pfCoeff_nonneg x
  · intro x hx
    have hxs : ∀ s, (0 : ℚ) ≤ x s := fun s => Finsupp.le_def.mp hx s
    have hsum : ∑ s ∈ x.support, single s (x s) = x := Finsupp.sum_single x
    have hmem : ∀ s ∈ x.support,
        single s (x s) ∈ AddMonoidHom.mrange (pfCoeffHom (Γ := Γ)) :=
      fun s _ => single_mem_range_pfCoeff hQ s (hxs s)
    have hall := AddSubmonoid.sum_mem (AddMonoidHom.mrange (pfCoeffHom (Γ := Γ))) hmem
    rw [hsum] at hall
    exact hall

/-! ## ★4. 同型として -/

variable (S) in
/-- ★`ℚ≥0[S] ⊆ ℚ[S]`。 -/
def nonnegRatFs : AddSubmonoid (S →₀ ℚ) where
  carrier := {x | 0 ≤ x}
  add_mem' ha hb := Set.mem_setOf.mpr (add_nonneg (Set.mem_setOf.mp ha) (Set.mem_setOf.mp hb))
  zero_mem' := Set.mem_setOf.mpr (le_refl 0)

theorem mem_nonnegRatFs {x : S →₀ ℚ} : x ∈ nonnegRatFs S ↔ 0 ≤ x := by
  simp [nonnegRatFs, AddSubmonoid.mem_mk]

/-- ★★★★★**`Φ^pf ≃ ℚ≥0[D_L]`** —— [FrdI] `Example 6.1` の
「`Φ(L)^pf = ℚ≥0[D_L]`(since `D_K` is `K`-`Q`-Cartier)」。 -/
noncomputable def pfEquivNonneg (hQ : IsQCartierSubgroup Γ) :
    Pf (effSub Γ) ≃+ nonnegRatFs S :=
  AddEquiv.ofBijective
    (pfCoeffHom.codRestrict (nonnegRatFs S) (fun x => mem_nonnegRatFs.mpr (pfCoeff_nonneg x)))
    ⟨fun x y h => pfCoeffHom_injective (congrArg Subtype.val h), by
      rintro ⟨y, hy⟩
      have hr : y ∈ Set.range (pfCoeffHom (Γ := Γ)) := by
        rw [pfCoeffHom_range hQ]
        exact mem_nonnegRatFs.mp hy
      obtain ⟨x, hx⟩ := hr
      exact ⟨x, Subtype.ext hx⟩⟩

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.1` の `Φ(L)^pf = ℚ≥0[D_L]`。 -/
def pfEquivNonneg.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Φ(L)^pf = ℚ≥0[D_L]",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.FrdI
