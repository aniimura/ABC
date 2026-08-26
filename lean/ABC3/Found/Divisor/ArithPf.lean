/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithMonoprime

/-!
# 実係数の `Φ^pf` —— `Φ^pf ≃ Γ^ℚ ∩ ℝ≥0[S]`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> as an effective arithmetic divisor on F, and to an element of the group

## ★★幾何版との違い

`Example 6.1`(幾何)では `Φ(L)^pf = ℚ≥0[D_L]` であり、
`CartierPf.lean` は `ℤ[S] → ℚ[S]` の係数拡大として書いた。

★算術では**係数がはじめから `ℝ`** なので係数拡大は要らない。
代わりに、完全化がやることは「`Γ` を **`ℚ` 倍で飽和させる**」ことである:

  `Γ^ℚ := {x | ∃ a : ℕ≥1, a·x ∈ Γ}`

非アルキメデス `v` では `ℤ·log(N v) ↦ ℚ·log(N v)`、
アルキメデス `v` では `ℝ ↦ ℝ`(何も増えない、`ℝ` はすでに可除)である。
★★これが原文の「`ord(O_v^▷) ≅ ℤ≥0` は完全化で `ℚ≥0` になるが、
`≅ ℝ≥0` の方は変わらない」の中身である。

## ★段取り

| 段 | 中身 |
|---|---|
| `satQ` | `Γ` の `ℚ` 飽和化(部分群) |
| `pfCoeffRHom` | `Φ^pf →+ (S →₀ ℝ)`、`m/a ↦ (1/a)·m` |
| 単射 | `ℝ` に捻れが無い |
| 像 | ちょうど `effR (satQ Γ)` |
-/

namespace ABC3.Found.Divisor

open Finsupp ABC3.Found.FrdI

variable {S : Type*} {Γ : AddSubgroup (S →₀ ℝ)}

/-! ## ★1. `Γ` の `ℚ` 飽和化 -/

/-- ★★**`Γ^ℚ := {x | ∃ a : ℕ≥1, a·x ∈ Γ}`** —— `Γ` を `ℚ` 倍で飽和させたもの。 -/
def satQ (Γ : AddSubgroup (S →₀ ℝ)) : AddSubgroup (S →₀ ℝ) where
  carrier := {x | ∃ a : ℕ+, ((a : ℕ)) • x ∈ Γ}
  add_mem' := by
    rintro x y ⟨a, ha⟩ ⟨b, hb⟩
    refine ⟨a * b, ?_⟩
    have h : (((a * b : ℕ+) : ℕ)) • (x + y)
        = (b : ℕ) • ((a : ℕ) • x) + (a : ℕ) • ((b : ℕ) • y) := by
      ext s
      simp only [PNat.mul_coe, Finsupp.add_apply, Finsupp.smul_apply, nsmul_eq_mul,
        smul_add, smul_smul]
      push_cast
      ring
    rw [h]
    exact Γ.add_mem (Γ.nsmul_mem ha _) (Γ.nsmul_mem hb _)
  zero_mem' := ⟨1, by simpa using Γ.zero_mem⟩
  neg_mem' := by
    rintro x ⟨a, ha⟩
    refine ⟨a, ?_⟩
    rw [smul_neg]
    exact Γ.neg_mem ha

theorem mem_satQ {x : S →₀ ℝ} : x ∈ satQ Γ ↔ ∃ a : ℕ+, ((a : ℕ)) • x ∈ Γ := Iff.rfl

theorem le_satQ : Γ ≤ satQ Γ := fun x hx => ⟨1, by simpa using hx⟩

/-! ## ★2. `Φ^pf → (S →₀ ℝ)` -/

/-- ★★**`Φ^pf → (S →₀ ℝ)`** —— `m/a ↦ (1/a)·m`。 -/
noncomputable def pfCoeffR : Pf (effR Γ) → (S →₀ ℝ) :=
  Quotient.lift (fun p : effR Γ × ℕ+ =>
      (((p.2 : ℕ) : ℝ))⁻¹ • ((p.1 : effR Γ) : S →₀ ℝ)) (by
    rintro ⟨m, a⟩ ⟨m', a'⟩ ⟨k, hk⟩
    ext s
    simp only [Finsupp.smul_apply, smul_eq_mul]
    have hkR : ((k : ℕ) * (a' : ℕ)) • ((m : effR Γ) : S →₀ ℝ)
        = ((k : ℕ) * (a : ℕ)) • ((m' : effR Γ) : S →₀ ℝ) :=
      congrArg (fun t : effR Γ => (t : S →₀ ℝ)) hk
    have h1 := congrArg (fun t : S →₀ ℝ => t s) hkR
    simp only [Finsupp.smul_apply, nsmul_eq_mul] at h1
    have ha : ((a : ℕ) : ℝ) ≠ 0 := by positivity
    have ha' : ((a' : ℕ) : ℝ) ≠ 0 := by positivity
    have hkne : ((k : ℕ) : ℝ) ≠ 0 := by positivity
    push_cast at h1
    have h2 : ((a' : ℕ) : ℝ) * ((m : effR Γ) : S →₀ ℝ) s
        = ((a : ℕ) : ℝ) * ((m' : effR Γ) : S →₀ ℝ) s := by
      refine mul_left_cancel₀ hkne ?_
      rw [← mul_assoc, ← mul_assoc]
      exact h1
    rw [inv_mul_eq_div, inv_mul_eq_div, div_eq_div_iff ha ha']
    rw [mul_comm] at h2
    rw [h2, mul_comm])

@[simp] theorem pfCoeffR_mk (m : effR Γ) (a : ℕ+) :
    pfCoeffR (Pf.mk m a) = (((a : ℕ) : ℝ))⁻¹ • ((m : effR Γ) : S →₀ ℝ) := rfl

/-- ★★`Φ^pf →+ (S →₀ ℝ)`。 -/
noncomputable def pfCoeffRHom : Pf (effR Γ) →+ (S →₀ ℝ) where
  toFun := pfCoeffR
  map_zero' := by
    show pfCoeffR (Pf.mk (0 : effR Γ) 1) = 0
    rw [pfCoeffR_mk]
    simp
  map_add' x y := by
    induction x using Pf.inductionOn with | _ m a =>
    induction y using Pf.inductionOn with | _ m' a' =>
    rw [Pf.mk_add_mk, pfCoeffR_mk, pfCoeffR_mk, pfCoeffR_mk]
    ext s
    have ha : ((a : ℕ) : ℝ) ≠ 0 := by positivity
    have ha' : ((a' : ℕ) : ℝ) ≠ 0 := by positivity
    have hcoe : (((((a' : ℕ) • m + (a : ℕ) • m' : effR Γ)) : S →₀ ℝ))
        = (a' : ℕ) • ((m : effR Γ) : S →₀ ℝ) + (a : ℕ) • ((m' : effR Γ) : S →₀ ℝ) := rfl
    simp only [Finsupp.add_apply, Finsupp.smul_apply, smul_eq_mul, hcoe,
      PNat.mul_coe, nsmul_eq_mul]
    push_cast
    field_simp

@[simp] theorem pfCoeffRHom_mk (m : effR Γ) (a : ℕ+) :
    (pfCoeffRHom (Pf.mk m a) : S →₀ ℝ)
      = (((a : ℕ) : ℝ))⁻¹ • ((m : effR Γ) : S →₀ ℝ) := rfl

/-- ★★`Φ^pf → (S →₀ ℝ)` は**単射** —— `ℝ` に捻れが無いから。 -/
theorem pfCoeffRHom_injective : Function.Injective (pfCoeffRHom (Γ := Γ)) := by
  intro x y h
  induction x using Pf.inductionOn with | _ m a =>
  induction y using Pf.inductionOn with | _ m' a' =>
  rw [pfCoeffRHom_mk, pfCoeffRHom_mk] at h
  have ha : ((a : ℕ) : ℝ) ≠ 0 := by positivity
  have ha' : ((a' : ℕ) : ℝ) ≠ 0 := by positivity
  refine Pf.sound 1 (Subtype.ext ?_)
  have hcoeL : ((((1 : ℕ+) : ℕ) * (a' : ℕ)) • m : effR Γ).1
      = (((1 : ℕ+) : ℕ) * (a' : ℕ)) • ((m : effR Γ) : S →₀ ℝ) := rfl
  have hcoeR : ((((1 : ℕ+) : ℕ) * (a : ℕ)) • m' : effR Γ).1
      = (((1 : ℕ+) : ℕ) * (a : ℕ)) • ((m' : effR Γ) : S →₀ ℝ) := rfl
  rw [hcoeL, hcoeR]
  ext s
  have h1 := congrArg (fun t : S →₀ ℝ => t s) h
  simp only [Finsupp.smul_apply, smul_eq_mul] at h1
  rw [inv_mul_eq_div, inv_mul_eq_div, div_eq_div_iff ha ha'] at h1
  simp only [Finsupp.smul_apply, nsmul_eq_mul, PNat.one_coe, one_mul]
  rw [mul_comm, mul_comm (((a : ℕ) : ℝ)) _]
  exact h1

theorem pfCoeffR_nonneg (x : Pf (effR Γ)) : 0 ≤ pfCoeffRHom x := by
  induction x using Pf.inductionOn with | _ m a =>
  rw [pfCoeffRHom_mk]
  refine Finsupp.le_def.mpr (fun s => ?_)
  have h1 : (0 : ℝ) ≤ ((m : effR Γ) : S →₀ ℝ) s := effR_nonneg m s
  have h3 : (0 : ℝ) ≤ (((a : ℕ)) : ℝ)⁻¹ := by positivity
  simp only [Finsupp.coe_zero, Pi.zero_apply, Finsupp.smul_apply, smul_eq_mul]
  exact mul_nonneg h3 h1

/-! ## ★3. 像はちょうど `Γ^ℚ ∩ ℝ≥0[S]` -/

theorem pfCoeffR_mem_satQ (x : Pf (effR Γ)) : (pfCoeffRHom x : S →₀ ℝ) ∈ satQ Γ := by
  induction x using Pf.inductionOn with | _ m a =>
  rw [pfCoeffRHom_mk]
  refine ⟨a, ?_⟩
  have ha : ((a : ℕ) : ℝ) ≠ 0 := by positivity
  have h : ((a : ℕ)) • ((((a : ℕ) : ℝ))⁻¹ • ((m : effR Γ) : S →₀ ℝ))
      = ((m : effR Γ) : S →₀ ℝ) := by
    ext s
    simp only [Finsupp.smul_apply, smul_eq_mul, nsmul_eq_mul]
    field_simp
  rw [h]
  exact (mem_effR.mp m.2).1

theorem pfCoeffR_mem (x : Pf (effR Γ)) : (pfCoeffRHom x : S →₀ ℝ) ∈ effR (satQ Γ) :=
  mem_effR.mpr ⟨pfCoeffR_mem_satQ x, fun s => Finsupp.le_def.mp (pfCoeffR_nonneg x) s⟩

/-- ★★★**像はちょうど `Γ^ℚ ∩ ℝ≥0[S]`**。

★`x ∈ effR (satQ Γ)` なら `a·x ∈ Γ` かつ `a·x ≥ 0`、すなわち `a·x ∈ Φ`。
そこで `pfCoeffR ((a·x)/a) = x` になる。 -/
theorem pfCoeffRHom_surjective {x : S →₀ ℝ} (hx : x ∈ effR (satQ Γ)) :
    ∃ y : Pf (effR Γ), pfCoeffRHom y = x := by
  obtain ⟨⟨a, ha⟩, hnn⟩ := mem_effR.mp hx
  have hmem : ((a : ℕ)) • x ∈ effR Γ := by
    refine mem_effR.mpr ⟨ha, fun s => ?_⟩
    have := hnn s
    simp only [Finsupp.smul_apply, nsmul_eq_mul]
    positivity
  refine ⟨Pf.mk ⟨_, hmem⟩ a, ?_⟩
  rw [pfCoeffRHom_mk]
  have haR : ((a : ℕ) : ℝ) ≠ 0 := by positivity
  ext s
  show (((a : ℕ) : ℝ))⁻¹ * (((a : ℕ)) • x) s = x s
  simp only [Finsupp.smul_apply, nsmul_eq_mul]
  field_simp

/-- ★★★★★**`Φ^pf ≃ Γ^ℚ ∩ ℝ≥0[S]`** —— 実係数版の
[FrdI] `Example 6.1` の `Φ(L)^pf = ℚ≥0[D_L]` に当たるもの。 -/
noncomputable def pfEquivEffRSatQ : Pf (effR Γ) ≃+ effR (satQ Γ) :=
  AddEquiv.ofBijective
    (pfCoeffRHom.codRestrict (effR (satQ Γ)) pfCoeffR_mem)
    ⟨fun x y h => pfCoeffRHom_injective (congrArg Subtype.val h), by
      rintro ⟨y, hy⟩
      obtain ⟨x, hx⟩ := pfCoeffRHom_surjective hy
      exact ⟨x, Subtype.ext hx⟩⟩

/-! ### ★出典の紐付け -/

/-- ★★locator —— `Example 6.3` の完全化 `Φ(L)^pf`。 -/
def pfEquivEffRSatQ.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(L)^pf は Γ の ℚ 飽和化の有効部分",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
