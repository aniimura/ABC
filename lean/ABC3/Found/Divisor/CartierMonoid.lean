/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.FreeDivisorial
import ABC3.Found.FrdI.MonoidTransport
import Mathlib.GroupTheory.Archimedean

/-!
# Cartier 有効因子の単系 `Φ(L) = Γ ∩ ℤ≥0[D_L]`

★★★[FrdI] `Example 6.1` の中核である。原文は

> Φ(L) ⊆ Z≥0[DL] ⊆ Z[DL]
> for the monoid of Cartier effective divisors D on V [L] with support in DL
> …
> Observe that Φ(L)gp ⊆ Z[DL] may be identified with the group of Cartier divisors
> on V [L] … one verifies immediately that Φ(L) is perf-factorial

と書く。**幾何を落として単系論だけ残すと、次の 2 つの定理になる。**

`S` を素因子の型、`Γ ≤ ℤ[S]` を **Cartier 因子の群**とし、`Φ := Γ ∩ ℤ≥0[S]` と置く。

1. **`Φ` は divisorial**(`isDivisorial_effSub`)。★**`Q`-Cartier 性は要らない** ——
   効くのは `Γ` が**部分群**であることだけ。saturated 性は
   「`n·(x−y) ≥ 0` かつ `n > 0` ⟹ `x−y ≥ 0`」、
   sharp 性は「成分ごとに `x s + y s = 0`、両方 `≥ 0`」から出る。
2. **`Q`-Cartier なら `Φ^gp ≃ Γ`**(`effSubGpEquiv`)。★ここで初めて
   「各素因子 `s` の正の倍数 `n·s` が `Γ` に入る」が効く ——
   `γ ∈ Γ` に対して `b := Σ_{s ∈ supp γ} (n_s · max 0 (−γ s))·s` を足せば
   `γ + b` が有効になり、`γ = (γ+b) − b` と書ける。

## ★幾何との対応(未接続)

`S = D_L`、`Γ =` (`V[L]` 上の Cartier 因子で台が `D_L` に入るもの) と読む。
`Γ` が部分群であること・`Q`-Cartier 性は、それぞれ
原文の「Cartier 因子は群をなす」「`D_K` は `K`-`Q`-Cartier」に当たる。
★**`V[L]`(相対正規化)と `D_L` の構成そのものは、まだスキーム層で書いていない。**

## ★副産物 —— `M^gp` の普遍性

`gpLift` / `gp_eq_sub` は一般の単系に使える語彙である
(`gpMap` は `Gp N` 行きしか作れず、**群 `G` 行き**が要ったので足した)。
-/

namespace ABC3.Found.FrdI

open Finsupp

/-! ## ★0. `M^gp` の普遍性 -/

/-- ★★**`M^gp` の普遍性** —— 群への準同型は一意に持ち上がる。

★既存の `gpMap` は `Gp M →+ Gp N` しか作れない。ここでは**群 `G` 行き**が要る。 -/
noncomputable def gpLift {M G : Type*} [AddCommMonoid M] [AddCommGroup G] (f : M →+ G) :
    Gp M →+ G :=
  Algebra.GrothendieckAddGroup.lift f

@[simp] theorem gpLift_toGp {M G : Type*} [AddCommMonoid M] [AddCommGroup G] (f : M →+ G)
    (m : M) : gpLift f (toGp M m) = f m := by
  have hb : toGp M m = (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) m := rfl
  rw [hb]
  simp only [gpLift, Algebra.GrothendieckAddGroup.lift, Equiv.coe_fn_mk]
  rw [AddSubmonoid.LocalizationMap.lift_eq]

/-- ★`M^gp` の元はつねに `toGp a − toGp b` の形。 -/
theorem gp_eq_sub {M : Type*} [AddCommMonoid M] (z : Gp M) :
    ∃ a b : M, z = toGp M a - toGp M b := by
  refine AddLocalization.induction_on z ?_
  rintro ⟨a, b⟩
  exact ⟨a, (b : M), by rw [eq_sub_iff_add_eq]; exact mk_add_toGp M a b⟩

variable {S : Type*}

/-! ## ★1. `Φ = Γ ∩ ℤ≥0[S]` -/

/-- ★★**`Γ` の有効な元のなす部分単系** —— [FrdI] `Example 6.1` の `Φ(L)`。 -/
def effSub (Γ : AddSubgroup (S →₀ ℤ)) : AddSubmonoid (S →₀ ℤ) where
  carrier := {x | x ∈ Γ ∧ 0 ≤ x}
  add_mem' ha hb := ⟨Γ.add_mem ha.1 hb.1, add_nonneg ha.2 hb.2⟩
  zero_mem' := ⟨Γ.zero_mem, le_refl 0⟩

theorem mem_effSub {Γ : AddSubgroup (S →₀ ℤ)} {x : S →₀ ℤ} :
    x ∈ effSub Γ ↔ x ∈ Γ ∧ ∀ s, 0 ≤ x s := by
  simp [effSub, AddSubmonoid.mem_mk, Finsupp.le_def]

/-! ## ★2. `Φ` は divisorial -/

/-- ★★`Γ ∩ ℤ≥0[S]` は sharp —— 成分ごとに `x s + y s = 0`、両方 `≥ 0`。 -/
theorem isSharp_effSub (Γ : AddSubgroup (S →₀ ℤ)) : IsSharp (effSub Γ) := by
  intro a ha
  obtain ⟨u, hu⟩ := ha
  have h0 : (a : S →₀ ℤ) + ((u.neg : effSub Γ) : S →₀ ℤ) = 0 := by
    have h := u.val_neg
    rw [hu] at h
    exact congrArg (fun t : effSub Γ => (t : S →₀ ℤ)) h
  have hane := (mem_effSub.mp a.2).2
  have hbne := (mem_effSub.mp (u.neg : effSub Γ).2).2
  refine Subtype.ext ?_
  ext s
  have h1 := congrArg (fun t : S →₀ ℤ => t s) h0
  simp only [Finsupp.add_apply, Finsupp.coe_zero, Pi.zero_apply] at h1
  have := hane s
  have := hbne s
  show (a : S →₀ ℤ) s = 0
  omega

/-- ★★`Γ ∩ ℤ≥0[S]` は saturated —— `n·(x−y) = w ≥ 0` かつ `n > 0` なら `x−y ≥ 0`、
しかも `Γ` は**群**だから `x−y ∈ Γ`。 -/
theorem isSaturated_effSub (Γ : AddSubgroup (S →₀ ℤ)) : IsSaturatedMonoid (effSub Γ) := by
  refine isSaturatedMonoid_of_cancel_of_nsmul_le _ ?_
  rintro x y n hn ⟨w, hw⟩
  have hw' : n • (y : S →₀ ℤ) + (w : S →₀ ℤ) = n • (x : S →₀ ℤ) :=
    congrArg (fun t : effSub Γ => (t : S →₀ ℤ)) hw
  have hmem : (x : S →₀ ℤ) - (y : S →₀ ℤ) ∈ effSub Γ := by
    refine mem_effSub.mpr ⟨Γ.sub_mem (mem_effSub.mp x.2).1 (mem_effSub.mp y.2).1, fun s => ?_⟩
    have h1 := congrArg (fun t : S →₀ ℤ => t s) hw'
    simp only [Finsupp.add_apply, Finsupp.smul_apply, nsmul_eq_mul] at h1
    have h2 := (mem_effSub.mp w.2).2 s
    have h3 : (0:ℤ) ≤ (n : ℤ) * ((x : S →₀ ℤ) s - (y : S →₀ ℤ) s) := by
      simp only [mul_sub]; omega
    have hn' : (0:ℤ) < (n : ℤ) := by exact_mod_cast hn
    simpa using nonneg_of_mul_nonneg_right h3 hn'
  exact ⟨⟨_, hmem⟩, Subtype.ext (by simp)⟩

/-- ★★★★**`Φ := Γ ∩ ℤ≥0[S]` は divisorial**。

★[FrdI] `Example 6.1` の `Φ(L) ⊆ ℤ≥0[D_L]`(`Γ` = Cartier 因子の群)がこれである。
★**`Q`-Cartier 性は要らない** —— 効くのは `Γ` が**部分群**であることだけ。 -/
theorem isDivisorial_effSub (Γ : AddSubgroup (S →₀ ℤ)) : IsDivisorial (effSub Γ) := by
  have hsharp := isSharp_effSub Γ
  exact ⟨⟨isIntegralMonoid_of_isCancelAdd _, isSaturated_effSub Γ,
    isOfCharacteristicType_of_isSharp _ hsharp⟩, hsharp⟩

/-! ## ★3. `Q`-Cartier 性と `Φ^gp ≃ Γ` -/

/-- ★★**`Γ` が `Q`-Cartier** —— 各素因子 `s` について、正の倍数 `n·s` が `Γ` に入る。

★[FrdI] `Example 6.1` の「`D_K` は `K`-`Q`-Cartier」を、`V[L]` を 1 つ固定して読んだもの。 -/
def IsQCartierSubgroup (Γ : AddSubgroup (S →₀ ℤ)) : Prop :=
  ∀ s : S, ∃ n : ℕ, 0 < n ∧ (single s (n : ℤ)) ∈ Γ

/-- ★★★**`Γ` が `Q`-Cartier なら、`Γ` の任意の元は有効な 2 元の差**。

★`γ` の台の各点 `s` で `n_s·s ∈ Γ` を取り、`b := Σ_s (n_s · max 0 (−γ s))·s` を足す。
`n_s ≥ 1` だから `γ + b ≥ 0` になる。 -/
theorem exists_effSub_sub (Γ : AddSubgroup (S →₀ ℤ)) (hQ : IsQCartierSubgroup Γ)
    (γ : S →₀ ℤ) (hγ : γ ∈ Γ) :
    ∃ a b : S →₀ ℤ, a ∈ effSub Γ ∧ b ∈ effSub Γ ∧ γ = a - b := by
  classical
  set N : S → ℕ := fun s => (hQ s).choose with hN
  have hNpos : ∀ s, 0 < N s := fun s => (hQ s).choose_spec.1
  have hNmem : ∀ s, (single s (N s : ℤ)) ∈ Γ := fun s => (hQ s).choose_spec.2
  set c : S → ℤ := fun s => (N s : ℤ) * max 0 (- γ s) with hc
  set b : S →₀ ℤ := γ.support.sum (fun s => single s (c s)) with hb
  have hbapp : ∀ t, b t = if t ∈ γ.support then c t else 0 := by
    intro t
    rw [hb, Finsupp.finsetSum_apply]
    simp only [Finsupp.single_apply]
    exact Finset.sum_ite_eq' γ.support t c
  have hbmem : b ∈ Γ := by
    rw [hb]
    refine AddSubgroup.sum_mem _ (fun s _ => ?_)
    have hsm : single s (c s) = (max 0 (- γ s)) • (single s (N s : ℤ)) := by
      rw [Finsupp.smul_single, hc]
      simp [mul_comm]
    rw [hsm]
    exact AddSubgroup.zsmul_mem _ (hNmem s) _
  have hcnn : ∀ t, 0 ≤ c t := fun t => mul_nonneg (by positivity) (le_max_left _ _)
  have hbnn : ∀ t, 0 ≤ b t := by
    intro t
    rw [hbapp t]
    split
    · exact hcnn t
    · exact le_refl 0
  have hgbnn : ∀ t, 0 ≤ γ t + b t := by
    intro t
    by_cases ht : t ∈ γ.support
    · have hbt : b t = c t := by simp [hbapp t, ht]
      rcases lt_or_ge (γ t) 0 with h | h
      · have hct : c t = (N t : ℤ) * (- γ t) := by
          simp [hc, max_eq_right (by omega : (0:ℤ) ≤ - γ t)]
        have h1 : (1 : ℤ) ≤ (N t : ℤ) := by exact_mod_cast hNpos t
        have key : (0:ℤ) ≤ ((N t : ℤ) - 1) * (- γ t) := mul_nonneg (by omega) (by omega)
        have expand : ((N t : ℤ) - 1) * (- γ t) = γ t + (N t : ℤ) * (- γ t) := by ring
        rw [hbt, hct, ← expand]
        exact key
      · have h2 := hcnn t
        rw [hbt]
        omega
    · have h0 : γ t = 0 := Finsupp.notMem_support_iff.mp ht
      have h2 := hbnn t
      omega
  refine ⟨γ + b, b, mem_effSub.mpr ⟨Γ.add_mem hγ hbmem, ?_⟩, mem_effSub.mpr ⟨hbmem, hbnn⟩,
    (add_sub_cancel_right γ b).symm⟩
  intro t
  simpa using hgbnn t

/-- ★`Φ = Γ ∩ ℤ≥0[S] ↪ Γ`。 -/
def effSubIncl (Γ : AddSubgroup (S →₀ ℤ)) : effSub Γ →+ Γ where
  toFun x := ⟨(x : S →₀ ℤ), (mem_effSub.mp x.2).1⟩
  map_zero' := rfl
  map_add' _ _ := rfl

/-- ★★★`Φ^gp →+ Γ`。 -/
noncomputable def effSubGpHom (Γ : AddSubgroup (S →₀ ℤ)) : Gp (effSub Γ) →+ Γ :=
  gpLift (effSubIncl Γ)

theorem effSubGpHom_injective (Γ : AddSubgroup (S →₀ ℤ)) :
    Function.Injective (effSubGpHom Γ) := by
  rw [injective_iff_map_eq_zero]
  intro z hz
  obtain ⟨a, b, rfl⟩ := gp_eq_sub z
  rw [effSubGpHom, map_sub, gpLift_toGp, gpLift_toGp, sub_eq_zero] at hz
  have h1 : (a : S →₀ ℤ) = (b : S →₀ ℤ) := congrArg (fun t : Γ => (t : S →₀ ℤ)) hz
  rw [Subtype.ext h1, sub_self]

/-- ★★★★**`Q`-Cartier なら `Φ^gp ≃ Γ`** ——
[FrdI] `Example 6.1` の「`Φ(L)^gp ⊆ ℤ[D_L]` は `V[L]` 上の Cartier 因子の群と
同一視できる」。 -/
noncomputable def effSubGpEquiv (Γ : AddSubgroup (S →₀ ℤ)) (hQ : IsQCartierSubgroup Γ) :
    Gp (effSub Γ) ≃+ Γ :=
  AddEquiv.ofBijective (effSubGpHom Γ) ⟨effSubGpHom_injective Γ, by
    rintro ⟨γ, hγ⟩
    obtain ⟨a, b, ha, hb, hab⟩ := exists_effSub_sub Γ hQ γ hγ
    refine ⟨toGp _ ⟨a, ha⟩ - toGp _ ⟨b, hb⟩, ?_⟩
    rw [effSubGpHom, map_sub, gpLift_toGp, gpLift_toGp]
    exact Subtype.ext hab.symm⟩

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.1` の `Φ(L)` は divisorial。 -/
def isDivisorial_effSub.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Φ(L) ⊆ ℤ≥0[D_L] は divisorial",
    sectionId := "frdi-example-6-1" }

/-- ★locator —— `Example 6.1` の `Φ(L)^gp = Cartier 因子の群`。 -/
def effSubGpEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Φ(L)^gp は V[L] 上の Cartier 因子の群",
    sectionId := "frdi-example-6-1" }

/-! ## ★★幾何側の non-dilating の骨組み(2026-08-21)

★★算術側(`Example 6.3`)は「primary 元で生成される」から出たが、
**幾何側は生成が成り立たない** —— `Φ(L)` は「台が `D_L` に入る **Cartier** 因子」で、
単独の素因子 `D` は Cartier とは限らない(`D_K` は `K`-`Q`-Cartier としか仮定していない)。

★★★そこで**置換版**の骨組みを置く: 引き戻しが**素因子の置換(係数 1)**で
誘導されるなら、non-dilating の仮定は置換が各素因子を固定することを言う。
★各素因子が primary 元を持つことは `K`-`Q`-Cartier(`qc`)から出る。 -/

/-- ★★台が 1 点の元は primary(ℤ 係数版)。 -/
theorem isPrimaryElt_effSub_single {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)}
    {s : S} {r : ℤ} (hr : 0 < r) (hmem : Finsupp.single s r ∈ effSub Γ) :
    IsPrimaryElt (⟨Finsupp.single s r, hmem⟩ : effSub Γ) := by
  refine ⟨?_, ?_⟩
  · intro h0
    have h : Finsupp.single s r = (0 : S →₀ ℤ) := congrArg Subtype.val h0
    rw [Finsupp.single_eq_zero] at h
    exact hr.ne' h
  · rintro b hb0 ⟨n, hn, c, hbc⟩
    have hcoe : (b : S →₀ ℤ) + (c : S →₀ ℤ) = n • (Finsupp.single s r : S →₀ ℤ) :=
      congrArg Subtype.val hbc
    have hbnn : ∀ t, 0 ≤ (b : S →₀ ℤ) t := fun t => by simpa using Finsupp.le_def.mp b.2.2 t
    have hcnn : ∀ t, 0 ≤ (c : S →₀ ℤ) t := fun t => by simpa using Finsupp.le_def.mp c.2.2 t
    have hbt : ∀ t, t ≠ s → (b : S →₀ ℤ) t = 0 := by
      intro t ht
      have h1 := congrArg (fun f : S →₀ ℤ => f t) hcoe
      have hz : (Finsupp.single s r : S →₀ ℤ) t = 0 := Finsupp.single_eq_of_ne ht
      simp only [Finsupp.add_apply, Finsupp.smul_apply] at h1
      rw [hz, smul_zero] at h1
      linarith [hbnn t, hcnn t]
    have hbsupp : (b : S →₀ ℤ) = Finsupp.single s ((b : S →₀ ℤ) s) := by
      refine Finsupp.ext fun t => ?_
      rcases eq_or_ne t s with rfl | ht
      · rw [Finsupp.single_eq_same]
      · rw [hbt t ht, Finsupp.single_eq_of_ne ht]
    have hbs : 0 < (b : S →₀ ℤ) s := by
      rcases lt_or_eq_of_le (hbnn s) with h | h
      · exact h
      · refine absurd (Subtype.ext ?_) hb0
        rw [hbsupp, ← h, Finsupp.single_zero]
        rfl
    obtain ⟨m, hm⟩ : ∃ m : ℕ, r ≤ (m : ℤ) * (b : S →₀ ℤ) s := ⟨r.toNat, by
      have h1 : (r.toNat : ℤ) = r := Int.toNat_of_nonneg hr.le
      nlinarith [hbs, h1]⟩
    have hm0 : 0 < m := by
      rcases Nat.eq_zero_or_pos m with h | h
      · rw [h] at hm; simp at hm; omega
      · exact h
    have hdmem : Finsupp.single s ((m : ℤ) * (b : S →₀ ℤ) s - r) ∈ effSub Γ := by
      refine ⟨?_, ?_⟩
      · have h1 : Finsupp.single s ((m : ℤ) * (b : S →₀ ℤ) s - r)
            = (m : ℕ) • (b : S →₀ ℤ) - Finsupp.single s r := by
          rw [hbsupp, Finsupp.smul_single, ← Finsupp.single_sub, nsmul_eq_mul]
          congr 1
          rw [hbsupp, Finsupp.single_eq_same]
        rw [h1]
        exact Γ.sub_mem (Γ.nsmul_mem b.2.1 m) (mem_effSub.mp hmem).1
      · refine Finsupp.le_def.mpr fun t => ?_
        rcases eq_or_ne t s with rfl | ht
        · simp only [Finsupp.coe_zero, Pi.zero_apply, Finsupp.single_eq_same]
          linarith
        · simp only [Finsupp.coe_zero, Pi.zero_apply, Finsupp.single_eq_of_ne ht]
          exact le_refl 0
    refine ⟨m, hm0, ⟨⟨_, hdmem⟩, Subtype.ext ?_⟩⟩
    show Finsupp.single s r + Finsupp.single s ((m : ℤ) * (b : S →₀ ℤ) s - r)
      = (m : ℕ) • (b : S →₀ ℤ)
    rw [← Finsupp.single_add, hbsupp, Finsupp.smul_single, Finsupp.single_eq_same, nsmul_eq_mul]
    congr 1
    ring

/-- ★★★★**素因子の置換で誘導される引き戻しは non-dilating**(幾何側の骨組み)。

★入力は 2 つ:
* `hprim` —— 各素因子 `s` に `n·s ∈ Φ(L)` となる `n > 0` がある(`K`-`Q`-Cartier から)
* `hα` —— 引き戻しが**素因子の置換(係数 1)**で誘導される

★★2 つ目は `CartierDatum` に条項が**無い**ので、
構造に 1 条足すか幾何側で作るかが残っている(節点 `thm62-iii-nondil`)。 -/
theorem isNonDilating_effSub_of_perm {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)} (σ : S ≃ S)
    (hprim : ∀ s : S, ∃ r : ℤ, 0 < r ∧ Finsupp.single s r ∈ effSub Γ)
    (α : effSub Γ →+ effSub Γ)
    (hα : ∀ x : effSub Γ, ((α x : effSub Γ) : S →₀ ℤ)
      = Finsupp.equivMapDomain σ ((x : effSub Γ) : S →₀ ℤ)) :
    IsNonDilating α := by
  refine isNonDilating_of_sharp (isSharp_effSub Γ) α (fun h => ?_)
  have hσ : ∀ s, σ s = s := by
    intro s
    obtain ⟨r, hr, hmem⟩ := hprim s
    obtain ⟨n, hn, c, hc⟩ := h _ (isPrimaryElt_effSub_single hr hmem)
    have hval : ((α ⟨Finsupp.single s r, hmem⟩ : effSub Γ) : S →₀ ℤ)
        = Finsupp.single (σ s) r := by
      rw [hα]
      exact Finsupp.equivMapDomain_single σ s r
    by_contra hne
    have hcoe : ((α ⟨Finsupp.single s r, hmem⟩ : effSub Γ) : S →₀ ℤ) + (c : S →₀ ℤ)
        = n • (Finsupp.single s r : S →₀ ℤ) := congrArg Subtype.val hc
    have h1 := congrArg (fun f : S →₀ ℤ => f (σ s)) hcoe
    simp only [Finsupp.add_apply, Finsupp.smul_apply] at h1
    rw [hval, Finsupp.single_eq_same, Finsupp.single_eq_of_ne hne, smul_zero] at h1
    have h2 : 0 ≤ (c : S →₀ ℤ) (σ s) := by simpa using Finsupp.le_def.mp c.2.2 (σ s)
    linarith
  have hid : σ = Equiv.refl S := Equiv.ext hσ
  refine AddMonoidHom.ext fun x => Subtype.ext ?_
  rw [hα, hid, Finsupp.equivMapDomain_refl]
  rfl

/-- ★★★locator —— 幾何側の non-dilating の骨組み。 -/
def isNonDilating_effSub_of_perm.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — 素因子の置換で誘導される引き戻しは non-dilating",
    sectionId := "frdi-thm-6-2" }

/-- ★★**`K`-`Q`-Cartier から「各素因子が primary 元を持つ」が出る**。

★これが `isNonDilating_effSub_of_perm` の入力 `hprim` である。 -/
theorem exists_single_mem_effSub_of_qc {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)}
    (hq : IsQCartierSubgroup Γ) (s : S) :
    ∃ r : ℤ, 0 < r ∧ Finsupp.single s r ∈ effSub Γ := by
  obtain ⟨n, hn, hmem⟩ := hq s
  refine ⟨(n : ℤ), by exact_mod_cast hn, hmem, ?_⟩
  refine Finsupp.le_def.mpr fun t => ?_
  rcases eq_or_ne t s with rfl | ht
  · simp only [Finsupp.coe_zero, Pi.zero_apply, Finsupp.single_eq_same]
    exact_mod_cast hn.le
  · simp only [Finsupp.coe_zero, Pi.zero_apply, Finsupp.single_eq_of_ne ht]
    exact le_refl 0

/-- ★★★★★**`K`-`Q`-Cartier の下で、素因子の置換で誘導される引き戻しは non-dilating**。

★★これで `Theorem 6.2, (iii)` の non-dilating に残る入力は
**「引き戻しが素因子の置換(係数 1)で誘導される」ただ 1 つ**になった。 -/
theorem isNonDilating_effSub_of_perm_of_qc {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)}
    (hq : IsQCartierSubgroup Γ) (σ : S ≃ S) (α : effSub Γ →+ effSub Γ)
    (hα : ∀ x : effSub Γ, ((α x : effSub Γ) : S →₀ ℤ)
      = Finsupp.equivMapDomain σ ((x : effSub Γ) : S →₀ ℤ)) :
    IsNonDilating α :=
  isNonDilating_effSub_of_perm σ (exists_single_mem_effSub_of_qc hq) α hα

/-! ## ★★★`ℤ` の部分群の非負部分 —— 幾何側 non-dilating の核(2026-08-21)

★★**幾何の引き戻しが「素因子の置換(係数 1)」であることを仮定しなくても済む道**が
これである —— 各素因子の成分は `ℤ` の部分群の非負部分なので**離散**であり、
★**加法的な全単射は恒等しかない**。
(算術側でこの道が使えないのは、アルキメデス素点の成分が `ℝ≥0` そのもので
`c ↦ 2c` という加法的全単射があるからである。) -/

/-- ★`H ≤ ℤ` の非負部分。 -/
def intNonneg (H : AddSubgroup ℤ) : AddSubmonoid ℤ where
  carrier := {n : ℤ | n ∈ H ∧ 0 ≤ n}
  add_mem' ha hb := ⟨H.add_mem ha.1 hb.1, add_nonneg ha.2 hb.2⟩
  zero_mem' := ⟨H.zero_mem, le_refl 0⟩

theorem mem_intNonneg {H : AddSubgroup ℤ} {n : ℤ} :
    n ∈ intNonneg H ↔ n ∈ H ∧ 0 ≤ n := Iff.rfl

/-- ★★非負部分は生成元の**自然数**倍全体(`ℤ` の部分群は巡回だから)。 -/
theorem exists_gen_intNonneg (H : AddSubgroup ℤ) :
    ∃ d : ℤ, 0 ≤ d ∧ d ∈ H ∧ ∀ n ∈ intNonneg H, ∃ k : ℕ, n = (k : ℤ) * d := by
  obtain ⟨a, ha⟩ := Int.subgroup_cyclic H
  refine ⟨|a|, abs_nonneg a, ?_, ?_⟩
  · rw [ha]
    rcases abs_choice a with h | h
    · rw [h]; exact AddSubgroup.subset_closure rfl
    · rw [h]; exact AddSubgroup.neg_mem _ (AddSubgroup.subset_closure rfl)
  · intro n hn
    have hnH : n ∈ H := hn.1
    rw [ha, AddSubgroup.mem_closure_singleton] at hnH
    obtain ⟨m, hm⟩ := hnH
    rw [smul_eq_mul] at hm
    rcases lt_trichotomy a 0 with hneg | rfl | hpos
    · rw [abs_of_neg hneg]
      refine ⟨(-m).toNat, ?_⟩
      have hm0 : 0 ≤ -m := by nlinarith [hn.2]
      rw [Int.toNat_of_nonneg hm0, ← hm]
      ring
    · refine ⟨0, ?_⟩
      simp at hm ⊢
      omega
    · rw [abs_of_pos hpos]
      refine ⟨m.toNat, ?_⟩
      have hm0 : 0 ≤ m := by nlinarith [hn.2]
      rw [Int.toNat_of_nonneg hm0, ← hm]

/-- ★★★★**離散だから、加法的な全単射は恒等しかない**。

★★これが幾何側 non-dilating の核である ——
`k` 倍の写像が全単射なら `k = 1`。 -/
theorem eq_id_of_bijective_intNonneg {H : AddSubgroup ℤ} (f : intNonneg H →+ intNonneg H)
    (hbij : Function.Bijective f) (x : intNonneg H) : f x = x := by
  obtain ⟨d, hd0, hdH, hgen⟩ := exists_gen_intNonneg H
  have hdmem : d ∈ intNonneg H := ⟨hdH, hd0⟩
  set D : intNonneg H := ⟨d, hdmem⟩ with hD
  have hcoe : ∀ (k : ℕ) (y : intNonneg H), (((k • y : intNonneg H)) : ℤ) = (k : ℤ) * (y : ℤ) := by
    intro k y; simp
  by_cases hd : d = 0
  · obtain ⟨k, hk⟩ := hgen (x : ℤ) x.2
    have hx : x = 0 := Subtype.ext (by rw [hk, hd, mul_zero]; rfl)
    rw [hx, map_zero]
  · obtain ⟨c, hc⟩ := hgen ((f D : intNonneg H) : ℤ) (f D).2
    obtain ⟨z, hz⟩ := hbij.2 D
    obtain ⟨k, hk⟩ := hgen (z : ℤ) z.2
    have hzD : z = k • D := Subtype.ext (by rw [hk, hcoe])
    have h2 : f (k • D) = D := by rw [← hzD]; exact hz
    have h3 := congrArg (fun y : intNonneg H => (y : ℤ)) h2
    simp only [map_nsmul, hcoe, hc] at h3
    have hkc : (k : ℤ) * (c : ℤ) = 1 := by
      have h4 : ((k : ℤ) * (c : ℤ)) * d = 1 * d := by
        rw [one_mul, mul_assoc]
        exact h3
      exact mul_right_cancel₀ hd h4
    have hnat : k * c = 1 := by exact_mod_cast hkc
    have hc1 : (c : ℤ) = 1 := by
      have hcc : c = 1 := Nat.eq_one_of_mul_eq_one_left hnat
      exact_mod_cast hcc
    have hfD : f D = D := Subtype.ext (by rw [hc, hc1, one_mul])
    obtain ⟨k', hk'⟩ := hgen (x : ℤ) x.2
    have hxD : x = k' • D := Subtype.ext (by rw [hk', hcoe])
    rw [hxD, map_nsmul, hfD]

/-- ★★★locator —— 幾何側 non-dilating の核(離散性)。 -/
def eq_id_of_bijective_intNonneg.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — 素因子の成分は離散なので加法的全単射は恒等",
    sectionId := "frdi-thm-6-2" }

/-! ## ★★素因子ごとの成分(幾何側 non-dilating の部品、2026-08-21) -/

/-- ★素因子 `s` の成分がなす `ℤ` の部分群。 -/
def compSubgroup {S : Type*} (Γ : AddSubgroup (S →₀ ℤ)) (s : S) : AddSubgroup ℤ where
  carrier := {n : ℤ | Finsupp.single s n ∈ Γ}
  add_mem' ha hb := by
    show Finsupp.single s _ ∈ Γ
    rw [Finsupp.single_add]
    exact Γ.add_mem ha hb
  zero_mem' := by
    show Finsupp.single s (0 : ℤ) ∈ Γ
    rw [Finsupp.single_zero]
    exact Γ.zero_mem
  neg_mem' ha := by
    show Finsupp.single s _ ∈ Γ
    rw [Finsupp.single_neg]
    exact Γ.neg_mem ha

theorem mem_compSubgroup {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)} {s : S} {n : ℤ} :
    n ∈ compSubgroup Γ s ↔ Finsupp.single s n ∈ Γ := Iff.rfl

/-- ★★**台が 1 点の元は成分の非負部分と対応する**。 -/
theorem single_mem_effSub_iff {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)} {s : S} {n : ℤ} :
    Finsupp.single s n ∈ effSub Γ ↔ n ∈ intNonneg (compSubgroup Γ s) := by
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨h1, ?_⟩
    have := Finsupp.le_def.mp h2 s
    simpa using this
  · rintro ⟨h1, h2⟩
    refine ⟨h1, Finsupp.le_def.mpr fun t => ?_⟩
    rcases eq_or_ne t s with rfl | ht
    · simpa using h2
    · simp [Finsupp.single_eq_of_ne ht]

/-- ★★**non-dilating の仮定は「台が動かない」と読める**。 -/
theorem alpha_single_supported {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)}
    (α : effSub Γ →+ effSub Γ)
    (h : ∀ a : effSub Γ, IsPrimaryElt a → MPrec (α a) a)
    {s : S} {n : ℤ} (hn : 0 < n) (hmem : Finsupp.single s n ∈ effSub Γ) {t : S} (ht : t ≠ s) :
    ((α ⟨Finsupp.single s n, hmem⟩ : effSub Γ) : S →₀ ℤ) t = 0 := by
  obtain ⟨m, hm, c, hc⟩ := h _ (isPrimaryElt_effSub_single hn hmem)
  have hcoe : ((α ⟨Finsupp.single s n, hmem⟩ : effSub Γ) : S →₀ ℤ) + (c : S →₀ ℤ)
      = m • (Finsupp.single s n : S →₀ ℤ) := congrArg Subtype.val hc
  have h1 := congrArg (fun f : S →₀ ℤ => f t) hcoe
  have hz : (Finsupp.single s n : S →₀ ℤ) t = 0 := Finsupp.single_eq_of_ne ht
  simp only [Finsupp.add_apply, Finsupp.smul_apply] at h1
  rw [hz, smul_zero] at h1
  have h2 : 0 ≤ ((α ⟨Finsupp.single s n, hmem⟩ : effSub Γ) : S →₀ ℤ) t := by
    simpa using Finsupp.le_def.mp (α ⟨Finsupp.single s n, hmem⟩).2.2 t
  have h3 : 0 ≤ (c : S →₀ ℤ) t := by simpa using Finsupp.le_def.mp c.2.2 t
  linarith

/-- ★★台が 1 点なら像も台が 1 点。 -/
theorem alpha_single_eq {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)} (α : effSub Γ →+ effSub Γ)
    (h : ∀ a : effSub Γ, IsPrimaryElt a → MPrec (α a) a)
    {s : S} {n : ℤ} (hn : 0 < n) (hmem : Finsupp.single s n ∈ effSub Γ) :
    ((α ⟨Finsupp.single s n, hmem⟩ : effSub Γ) : S →₀ ℤ)
      = Finsupp.single s (((α ⟨Finsupp.single s n, hmem⟩ : effSub Γ) : S →₀ ℤ) s) := by
  refine Finsupp.ext fun t => ?_
  rcases eq_or_ne t s with rfl | ht
  · rw [Finsupp.single_eq_same]
  · rw [alpha_single_supported α h hn hmem ht, Finsupp.single_eq_of_ne ht]

/-! ## ★★★★★幾何側 non-dilating —— 幾何の入力なしで閉じる(2026-08-21)

★★★**「引き戻しが素因子の置換(係数 1)で誘導される」を仮定しなくてよい**ことが分かった。
必要なのは **`K`-`Q`-Cartier** と **引き戻しが全単射**(自己射は `𝒟` で同型なので
`pull_id` / `pull_comp` から出る)だけである。

★筋:
1. non-dilating の仮定 ⟹ 台が 1 点の元の像も台が 1 点(`alpha_single_supported`)
2. 各素因子の成分は `ℤ` の部分群の非負部分 ⟹ **離散**
   ⟹ そこでの加法的全単射は恒等(`eq_id_of_bijective_intNonneg`)
3. `K`-`Q`-Cartier ⟹ どの元も `N` 倍すれば台が 1 点の元の**有限和**(`exists_common_mult`)
   ⟹ 引き戻しは**座標ごと**(`alpha_coord`)
4. 2 と 3 で `α = id`

★★算術側(`Example 6.3`)でこの道が使えないのは、アルキメデス素点の成分が
`ℝ≥0` そのもので `c ↦ 2c` という加法的全単射があるからである。 -/

section GeomND
variable {S : Type*} [DecidableEq S] {Γ : AddSubgroup (S →₀ ℤ)}

/-- ★★`K`-`Q`-Cartier ⟹ どの元も `N` 倍すれば各座標が単独で `Φ(L)` に入る。 -/
theorem exists_common_mult (hq : IsQCartierSubgroup Γ) (x : effSub Γ) :
    ∃ N : ℕ, 0 < N ∧ ∀ t : S, Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)) ∈ effSub Γ := by
  classical
  set n : S → ℕ := fun t => (hq t).choose with hn
  have hnpos : ∀ t, 0 < n t := fun t => (hq t).choose_spec.1
  have hnmem : ∀ t, Finsupp.single t ((n t : ℤ)) ∈ Γ := fun t => (hq t).choose_spec.2
  set N : ℕ := ∏ t ∈ ((x : S →₀ ℤ)).support, n t with hN
  have hNpos : 0 < N := Finset.prod_pos (fun t _ => hnpos t)
  refine ⟨N, hNpos, fun t => ?_⟩
  by_cases ht : t ∈ ((x : S →₀ ℤ)).support
  · have hdvd : (n t : ℤ) ∣ (N : ℤ) :=
      Int.natCast_dvd_natCast.mpr (Finset.dvd_prod_of_mem n ht)
    obtain ⟨q, hqq⟩ := hdvd
    refine ⟨?_, ?_⟩
    · have heq : Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t))
          = (q * ((x : S →₀ ℤ) t)) • Finsupp.single t ((n t : ℤ)) := by
        rw [Finsupp.smul_single, hqq]
        congr 1
        ring
      rw [heq]
      exact Γ.zsmul_mem (hnmem t) _
    · refine Finsupp.le_def.mpr fun u => ?_
      rcases eq_or_ne u t with rfl | hu
      · simp only [Finsupp.coe_zero, Pi.zero_apply, Finsupp.single_eq_same]
        have hx0 : 0 ≤ ((x : S →₀ ℤ)) u := by simpa using Finsupp.le_def.mp x.2.2 u
        positivity
      · have hz : (Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)) : S →₀ ℤ) u = 0 :=
          Finsupp.single_eq_of_ne hu
        simp only [Finsupp.coe_zero, Pi.zero_apply, hz]
        exact le_refl 0
  · have hx0 : ((x : S →₀ ℤ)) t = 0 := by simpa using Finsupp.notMem_support_iff.mp ht
    rw [hx0, mul_zero, Finsupp.single_zero]
    exact (effSub Γ).zero_mem

/-- ★`N` 倍は座標ごとの有限和。 -/
theorem nsmul_eq_sum_single (x : effSub Γ) (N : ℕ)
    (hmem : ∀ t : S, Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)) ∈ effSub Γ) :
    (N • x : effSub Γ)
      = ∑ t ∈ ((x : S →₀ ℤ)).support,
          (⟨Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)), hmem t⟩ : effSub Γ) := by
  classical
  refine Subtype.ext ?_
  rw [AddSubmonoid.coe_finsetSum]
  refine Finsupp.ext fun u => ?_
  rw [Finsupp.finsetSum_apply]
  by_cases hu : u ∈ ((x : S →₀ ℤ)).support
  · rw [Finset.sum_eq_single u]
    · simp
    · intro t _ htu
      exact Finsupp.single_eq_of_ne (Ne.symm htu)
    · intro hcon
      exact absurd hu hcon
  · have hx0 : ((x : S →₀ ℤ)) u = 0 := by simpa using Finsupp.notMem_support_iff.mp hu
    have hz : ∀ t ∈ ((x : S →₀ ℤ)).support,
        (Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)) : S →₀ ℤ) u = 0 := by
      intro t ht
      refine Finsupp.single_eq_of_ne ?_
      intro hc
      exact hu (hc ▸ ht)
    rw [Finset.sum_eq_zero hz]
    simp [hx0]

/-- ★★★**引き戻しは座標ごと**(`N` 倍して見れば)。 -/
theorem alpha_coord (α : effSub Γ →+ effSub Γ)
    (h : ∀ a : effSub Γ, IsPrimaryElt a → MPrec (α a) a)
    (x : effSub Γ) {N : ℕ} (hN : 0 < N)
    (hmem : ∀ t : S, Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)) ∈ effSub Γ) (u : S) :
    (N : ℤ) * (((α x : effSub Γ) : S →₀ ℤ) u)
      = ((α ⟨Finsupp.single u ((N : ℤ) * ((x : S →₀ ℤ) u)), hmem u⟩ : effSub Γ) : S →₀ ℤ) u := by
  classical
  have hxnn : ∀ t, 0 ≤ ((x : S →₀ ℤ)) t := fun t => by simpa using Finsupp.le_def.mp x.2.2 t
  have hstep : α (N • x)
      = ∑ t ∈ ((x : S →₀ ℤ)).support,
          α (⟨Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)), hmem t⟩ : effSub Γ) := by
    rw [nsmul_eq_sum_single x N hmem, map_sum]
  have hL : (((α (N • x) : effSub Γ)) : S →₀ ℤ) u
      = (N : ℤ) * (((α x : effSub Γ) : S →₀ ℤ) u) := by
    rw [map_nsmul]
    simp
  have hR := congrArg (fun y : effSub Γ => ((y : S →₀ ℤ)) u) hstep
  rw [hL] at hR
  rw [hR, AddSubmonoid.coe_finsetSum, Finsupp.finsetSum_apply]
  by_cases hu : u ∈ ((x : S →₀ ℤ)).support
  · refine Finset.sum_eq_single u ?_ ?_
    · intro t ht htu
      have hxt : 0 < ((x : S →₀ ℤ)) t :=
        lt_of_le_of_ne (hxnn t) (Ne.symm (Finsupp.mem_support_iff.mp ht))
      have hpos : 0 < (N : ℤ) * ((x : S →₀ ℤ)) t := by positivity
      exact alpha_single_supported α h hpos (hmem t) (Ne.symm htu)
    · intro hcon; exact absurd hu hcon
  · have hx0 : ((x : S →₀ ℤ)) u = 0 := by simpa using Finsupp.notMem_support_iff.mp hu
    have hz : ∀ t ∈ ((x : S →₀ ℤ)).support,
        ((α (⟨Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)), hmem t⟩ : effSub Γ)
          : effSub Γ) : S →₀ ℤ) u = 0 := by
      intro t ht
      have hxt : 0 < ((x : S →₀ ℤ)) t :=
        lt_of_le_of_ne (hxnn t) (Ne.symm (Finsupp.mem_support_iff.mp ht))
      have hpos : 0 < (N : ℤ) * ((x : S →₀ ℤ)) t := by positivity
      refine alpha_single_supported α h hpos (hmem t) ?_
      intro hc
      exact hu (hc ▸ ht)
    have hzero : (⟨Finsupp.single u ((N : ℤ) * ((x : S →₀ ℤ) u)), hmem u⟩ : effSub Γ) = 0 := by
      refine Subtype.ext ?_
      show Finsupp.single u ((N : ℤ) * ((x : S →₀ ℤ) u)) = (0 : S →₀ ℤ)
      rw [Finsupp.single_eq_zero, hx0, mul_zero]
    rw [Finset.sum_eq_zero hz, hzero, map_zero]
    rfl

theorem alpha_single_val_mem (α : effSub Γ →+ effSub Γ)
    (h : ∀ a : effSub Γ, IsPrimaryElt a → MPrec (α a) a)
    {u : S} {n : ℤ} (hn : n ∈ intNonneg (compSubgroup Γ u)) :
    ((α ⟨Finsupp.single u n, single_mem_effSub_iff.mpr hn⟩ : effSub Γ) : S →₀ ℤ) u
      ∈ intNonneg (compSubgroup Γ u) := by
  set y : effSub Γ := α ⟨Finsupp.single u n, single_mem_effSub_iff.mpr hn⟩ with hy
  rcases eq_or_lt_of_le hn.2 with hn0 | hpos
  · have hz : (⟨Finsupp.single u n, single_mem_effSub_iff.mpr hn⟩ : effSub Γ) = 0 := by
      refine Subtype.ext ?_
      show Finsupp.single u n = (0 : S →₀ ℤ)
      rw [Finsupp.single_eq_zero, ← hn0]
    have hy0 : y = 0 := by rw [hy, hz, map_zero]
    rw [hy0]
    exact (intNonneg (compSubgroup Γ u)).zero_mem
  · have heq : ((y : effSub Γ) : S →₀ ℤ) = Finsupp.single u (((y : effSub Γ) : S →₀ ℤ) u) :=
      alpha_single_eq α h hpos (single_mem_effSub_iff.mpr hn)
    refine single_mem_effSub_iff.mp ?_
    rw [← heq]
    exact y.2

/-- ★★成分ごとの写像。 -/
noncomputable def compMap (α : effSub Γ →+ effSub Γ)
    (h : ∀ a : effSub Γ, IsPrimaryElt a → MPrec (α a) a) (u : S) :
    intNonneg (compSubgroup Γ u) →+ intNonneg (compSubgroup Γ u) where
  toFun n := ⟨((α ⟨Finsupp.single u (n : ℤ), single_mem_effSub_iff.mpr n.2⟩ : effSub Γ)
      : S →₀ ℤ) u, alpha_single_val_mem α h n.2⟩
  map_zero' := by
    refine Subtype.ext ?_
    have hz : (⟨Finsupp.single u ((0 : intNonneg (compSubgroup Γ u)) : ℤ),
        single_mem_effSub_iff.mpr (0 : intNonneg (compSubgroup Γ u)).2⟩ : effSub Γ) = 0 := by
      refine Subtype.ext ?_
      show Finsupp.single u ((0 : intNonneg (compSubgroup Γ u)) : ℤ) = (0 : S →₀ ℤ)
      rw [Finsupp.single_eq_zero]
      rfl
    show ((α _ : effSub Γ) : S →₀ ℤ) u = 0
    rw [hz, map_zero]
    rfl
  map_add' m n := by
    refine Subtype.ext ?_
    have hadd : (⟨Finsupp.single u ((m + n : intNonneg (compSubgroup Γ u)) : ℤ),
        single_mem_effSub_iff.mpr (m + n).2⟩ : effSub Γ)
        = ⟨Finsupp.single u (m : ℤ), single_mem_effSub_iff.mpr m.2⟩
          + ⟨Finsupp.single u (n : ℤ), single_mem_effSub_iff.mpr n.2⟩ := by
      refine Subtype.ext ?_
      show Finsupp.single u ((m : ℤ) + (n : ℤ)) = _
      rw [Finsupp.single_add]
      rfl
    show ((α _ : effSub Γ) : S →₀ ℤ) u = _
    rw [hadd, map_add]
    rfl

@[simp] theorem compMap_val (α : effSub Γ →+ effSub Γ)
    (h : ∀ a : effSub Γ, IsPrimaryElt a → MPrec (α a) a) (u : S)
    (n : intNonneg (compSubgroup Γ u)) :
    ((compMap α h u n : intNonneg (compSubgroup Γ u)) : ℤ)
      = ((α ⟨Finsupp.single u (n : ℤ), single_mem_effSub_iff.mpr n.2⟩ : effSub Γ) : S →₀ ℤ) u :=
  rfl

theorem alpha_single_coe (α : effSub Γ →+ effSub Γ)
    (h : ∀ a : effSub Γ, IsPrimaryElt a → MPrec (α a) a)
    {u : S} (n : intNonneg (compSubgroup Γ u)) :
    ((α ⟨Finsupp.single u (n : ℤ), single_mem_effSub_iff.mpr n.2⟩ : effSub Γ) : S →₀ ℤ)
      = Finsupp.single u
        (((α ⟨Finsupp.single u (n : ℤ), single_mem_effSub_iff.mpr n.2⟩ : effSub Γ)
          : S →₀ ℤ) u) := by
  rcases eq_or_lt_of_le n.2.2 with hn0 | hpos
  · have hz : (⟨Finsupp.single u (n : ℤ), single_mem_effSub_iff.mpr n.2⟩ : effSub Γ) = 0 := by
      refine Subtype.ext ?_
      show Finsupp.single u (n : ℤ) = (0 : S →₀ ℤ)
      rw [Finsupp.single_eq_zero, ← hn0]
    rw [hz, map_zero]
    show (0 : S →₀ ℤ) = Finsupp.single u ((0 : S →₀ ℤ) u)
    simp
  · exact alpha_single_eq α h hpos (single_mem_effSub_iff.mpr n.2)

theorem compMap_injective (α : effSub Γ →+ effSub Γ)
    (h : ∀ a : effSub Γ, IsPrimaryElt a → MPrec (α a) a)
    (hinj : Function.Injective α) (u : S) : Function.Injective (compMap α h u) := by
  intro m n hmn
  have h1 : ((α ⟨Finsupp.single u (m : ℤ), single_mem_effSub_iff.mpr m.2⟩ : effSub Γ)
        : S →₀ ℤ) u
      = ((α ⟨Finsupp.single u (n : ℤ), single_mem_effSub_iff.mpr n.2⟩ : effSub Γ) : S →₀ ℤ) u :=
    congrArg Subtype.val hmn
  have h2 : (α ⟨Finsupp.single u (m : ℤ), single_mem_effSub_iff.mpr m.2⟩ : effSub Γ)
      = α ⟨Finsupp.single u (n : ℤ), single_mem_effSub_iff.mpr n.2⟩ := by
    refine Subtype.ext ?_
    rw [alpha_single_coe α h m, alpha_single_coe α h n, h1]
  have h3 := hinj h2
  have h4 : Finsupp.single u (m : ℤ) = Finsupp.single u (n : ℤ) := congrArg Subtype.val h3
  have h5 := congrArg (fun f : S →₀ ℤ => f u) h4
  simp only [Finsupp.single_eq_same] at h5
  exact Subtype.ext h5

theorem compMap_surjective (hq : IsQCartierSubgroup Γ) (α : effSub Γ →+ effSub Γ)
    (h : ∀ a : effSub Γ, IsPrimaryElt a → MPrec (α a) a)
    (hbij : Function.Bijective α) (u : S) : Function.Surjective (compMap α h u) := by
  intro n
  obtain ⟨x, hx⟩ := hbij.2 ⟨Finsupp.single u (n : ℤ), single_mem_effSub_iff.mpr n.2⟩
  obtain ⟨N, hN, hmem⟩ := exists_common_mult hq x
  have hxnn : ∀ t, 0 ≤ ((x : S →₀ ℤ)) t := fun t => by simpa using Finsupp.le_def.mp x.2.2 t
  have hzero : ∀ t : S, t ≠ u → ((x : S →₀ ℤ)) t = 0 := by
    intro t ht
    by_contra hne
    have hpos : 0 < ((x : S →₀ ℤ)) t := lt_of_le_of_ne (hxnn t) (Ne.symm hne)
    have hNpos : 0 < (N : ℤ) * ((x : S →₀ ℤ)) t := by positivity
    have hc := alpha_coord α h x hN hmem t
    have hax : ((α x : effSub Γ) : S →₀ ℤ) t = 0 := by
      rw [hx]
      show (Finsupp.single u (n : ℤ) : S →₀ ℤ) t = 0
      exact Finsupp.single_eq_of_ne ht
    rw [hax, mul_zero] at hc
    have hsupp := alpha_single_eq α h hNpos (hmem t)
    have hz : (α ⟨Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)), hmem t⟩ : effSub Γ) = 0 := by
      refine Subtype.ext ?_
      rw [hsupp, ← hc, Finsupp.single_zero]
      rfl
    have hne0 : (⟨Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)), hmem t⟩ : effSub Γ) ≠ 0 := by
      intro hc2
      have hcc : Finsupp.single t ((N : ℤ) * ((x : S →₀ ℤ) t)) = (0 : S →₀ ℤ) :=
        congrArg Subtype.val hc2
      rw [Finsupp.single_eq_zero] at hcc
      exact hNpos.ne' hcc
    exact hne0 (hbij.1 (by rw [hz, map_zero]))
  have hxsingle : ((x : S →₀ ℤ)) = Finsupp.single u (((x : S →₀ ℤ)) u) := by
    refine Finsupp.ext fun t => ?_
    rcases eq_or_ne t u with rfl | ht
    · rw [Finsupp.single_eq_same]
    · rw [hzero t ht, Finsupp.single_eq_of_ne ht]
  have hmemu : ((x : S →₀ ℤ)) u ∈ intNonneg (compSubgroup Γ u) :=
    single_mem_effSub_iff.mp (by rw [← hxsingle]; exact x.2)
  refine ⟨⟨((x : S →₀ ℤ)) u, hmemu⟩, ?_⟩
  refine Subtype.ext ?_
  rw [compMap_val]
  have hxeq : (⟨Finsupp.single u (((x : S →₀ ℤ)) u), single_mem_effSub_iff.mpr hmemu⟩
      : effSub Γ) = x := Subtype.ext hxsingle.symm
  rw [hxeq, hx]
  show (Finsupp.single u (n : ℤ) : S →₀ ℤ) u = (n : ℤ)
  rw [Finsupp.single_eq_same]

/-- ★★★★★★**幾何側 non-dilating** ——
`K`-`Q`-Cartier の下で、**全単射な引き戻しは non-dilating**。

★★幾何の入力(素因子の置換であること)は**要らない** ——
各素因子の成分が離散であることだけで出る。 -/
theorem isNonDilating_effSub_of_bijective (hq : IsQCartierSubgroup Γ)
    (α : effSub Γ →+ effSub Γ) (hbij : Function.Bijective α) : IsNonDilating α := by
  refine isNonDilating_of_sharp (isSharp_effSub Γ) α (fun h => ?_)
  refine AddMonoidHom.ext fun x => Subtype.ext (Finsupp.ext fun u => ?_)
  obtain ⟨N, hN, hmem⟩ := exists_common_mult hq x
  have hc := alpha_coord α h x hN hmem u
  have hmemu : ((N : ℤ) * ((x : S →₀ ℤ) u)) ∈ intNonneg (compSubgroup Γ u) :=
    single_mem_effSub_iff.mp (hmem u)
  have hid := eq_id_of_bijective_intNonneg (compMap α h u)
    ⟨compMap_injective α h hbij.1 u, compMap_surjective hq α h hbij u⟩ ⟨_, hmemu⟩
  have hval := congrArg Subtype.val hid
  rw [compMap_val] at hval
  rw [hval] at hc
  have hNne : (N : ℤ) ≠ 0 := by exact_mod_cast hN.ne'
  exact mul_left_cancel₀ hNne hc

/-- ★★★★★locator —— 幾何側 non-dilating(幾何の入力なし)。 -/
def isNonDilating_effSub_of_bijective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — K-Q-Cartier と全単射だけで non-dilating",
    sectionId := "frdi-thm-6-2" }

end GeomND

/-! ## ★★★`Γ` の元を「正・負」に割る —— Cartier の壁は壁でなかった(2026-08-21)

★`IsStrictlyRational` は `y = toGp a − toGp b` を
「`p` で正・`p` を含まない」形に取り替えることを要求する。
★★**Cartier 因子の正部分は Cartier とは限らない**ので素朴には割れないが、
★★★**`K`-`Q`-Cartier があれば割れる** —— 負の座標を
`n_t·|y_t|`(`n_t·[t] ∈ Γ`)で埋めればよい。`p` の座標は `y_p > 0` なので埋めなくてよい。 -/

section Split
variable {S : Type*} [DecidableEq S] {Γ : AddSubgroup (S →₀ ℤ)}

/-- ★★★★**`K`-`Q`-Cartier の下では、`Γ` の元は「`p` で正・`p` を含まない」差に書ける**。 -/
theorem exists_split_of_qc (hq : IsQCartierSubgroup Γ) {y : S →₀ ℤ} (hy : y ∈ Γ)
    {p : S} (hp : 0 < y p) :
    ∃ a b : S →₀ ℤ, a ∈ effSub Γ ∧ b ∈ effSub Γ ∧ a = b + y ∧ 0 < a p ∧ b p = 0 := by
  classical
  set n : S → ℕ := fun t => (hq t).choose with hn
  have hnpos : ∀ t, 0 < n t := fun t => (hq t).choose_spec.1
  have hnmem : ∀ t, Finsupp.single t ((n t : ℤ)) ∈ Γ := fun t => (hq t).choose_spec.2
  set b : S →₀ ℤ :=
    ∑ t ∈ y.support.erase p, Finsupp.single t ((n t : ℤ) * |y t|) with hb
  have hbmem : b ∈ Γ := by
    refine AddSubgroup.sum_mem _ fun t _ => ?_
    have hsm : Finsupp.single t ((n t : ℤ) * |y t|) = |y t| • Finsupp.single t ((n t : ℤ)) := by
      rw [Finsupp.smul_single]
      congr 1
      ring
    rw [hsm]
    exact Γ.zsmul_mem (hnmem t) _
  have hbapp : ∀ u : S, b u = if u ∈ y.support.erase p then (n u : ℤ) * |y u| else 0 := by
    intro u
    rw [hb, Finsupp.finsetSum_apply]
    by_cases hu : u ∈ y.support.erase p
    · rw [if_pos hu, Finset.sum_eq_single u]
      · rw [Finsupp.single_eq_same]
      · intro t _ htu
        exact Finsupp.single_eq_of_ne (Ne.symm htu)
      · intro hcon; exact absurd hu hcon
    · rw [if_neg hu, Finset.sum_eq_zero]
      intro t ht
      refine Finsupp.single_eq_of_ne ?_
      intro hc
      exact hu (hc ▸ ht)
  have hbnn : ∀ u, 0 ≤ b u := by
    intro u
    rw [hbapp u]
    by_cases hu : u ∈ y.support.erase p
    · rw [if_pos hu]
      have hn1 : (0 : ℤ) ≤ (n u : ℤ) := Int.natCast_nonneg _
      positivity
    · rw [if_neg hu]
  have hbp : b p = 0 := by
    rw [hbapp p, if_neg (by simp)]
  refine ⟨b + y, b, ⟨Γ.add_mem hbmem hy, ?_⟩, ⟨hbmem, Finsupp.le_def.mpr fun u => by
      simpa using hbnn u⟩, rfl, ?_, hbp⟩
  · refine Finsupp.le_def.mpr fun u => ?_
    simp only [Finsupp.coe_zero, Pi.zero_apply, Finsupp.add_apply]
    by_cases hu : u ∈ y.support.erase p
    · rw [hbapp u, if_pos hu]
      have h1 : (1 : ℤ) ≤ (n u : ℤ) := by exact_mod_cast hnpos u
      nlinarith [abs_nonneg (y u), neg_abs_le (y u), le_abs_self (y u)]
    · rw [hbapp u, if_neg hu]
      rcases eq_or_ne u p with rfl | hup
      · linarith
      · have hy0 : y u = 0 := by
          by_contra hc
          exact hu (Finset.mem_erase.mpr ⟨hup, Finsupp.mem_support_iff.mpr hc⟩)
        rw [hy0]
        exact le_refl 0
  · simp only [Finsupp.add_apply, hbp]
    linarith

/-- ★★★locator —— `Theorem 6.2, (iii)` の rationally standard に要る分割。 -/
def exists_split_of_qc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — K-Q-Cartier なら Γ の元は p で正・p を含まない差に割れる",
    sectionId := "frdi-thm-6-2" }

end Split

end ABC3.Found.FrdI
