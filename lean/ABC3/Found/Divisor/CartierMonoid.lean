/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.FreeDivisorial
import ABC3.Found.FrdI.MonoidTransport

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

end ABC3.Found.FrdI
