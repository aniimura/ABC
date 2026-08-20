/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.CartierPrime

/-!
# 実係数の因子単系 `Φ = Γ ∩ ℝ≥0[S]` —— `Prime(Φ) ≃ S` まで

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> Prime(Φ(L))

## ★★なぜ実係数版が要るか

[FrdI] `Example 6.1`(幾何)の `Φ(L)` は `ℤ≥0[D_L]` の部分単系であり、
`Found/Divisor/Cartier*.lean` がそれを扱う。
★ところが `Example 6.3`(算術)の `Φ(L) = ⊕_v ord(O_v^▷)` は
**アルキメデス成分が `ℝ≥0`** なので、係数が `ℤ` では収まらない。

★★そこで本ファイルは `Cartier*.lean` の議論を **`S →₀ ℝ` の上で書き直す**。
`Example 6.3` の `Φ(L)` は `Γ = arithDivGroup L`(`ArithDivisor.lean`)を取った場合である。

## ★離散性が消える所

幾何版の `mprec_effSub_iff`(`⪯` は台の包含)の `⟸` は
「`b s ≠ 0` かつ `b s ≥ 0` なら `b s ≥ 1`」という **`ℤ` の離散性**を使い、
`n = 1 + Σ a s` で足りた。★実係数ではこれが使えないので、
`n = 1 + Σ_{s ∈ supp a} ⌈a s / b s⌉` と**台の上の有限和**で取る
(`mprec_effR_iff`)。★台が有限であることだけが効いている。

## ★段取り(幾何版と 1 対 1)

| 段 | 幾何版 | 本ファイル |
|---|---|---|
| `Φ` の定義 | `effSub` | `effR` |
| divisorial | `isDivisorial_effSub` | `isDivisorial_effR` |
| `≤` の翻訳 | `mle_effSub_iff` | `mle_effR_iff` |
| `⪯` の翻訳 | `mprec_effSub_iff` | `mprec_effR_iff` |
| 生成元 | `qcGen` | `genR` |
| `Prime ≃ S` | `effSubPrimeEquiv` | `effRPrimeEquiv` |
| 台 | `exists_effSub_support_eq` | `exists_effR_support_eq` |
-/

namespace ABC3.Found.Divisor

open Finsupp ABC3.Found.FrdI

variable {S : Type*} {Γ : AddSubgroup (S →₀ ℝ)}

/-! ## ★1. `Φ = Γ ∩ ℝ≥0[S]` -/

/-- ★★**`Γ` の有効な元のなす部分単系**(実係数版)。 -/
def effR (Γ : AddSubgroup (S →₀ ℝ)) : AddSubmonoid (S →₀ ℝ) where
  carrier := {x | x ∈ Γ ∧ 0 ≤ x}
  add_mem' ha hb := ⟨Γ.add_mem ha.1 hb.1, add_nonneg ha.2 hb.2⟩
  zero_mem' := ⟨Γ.zero_mem, le_refl 0⟩

theorem mem_effR {x : S →₀ ℝ} : x ∈ effR Γ ↔ x ∈ Γ ∧ ∀ s, 0 ≤ x s := by
  simp [effR, AddSubmonoid.mem_mk, Finsupp.le_def]

theorem effR_nonneg (a : effR Γ) (s : S) : 0 ≤ (a : S →₀ ℝ) s :=
  (mem_effR.mp a.2).2 s

theorem effR_eq_zero_iff (a : effR Γ) : a = 0 ↔ (a : S →₀ ℝ) = 0 :=
  ⟨fun h => by rw [h]; rfl, fun h => Subtype.ext h⟩

/-! ## ★2. `Φ` は divisorial -/

/-- ★`Γ ∩ ℝ≥0[S]` は sharp。 -/
theorem isSharp_effR (Γ : AddSubgroup (S →₀ ℝ)) : IsSharp (effR Γ) := by
  intro a ha
  obtain ⟨u, hu⟩ := ha
  have h0 : (a : S →₀ ℝ) + ((u.neg : effR Γ) : S →₀ ℝ) = 0 := by
    have h := u.val_neg
    rw [hu] at h
    exact congrArg (fun t : effR Γ => (t : S →₀ ℝ)) h
  have hane := (mem_effR.mp a.2).2
  have hbne := (mem_effR.mp (u.neg : effR Γ).2).2
  refine Subtype.ext ?_
  ext s
  have h1 := congrArg (fun t : S →₀ ℝ => t s) h0
  simp only [Finsupp.add_apply, Finsupp.coe_zero, Pi.zero_apply] at h1
  have h2 := hane s
  have h3 := hbne s
  show (a : S →₀ ℝ) s = 0
  linarith

/-- ★`Γ ∩ ℝ≥0[S]` は saturated。 -/
theorem isSaturated_effR (Γ : AddSubgroup (S →₀ ℝ)) : IsSaturatedMonoid (effR Γ) := by
  refine isSaturatedMonoid_of_cancel_of_nsmul_le _ ?_
  rintro x y n hn ⟨w, hw⟩
  have hw' : n • (y : S →₀ ℝ) + (w : S →₀ ℝ) = n • (x : S →₀ ℝ) :=
    congrArg (fun t : effR Γ => (t : S →₀ ℝ)) hw
  have hmem : (x : S →₀ ℝ) - (y : S →₀ ℝ) ∈ effR Γ := by
    refine mem_effR.mpr ⟨Γ.sub_mem (mem_effR.mp x.2).1 (mem_effR.mp y.2).1, fun s => ?_⟩
    have h1 := congrArg (fun t : S →₀ ℝ => t s) hw'
    simp only [Finsupp.add_apply, Finsupp.smul_apply, nsmul_eq_mul] at h1
    have h2 := (mem_effR.mp w.2).2 s
    have hn' : (0:ℝ) < (n : ℝ) := by exact_mod_cast hn
    simp only [Finsupp.sub_apply]
    nlinarith
  exact ⟨⟨_, hmem⟩, Subtype.ext (by simp)⟩

/-- ★★★**`Φ := Γ ∩ ℝ≥0[S]` は divisorial**。 -/
theorem isDivisorial_effR (Γ : AddSubgroup (S →₀ ℝ)) : IsDivisorial (effR Γ) := by
  have hsharp := isSharp_effR Γ
  exact ⟨⟨isIntegralMonoid_of_isCancelAdd _, isSaturated_effR Γ,
    isOfCharacteristicType_of_isSharp _ hsharp⟩, hsharp⟩

/-! ## ★3. `≤` と `⪯` の係数による翻訳 -/

/-- ★`Φ` の `≤` は**係数ごとの `≤`**。 -/
theorem mle_effR_iff (a b : effR Γ) :
    MLe a b ↔ ∀ s, (a : S →₀ ℝ) s ≤ (b : S →₀ ℝ) s := by
  constructor
  · rintro ⟨c, rfl⟩ s
    have := effR_nonneg c s
    show (a : S →₀ ℝ) s ≤ ((a : S →₀ ℝ) + (c : S →₀ ℝ)) s
    simp only [Finsupp.add_apply]
    linarith
  · intro h
    have hmem : (b : S →₀ ℝ) - (a : S →₀ ℝ) ∈ effR Γ := by
      refine mem_effR.mpr ⟨Γ.sub_mem (mem_effR.mp b.2).1 (mem_effR.mp a.2).1, fun s => ?_⟩
      have := h s
      simp only [Finsupp.sub_apply]
      linarith
    exact ⟨⟨_, hmem⟩, Subtype.ext (by simp)⟩

/-- ★★**`Φ` の `⪯` は台の包含**(実係数版)。

★★`⟸` の `n` は幾何版の `1 + Σ a s` では**足りない**(実係数には離散性が無い)。
台の上の `⌈a s / b s⌉` の和を取る —— **台が有限であること**だけが効く。 -/
theorem mprec_effR_iff (a b : effR Γ) :
    MPrec a b ↔ (a : S →₀ ℝ).support ⊆ (b : S →₀ ℝ).support := by
  classical
  constructor
  · rintro ⟨n, hn, hle⟩ s hs
    rw [mle_effR_iff] at hle
    have h1 := hle s
    have h2 : ((n • b : effR Γ) : S →₀ ℝ) s = (n : ℝ) * (b : S →₀ ℝ) s := by
      show ((n • (b : S →₀ ℝ)) : S →₀ ℝ) s = _
      simp
    rw [h2] at h1
    have h3 : (a : S →₀ ℝ) s ≠ 0 := Finsupp.mem_support_iff.mp hs
    have h4 := effR_nonneg a s
    refine Finsupp.mem_support_iff.mpr (fun hz => ?_)
    rw [hz, mul_zero] at h1
    exact h3 (le_antisymm h1 h4)
  · intro hsub
    set T := (a : S →₀ ℝ).support with hT
    set n : ℕ := 1 + T.sum (fun s => ⌈(a : S →₀ ℝ) s / (b : S →₀ ℝ) s⌉₊) with hn
    have hnpos : 0 < n := by omega
    refine ⟨n, hnpos, ?_⟩
    rw [mle_effR_iff]
    intro s
    have h2 : ((n • b : effR Γ) : S →₀ ℝ) s = (n : ℝ) * (b : S →₀ ℝ) s := by
      show ((n • (b : S →₀ ℝ)) : S →₀ ℝ) s = _
      simp
    rw [h2]
    by_cases hs : s ∈ T
    · have hbs : (b : S →₀ ℝ) s ≠ 0 := Finsupp.mem_support_iff.mp (hsub hs)
      have hbpos : 0 < (b : S →₀ ℝ) s := lt_of_le_of_ne (effR_nonneg b s) (Ne.symm hbs)
      have hle : ⌈(a : S →₀ ℝ) s / (b : S →₀ ℝ) s⌉₊
          ≤ T.sum (fun t => ⌈(a : S →₀ ℝ) t / (b : S →₀ ℝ) t⌉₊) :=
        Finset.single_le_sum
          (f := fun t => ⌈(a : S →₀ ℝ) t / (b : S →₀ ℝ) t⌉₊) (fun _ _ => Nat.zero_le _) hs
      have hcast : (a : S →₀ ℝ) s / (b : S →₀ ℝ) s ≤ (n : ℝ) := by
        refine le_trans (Nat.le_ceil _) ?_
        have : (⌈(a : S →₀ ℝ) s / (b : S →₀ ℝ) s⌉₊ : ℝ) ≤ (n : ℝ) := by
          exact_mod_cast Nat.le_of_lt_succ (by omega)
        exact this
      calc (a : S →₀ ℝ) s = ((a : S →₀ ℝ) s / (b : S →₀ ℝ) s) * (b : S →₀ ℝ) s := by
            field_simp
        _ ≤ (n : ℝ) * (b : S →₀ ℝ) s := mul_le_mul_of_nonneg_right hcast hbpos.le
    · have h0 : (a : S →₀ ℝ) s = 0 := Finsupp.notMem_support_iff.mp hs
      have hbnn := effR_nonneg b s
      rw [h0]
      positivity

/-! ## ★4. 各素点の生成元 -/

/-- ★★**`Γ` が各素点で正の元を持つ**(幾何版の `Q`-Cartier の実係数版)。

★`Example 6.3` では非アルキメデス素点で `log(N v)`、アルキメデス素点で任意の正数が取れる。 -/
def IsGenSubgroupR (Γ : AddSubgroup (S →₀ ℝ)) : Prop :=
  ∀ s : S, ∃ r : ℝ, 0 < r ∧ (single s r) ∈ Γ

/-- ★台が `{s}` の `Φ` の元。 -/
noncomputable def genR (hG : IsGenSubgroupR Γ) (s : S) : effR Γ :=
  ⟨single s (hG s).choose, mem_effR.mpr ⟨(hG s).choose_spec.2, fun t => by
    rcases eq_or_ne s t with rfl | hst
    · simpa using (hG s).choose_spec.1.le
    · simp [hst]⟩⟩

theorem genR_support (hG : IsGenSubgroupR Γ) (s : S) :
    ((genR hG s : effR Γ) : S →₀ ℝ).support = {s} :=
  Finsupp.support_single_ne_zero s (hG s).choose_spec.1.ne'

/-! ## ★5. primary ⟺ 台が 1 点 -/

/-- ★台が 1 点なら primary。 -/
theorem isPrimaryElt_of_support_singleton_R {a : effR Γ} {s : S}
    (h : (a : S →₀ ℝ).support = {s}) : IsPrimaryElt a := by
  classical
  refine ⟨?_, ?_⟩
  · intro hz
    rw [effR_eq_zero_iff] at hz
    rw [hz] at h
    simp at h
  · intro b hb hba
    rw [mprec_effR_iff] at hba ⊢
    rw [h] at hba ⊢
    have hbne : (b : S →₀ ℝ) ≠ 0 := fun hz => hb ((effR_eq_zero_iff b).mpr hz)
    obtain ⟨t, ht⟩ := Finsupp.support_nonempty_iff.mpr hbne
    have : t = s := Finset.mem_singleton.mp (hba ht)
    subst this
    exact Finset.singleton_subset_iff.mpr ht

/-- ★★primary なら台は 1 点。 -/
theorem support_singleton_of_isPrimaryElt_R (hG : IsGenSubgroupR Γ) {a : effR Γ}
    (ha : IsPrimaryElt a) : ∃ s : S, (a : S →₀ ℝ).support = {s} := by
  classical
  have hane : (a : S →₀ ℝ) ≠ 0 := fun hz => ha.1 ((effR_eq_zero_iff a).mpr hz)
  obtain ⟨s, hs⟩ := Finsupp.support_nonempty_iff.mpr hane
  refine ⟨s, ?_⟩
  have hgne : genR hG s ≠ 0 := by
    intro hz
    have hsup := genR_support hG s
    rw [(effR_eq_zero_iff _).mp hz] at hsup
    simp at hsup
  have hprec : MPrec (genR hG s) a := by
    rw [mprec_effR_iff, genR_support hG s]
    exact Finset.singleton_subset_iff.mpr hs
  have hback := ha.2 _ hgne hprec
  rw [mprec_effR_iff, genR_support hG s] at hback
  exact Finset.Subset.antisymm hback (Finset.singleton_subset_iff.mpr hs)

theorem genR_isPrimaryElt (hG : IsGenSubgroupR Γ) (s : S) : IsPrimaryElt (genR hG s) :=
  isPrimaryElt_of_support_singleton_R (genR_support hG s)

/-! ## ★6. `Prime(Φ) ≃ S` -/

/-- ★素元が定める素点。 -/
noncomputable def primePtR (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ)) : S :=
  Quotient.liftOn p (fun x => (support_singleton_of_isPrimaryElt_R hG x.2).choose) (by
    rintro ⟨x, hx⟩ ⟨y, hy⟩ (h : MPrec x y)
    have h1 : (x : S →₀ ℝ).support ⊆ (y : S →₀ ℝ).support := (mprec_effR_iff x y).mp h
    have h2 : (y : S →₀ ℝ).support ⊆ (x : S →₀ ℝ).support :=
      (mprec_effR_iff y x).mp (mprec_symm_of_primary hy hx.1 h)
    have heq : (x : S →₀ ℝ).support = (y : S →₀ ℝ).support := Finset.Subset.antisymm h1 h2
    have hxs := (support_singleton_of_isPrimaryElt_R hG hx).choose_spec
    have hys := (support_singleton_of_isPrimaryElt_R hG hy).choose_spec
    have hsing : ({(support_singleton_of_isPrimaryElt_R hG hx).choose} : Finset S)
        = {(support_singleton_of_isPrimaryElt_R hG hy).choose} := by rw [← hxs, ← hys, heq]
    exact Finset.singleton_injective hsing)

theorem primePtR_toPrime (hG : IsGenSubgroupR Γ) {a : effR Γ} (ha : IsPrimaryElt a)
    {s : S} (hs : (a : S →₀ ℝ).support = {s}) : primePtR hG (toPrime _ a ha) = s := by
  have hspec := (support_singleton_of_isPrimaryElt_R hG ha).choose_spec
  have hsing : ({(support_singleton_of_isPrimaryElt_R hG ha).choose} : Finset S) = {s} := by
    rw [← hspec, hs]
  exact Finset.singleton_injective hsing

/-- ★★★★**`Prime(Φ) ≃ S`**(実係数版)——
[FrdI] `Example 6.3` の「`Prime(Φ(L)) ≃ V(L)`」の骨格。 -/
noncomputable def effRPrimeEquiv (hG : IsGenSubgroupR Γ) : Prime (effR Γ) ≃ S where
  toFun := primePtR hG
  invFun s := toPrime _ (genR hG s) (genR_isPrimaryElt hG s)
  left_inv p := by
    refine Quotient.inductionOn p ?_
    rintro ⟨x, hx⟩
    obtain ⟨s, hs⟩ := support_singleton_of_isPrimaryElt_R hG hx
    have hpt : primePtR hG (toPrime _ x hx) = s := primePtR_toPrime hG hx hs
    show toPrime _ (genR hG (primePtR hG (toPrime _ x hx))) _ = toPrime _ x hx
    refine (toPrime_eq_iff _ hx).mpr ?_
    rw [mprec_effR_iff, genR_support hG (primePtR hG (toPrime _ x hx)), hpt, hs]
  right_inv s := primePtR_toPrime hG (genR_isPrimaryElt hG s) (genR_support hG s)

/-! ## ★7. 台は「`S` の有限部分集合ちょうど」 -/

/-- ★★★**任意の有限集合が `Φ` の元の台になる**(実係数版)。 -/
theorem exists_effR_support_eq (hG : IsGenSubgroupR Γ) (T : Finset S) :
    ∃ a : effR Γ, (a : S →₀ ℝ).support = T := by
  classical
  refine ⟨∑ s ∈ T, genR hG s, ?_⟩
  have hcoe : ((∑ s ∈ T, genR hG s : effR Γ) : S →₀ ℝ)
      = ∑ s ∈ T, ((genR hG s : effR Γ) : S →₀ ℝ) :=
    AddSubmonoid.coe_finsetSum _ _ _
  rw [hcoe, Finsupp.support_sum_eq_biUnion T (fun i₁ i₂ h => by
    rw [genR_support hG i₁, genR_support hG i₂]
    simpa using h)]
  simp [genR_support hG]

/-! ### ★出典の紐付け -/

/-- ★★locator —— `Example 6.3` の `Prime(Φ(L)) ≃ V(L)`(実係数の骨格)。 -/
def effRPrimeEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Prime(Φ(L)) ≃ V(L)(実係数の骨格)",
    sectionId := "frdi-example-6-3" }

/-- ★★locator —— `Example 6.3` の「台は `V(L)` の有限部分集合ちょうど」。 -/
def exists_effR_support_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(L) の元の台は V(L) の有限部分集合ちょうど",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
