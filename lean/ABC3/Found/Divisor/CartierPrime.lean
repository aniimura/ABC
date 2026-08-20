/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.CartierMonoid
import ABC3.Found.FrdI.MonoidPrime

/-!
# `Prime(Φ(L)) ≃ D_L` と、台の同定

★★★[FrdI] `Example 6.1` は `Φ(L) ⊆ ℤ≥0[D_L]` について

> one verifies immediately that Φ(L) is perf-factorial, that there is a natural
> bijection Prime(Φ(L)) → DL, and that the supports of elements of Φ(L) are
> precisely the finite subsets of DL

と書く。★**後半 2 つをここで閉じる**(前半 perf-factorial は 11 条のうち
divisorial(`CartierMonoid.lean`)とこの 2 つが済んだ形になる)。

## ★段取り

`Φ := Γ ∩ ℤ≥0[S]` の上で §0 の 2 つの順序を**係数の言葉に翻訳する**のが要点である。

* `mle_effSub_iff` —— `a ≤ b ⟺ 係数ごとに a s ≤ b s`
  (`b − a` は `Γ` が群だから自動で `Γ` に入る)
* `mprec_effSub_iff` —— **`a ⪯ b ⟺ supp a ⊆ supp b`**。
  ★`⟸` の `n` は `1 + Σ_{s ∈ supp a} a s` で足りる —— `supp b ⊇ supp a` の上で
  `b s ≥ 1` だから `n · b s ≥ n > a s`。

これで `primary` が**台が 1 点**と同値になり(`Q`-Cartier 性は `⟹` にだけ効く ——
各 `s` に台 `{s}` の元 `qcGen` が要る)、`Prime(Φ) ≃ S` が出る。

★台の同定のほうは `qcGen` を足し合わせるだけである
(`Finsupp.support_sum_eq_biUnion`)。
-/

namespace ABC3.Found.FrdI

open Finsupp

variable {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)}

/-! ## ★1. `≤` と `⪯` の係数による翻訳 -/

theorem effSub_nonneg (a : effSub Γ) (s : S) : 0 ≤ (a : S →₀ ℤ) s :=
  (mem_effSub.mp a.2).2 s

theorem effSub_eq_zero_iff (a : effSub Γ) : a = 0 ↔ (a : S →₀ ℤ) = 0 :=
  ⟨fun h => by rw [h]; rfl, fun h => Subtype.ext h⟩

/-- ★`Φ` の `≤` は**係数ごとの `≤`**。

★`b − a` が `Γ` に入るのは `Γ` が**群**だから。 -/
theorem mle_effSub_iff (a b : effSub Γ) :
    MLe a b ↔ ∀ s, (a : S →₀ ℤ) s ≤ (b : S →₀ ℤ) s := by
  constructor
  · rintro ⟨c, rfl⟩ s
    have := effSub_nonneg c s
    show (a : S →₀ ℤ) s ≤ ((a : S →₀ ℤ) + (c : S →₀ ℤ)) s
    simp only [Finsupp.add_apply]
    omega
  · intro h
    have hmem : (b : S →₀ ℤ) - (a : S →₀ ℤ) ∈ effSub Γ := by
      refine mem_effSub.mpr ⟨Γ.sub_mem (mem_effSub.mp b.2).1 (mem_effSub.mp a.2).1, fun s => ?_⟩
      have := h s
      simp only [Finsupp.sub_apply]
      omega
    exact ⟨⟨_, hmem⟩, Subtype.ext (by simp)⟩

/-- ★★**`Φ` の `⪯` は台の包含**。

★`⟸` の `n` は `1 + Σ_{s ∈ supp a} a s` で足りる。 -/
theorem mprec_effSub_iff (a b : effSub Γ) :
    MPrec a b ↔ (a : S →₀ ℤ).support ⊆ (b : S →₀ ℤ).support := by
  classical
  constructor
  · rintro ⟨n, hn, hle⟩ s hs
    rw [mle_effSub_iff] at hle
    have h1 := hle s
    have h2 : ((n • b : effSub Γ) : S →₀ ℤ) s = (n : ℤ) * (b : S →₀ ℤ) s := by
      show ((n • (b : S →₀ ℤ)) : S →₀ ℤ) s = _
      simp
    rw [h2] at h1
    have h3 : (a : S →₀ ℤ) s ≠ 0 := Finsupp.mem_support_iff.mp hs
    have h4 := effSub_nonneg a s
    refine Finsupp.mem_support_iff.mpr (fun hz => ?_)
    rw [hz, mul_zero] at h1
    omega
  · intro hsub
    set n : ℕ := 1 + ((a : S →₀ ℤ).support.sum (fun s => ((a : S →₀ ℤ) s).toNat)) with hn
    have hnpos : 0 < n := by omega
    refine ⟨n, hnpos, ?_⟩
    rw [mle_effSub_iff]
    intro s
    have h2 : ((n • b : effSub Γ) : S →₀ ℤ) s = (n : ℤ) * (b : S →₀ ℤ) s := by
      show ((n • (b : S →₀ ℤ)) : S →₀ ℤ) s = _
      simp
    rw [h2]
    by_cases hs : s ∈ (a : S →₀ ℤ).support
    · have hbs : (b : S →₀ ℤ) s ≠ 0 := Finsupp.mem_support_iff.mp (hsub hs)
      have hbnn := effSub_nonneg b s
      have hb1 : (1 : ℤ) ≤ (b : S →₀ ℤ) s := by omega
      have hle : ((a : S →₀ ℤ) s).toNat
          ≤ (a : S →₀ ℤ).support.sum (fun t => ((a : S →₀ ℤ) t).toNat) :=
        Finset.single_le_sum (f := fun t => ((a : S →₀ ℤ) t).toNat) (fun _ _ => Nat.zero_le _) hs
      have hann := effSub_nonneg a s
      have hcast : ((a : S →₀ ℤ) s) ≤ (n : ℤ) := by
        have : ((a : S →₀ ℤ) s).toNat ≤ n := by omega
        omega
      calc (a : S →₀ ℤ) s ≤ (n : ℤ) := hcast
        _ = (n : ℤ) * 1 := (mul_one _).symm
        _ ≤ (n : ℤ) * (b : S →₀ ℤ) s := mul_le_mul_of_nonneg_left hb1 (by positivity)
    · have h0 : (a : S →₀ ℤ) s = 0 := Finsupp.notMem_support_iff.mp hs
      have hbnn := effSub_nonneg b s
      rw [h0]
      positivity

/-! ## ★2. `Q`-Cartier が与える「素因子そのもの」 -/

/-- ★`Q`-Cartier が与える、台が `{s}` の `Φ` の元。 -/
noncomputable def qcGen (hQ : IsQCartierSubgroup Γ) (s : S) : effSub Γ :=
  ⟨single s ((hQ s).choose : ℤ), mem_effSub.mpr ⟨(hQ s).choose_spec.2, fun t => by
    rcases eq_or_ne s t with rfl | hst
    · simp
    · simp [hst]⟩⟩

theorem qcGen_support (hQ : IsQCartierSubgroup Γ) (s : S) :
    ((qcGen hQ s : effSub Γ) : S →₀ ℤ).support = {s} := by
  have h0 : (hQ s).choose ≠ 0 := Nat.pos_iff_ne_zero.mp (hQ s).choose_spec.1
  have hne : (((hQ s).choose : ℤ)) ≠ 0 := by exact_mod_cast h0
  exact Finsupp.support_single_ne_zero s hne

/-! ## ★3. primary ⟺ 台が 1 点 -/

/-- ★台が 1 点なら primary。 -/
theorem isPrimaryElt_of_support_singleton {a : effSub Γ} {s : S}
    (h : (a : S →₀ ℤ).support = {s}) : IsPrimaryElt a := by
  classical
  refine ⟨?_, ?_⟩
  · intro hz
    rw [effSub_eq_zero_iff] at hz
    rw [hz] at h
    simp at h
  · intro b hb hba
    rw [mprec_effSub_iff] at hba ⊢
    rw [h] at hba ⊢
    have hbne : (b : S →₀ ℤ) ≠ 0 := fun hz => hb ((effSub_eq_zero_iff b).mpr hz)
    obtain ⟨t, ht⟩ := Finsupp.support_nonempty_iff.mpr hbne
    have : t = s := Finset.mem_singleton.mp (hba ht)
    subst this
    exact Finset.singleton_subset_iff.mpr ht

/-- ★★`Q`-Cartier のもとで、primary なら台は 1 点。

★ここで初めて `Q`-Cartier 性が効く —— `s ∈ supp a` に対し台 `{s}` の元 `qcGen hQ s`
を当てて、primary 性から `supp a ⊆ {s}` を得る。 -/
theorem support_singleton_of_isPrimaryElt (hQ : IsQCartierSubgroup Γ) {a : effSub Γ}
    (ha : IsPrimaryElt a) : ∃ s : S, (a : S →₀ ℤ).support = {s} := by
  classical
  have hane : (a : S →₀ ℤ) ≠ 0 := fun hz => ha.1 ((effSub_eq_zero_iff a).mpr hz)
  obtain ⟨s, hs⟩ := Finsupp.support_nonempty_iff.mpr hane
  refine ⟨s, ?_⟩
  have hgne : qcGen hQ s ≠ 0 := by
    intro hz
    have hsup := qcGen_support hQ s
    rw [(effSub_eq_zero_iff _).mp hz] at hsup
    simp at hsup
  have hprec : MPrec (qcGen hQ s) a := by
    rw [mprec_effSub_iff, qcGen_support hQ s]
    exact Finset.singleton_subset_iff.mpr hs
  have hback := ha.2 _ hgne hprec
  rw [mprec_effSub_iff, qcGen_support hQ s] at hback
  exact Finset.Subset.antisymm hback (Finset.singleton_subset_iff.mpr hs)

theorem qcGen_isPrimaryElt (hQ : IsQCartierSubgroup Γ) (s : S) :
    IsPrimaryElt (qcGen hQ s) :=
  isPrimaryElt_of_support_singleton (qcGen_support hQ s)

/-! ## ★4. `Prime(Φ) ≃ D_L` -/

/-- ★素元が定める素因子。 -/
noncomputable def primePt (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) : S :=
  Quotient.liftOn p (fun x => (support_singleton_of_isPrimaryElt hQ x.2).choose) (by
    rintro ⟨x, hx⟩ ⟨y, hy⟩ (h : MPrec x y)
    have h1 : (x : S →₀ ℤ).support ⊆ (y : S →₀ ℤ).support := (mprec_effSub_iff x y).mp h
    have h2 : (y : S →₀ ℤ).support ⊆ (x : S →₀ ℤ).support :=
      (mprec_effSub_iff y x).mp (mprec_symm_of_primary hy hx.1 h)
    have heq : (x : S →₀ ℤ).support = (y : S →₀ ℤ).support := Finset.Subset.antisymm h1 h2
    have hxs := (support_singleton_of_isPrimaryElt hQ hx).choose_spec
    have hys := (support_singleton_of_isPrimaryElt hQ hy).choose_spec
    have hsing : ({(support_singleton_of_isPrimaryElt hQ hx).choose} : Finset S)
        = {(support_singleton_of_isPrimaryElt hQ hy).choose} := by rw [← hxs, ← hys, heq]
    exact Finset.singleton_injective hsing)

theorem primePt_toPrime (hQ : IsQCartierSubgroup Γ) {a : effSub Γ} (ha : IsPrimaryElt a)
    {s : S} (hs : (a : S →₀ ℤ).support = {s}) : primePt hQ (toPrime _ a ha) = s := by
  have hspec := (support_singleton_of_isPrimaryElt hQ ha).choose_spec
  have hsing : ({(support_singleton_of_isPrimaryElt hQ ha).choose} : Finset S) = {s} := by
    rw [← hspec, hs]
  exact Finset.singleton_injective hsing

/-- ★★★★**`Prime(Φ) ≃ D_L`** —— [FrdI] `Example 6.1` の
「there is a natural bijection `Prime(Φ(L))` `→` `D_L`」。 -/
noncomputable def effSubPrimeEquiv (hQ : IsQCartierSubgroup Γ) : Prime (effSub Γ) ≃ S where
  toFun := primePt hQ
  invFun s := toPrime _ (qcGen hQ s) (qcGen_isPrimaryElt hQ s)
  left_inv p := by
    refine Quotient.inductionOn p ?_
    rintro ⟨x, hx⟩
    obtain ⟨s, hs⟩ := support_singleton_of_isPrimaryElt hQ hx
    have hpt : primePt hQ (toPrime _ x hx) = s := primePt_toPrime hQ hx hs
    show toPrime _ (qcGen hQ (primePt hQ (toPrime _ x hx))) _ = toPrime _ x hx
    refine (toPrime_eq_iff _ hx).mpr ?_
    rw [mprec_effSub_iff, qcGen_support hQ (primePt hQ (toPrime _ x hx)), hpt, hs]
  right_inv s := primePt_toPrime hQ (qcGen_isPrimaryElt hQ s) (qcGen_support hQ s)

/-! ## ★5. 台は「`D_L` の有限部分集合ちょうど」 -/

/-- ★★★**任意の有限集合が `Φ` の元の台になる** —— [FrdI] `Example 6.1` の
「the supports of elements of `Φ(L)` are precisely the finite subsets of `D_L`」。

★逆向き(台が有限集合であること)は `Finsupp.support` が `Finset` であることそのもの。 -/
theorem exists_effSub_support_eq (hQ : IsQCartierSubgroup Γ) (T : Finset S) :
    ∃ a : effSub Γ, (a : S →₀ ℤ).support = T := by
  classical
  refine ⟨∑ s ∈ T, qcGen hQ s, ?_⟩
  have hcoe : ((∑ s ∈ T, qcGen hQ s : effSub Γ) : S →₀ ℤ)
      = ∑ s ∈ T, ((qcGen hQ s : effSub Γ) : S →₀ ℤ) :=
    AddSubmonoid.coe_finsetSum _ _ _
  rw [hcoe, Finsupp.support_sum_eq_biUnion T (fun i₁ i₂ h => by
    rw [qcGen_support hQ i₁, qcGen_support hQ i₂]
    simpa using h)]
  simp [qcGen_support hQ]

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.1` の `Prime(Φ(L)) ≃ D_L`。 -/
def effSubPrimeEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Prime(Φ(L)) ≃ D_L",
    sectionId := "frdi-example-6-1" }

/-- ★locator —— `Example 6.1` の「台は `D_L` の有限部分集合ちょうど」。 -/
def exists_effSub_support_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Φ(L) の元の台は D_L の有限部分集合ちょうど",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.FrdI
