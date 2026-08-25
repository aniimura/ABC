/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55RlfRefl

/-!
# [FrdI] §0 —— `Prime(M) → Prime(M^rlf)` は全単射

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.12 / p.48。

原文 (FrdI p.12):
> the set of primes of M. If p

## ★★何のために要るか

`Theorem 6.4, (i)` は `𝒞, 𝒞^pf, 𝒞^rlf, 𝒞^un-tr, (𝒞^pf)^un-tr` の 5 圏が
rationally standard だと言う。★そのうち **`𝒞^rlf`** を出すには
`Definition 4.5, (ii)` の rational 型が要り、それは

    ∀ p : Prime(Φ^rlf(A)), ∃ a b, toGp a − toGp b ∈ Φ^birat ∧ p ∈ Supp a ∧ p ∉ Supp b

と**すべての素点について**量化する。★したがって
`Prime(Φ^rlf)` と `Prime(Φ)` の対応が要る。

★`𝒞^pf` のときは在庫の `primeToPf_bijective`(`MonoidPrime.lean`)が使えたが、
**実化の側には無かった**(依存グラフの `rlf-agree` の隣)。本ファイルがそれを埋める。

## ★★★中身 —— 台が 1 点であることだけ

`M^rlf = ℝ≥0 ⊗_ℕ M` は `M` より**真に大きい**
(ℕ 上のテンソルなので `r ⊗ a ≠ 1 ⊗ (r·a)` が一般に起きる)。
★★それでも **`Prime` は変わらない**。理由は 3 段:

1. **係数は均せる** —— `r ≠ 0` なら `r ⊗ a` と `1 ⊗ a` は `⪯` で同値
   (`mprec_tmul_toSc` / `mprec_toSc_tmul`。`n ≥ r` あるいは `n·r ≥ 1` を取るだけ)。
2. **primary な元は単項式に落ちる** —— `exists_tmul_mle`(在庫)で下に `r ⊗ a` を取り、
   primary 性で戻せば `w ~ r ⊗ a ~ 1 ⊗ a`(`exists_toSc_mprec_of_primary`)。
3. **`1 ⊗ a` が primary ⟺ `a` が primary** —— どちらも
   **`SuppElt ι a` が 1 点**であることに帰着する
   (`isPrimaryElt_iff_suppElt_singleton` ＋ `suppElt_subset_of_mprec_sc`)。

★★★鍵は在庫の **`suppElt_subset_of_mprec_sc`**(`Prop55RlfRefl.lean`)——
「`ℝ≥0 ⊗` の上の `⪯` は台の包含を与える」。
★これは**片道の橋** `scToFactor` の核が自明であることから出ており、
2 模型の一致(`rlf-agree`)は要らない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal TensorProduct

universe w

variable {M : Type w} [AddCommMonoid M] {ι : Prime M → Pf M → ℝ≥0}

/-! ## ★1. 係数は `1` に均せる -/

/-- ★`n • (x ⊗ a) = (n·x) ⊗ a`。 -/
theorem nsmul_tmul_nnreal (n : ℕ) (x : ℝ≥0) (a : M) :
    n • (x ⊗ₜ[ℕ] a) = (((n : ℝ≥0)) * x) ⊗ₜ[ℕ] a := by
  induction n with
  | zero => simp
  | succ k ih =>
      rw [succ_nsmul, ih, ← TensorProduct.add_tmul]
      congr 1
      push_cast
      ring

theorem toSc_eq_one_tmul (a : M) : toSc (S := ℝ≥0) a = (1 : ℝ≥0) ⊗ₜ[ℕ] a := rfl

/-- ★★`r ⊗ a ⪯ 1 ⊗ a` —— `n ≥ r` を取ればよい。 -/
theorem mprec_tmul_toSc (r : ℝ≥0) (a : M) :
    MPrec (r ⊗ₜ[ℕ] a) (toSc (S := ℝ≥0) a) := by
  obtain ⟨n, hn⟩ := exists_nat_ge (r : ℝ)
  have hle : r ≤ ((n : ℝ≥0)) := by exact_mod_cast hn
  refine ⟨n + 1, Nat.succ_pos n, (((n + 1 : ℕ) : ℝ≥0) - r) ⊗ₜ[ℕ] a, ?_⟩
  rw [toSc_eq_one_tmul, nsmul_tmul_nnreal, mul_one, ← TensorProduct.add_tmul,
    add_tsub_cancel_of_le (le_trans hle (by push_cast; exact le_add_of_nonneg_right zero_le_one))]

/-- ★★`1 ⊗ a ⪯ r ⊗ a`(`r ≠ 0` なら `n·r ≥ 1` を取ればよい)。 -/
theorem mprec_toSc_tmul {r : ℝ≥0} (hr : r ≠ 0) (a : M) :
    MPrec (toSc (S := ℝ≥0) a) (r ⊗ₜ[ℕ] a) := by
  obtain ⟨n, hn⟩ := exists_nat_ge ((r : ℝ)⁻¹)
  have hrpos : (0 : ℝ) < r := lt_of_le_of_ne r.2 (fun h => hr (NNReal.coe_eq_zero.mp h.symm))
  have h1 : (1 : ℝ≥0) ≤ (((n + 1 : ℕ) : ℝ≥0)) * r := by
    have hR : (1 : ℝ) ≤ (((n : ℝ)) + 1) * r := by
      have hn' : (r : ℝ)⁻¹ ≤ (n : ℝ) + 1 := le_trans hn (by linarith)
      have hinv : (r : ℝ)⁻¹ * r = 1 := by field_simp
      nlinarith [hrpos, hn']
    have hc : ((1 : ℝ≥0) : ℝ) ≤ ((((n + 1 : ℕ) : ℝ≥0)) * r : ℝ≥0) := by push_cast; exact hR
    exact_mod_cast hc
  refine ⟨n + 1, Nat.succ_pos n, ((((n + 1 : ℕ) : ℝ≥0) * r) - 1) ⊗ₜ[ℕ] a, ?_⟩
  rw [toSc_eq_one_tmul, nsmul_tmul_nnreal, ← TensorProduct.add_tmul,
    add_tsub_cancel_of_le h1]

/-! ## ★2. primary な元は `1 ⊗ a` に落ちる -/

/-- ★`a ≠ 0` かつ `r ≠ 0` なら `r ⊗ a ≠ 0`(片道の橋の核が自明だから)。 -/
theorem tmul_ne_zero (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {r : ℝ≥0} (hr : r ≠ 0) {a : M} (ha : a ≠ 0) : (r ⊗ₜ[ℕ] a : ScT ℝ≥0 M) ≠ 0 := by
  intro h0
  have hkey : scToFactor H (r ⊗ₜ[ℕ] a) = 0 := by rw [h0, map_zero]
  rw [scToFactor_tmul] at hkey
  have hμ : factorMap ι (Pf.of a) ≠ 0 := by
    intro hz
    refine ha (suppFactorHom_injective H hdiv ?_)
    rw [suppFactorHom_apply, suppFactorHom_apply, hz, map_zero]
    exact (factorMap_zero H).symm
  refine hμ (funext fun p => ?_)
  have hp := congrFun hkey p
  rw [Pi.smul_apply, Pi.zero_apply, smul_eq_mul] at hp
  exact (mul_eq_zero.mp hp).resolve_left hr

/-- ★★★★**`M^rlf` の primary な元はすべて `1 ⊗ a` と同値**。 -/
theorem exists_toSc_mprec_of_primary (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {w : ScT ℝ≥0 M} (hw : IsPrimaryElt w) :
    ∃ a : M, a ≠ 0 ∧ MPrec w (toSc (S := ℝ≥0) a) ∧ MPrec (toSc (S := ℝ≥0) a) w := by
  obtain ⟨r, a, hr, ha, hle⟩ := exists_tmul_mle hw.1
  have hne : (r ⊗ₜ[ℕ] a : ScT ℝ≥0 M) ≠ 0 := tmul_ne_zero H hdiv hr ha
  have h1 : MPrec (r ⊗ₜ[ℕ] a : ScT ℝ≥0 M) w := mprec_of_mle hle
  have h2 : MPrec w (r ⊗ₜ[ℕ] a : ScT ℝ≥0 M) := hw.2 _ hne h1
  exact ⟨a, ha, mprec_trans h2 (mprec_tmul_toSc r a),
    mprec_trans (mprec_toSc_tmul hr a) h1⟩

/-! ## ★3. `1 ⊗ a` が primary ⟺ `a` が primary -/

/-- ★primary 性は `⪯` の同値で移る。 -/
theorem isPrimaryElt_of_mprec_equiv {N : Type w} [AddCommMonoid N] {a b : N}
    (ha : IsPrimaryElt a) (hb : b ≠ 0) (_h1 : MPrec a b) (h2 : MPrec b a) :
    IsPrimaryElt b :=
  ⟨hb, fun c hc hcb => mprec_trans h2 (ha.2 c hc (mprec_trans hcb h2))⟩

/-- ★★★★★**`1 ⊗ a` が primary なら `a` も primary**。

★台を 1 点に切り出して(`exists_split_suppElt`)、
`toSc` を通した primary 性で戻し(`suppElt_subset_of_mprec_sc`)、
`SuppElt ι a` が 1 点だと結論する。 -/
theorem isPrimaryElt_of_toSc (H : IsPerfFactorialWith M ι) (hperf : IsPerfectMonoid M)
    (hdiv : IsDivisorial M) {a : M} (h : IsPrimaryElt (toSc (S := ℝ≥0) a)) :
    IsPrimaryElt a := by
  classical
  have ha : a ≠ 0 := by
    intro h0
    exact h.1 (by rw [h0, map_zero])
  obtain ⟨p, hp⟩ := Set.nonempty_iff_ne_empty.mpr (suppElt_ne_empty H hdiv ha)
  obtain ⟨y, z, hsum, hy, -⟩ := exists_split_suppElt H hperf hdiv a {p}
  have hyS : SuppElt ι y = {p} := by
    rw [hy]
    exact Set.inter_eq_left.mpr (Set.singleton_subset_iff.mpr hp)
  have hy0 : y ≠ 0 := by
    intro h0
    rw [h0, suppElt_zero H] at hyS
    exact absurd hyS.symm (Set.singleton_ne_empty p)
  have hle : MPrec (toSc (S := ℝ≥0) y) (toSc (S := ℝ≥0) a) :=
    mprec_map _ (mprec_of_mle ⟨z, hsum.symm⟩)
  have hback : MPrec (toSc (S := ℝ≥0) a) (toSc (S := ℝ≥0) y) :=
    h.2 _ (toSc_ne_zero H hdiv hy0) hle
  have hsub : SuppElt ι a ⊆ SuppElt ι y := suppElt_subset_of_mprec_sc H hback
  refine isPrimaryElt_of_suppElt_singleton H hdiv ha (P := p) ?_
  rw [hyS] at hsub
  exact Set.Subset.antisymm hsub (Set.singleton_subset_iff.mpr hp)

/-- ★★**`a` が primary なら `1 ⊗ a` も primary**。

★★在庫の `isPrimaryElt_toSc`(`Prop55RlfRefl.lean`)は台が 1 点であることを
直接受け取る形なので、`suppElt_singleton_of_primary` を挟むだけ。 -/
theorem isPrimaryElt_toSc_of_primary (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M) {a : M} (h : IsPrimaryElt a) :
    IsPrimaryElt (toSc (S := ℝ≥0) a) :=
  let ⟨_, hp⟩ := suppElt_singleton_of_primary H hperf hdiv h
  isPrimaryElt_toSc H hdiv hp h.1

/-! ## ★4. `Prime(M) ≃ Prime(M^rlf)` -/

variable (H : IsPerfFactorialWith M ι) (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M)

/-- ★★★★**`Prime(M) → Prime(M^rlf)`**(`primeToPf` の実化版)。 -/
noncomputable def primeToSc : Prime M → Prime (ScT ℝ≥0 M) :=
  Quotient.map (fun x => ⟨toSc (S := ℝ≥0) x.1, isPrimaryElt_toSc_of_primary H hperf hdiv x.2⟩)
    (fun _ _ h => mprec_map _ h)

@[simp] theorem primeToSc_mk {a : M} (ha : IsPrimaryElt a) :
    primeToSc H hperf hdiv (toPrime M a ha)
      = toPrime (ScT ℝ≥0 M) (toSc a) (isPrimaryElt_toSc_of_primary H hperf hdiv ha) := rfl

/-- ★★★★★★**`Prime(M) → Prime(M^rlf)` は全単射**。

★★これが `Theorem 6.4, (i)` の 5 圏のうち `𝒞^rlf` を出すための土台である
(`𝒞^pf` のときの `primeToPf_bijective` に当たるもの)。

★単射性は台が 1 点であること(`suppElt_subset_of_mprec_sc` ＋
`mprec_of_suppElt_eq_singleton`)、
全射性は「primary な元はすべて `1 ⊗ a` と同値」(`exists_toSc_mprec_of_primary`)。 -/
theorem primeToSc_bijective : Function.Bijective (primeToSc H hperf hdiv) := by
  classical
  constructor
  · refine Quotient.ind fun x => Quotient.ind fun y => ?_
    intro hxy
    have hsc : MPrec (toSc (S := ℝ≥0) x.1) (toSc (S := ℝ≥0) y.1) := Quotient.exact hxy
    have hsub : SuppElt ι x.1 ⊆ SuppElt ι y.1 := suppElt_subset_of_mprec_sc H hsc
    obtain ⟨p, hp⟩ := suppElt_singleton_of_primary H hperf hdiv x.2
    obtain ⟨q, hq⟩ := suppElt_singleton_of_primary H hperf hdiv y.2
    have hpq : p = q := by
      have h1 : p ∈ SuppElt ι x.1 := by rw [hp]; exact rfl
      have h2 := hsub h1
      rwa [hq, Set.mem_singleton_iff] at h2
    subst hpq
    exact Quotient.sound (mprec_of_suppElt_eq_singleton H hdiv hp hq)
  · refine Quotient.ind fun x => ?_
    obtain ⟨a, ha, h1, h2⟩ := exists_toSc_mprec_of_primary H hdiv x.2
    have hsc : IsPrimaryElt (toSc (S := ℝ≥0) a) :=
      isPrimaryElt_of_mprec_equiv x.2 (toSc_ne_zero H hdiv ha) h1 h2
    have haP : IsPrimaryElt a := isPrimaryElt_of_toSc H hperf hdiv hsc
    exact ⟨toPrime M a haP, Quotient.sound h2⟩

/-! ### ★出典の紐付け -/

def primeToSc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 12,
    item := "§0 Monoids — Prime(M) → Prime(M^rlf)",
    sectionId := "frdi-s0-prime" }

def primeToSc_bijective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 12,
    item := "§0 Monoids — Prime(M) → Prime(M^rlf) は全単射",
    sectionId := "frdi-s0-prime" }

def primeToSc_bijective.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "suppElt_subset_of_mprec_sc(ℝ≥0 ⊗ の上の ⪯ は台の包含を与える)"
      (.inProject "ABC3" "ABC3.Found.FrdI.suppElt_subset_of_mprec_sc") 48,
    .citation "[ABC3]" "exists_tmul_mle(0 でないテンソルには 0 でない単項式が下にある)"
      (.inProject "ABC3" "ABC3.Found.FrdI.exists_tmul_mle") 48,
    .citation "[ABC3]" "isPrimaryElt_of_suppElt_singleton / mprec_of_suppElt_eq_singleton"
      (.inProject "ABC3" "ABC3.Found.FrdI.isPrimaryElt_of_suppElt_singleton") 12,
    .implicitStep
      ("★`M^rlf = ℝ≥0 ⊗_ℕ M` は `M` より真に大きい(ℕ 上のテンソルなので " ++
       "r ⊗ a ≠ 1 ⊗ (r·a) が一般に起きる)が、**素点は変わらない**。" ++
       "係数を 1 に均せる(mprec_tmul_toSc / mprec_toSc_tmul)のが要点で、" ++
       "2 模型の一致(rlf-agree)は要らない") 48 ]

end ABC3.Found.FrdI
