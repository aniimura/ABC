/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24SuppIntr

/-!
# 単系の同型に沿った `M^pf` / `Prime(M)` / `Supp` の移送

原文 (FrdI p.48):
> Prime(M ), as a ranges over the elements of M pf; if a, b

★★`Remark 4.5.1` の `rational` 型の移送に要る ——
`B` と `B^istr` の底は**同型だが等しくない**ので、
`Prime` と `Supp` を単系同型で運ぶ必要がある。

★`Supp` が `ι` を含まない形に書けること(`mem_supp_factorMap_iff`)と合わせて、
**`Supp` は単系の同型だけで移る**ことが言える。
-/

namespace ABC3.Found.FrdI

open scoped NNReal

universe w w'

variable {M : Type w} [AddCommMonoid M] {N : Type w'} [AddCommMonoid N]

/-! ## ★1. `M^pf` の移送 -/

/-- ★`M^pf → N^pf` —— 分母はそのまま、分子に `e` を当てる。 -/
def pfMap (e : M →+ N) : Pf M → Pf N :=
  Quotient.map (fun x : M × ℕ+ => (e x.1, x.2)) (by
    rintro ⟨a, n⟩ ⟨b, m⟩ ⟨k, hk⟩
    refine ⟨k, ?_⟩
    show ((k : ℕ) * (m : ℕ)) • e a = ((k : ℕ) * (n : ℕ)) • e b
    rw [← map_nsmul, ← map_nsmul, hk])

@[simp] theorem pfMap_mk (e : M →+ N) (a : M) (n : ℕ+) :
    pfMap e (Pf.mk a n) = Pf.mk (e a) n := rfl

theorem pfMap_comp {O : Type w'} [AddCommMonoid O] (e : M →+ N) (f : N →+ O) (x : Pf M) :
    pfMap f (pfMap e x) = pfMap (f.comp e) x := by
  refine Pf.inductionOn x (fun a n => ?_)
  simp [pfMap_mk]

theorem pfMap_id (x : Pf M) : pfMap (AddMonoidHom.id M) x = x := by
  refine Pf.inductionOn x (fun a n => ?_)
  rfl

/-- ★★`M^pf ≃ N^pf` —— 単系同型から。 -/
def pfEquiv (e : M ≃+ N) : Pf M ≃ Pf N where
  toFun := pfMap e.toAddMonoidHom
  invFun := pfMap e.symm.toAddMonoidHom
  left_inv x := by
    refine Pf.inductionOn x (fun a n => ?_)
    show Pf.mk (e.symm (e a)) n = Pf.mk a n
    rw [e.symm_apply_apply]
  right_inv y := by
    refine Pf.inductionOn y (fun b n => ?_)
    show Pf.mk (e (e.symm b)) n = Pf.mk b n
    rw [e.apply_symm_apply]

/-! ## ★2. `Prime(M)` の移送 -/

/-- ★加法準同型は `≼` を保つ(`Prop110.lean` の `MPrec.map` と同じ内容を、
重い import を避けてここに置く)。 -/
theorem mprec_transport (f : M →+ N) {a b : M} (h : MPrec a b) : MPrec (f a) (f b) := by
  obtain ⟨n, hn, c, hc⟩ := h
  exact ⟨n, hn, f c, by rw [← map_add, hc, map_nsmul]⟩


/-- ★★単系同型は primary 元を保つ。 -/
theorem isPrimaryElt_map (e : M ≃+ N) {a : M} (h : IsPrimaryElt a) : IsPrimaryElt (e a) := by
  refine ⟨fun hz => h.1 ?_, ?_⟩
  · have := congrArg e.symm hz
    rwa [e.symm_apply_apply, map_zero] at this
  · intro b hb hprec
    have hb' : e.symm b ≠ 0 := by
      intro hz
      exact hb (by rw [← e.apply_symm_apply b, hz, map_zero])
    have hprec' : MPrec (e.symm b) a := by
      have h2 := mprec_transport (e.symm.toAddMonoidHom) hprec
      simpa using h2
    have h3 := mprec_transport (e.toAddMonoidHom) (h.2 _ hb' hprec')
    simpa using h3

/-- ★`Prime(M) → Prime(N)`。 -/
def primeMap (e : M ≃+ N) : Prime M → Prime N :=
  Quotient.map (fun x : {a : M // IsPrimaryElt a} => ⟨e x.1, isPrimaryElt_map e x.2⟩)
    (by
      rintro ⟨a, ha⟩ ⟨b, hb⟩ h
      exact mprec_transport (e.toAddMonoidHom) h)

/-- ★★`Prime(M) ≃ Prime(N)`。 -/
def primeEquiv (e : M ≃+ N) : Prime M ≃ Prime N where
  toFun := primeMap e
  invFun := primeMap e.symm
  left_inv p := by
    refine Quotient.inductionOn p (fun x => ?_)
    refine Quotient.sound ?_
    show MPrec (e.symm (e x.1)) x.1
    rw [e.symm_apply_apply]
    exact mprec_refl _
  right_inv q := by
    refine Quotient.inductionOn q (fun y => ?_)
    refine Quotient.sound ?_
    show MPrec (e (e.symm y.1)) y.1
    rw [e.apply_symm_apply]
    exact mprec_refl _

def primeEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 12,
    item := "§0 Monoids — Prime(M) は単系同型で移る",
    sectionId := "frdi-s0-prime" }

/-! ## ★3. 加法性と `Pf.of` -/

theorem pfMap_add (e : M →+ N) (x y : Pf M) :
    pfMap e (x + y) = pfMap e x + pfMap e y := by
  induction x using Pf.inductionOn with | _ a n =>
  induction y using Pf.inductionOn with | _ b m =>
  show pfMap e (Pf.mk ((m : ℕ) • a + (n : ℕ) • b) (n * m))
    = Pf.mk ((m : ℕ) • e a + (n : ℕ) • e b) (n * m)
  rw [pfMap_mk, map_add, map_nsmul, map_nsmul]

theorem pfMap_zero (e : M →+ N) : pfMap e (0 : Pf M) = 0 := by
  show pfMap e (Pf.mk 0 1) = Pf.mk 0 1
  rw [pfMap_mk, map_zero]

theorem pfMap_of (e : M →+ N) (a : M) : pfMap e (Pf.of a) = Pf.of (e a) := rfl

/-- ★★`M^pf ≃+ N^pf`。 -/
def pfAddEquiv (e : M ≃+ N) : Pf M ≃+ Pf N where
  toEquiv := pfEquiv e
  map_add' := pfMap_add e.toAddMonoidHom

@[simp] theorem pfAddEquiv_mk (e : M ≃+ N) (a : M) (n : ℕ+) :
    pfAddEquiv e (Pf.mk a n) = Pf.mk (e a) n := rfl

@[simp] theorem pfAddEquiv_of (e : M ≃+ N) (a : M) :
    pfAddEquiv e (Pf.of a) = Pf.of (e a) := rfl

/-! ## ★4. `primeCarrier` と `pCarrierPf` の対応 -/

theorem mem_primeCarrier_map (e : M ≃+ N) (p : Prime M) (a : M) :
    a ∈ primeCarrier M p ↔ e a ∈ primeCarrier N (primeEquiv e p) := by
  constructor
  · rintro ⟨ha, rfl⟩
    exact ⟨isPrimaryElt_map e ha, rfl⟩
  · rintro ⟨hea, hq⟩
    have ha : IsPrimaryElt a := by
      have h1 := isPrimaryElt_map e.symm hea
      rwa [e.symm_apply_apply] at h1
    refine ⟨ha, ?_⟩
    refine (primeEquiv e).injective ?_
    show toPrime N (e a) (isPrimaryElt_map e ha) = primeEquiv e p
    rw [← hq]

theorem mem_pCarrierPf_map (e : M ≃+ N) (p : Prime M) (x : Pf M) :
    x ∈ pCarrierPf M p ↔ pfAddEquiv e x ∈ pCarrierPf N (primeEquiv e p) := by
  constructor
  · rintro (⟨n, b, hb, hnx⟩ | hz)
    · refine Or.inl ⟨n, e b, (mem_primeCarrier_map e p b).mp hb, ?_⟩
      have h1 : pfAddEquiv e (((n : ℕ+) : ℕ) • x) = pfAddEquiv e (Pf.of b) :=
        congrArg (pfAddEquiv e) hnx
      rw [map_nsmul, pfAddEquiv_of] at h1
      exact h1
    · rw [show x = 0 from hz, map_zero]
      exact Or.inr rfl
  · rintro (⟨n, c, hc, hnx⟩ | hz)
    · have hiff := mem_primeCarrier_map e p (e.symm c)
      rw [e.apply_symm_apply] at hiff
      refine Or.inl ⟨n, e.symm c, hiff.mpr hc, ?_⟩
      refine (pfAddEquiv e).injective ?_
      rw [map_nsmul, hnx, pfAddEquiv_of, e.apply_symm_apply]
    · refine Or.inr ?_
      refine (pfAddEquiv e).injective ?_
      rw [map_zero]
      exact hz

/-! ## ★5. `Bound` と `Supp` の対応 -/

theorem mem_bound_map (e : M ≃+ N) (p : Prime M) (a x : Pf M) :
    x ∈ Bound (Pf M) (pCarrierPf M p) a
      ↔ pfAddEquiv e x ∈ Bound (Pf N) (pCarrierPf N (primeEquiv e p)) (pfAddEquiv e a) := by
  constructor
  · rintro ⟨hx, c, hc⟩
    exact ⟨(mem_pCarrierPf_map e p x).mp hx, pfAddEquiv e c, by rw [← map_add, hc]⟩
  · rintro ⟨hx, c, hc⟩
    refine ⟨(mem_pCarrierPf_map e p x).mpr hx, (pfAddEquiv e).symm c, ?_⟩
    refine (pfAddEquiv e).injective ?_
    rw [map_add, (pfAddEquiv e).apply_symm_apply, hc]

/-- ★★★★★**`Supp` は単系同型で移る**。

★`mem_supp_factorMap_iff` で `ι` を消してあるので、
**両側で違う `ι` を使っていてよい**。 -/
theorem mem_supp_map {ιM : Prime M → Pf M → ℝ≥0} {ιN : Prime N → Pf N → ℝ≥0}
    (HM : IsPerfFactorialWith M ιM) (HN : IsPerfFactorialWith N ιN)
    (e : M ≃+ N) (a : Pf M) (p : Prime M) :
    p ∈ Supp (factorMap ιM a)
      ↔ primeEquiv e p ∈ Supp (factorMap ιN (pfAddEquiv e a)) := by
  rw [mem_supp_factorMap_iff HM, mem_supp_factorMap_iff HN]
  constructor
  · rintro ⟨x, hx, hx0⟩
    refine ⟨pfAddEquiv e x, (mem_bound_map e p a x).mp hx, ?_⟩
    intro hz
    exact hx0 ((pfAddEquiv e).injective (by rw [hz, map_zero]))
  · rintro ⟨y, hy, hy0⟩
    refine ⟨(pfAddEquiv e).symm y, ?_, ?_⟩
    · refine (mem_bound_map e p a _).mpr ?_
      rwa [(pfAddEquiv e).apply_symm_apply]
    · intro hz
      exact hy0 (by rw [← (pfAddEquiv e).apply_symm_apply y, hz, map_zero])

def mem_supp_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 48,
    item := "Definition 2.4, (i), (d) — Supp は単系同型で移る",
    sectionId := "frdi-def-2-4" }

/-! ## ★sharp なモノイドでの non-dilating —— `M^char` を `M` に置き換える

★`IsNonDilating` の言明は `M^char` の上に書かれているが、**sharp なら `M ≃+ M^char`** なので
`M` の上の条件に置き換えられる。★`Theorem 6.2, (iii)` / `Theorem 6.4, (i)` はどちらも
sharp な `Φ(L)`(有効因子の単系)を相手にするので、この形が要る。 -/

/-- ★★**sharp なら `M ≃+ M^char`**。 -/
noncomputable def mCharEquivOfSharp (hs : IsSharp M) : M ≃+ MChar M :=
  AddEquiv.ofBijective toChar ⟨toChar_injective_of_isSharp hs, toChar_surjective M⟩

@[simp] theorem mCharEquivOfSharp_apply (hs : IsSharp M) (a : M) :
    mCharEquivOfSharp hs a = toChar a := rfl

/-- ★★★★**non-dilating の十分条件(sharp 版)** ——
`M` が primary 元で生成され、primary 元の上で `α` が恒等なら `α` は non-dilating。

★★これが `Theorem 6.2, (iii)` / `Theorem 6.4, (i)` が「immediately / clearly」と
畳んだ段の骨組みである ——
`hfix` が「自己同型は素因子を素因子へ(係数 1 で)移す」に、
`hgen` が「`Φ(L)` は素因子で生成される」に対応する。 -/
theorem isNonDilating_of_primary_sharp (hs : IsSharp M) (α : M →+ M)
    (hgen : AddSubmonoid.closure {a : M | IsPrimaryElt a} = ⊤)
    (hfix : ∀ a : M, IsPrimaryElt a → MPrec (α a) a → α a = a) :
    IsNonDilating α := by
  refine isNonDilating_of_primary α ?_ ?_
  · refine eq_top_iff.mpr fun x _ => ?_
    obtain ⟨y, rfl⟩ := toChar_surjective M x
    have hy : y ∈ AddSubmonoid.closure {a : M | IsPrimaryElt a} := by rw [hgen]; trivial
    refine AddSubmonoid.closure_induction ?_ ?_ ?_ hy
    · intro a ha
      exact AddSubmonoid.subset_closure (isPrimaryElt_map (mCharEquivOfSharp hs) ha)
    · rw [map_zero]; exact AddSubmonoid.zero_mem _
    · intro u v _ _ hu hv
      rw [map_add]; exact AddSubmonoid.add_mem _ hu hv
  · intro a ha hprec
    obtain ⟨b, rfl⟩ := toChar_surjective M a
    have hb : IsPrimaryElt b :=
      (mCharEquivOfSharp hs).symm_apply_apply b ▸ isPrimaryElt_map (mCharEquivOfSharp hs).symm ha
    rw [charMap_toChar] at hprec ⊢
    have hp : MPrec (α b) b := by
      have h2 := mprec_map ((mCharEquivOfSharp hs).symm : MChar M →+ M) hprec
      have h3 : ((mCharEquivOfSharp hs).symm : MChar M →+ M) (toChar (α b)) = α b :=
        (mCharEquivOfSharp hs).symm_apply_apply (α b)
      have h4 : ((mCharEquivOfSharp hs).symm : MChar M →+ M) (toChar b) = b :=
        (mCharEquivOfSharp hs).symm_apply_apply b
      rwa [h3, h4] at h2
    rw [hfix b hb hp]

/-- ★★★locator —— sharp なモノイドでの non-dilating の十分条件。 -/
def isNonDilating_of_primary_sharp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 19,
    item := "Definition 1.1, (i) — non-dilating の十分条件(sharp 版)",
    sectionId := "frdi-def-1-1-i" }

end ABC3.Found.FrdI
