/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeCartierPull
import ABC3.Found.Divisor.NormCodim
import ABC3.Found.Divisor.NormKQC

/-!
# `Φ(L)` の関手性(鎖 `normalize` の `phi-functor`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> If for every Spec(L) ∈Ob(B(G)0) [cf. §0], every prime divisor of DL is Q-Cartier,

## ★★台が `D_L` に入ることは引き戻しでも保たれる

`w ∉ D_M` なら `g(w) ∉ D_L`(`not_mem_DLSet_of_not_mem`)なので、
`D` の局所方程式 `f` は `g(w)` で `ord = 0`、したがって茎の単元、
したがって `ord_w(g^*f) = 0` である。
★`g(w)` が余次元 1 でないときは茎が体になるので、やはり `0`。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `pullCoeff_eq_zero_of_notMem_DL` | ★★★**台は `D_M` に入る** |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

/-- ★★★**引き戻しの台は `D_M` に入る** —— `D` の台が `D_L` に入るとき。 -/
theorem pullCoeff_eq_zero_of_notMem_DL (DK : Set (PrimeDivisorPt V))
    {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    [IsLocallyNoetherian (normObj V L)] [IsLocallyNoetherian (normObj V M)]
    {D : WeilDiv (normObj V L)} (hD : IsCartierDiv (normObj_isNormalScheme V L) D)
    (hsupp : ∀ z : PrimeDivisorPt (normObj V L), z ∉ DLSet V DK L → D z = 0)
    (w : PrimeDivisorPt (normObj V M)) (hw : w ∉ DLSet V DK M) :
    pullCoeff (normMap V f) (normObj_isNormalScheme V M)
      (normObj_isNormalScheme V L) hD w = 0 := by
  have hdim : ringKrullDim ((normObj V L).presheaf.stalk ((normMap V f).base w.1)) ≤ 1 :=
    w.2 ▸ ringKrullDim_stalk_normMap_le V f w.1
  obtain ⟨U, hwU, q, hq, hDU⟩ := hD ((normMap V f).base w.1)
  rw [pullCoeff_eq (normMap V f) (normObj_isNormalScheme V M) (normObj_isNormalScheme V L)
    hD w hdim hq hwU hDU]
  obtain ⟨t, ht⟩ := exists_unit_of_ordPt_of_dim_le_one (normObj_isNormalScheme V L)
    ((normMap V f).base w.1) hdim (Units.mk0 q hq) (by
      intro hc
      have hznot : (⟨(normMap V f).base w.1, hc⟩ : PrimeDivisorPt (normObj V L))
          ∉ DLSet V DK L := not_mem_DLSet_of_not_mem V DK f w hw _ rfl
      rw [show ((Units.mk0 q hq : ((normObj V L).functionField)ˣ)
        : (normObj V L).functionField) = q from rfl,
        ← hDU ⟨(normMap V f).base w.1, hc⟩ hwU]
      exact hsupp _ hznot)
  rw [show q = ((Units.mk0 q hq : ((normObj V L).functionField)ˣ)
    : (normObj V L).functionField) from rfl, ← ht]
  exact ordPt_ffMap_eq_zero_of_isUnit (normMap V f) (normObj_isNormalScheme V M) w t

/-! ## ★2. `Φ(L)^gp` の引き戻し -/

variable {V} in
/-- ★`codim-preserved` が与える次元の上界(引き戻しに渡す形)。 -/
theorem normMap_hdim {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    (w : PrimeDivisorPt (normObj V M)) :
    ringKrullDim ((normObj V L).presheaf.stalk ((normMap V f).base w.1)) ≤ 1 :=
  w.2 ▸ ringKrullDim_stalk_normMap_le V f w.1

/-- ★台が `D_M` に入る Weil 因子は `D_M` 添字の `Finsupp` から来る。 -/
theorem embDomain_subtypeDomain (DK : Set (PrimeDivisorPt V))
    (M : FinSub V.functionField Kbar) (E : WeilDiv (normObj V M))
    (hE : ∀ w : PrimeDivisorPt (normObj V M), w ∉ DLSet V DK M → E w = 0) :
    toWeilOnDL V DK M (Finsupp.subtypeDomain (· ∈ DLSet V DK M) E) = E := by
  refine Finsupp.ext fun x => ?_
  show Finsupp.embDomain (DLEmb V DK M) (Finsupp.subtypeDomain _ E) x = E x
  by_cases hx : x ∈ DLSet V DK M
  · show Finsupp.embDomain (DLEmb V DK M) (Finsupp.subtypeDomain _ E)
        (DLEmb V DK M ⟨x, hx⟩) = E x
    rw [Finsupp.embDomain_apply_self]
    rfl
  · rw [Finsupp.embDomain_notin_range, (hE x hx).symm]
    rintro ⟨y, rfl⟩
    exact hx y.2

open scoped Classical in
/-- ★★★★★★**`Φ(L)^gp → Φ(M)^gp`** —— Cartier 因子の引き戻し。 -/
noncomputable def phiPull (DK : Set (PrimeDivisorPt V))
    {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    [IsLocallyNoetherian (normObj V L)] [IsLocallyNoetherian (normObj V M)]
    [CompactSpace (normObj V M)] :
    cartierOnDL V DK L (normObj_isNormalScheme V L)
      →+ cartierOnDL V DK M (normObj_isNormalScheme V M) where
  toFun d := ⟨Finsupp.subtypeDomain (· ∈ DLSet V DK M)
      (cartierPullback (normMap V f) (normObj_isNormalScheme V M)
        (normObj_isNormalScheme V L) d.2 (normMap_hdim f)), by
    show toWeilOnDL V DK M _ ∈ cartierSubgroup (normObj_isNormalScheme V M)
    rw [embDomain_subtypeDomain V DK M _ (fun w hw => by
      rw [cartierPullback_apply]
      exact pullCoeff_eq_zero_of_notMem_DL V DK f d.2
        (fun z hz => by
          show Finsupp.embDomain (DLEmb V DK L) (d : (DLSet V DK L) →₀ ℤ) z = 0
          refine Finsupp.embDomain_notin_range _ _ _ ?_
          rintro ⟨y, rfl⟩
          exact hz y.2) w hw)]
    exact isCartierDiv_cartierPullback (normMap V f) (normObj_isNormalScheme V M)
      (normObj_isNormalScheme V L) d.2 (normMap_hdim f)⟩
  map_zero' := by
    refine Subtype.ext (Finsupp.ext fun x => ?_)
    show pullCoeff (normMap V f) (normObj_isNormalScheme V M)
      (normObj_isNormalScheme V L) (0 : cartierOnDL V DK L _).2 x.1 = 0
    rw [pullCoeff_eq (normMap V f) (normObj_isNormalScheme V M) (normObj_isNormalScheme V L)
      (0 : cartierOnDL V DK L _).2 x.1 (normMap_hdim f x.1) (U := ⊤) (f := 1) (one_ne_zero)
      (Set.mem_univ _) (fun v _ => by simp [toWeilOnDL, ordPt_one])]
    simp [ordPt_one]
  map_add' d e := by
    refine Subtype.ext (Finsupp.ext fun x => ?_)
    show pullCoeff (normMap V f) (normObj_isNormalScheme V M)
        (normObj_isNormalScheme V L) (d + e).2 x.1
      = pullCoeff (normMap V f) (normObj_isNormalScheme V M)
          (normObj_isNormalScheme V L) d.2 x.1
        + pullCoeff (normMap V f) (normObj_isNormalScheme V M)
          (normObj_isNormalScheme V L) e.2 x.1
    have hDeq : toWeilOnDL V DK L ((d : (DLSet V DK L) →₀ ℤ) + (e : (DLSet V DK L) →₀ ℤ))
        = toWeilOnDL V DK L (d : (DLSet V DK L) →₀ ℤ)
          + toWeilOnDL V DK L (e : (DLSet V DK L) →₀ ℤ) := map_add _ _ _
    have hde : IsCartierDiv (normObj_isNormalScheme V L)
        (toWeilOnDL V DK L (d : (DLSet V DK L) →₀ ℤ)
          + toWeilOnDL V DK L (e : (DLSet V DK L) →₀ ℤ)) := hDeq ▸ (d + e).2
    rw [pullCoeff_congr (normMap V f) (normObj_isNormalScheme V M)
      (normObj_isNormalScheme V L) (d + e).2 hde hDeq x.1]
    exact pullCoeff_add (normMap V f) (normObj_isNormalScheme V M)
      (normObj_isNormalScheme V L) d.2 e.2 hde x.1 (normMap_hdim f x.1)

/-! ## ★3. 関手則 -/

/-- ★`Φ(L)^gp` の元としての引き戻しは、Weil 因子として見ても係数そのもの。 -/
theorem toWeilOnDL_phiPull (DK : Set (PrimeDivisorPt V))
    {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    [IsLocallyNoetherian (normObj V L)] [IsLocallyNoetherian (normObj V M)]
    [CompactSpace (normObj V M)]
    (d : cartierOnDL V DK L (normObj_isNormalScheme V L))
    (v : PrimeDivisorPt (normObj V M)) :
    toWeilOnDL V DK M ((phiPull V DK f d : cartierOnDL V DK M _)
        : (DLSet V DK M) →₀ ℤ) v
      = pullCoeff (normMap V f) (normObj_isNormalScheme V M)
        (normObj_isNormalScheme V L) d.2 v := by
  have hsupp : ∀ z : PrimeDivisorPt (normObj V L), z ∉ DLSet V DK L →
      toWeilOnDL V DK L (d : (DLSet V DK L) →₀ ℤ) z = 0 := by
    intro z hz
    show Finsupp.embDomain (DLEmb V DK L) (d : (DLSet V DK L) →₀ ℤ) z = 0
    refine Finsupp.embDomain_notin_range _ _ _ ?_
    rintro ⟨y, rfl⟩
    exact hz y.2
  have := embDomain_subtypeDomain V DK M
    (cartierPullback (normMap V f) (normObj_isNormalScheme V M)
      (normObj_isNormalScheme V L) d.2 (normMap_hdim f))
    (fun w hw => by
      rw [cartierPullback_apply]
      exact pullCoeff_eq_zero_of_notMem_DL V DK f d.2 hsupp w hw)
  rw [show ((phiPull V DK f d : cartierOnDL V DK M _) : (DLSet V DK M) →₀ ℤ)
      = Finsupp.subtypeDomain (· ∈ DLSet V DK M)
        (cartierPullback (normMap V f) (normObj_isNormalScheme V M)
          (normObj_isNormalScheme V L) d.2 (normMap_hdim f)) from rfl, this]
  rfl

/-- ★★★**`pull_id`**。 -/
theorem phiPull_id (DK : Set (PrimeDivisorPt V)) (L : FinSub V.functionField Kbar)
    [IsLocallyNoetherian (normObj V L)] [CompactSpace (normObj V L)]
    (d : cartierOnDL V DK L (normObj_isNormalScheme V L)) :
    phiPull V DK (𝟙 L) d = d := by
  refine Subtype.ext (Finsupp.ext fun x => ?_)
  show pullCoeff (normMap V (𝟙 L)) (normObj_isNormalScheme V L)
    (normObj_isNormalScheme V L) d.2 x.1 = (d : (DLSet V DK L) →₀ ℤ) x
  have hpt : (normMap V (𝟙 L)).base x.1.1 = x.1.1 := by
    rw [normMap_id]
    rfl
  obtain ⟨U, hwU, q, hq, hDU⟩ := d.2 ((normMap V (𝟙 L)).base x.1.1)
  rw [pullCoeff_eq (normMap V (𝟙 L)) (normObj_isNormalScheme V L)
    (normObj_isNormalScheme V L) d.2 x.1 (normMap_hdim (𝟙 L) x.1) hq hwU hDU]
  have hff : ffMap (normMap V (𝟙 L)) q = q := by
    have h1 := normFF_id V L
    rw [normFF_eq_ffMap] at h1
    rw [h1]
    rfl
  rw [hff, ← hDU x.1 (hpt ▸ hwU)]
  exact Finsupp.embDomain_apply_self (DLEmb V DK L) (d : (DLSet V DK L) →₀ ℤ) x

/-- ★★★★**`pull_comp`** —— `g^*f^* = (f≫g)^*`。 -/
theorem phiPull_comp (DK : Set (PrimeDivisorPt V))
    {L M N : FinSub V.functionField Kbar} (f : L ⟶ M) (g : M ⟶ N)
    [IsLocallyNoetherian (normObj V L)] [IsLocallyNoetherian (normObj V M)]
    [IsLocallyNoetherian (normObj V N)]
    [CompactSpace (normObj V M)] [CompactSpace (normObj V N)]
    (d : cartierOnDL V DK L (normObj_isNormalScheme V L)) :
    phiPull V DK (f ≫ g) d = phiPull V DK g (phiPull V DK f d) := by
  refine Subtype.ext (Finsupp.ext fun x => ?_)
  show pullCoeff (normMap V (f ≫ g)) (normObj_isNormalScheme V N)
      (normObj_isNormalScheme V L) d.2 x.1
    = pullCoeff (normMap V g) (normObj_isNormalScheme V N)
      (normObj_isNormalScheme V M) (phiPull V DK f d).2 x.1
  obtain ⟨U, hwU, q, hq, hDU⟩ := d.2 ((normMap V (f ≫ g)).base x.1.1)
  have hpt : (normMap V (f ≫ g)).base x.1.1
      = (normMap V f).base ((normMap V g).base x.1.1) := by
    rw [normMap_comp]
    rfl
  rw [pullCoeff_eq (normMap V (f ≫ g)) (normObj_isNormalScheme V N)
    (normObj_isNormalScheme V L) d.2 x.1 (normMap_hdim (f ≫ g) x.1) hq hwU hDU]
  have hmem : (normMap V g).base x.1.1 ∈ (normMap V f) ⁻¹ᵁ U := by
    show (normMap V f).base ((normMap V g).base x.1.1) ∈ U
    rw [← hpt]
    exact hwU
  have hffq : ffMap (normMap V f) q ≠ 0 := by
    intro h
    exact hq ((ffMap (normMap V f)).hom.injective (by simpa using h))
  have hloc : ∀ v : PrimeDivisorPt (normObj V M), (v : normObj V M) ∈ (normMap V f) ⁻¹ᵁ U →
      toWeilOnDL V DK M ((phiPull V DK f d : cartierOnDL V DK M _)
          : (DLSet V DK M) →₀ ℤ) v
        = ordPt (normObj V M) (normObj_isNormalScheme V M) v (ffMap (normMap V f) q) := by
    intro v hv
    rw [toWeilOnDL_phiPull V DK f d v]
    exact pullCoeff_eq (normMap V f) (normObj_isNormalScheme V M)
      (normObj_isNormalScheme V L) d.2 v (normMap_hdim f v) hq hv hDU
  rw [pullCoeff_eq (normMap V g) (normObj_isNormalScheme V N)
    (normObj_isNormalScheme V M) (phiPull V DK f d).2 x.1 (normMap_hdim g x.1)
    hffq hmem hloc]
  have hcomp : ffMap (normMap V (f ≫ g)) q
      = ffMap (normMap V g) (ffMap (normMap V f) q) := by
    have h1 := normFF_comp V f g
    rw [normFF_eq_ffMap, normFF_eq_ffMap, normFF_eq_ffMap] at h1
    rw [h1]
    rfl
  rw [hcomp]

/-- ★★★**`pull_nonneg`** —— 有効因子の引き戻しは有効。 -/
theorem phiPull_nonneg (DK : Set (PrimeDivisorPt V))
    {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    [IsLocallyNoetherian (normObj V L)] [IsLocallyNoetherian (normObj V M)]
    [CompactSpace (normObj V M)]
    (d : cartierOnDL V DK L (normObj_isNormalScheme V L))
    (hd : (0 : (DLSet V DK L) →₀ ℤ) ≤ (d : (DLSet V DK L) →₀ ℤ)) :
    (0 : (DLSet V DK M) →₀ ℤ)
      ≤ ((phiPull V DK f d : cartierOnDL V DK M _) : (DLSet V DK M) →₀ ℤ) := by
  refine Finsupp.le_def.mpr fun x => ?_
  show 0 ≤ pullCoeff (normMap V f) (normObj_isNormalScheme V M)
    (normObj_isNormalScheme V L) d.2 x.1
  refine pullCoeff_nonneg (normMap V f) (normObj_isNormalScheme V M)
    (normObj_isNormalScheme V L) d.2 x.1 (normMap_hdim f x.1) (fun v => ?_)
  show 0 ≤ Finsupp.embDomain (DLEmb V DK L) (d : (DLSet V DK L) →₀ ℤ) v
  by_cases hv : v ∈ DLSet V DK L
  · rw [show v = DLEmb V DK L ⟨v, hv⟩ from rfl, Finsupp.embDomain_apply_self]
    exact Finsupp.le_def.mp hd ⟨v, hv⟩
  · rw [Finsupp.embDomain_notin_range]
    rintro ⟨y, rfl⟩
    exact hv y.2

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の `Φ(L)` の関手性(台が `D_L` に入ること)。 -/
def pullCoeff_eq_zero_of_notMem_DL.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 引き戻しの台は D_M に入る",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
