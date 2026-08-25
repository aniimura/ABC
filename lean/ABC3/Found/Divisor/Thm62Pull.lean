/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormPhiFunctor
import ABC3.Found.Divisor.NormMapBase

/-!
# [FrdI] Theorem 6.2, (i) —— 底が動くときの `Φ(L)^gp` の引き戻し

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> in CV,K,DK may be thought of as consisting of the following data: (a) a morphism

## ★★何をするか

`NormPhiFunctor.lean` の `phiPull` は**同じ `V` の中**の
`V[M] ⟶ V[L]`(`normMap`)に沿った引き戻しである。
★`Theorem 6.2, (i)` が要求するのは**底が動く**版

```
V₂[L₂] ⟶ V₁[L]      (normMapOfBase、2026-08-25 に閉じた)
```

に沿った `Φ₁(L)^gp ⟶ Φ₂(L₂)^gp` であり、本ファイルはそれを置く。

## ★★★入力(原文の仮定 (a) に対応)

| 仮引数 | 原文 |
|---|---|
| `π` | `ψ : V₂ → V₁` が誘導する `V₂[L₂] ⟶ V₁[L]` |
| `hdim` | 余次元 1 の点が余次元 ≤ 1 の点へ行くこと(`codim-preserved`) |
| `hpull` | **`D_{K₁}` の引き戻しが `D_{K₂}` に入る**こと(原文の仮定 (a)) |

★★`hpull` の対偶を取ると「`D_{L₂}` の外の点は `D_L` の外へ行く」で、
それがちょうど引き戻しの台が `D_{L₂}` に収まる理由になる。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI ABC3.Meta

universe u

/-! ## ★1. 台の条件を抽象化した引き戻しの消滅 -/

/-- ★★★**引き戻しの台は `S₂` に入る**(一般形)——
`D` の台が `S₁` に入り、`w ∉ S₂` の像が `S₁` の外なら係数は `0`。

★`NormPhiFunctor.lean` の `pullCoeff_eq_zero_of_notMem_DL` は
`normMap` に固定された形なので、底が動く場合のために抽象化する。 -/
theorem pullCoeff_eq_zero_of_notMem {W₁ W₂ : Scheme.{u}} (π : W₂ ⟶ W₁) [IsDominant π]
    [IsIntegral W₁] [IsIntegral W₂] [IsLocallyNoetherian W₁] [IsLocallyNoetherian W₂]
    (hnorm₂ : IsNormalScheme W₂) (hnorm₁ : IsNormalScheme W₁)
    {D : WeilDiv W₁} (hD : IsCartierDiv hnorm₁ D)
    (S₁ : Set (PrimeDivisorPt W₁))
    (hsupp : ∀ z : PrimeDivisorPt W₁, z ∉ S₁ → D z = 0)
    (w : PrimeDivisorPt W₂)
    (hdim : ringKrullDim (W₁.presheaf.stalk (π.base w.1)) ≤ 1)
    (hpull : ∀ hc : IsCodimOnePt W₁ (π.base w.1),
      (⟨π.base w.1, hc⟩ : PrimeDivisorPt W₁) ∉ S₁) :
    pullCoeff π hnorm₂ hnorm₁ hD w = 0 := by
  obtain ⟨U, hwU, q, hq, hDU⟩ := hD (π.base w.1)
  rw [pullCoeff_eq π hnorm₂ hnorm₁ hD w hdim hq hwU hDU]
  obtain ⟨t, ht⟩ := exists_unit_of_ordPt_of_dim_le_one hnorm₁ (π.base w.1) hdim
    (Units.mk0 q hq) (by
      intro hc
      rw [show ((Units.mk0 q hq : (W₁.functionField)ˣ) : W₁.functionField) = q from rfl,
        ← hDU ⟨π.base w.1, hc⟩ hwU]
      exact hsupp _ (hpull hc))
  rw [show q = ((Units.mk0 q hq : (W₁.functionField)ˣ) : W₁.functionField) from rfl, ← ht]
  exact ordPt_ffMap_eq_zero_of_isUnit π hnorm₂ w t

/-! ## ★2. 底が動くときの `Φ(L)^gp` の引き戻し -/

variable {V₁ V₂ : Scheme.{u}} [IsIntegral V₁] [IsIntegral V₂]
  {Kbar₁ Kbar₂ : Type u} [Field Kbar₁] [Field Kbar₂]
  [Algebra V₁.functionField Kbar₁] [Algebra V₂.functionField Kbar₂]

variable (DK₁ : Set (PrimeDivisorPt V₁)) (DK₂ : Set (PrimeDivisorPt V₂))
  (L : FinSub V₁.functionField Kbar₁) (FL : FinSub V₂.functionField Kbar₂)
  [IsLocallyNoetherian (normObj V₁ L)] [IsLocallyNoetherian (normObj V₂ FL)]
  [CompactSpace (normObj V₂ FL)]
  (π : normObj V₂ FL ⟶ normObj V₁ L) [IsDominant π]
  (hdim : ∀ w : PrimeDivisorPt (normObj V₂ FL),
    ringKrullDim ((normObj V₁ L).presheaf.stalk (π.base w.1)) ≤ 1)
  (hpull : ∀ w : PrimeDivisorPt (normObj V₂ FL), w ∉ DLSet V₂ DK₂ FL →
    ∀ hc : IsCodimOnePt (normObj V₁ L) (π.base w.1),
      (⟨π.base w.1, hc⟩ : PrimeDivisorPt (normObj V₁ L)) ∉ DLSet V₁ DK₁ L)

/-- ★`d ∈ Φ(L)^gp` の台は `D_L` に入る。 -/
theorem toWeilOnDL_eq_zero_of_notMem
    (d : cartierOnDL V₁ DK₁ L (normObj_isNormalScheme V₁ L))
    (z : PrimeDivisorPt (normObj V₁ L)) (hz : z ∉ DLSet V₁ DK₁ L) :
    toWeilOnDL V₁ DK₁ L (d : (DLSet V₁ DK₁ L) →₀ ℤ) z = 0 := by
  show Finsupp.embDomain (DLEmb V₁ DK₁ L) (d : (DLSet V₁ DK₁ L) →₀ ℤ) z = 0
  refine Finsupp.embDomain_notin_range _ _ _ ?_
  rintro ⟨y, rfl⟩
  exact hz y.2

open scoped Classical in
/-- ★★★★★★**`Φ₁(L)^gp → Φ₂(L₂)^gp`**(底が動く版)——
`Theorem 6.2, (i)` の `Φ₁ → Φ₂|𝒟₁` の中身。

★`NormPhiFunctor.lean` の `phiPull` を `normMap` から一般の支配射 `π` へ移したもの。 -/
noncomputable def phiPullBase :
    cartierOnDL V₁ DK₁ L (normObj_isNormalScheme V₁ L)
      →+ cartierOnDL V₂ DK₂ FL (normObj_isNormalScheme V₂ FL) where
  toFun d := ⟨Finsupp.subtypeDomain (· ∈ DLSet V₂ DK₂ FL)
      (cartierPullback π (normObj_isNormalScheme V₂ FL)
        (normObj_isNormalScheme V₁ L) d.2 hdim), by
    show toWeilOnDL V₂ DK₂ FL _ ∈ cartierSubgroup (normObj_isNormalScheme V₂ FL)
    rw [embDomain_subtypeDomain V₂ DK₂ FL _ (fun w hw => by
      rw [cartierPullback_apply]
      exact pullCoeff_eq_zero_of_notMem π (normObj_isNormalScheme V₂ FL)
        (normObj_isNormalScheme V₁ L) d.2 (DLSet V₁ DK₁ L)
        (toWeilOnDL_eq_zero_of_notMem DK₁ L d) w (hdim w) (hpull w hw))]
    exact isCartierDiv_cartierPullback π (normObj_isNormalScheme V₂ FL)
      (normObj_isNormalScheme V₁ L) d.2 hdim⟩
  map_zero' := by
    refine Subtype.ext (Finsupp.ext fun x => ?_)
    show pullCoeff π (normObj_isNormalScheme V₂ FL)
      (normObj_isNormalScheme V₁ L) (0 : cartierOnDL V₁ DK₁ L _).2 x.1 = 0
    rw [pullCoeff_eq π (normObj_isNormalScheme V₂ FL) (normObj_isNormalScheme V₁ L)
      (0 : cartierOnDL V₁ DK₁ L _).2 x.1 (hdim x.1) (U := ⊤) (f := 1) (one_ne_zero)
      (Set.mem_univ _) (fun v _ => by simp [toWeilOnDL, ordPt_one])]
    simp [ordPt_one]
  map_add' d e := by
    refine Subtype.ext (Finsupp.ext fun x => ?_)
    show pullCoeff π (normObj_isNormalScheme V₂ FL)
        (normObj_isNormalScheme V₁ L) (d + e).2 x.1
      = pullCoeff π (normObj_isNormalScheme V₂ FL)
          (normObj_isNormalScheme V₁ L) d.2 x.1
        + pullCoeff π (normObj_isNormalScheme V₂ FL)
          (normObj_isNormalScheme V₁ L) e.2 x.1
    have hDeq : toWeilOnDL V₁ DK₁ L
        ((d : (DLSet V₁ DK₁ L) →₀ ℤ) + (e : (DLSet V₁ DK₁ L) →₀ ℤ))
        = toWeilOnDL V₁ DK₁ L (d : (DLSet V₁ DK₁ L) →₀ ℤ)
          + toWeilOnDL V₁ DK₁ L (e : (DLSet V₁ DK₁ L) →₀ ℤ) := map_add _ _ _
    have hde : IsCartierDiv (normObj_isNormalScheme V₁ L)
        (toWeilOnDL V₁ DK₁ L (d : (DLSet V₁ DK₁ L) →₀ ℤ)
          + toWeilOnDL V₁ DK₁ L (e : (DLSet V₁ DK₁ L) →₀ ℤ)) := hDeq ▸ (d + e).2
    rw [pullCoeff_congr π (normObj_isNormalScheme V₂ FL)
      (normObj_isNormalScheme V₁ L) (d + e).2 hde hDeq x.1]
    exact pullCoeff_add π (normObj_isNormalScheme V₂ FL)
      (normObj_isNormalScheme V₁ L) d.2 e.2 hde x.1 (hdim x.1)

/-- ★引き戻しを Weil 因子として見ると係数そのもの。 -/
theorem toWeilOnDL_phiPullBase
    (d : cartierOnDL V₁ DK₁ L (normObj_isNormalScheme V₁ L))
    (v : PrimeDivisorPt (normObj V₂ FL)) :
    toWeilOnDL V₂ DK₂ FL
        ((phiPullBase DK₁ DK₂ L FL π hdim hpull d : cartierOnDL V₂ DK₂ FL _)
          : (DLSet V₂ DK₂ FL) →₀ ℤ) v
      = pullCoeff π (normObj_isNormalScheme V₂ FL)
        (normObj_isNormalScheme V₁ L) d.2 v := by
  have hemb := embDomain_subtypeDomain V₂ DK₂ FL
    (cartierPullback π (normObj_isNormalScheme V₂ FL)
      (normObj_isNormalScheme V₁ L) d.2 hdim)
    (fun w hw => by
      rw [cartierPullback_apply]
      exact pullCoeff_eq_zero_of_notMem π (normObj_isNormalScheme V₂ FL)
        (normObj_isNormalScheme V₁ L) d.2 (DLSet V₁ DK₁ L)
        (toWeilOnDL_eq_zero_of_notMem DK₁ L d) w (hdim w) (hpull w hw))
  rw [show ((phiPullBase DK₁ DK₂ L FL π hdim hpull d : cartierOnDL V₂ DK₂ FL _)
      : (DLSet V₂ DK₂ FL) →₀ ℤ)
      = Finsupp.subtypeDomain (· ∈ DLSet V₂ DK₂ FL)
        (cartierPullback π (normObj_isNormalScheme V₂ FL)
          (normObj_isNormalScheme V₁ L) d.2 hdim) from rfl, hemb]
  rfl

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def phiPullBase.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — Φ₁ → Φ₂|𝒟₁(底が動く版の因子の引き戻し)",
    sectionId := "frdi-thm-6-2" }

def phiPullBase.needs : List ProofObligation :=
  [ .citation "[ABC3]" "cartierPullback(Cartier 因子の引き戻し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.cartierPullback") 110,
    .citation "[ABC3]" "normMapOfBase(底が動くときの V[L] の射)"
      (.inProject "ABC3" "ABC3.Found.Divisor.normMapOfBase") 110,
    .citation "[ABC3]" "pullCoeff_comp(引き戻しの関手性、自然性の中身)"
      (.inProject "ABC3" "ABC3.Found.Divisor.pullCoeff_comp") 110,
    .derivation
      "hpull(D_{K₁} の引き戻しが D_{K₂} に入る)の対偶で、引き戻しの台が D_{L₂} に収まる" 110,
    .implicitStep
      "★原文は仮定 (a) から因子の引き戻しを「(i)」の一言で置いている" 110 ]

def pullCoeff_eq_zero_of_notMem.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 引き戻しの台は D_L に入る(一般形)",
    sectionId := "frdi-thm-6-2" }

def pullCoeff_eq_zero_of_notMem.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_unit_of_ordPt_of_dim_le_one"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_unit_of_ordPt_of_dim_le_one") 110,
    .citation "[ABC3]" "ordPt_ffMap_eq_zero_of_isUnit"
      (.inProject "ABC3" "ABC3.Found.Divisor.ordPt_ffMap_eq_zero_of_isUnit") 110 ]

end ABC3.Found.Divisor
