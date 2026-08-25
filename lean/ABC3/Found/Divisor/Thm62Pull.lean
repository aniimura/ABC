/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormPhiFunctor
import ABC3.Found.Divisor.NormMapBase
import ABC3.Found.Divisor.Ex61Model

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

/-! ## ★1.5. 射が等しければ引き戻しも等しい -/

/-- ★**射が等しければ関数体の射も等しい**(`IsDominant` は `Prop` なので証明無関係)。 -/
theorem ffMap_congr {X Y : Scheme.{u}} {g g' : X ⟶ Y} [IsDominant g] [IsDominant g']
    [IrreducibleSpace X] [IrreducibleSpace Y] (h : g = g') : ffMap g = ffMap g' := by
  subst h; rfl

/-- ★**射が等しければ引き戻しの係数も等しい**。 -/
theorem pullCoeff_congr_hom {X Y : Scheme.{u}} {g g' : X ⟶ Y} [IsDominant g] [IsDominant g']
    [IsIntegral X] [IsIntegral Y] [IsLocallyNoetherian X] [IsLocallyNoetherian Y]
    (hnormX : IsNormalScheme X) (hnormY : IsNormalScheme Y)
    {D : WeilDiv Y} (hD : IsCartierDiv hnormY D) (h : g = g') (w : PrimeDivisorPt X) :
    pullCoeff g hnormX hnormY hD w = pullCoeff g' hnormX hnormY hD w := by
  subst h; rfl

/-- ★★★★★**可換な四角形では引き戻しが一致する** ——
`thm62-i-pull` の `Div` の**自然性**の中身。

★`pullCoeff_comp`(引き戻しの関手性)を 2 回使い、四角形の可換性で繋ぐ。 -/
theorem pullCoeff_square {X Y Z W : Scheme.{u}}
    (a : X ⟶ Y) (b : Y ⟶ W) (c : X ⟶ Z) (e : Z ⟶ W)
    [IsDominant a] [IsDominant b] [IsDominant c] [IsDominant e]
    [IsDominant (a ≫ b)] [IsDominant (c ≫ e)]
    [IsIntegral X] [IsIntegral Y] [IsIntegral Z] [IsIntegral W]
    [IsLocallyNoetherian X] [IsLocallyNoetherian Y] [IsLocallyNoetherian Z]
    [IsLocallyNoetherian W]
    (hnX : IsNormalScheme X) (hnY : IsNormalScheme Y) (hnZ : IsNormalScheme Z)
    (hnW : IsNormalScheme W) {D : WeilDiv W} (hD : IsCartierDiv hnW D)
    [CompactSpace Y] [CompactSpace Z]
    (hdimb : ∀ v : PrimeDivisorPt Y, ringKrullDim (W.presheaf.stalk (b.base v.1)) ≤ 1)
    (hdime : ∀ v : PrimeDivisorPt Z, ringKrullDim (W.presheaf.stalk (e.base v.1)) ≤ 1)
    (x : PrimeDivisorPt X)
    (hdima : ringKrullDim (Y.presheaf.stalk (a.base x.1)) ≤ 1)
    (hdimc : ringKrullDim (Z.presheaf.stalk (c.base x.1)) ≤ 1)
    (hdimab : ringKrullDim (W.presheaf.stalk ((a ≫ b).base x.1)) ≤ 1)
    (hdimce : ringKrullDim (W.presheaf.stalk ((c ≫ e).base x.1)) ≤ 1)
    (hsq : a ≫ b = c ≫ e) :
    pullCoeff a hnX hnY (isCartierDiv_cartierPullback b hnY hnW hD hdimb) x
      = pullCoeff c hnX hnZ (isCartierDiv_cartierPullback e hnZ hnW hD hdime) x := by
  rw [← pullCoeff_comp a b hnX hnY hnW hD hdimb x hdima hdimab,
      ← pullCoeff_comp c e hnX hnZ hnW hD hdime x hdimc hdimce,
      pullCoeff_congr_hom hnX hnW hD hsq x]

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

/-! ## ★3. 底が動くときの `B(L)` の引き戻し -/

/-- ★★★**`S₂` の外では `ord = 0` が保たれる**(一般形)——
`B(L) → B(L₂)` の中身。 -/
theorem ordPt_ffMap_eq_zero_of_notMem {W₁ W₂ : Scheme.{u}} (π : W₂ ⟶ W₁) [IsDominant π]
    [IsIntegral W₁] [IsIntegral W₂] [IsLocallyNoetherian W₁] [IsLocallyNoetherian W₂]
    (hn₂ : IsNormalScheme W₂) (hn₁ : IsNormalScheme W₁)
    (S₁ : Set (PrimeDivisorPt W₁)) (S₂ : Set (PrimeDivisorPt W₂))
    (hdim : ∀ x : PrimeDivisorPt W₂, ringKrullDim (W₁.presheaf.stalk (π.base x.1)) ≤ 1)
    (hpull : ∀ x : PrimeDivisorPt W₂, x ∉ S₂ → ∀ hc : IsCodimOnePt W₁ (π.base x.1),
      (⟨π.base x.1, hc⟩ : PrimeDivisorPt W₁) ∉ S₁)
    (u : (W₁.functionField)ˣ)
    (hu : ∀ z : PrimeDivisorPt W₁, z ∉ S₁ → ordPt W₁ hn₁ z (u : W₁.functionField) = 0)
    (x : PrimeDivisorPt W₂) (hx : x ∉ S₂) :
    ordPt W₂ hn₂ x (ffMap π (u : W₁.functionField)) = 0 := by
  obtain ⟨t, ht⟩ : ∃ t : (W₁.presheaf.stalk (π.base x.1))ˣ,
      algebraMap (W₁.presheaf.stalk (π.base x.1)) W₁.functionField
        (t : W₁.presheaf.stalk (π.base x.1)) = (u : W₁.functionField) := by
    by_cases hcod : IsCodimOnePt W₁ (π.base x.1)
    · exact exists_unit_of_ordPt_eq_zero hn₁ _ u.ne_zero (hu _ (hpull x hx hcod))
    · exact exists_unit_of_ringKrullDim_le_zero
        (withBot_le_zero_of_le_one_of_ne (hdim x) hcod) u
  rw [← ht]
  exact ordPt_ffMap_eq_zero_of_isUnit π hn₂ x t

/-- ★★★★★★**`B₁(L) → B₂(L₂)`**(底が動く版)——
`Theorem 6.2, (i)` の `B₁ → B₂|𝒟₁` の中身。 -/
noncomputable def bPullBase :
    BSubgroup V₁ DK₁ L (normObj_isNormalScheme V₁ L)
      →* BSubgroup V₂ DK₂ FL (normObj_isNormalScheme V₂ FL) where
  toFun u := ⟨Units.map (ffMap π).hom.toMonoidHom (u : ((normObj V₁ L).functionField)ˣ), by
    intro x hx
    exact ordPt_ffMap_eq_zero_of_notMem π (normObj_isNormalScheme V₂ FL)
      (normObj_isNormalScheme V₁ L) (DLSet V₁ DK₁ L) (DLSet V₂ DK₂ FL) hdim hpull
      (u : ((normObj V₁ L).functionField)ˣ) u.2 x hx⟩
  map_one' := Subtype.ext (map_one _)
  map_mul' a b := Subtype.ext (map_mul _ _ _)

/-- ★★★★**可換な四角形では関数体の引き戻しが一致する** ——
`B` の**自然性**の中身。 -/
theorem ffMap_square {X Y Z W : Scheme.{u}}
    (a : X ⟶ Y) (b : Y ⟶ W) (c : X ⟶ Z) (e : Z ⟶ W)
    [IsDominant a] [IsDominant b] [IsDominant c] [IsDominant e]
    [IsDominant (a ≫ b)] [IsDominant (c ≫ e)]
    [IrreducibleSpace X] [IrreducibleSpace Y] [IrreducibleSpace Z] [IrreducibleSpace W]
    (hsq : a ≫ b = c ≫ e) (t : W.functionField) :
    ffMap a (ffMap b t) = ffMap c (ffMap e t) := by
  have h1 : ffMap (a ≫ b) t = ffMap a (ffMap b t) := by rw [ffMap_comp]; rfl
  have h2 : ffMap (c ≫ e) t = ffMap c (ffMap e t) := by rw [ffMap_comp]; rfl
  rw [← h1, ← h2, ffMap_congr hsq]

/-! ## ★3.5. 有効因子は有効因子へ —— `Φ` は `effSub` なので必要 -/

omit [IsLocallyNoetherian (normObj V₁ L)] in
/-- ★**非負の `Finsupp` を `D_L` 添字から広げても非負**。 -/
theorem toWeilOnDL_nonneg (d : (DLSet V₁ DK₁ L) →₀ ℤ) (hd : 0 ≤ d)
    (v : PrimeDivisorPt (normObj V₁ L)) :
    0 ≤ toWeilOnDL V₁ DK₁ L d v := by
  show 0 ≤ Finsupp.embDomain (DLEmb V₁ DK₁ L) d v
  by_cases h : ∃ y, DLEmb V₁ DK₁ L y = v
  · obtain ⟨y, rfl⟩ := h
    rw [Finsupp.embDomain_apply]
    simpa using (Finsupp.le_def.mp hd) y
  · rw [Finsupp.embDomain_notin_range _ _ _ (by rintro ⟨y, hy⟩; exact h ⟨y, hy⟩)]

/-- ★★★**有効因子の引き戻しは有効** ——
`Φ(L) = effSub (Φ(L)^gp)` なので `phiHom` を作るのに要る。 -/
theorem phiPullBase_nonneg (d : cartierOnDL V₁ DK₁ L (normObj_isNormalScheme V₁ L))
    (hd : 0 ≤ (d : (DLSet V₁ DK₁ L) →₀ ℤ)) :
    0 ≤ ((phiPullBase DK₁ DK₂ L FL π hdim hpull d
      : cartierOnDL V₂ DK₂ FL (normObj_isNormalScheme V₂ FL))
        : (DLSet V₂ DK₂ FL) →₀ ℤ) := by
  refine Finsupp.le_def.mpr fun x => ?_
  have hval : ((phiPullBase DK₁ DK₂ L FL π hdim hpull d
      : cartierOnDL V₂ DK₂ FL (normObj_isNormalScheme V₂ FL))
        : (DLSet V₂ DK₂ FL) →₀ ℤ) x
      = pullCoeff π (normObj_isNormalScheme V₂ FL)
        (normObj_isNormalScheme V₁ L) d.2 x.1 := rfl
  rw [show ((0 : (DLSet V₂ DK₂ FL) →₀ ℤ)) x = 0 from rfl, hval]
  exact pullCoeff_nonneg π (normObj_isNormalScheme V₂ FL) (normObj_isNormalScheme V₁ L)
    d.2 x.1 (hdim x.1) (fun v => toWeilOnDL_nonneg DK₁ L (d : (DLSet V₁ DK₁ L) →₀ ℤ) hd v)

/-- ★★**有効な部分へ制限する**(一般形)—— `Φ = effSub (Φ^gp)` なので要る。 -/
noncomputable def restrictEff {S T : Type u} {Γ₁ : AddSubgroup (S →₀ ℤ)}
    {Γ₂ : AddSubgroup (T →₀ ℤ)} (φ : Γ₁ →+ Γ₂)
    (hnn : ∀ a : Γ₁, 0 ≤ (a : S →₀ ℤ) → 0 ≤ ((φ a : Γ₂) : T →₀ ℤ)) :
    effSub Γ₁ →+ effSub Γ₂ where
  toFun a := ⟨((φ ⟨a.1, a.2.1⟩ : Γ₂) : T →₀ ℤ), (φ ⟨a.1, a.2.1⟩).2, hnn _ a.2.2⟩
  map_zero' := by
    refine Subtype.ext ?_
    have h0 : (⟨(0 : effSub Γ₁).1, (0 : effSub Γ₁).2.1⟩ : Γ₁) = 0 := Subtype.ext rfl
    show ((φ ⟨(0 : effSub Γ₁).1, (0 : effSub Γ₁).2.1⟩ : Γ₂) : T →₀ ℤ) = 0
    rw [h0, map_zero]
    rfl
  map_add' a b := by
    refine Subtype.ext ?_
    have hab : (⟨(a + b).1, (a + b).2.1⟩ : Γ₁) = ⟨a.1, a.2.1⟩ + ⟨b.1, b.2.1⟩ := Subtype.ext rfl
    show ((φ ⟨(a + b).1, (a + b).2.1⟩ : Γ₂) : T →₀ ℤ) = _
    rw [hab, map_add]
    rfl

/-- ★★★★★**`Φ₁(L) → Φ₂(L₂)`**(有効因子の側)——
`ModelDataHomOver.phiHom` がそのまま要求する形。 -/
noncomputable def phiPullEff :
    effSub (cartierOnDL V₁ DK₁ L (normObj_isNormalScheme V₁ L))
      →+ effSub (cartierOnDL V₂ DK₂ FL (normObj_isNormalScheme V₂ FL)) :=
  restrictEff (phiPullBase DK₁ DK₂ L FL π hdim hpull)
    (fun a ha => phiPullBase_nonneg DK₁ DK₂ L FL π hdim hpull a ha)

/-! ## ★3.9. `π_L` を作るための可換性

原文 (FrdI p.110):
> in CV,K,DK may be thought of as consisting of the following data: (a) a morphism
-/

/-- ★★★★★**`normMapOfBase` に渡す可換性** ——
体の埋め込み `φ : L ↪ L₂` が `K₁ → K₂` と両立すれば
`Spec L₂ ⟶ V₂ ⟶ V₁` と `Spec L₂ ⟶ Spec L ⟶ V₁` が一致する。

★★これで `normMapOfBase` が `V₂[L₂] ⟶ V₁[L]` を与える。
★中身は `Scheme.SpecMap_stalkMap_fromSpecStalk`(生成点からの射の自然性)と
`Scheme.SpecMap_stalkSpecializes_fromSpecStalk`(支配射は生成点を生成点へ送る)。 -/
theorem specToV_compat {V₁ V₂ : Scheme.{u}} [IsIntegral V₁] [IsIntegral V₂]
    (ψ : V₂ ⟶ V₁) [IsDominant ψ]
    {Kbar₁ Kbar₂ : Type u} [Field Kbar₁] [Field Kbar₂]
    [Algebra V₁.functionField Kbar₁] [Algebra V₂.functionField Kbar₂]
    (L : FinSub V₁.functionField Kbar₁) (FL : FinSub V₂.functionField Kbar₂)
    (φ : CommRingCat.of L.toIF ⟶ CommRingCat.of FL.toIF)
    (hφ : CommRingCat.ofHom (algebraMap V₁.functionField L.toIF) ≫ φ
        = ffMap ψ ≫ CommRingCat.ofHom (algebraMap V₂.functionField FL.toIF)) :
    specToV V₂ FL ≫ ψ = Spec.map φ ≫ specToV V₁ L := by
  have hgen : ψ.base (genericPoint (V₂ : Type u)) = genericPoint (V₁ : Type u) :=
    genericPoint_eq_of_isDominant ψ
  have hsp : ψ.base (genericPoint (V₂ : Type u)) ⤳ genericPoint (V₁ : Type u) :=
    (Inseparable.of_eq hgen).specializes
  have hfs : Spec.map (V₁.presheaf.stalkSpecializes hsp) ≫ V₁.fromSpecStalk
      (genericPoint (V₁ : Type u)) = V₁.fromSpecStalk (ψ.base (genericPoint (V₂ : Type u))) :=
    Scheme.SpecMap_stalkSpecializes_fromSpecStalk hsp
  have hff : V₁.presheaf.stalkSpecializes hsp ≫ ψ.stalkMap (genericPoint (V₂ : Type u))
      = ffMap ψ := rfl
  rw [specToV, specToV, Category.assoc, ← Scheme.SpecMap_stalkMap_fromSpecStalk ψ,
    ← hfs, ← Category.assoc, ← Spec.map_comp, ← Category.assoc, ← Spec.map_comp,
    ← Category.assoc (V₁.presheaf.stalkSpecializes hsp), hff,
    ← Category.assoc (Spec.map φ), ← Spec.map_comp, hφ]
  rfl

/-- ★★★★★★**幾何のデータの射が `V₂[L₂] ⟶ V₁[L]` を与える** ——
`Theorem 6.2, (i)` の因子・有理函数の引き戻しの土台。 -/
noncomputable def normMapGeom {V₁ V₂ : Scheme.{u}} [IsIntegral V₁] [IsIntegral V₂]
    (ψ : V₂ ⟶ V₁) [IsDominant ψ]
    {Kbar₁ Kbar₂ : Type u} [Field Kbar₁] [Field Kbar₂]
    [Algebra V₁.functionField Kbar₁] [Algebra V₂.functionField Kbar₂]
    (L : FinSub V₁.functionField Kbar₁) (FL : FinSub V₂.functionField Kbar₂)
    (φ : CommRingCat.of L.toIF ⟶ CommRingCat.of FL.toIF)
    (hφ : CommRingCat.ofHom (algebraMap V₁.functionField L.toIF) ≫ φ
        = ffMap ψ ≫ CommRingCat.ofHom (algebraMap V₂.functionField FL.toIF)) :
    normObj V₂ FL ⟶ normObj V₁ L :=
  normMapOfBase (specToV V₁ L) (specToV V₂ FL) ψ (Spec.map φ)
    (specToV_compat ψ L FL φ hφ)

@[reassoc] theorem normMapGeom_normDown {V₁ V₂ : Scheme.{u}} [IsIntegral V₁] [IsIntegral V₂]
    (ψ : V₂ ⟶ V₁) [IsDominant ψ]
    {Kbar₁ Kbar₂ : Type u} [Field Kbar₁] [Field Kbar₂]
    [Algebra V₁.functionField Kbar₁] [Algebra V₂.functionField Kbar₂]
    (L : FinSub V₁.functionField Kbar₁) (FL : FinSub V₂.functionField Kbar₂)
    (φ : CommRingCat.of L.toIF ⟶ CommRingCat.of FL.toIF)
    (hφ : CommRingCat.ofHom (algebraMap V₁.functionField L.toIF) ≫ φ
        = ffMap ψ ≫ CommRingCat.ofHom (algebraMap V₂.functionField FL.toIF)) :
    normMapGeom ψ L FL φ hφ ≫ normDown V₁ L = normDown V₂ FL ≫ ψ :=
  normMapOfBase_fromNormalization _ _ _ _ _

/-! ## ★4. `divCompat` の中身 —— 主因子の引き戻しは引き戻しの主因子 -/

/-- ★★★★★★**主因子の引き戻しは引き戻しの主因子** ——
`div(g^*f) = g^*(div f)`。

★★これが `ModelDataHomOver` の `divCompat`(`Div_B` との両立)の中身である。
★中身は 1 行 —— **主因子は `⊤` の上で `f` 自身を局所方程式に持つ**ので、
`pullCoeff_eq` をそのまま当てればよい。 -/
theorem pullCoeff_weilDivOfFn {W₁ W₂ : Scheme.{u}} (π : W₂ ⟶ W₁) [IsDominant π]
    [IsIntegral W₁] [IsIntegral W₂] [IsLocallyNoetherian W₁] [IsLocallyNoetherian W₂]
    [CompactSpace W₁]
    (hn₂ : IsNormalScheme W₂) (hn₁ : IsNormalScheme W₁)
    {f : W₁.functionField} (hf : f ≠ 0)
    (x : PrimeDivisorPt W₂)
    (hdimx : ringKrullDim (W₁.presheaf.stalk (π.base x.1)) ≤ 1) :
    pullCoeff π hn₂ hn₁ (isCartierDiv_weilDivOfFn hn₁ hf) x = ordPt W₂ hn₂ x (ffMap π f) :=
  pullCoeff_eq π hn₂ hn₁ _ x hdimx hf (U := ⊤) trivial (fun _ _ => rfl)

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

def pullCoeff_square.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 可換な四角形では引き戻しが一致する(Div の自然性)",
    sectionId := "frdi-thm-6-2" }

def pullCoeff_square.needs : List ProofObligation :=
  [ .citation "[ABC3]" "pullCoeff_comp(引き戻しの関手性)"
      (.inProject "ABC3" "ABC3.Found.Divisor.pullCoeff_comp") 110,
    .citation "[ABC3]" "normMapOfBase_naturality(四角形の可換性を与える側)"
      (.inProject "ABC3" "ABC3.Found.Divisor.normMapOfBase_naturality") 110,
    .derivation "関手性を 2 回使い、四角形の可換性で繋ぐ。IsDominant は Prop なので射の等式で移送できる" 110 ]

def bPullBase.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — B₁ → B₂|𝒟₁(底が動く版の有理函数の引き戻し)",
    sectionId := "frdi-thm-6-2" }

def bPullBase.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ffMap(支配射に沿った関数体の射)"
      (.inProject "ABC3" "ABC3.Found.Divisor.ffMap") 110,
    .citation "[ABC3]" "exists_unit_of_ordPt_eq_zero / exists_unit_of_ringKrullDim_le_zero"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_unit_of_ordPt_eq_zero") 110,
    .derivation
      "D_{L₂} の外の点は D_L の外へ行くので、そこでは u が茎の単元になり ord = 0 が保たれる" 110,
    .implicitStep
      "★原文は仮定 (a) から有理函数の引き戻しを「(i)」の一言で置いている" 110 ]

def ffMap_square.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 可換な四角形では関数体の引き戻しが一致する(B の自然性)",
    sectionId := "frdi-thm-6-2" }

def ffMap_square.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ffMap_comp"
      (.inProject "ABC3" "ABC3.Found.Divisor.ffMap_comp") 110 ]

def pullCoeff_weilDivOfFn.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 主因子の引き戻しは引き戻しの主因子(divCompat)",
    sectionId := "frdi-thm-6-2" }

def pullCoeff_weilDivOfFn.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isCartierDiv_weilDivOfFn(主因子は ⊤ の上で f 自身を局所方程式に持つ)"
      (.inProject "ABC3" "ABC3.Found.Divisor.isCartierDiv_weilDivOfFn") 110,
    .citation "[ABC3]" "pullCoeff_eq"
      (.inProject "ABC3" "ABC3.Found.Divisor.pullCoeff_eq") 110 ]


def specToV_compat.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 体の埋め込みが Spec L 上の四角形を可換にする",
    sectionId := "frdi-thm-6-2" }

def specToV_compat.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Scheme.SpecMap_stalkMap_fromSpecStalk"
      (.inMathlib "AlgebraicGeometry.Scheme.SpecMap_stalkMap_fromSpecStalk") 110,
    .citation "[mathlib]" "Scheme.SpecMap_stalkSpecializes_fromSpecStalk"
      (.inMathlib "AlgebraicGeometry.Scheme.SpecMap_stalkSpecializes_fromSpecStalk") 110,
    .derivation "支配射は生成点を生成点へ送るので、生成点からの射の四角形が閉じる" 110 ]

def normMapGeom.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 幾何のデータの射が V₂[L₂] ⟶ V₁[L] を与える",
    sectionId := "frdi-thm-6-2" }

def normMapGeom.needs : List ProofObligation :=
  [ .citation "[ABC3]" "normMapOfBase(底変換 ＋ normalizationDesc)"
      (.inProject "ABC3" "ABC3.Found.Divisor.normMapOfBase") 110,
    .citation "[ABC3]" "specToV_compat"
      (.inProject "ABC3" "ABC3.Found.Divisor.specToV_compat") 110 ]


/-- ★★**`normMapGeom` は支配的** —— `normUp` が支配的で、それを経由するから。 -/
instance normMapGeom_isDominant {V₁ V₂ : Scheme.{u}} [IsIntegral V₁] [IsIntegral V₂]
    (ψ : V₂ ⟶ V₁) [IsDominant ψ]
    {Kbar₁ Kbar₂ : Type u} [Field Kbar₁] [Field Kbar₂]
    [Algebra V₁.functionField Kbar₁] [Algebra V₂.functionField Kbar₂]
    (L : FinSub V₁.functionField Kbar₁) (FL : FinSub V₂.functionField Kbar₂)
    (φ : CommRingCat.of L.toIF ⟶ CommRingCat.of FL.toIF)
    (hφ : CommRingCat.ofHom (algebraMap V₁.functionField L.toIF) ≫ φ
        = ffMap ψ ≫ CommRingCat.ofHom (algebraMap V₂.functionField FL.toIF)) :
    IsDominant (normMapGeom ψ L FL φ hφ) := by
  haveI : Surjective (Spec.map φ) := ⟨fun _ => ⟨default, Subsingleton.elim _ _⟩⟩
  have heq : normUp V₂ FL ≫ normMapGeom ψ L FL φ hφ = Spec.map φ ≫ normUp V₁ L :=
    toNormalization_normMapOfBase (specToV V₁ L) (specToV V₂ FL) ψ (Spec.map φ) _
  have h : IsDominant (normUp V₂ FL ≫ normMapGeom ψ L FL φ hφ) := by
    rw [heq]; infer_instance
  exact IsDominant.of_comp (normUp V₂ FL) _

end ABC3.Found.Divisor
