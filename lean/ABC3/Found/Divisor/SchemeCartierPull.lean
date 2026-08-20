/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeDivFinite
import ABC3.Found.Divisor.SchemeFFMap
import ABC3.Found.Divisor.NormBFunctor

/-!
# Cartier 因子の引き戻し —— 局所方程式の取り替えに依らないこと(鎖 `cartier` の `cartier-pullback`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> If for every Spec(L) ∈Ob(B(G)0) [cf. §0], every prime divisor of DL is Q-Cartier,

## ★★引き戻しの定義と、その要

`IsCartierDiv` は**定義に局所方程式を持っている** ——
各点 `y` に開 `U ∋ y` と `f ∈ K(Y)^×` があって `D v = ord_v f`(`v ∈ U`)。
★したがって引き戻しの係数は `ord_w(g^* f)` で定めればよい。

★★★**要は「局所方程式の取り替えに依らないこと」**である。
2 つの局所方程式 `f₁, f₂` が `g(w)` の近くで両方 `D` を与えるなら、
`u := f₁/f₂` は `ord_{g(w)}(u) = 0` を満たすので
**`𝒪_{Y,g(w)}` の単元**(`SchemeOrdUnit.lean`)、したがって
`g^*u` は `𝒪_{X,w}` の単元(`SchemeFFMap.lean` の四角形)、
すなわち `ord_w(g^*f₁) = ord_w(g^*f₂)` である。

★`g(w)` が余次元 1 でない場合は `dim 𝒪_{Y,g(w)} ≤ 1` から茎が**体**になるので、
やはり単元である(`codim-preserved` の系)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `exists_unit_of_ordPt_of_dim_le_one` | ★`ord = 0`(余次元 1 のとき)＋ `dim ≤ 1` ⟹ 茎の単元 |
| `ordPt_ffMap_congr` | ★★★**局所方程式の取り替えに依らない** |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory

universe u

/-! ## ★1. 茎の単元を取り出す(まとめた形) -/

/-- ★★**`dim 𝒪_{Y,z} ≤ 1` のとき、`z` が余次元 1 なら `ord_z = 0`、
そうでなければ無条件に、`u` は茎の単元から来る**。 -/
theorem exists_unit_of_ordPt_of_dim_le_one {Y : Scheme.{u}} [IsIntegral Y]
    [IsLocallyNoetherian Y] (hnorm : IsNormalScheme Y) (z : Y)
    (hdim : ringKrullDim (Y.presheaf.stalk z) ≤ 1) (u : (Y.functionField)ˣ)
    (hord : ∀ hc : IsCodimOnePt Y z,
      ordPt Y hnorm ⟨z, hc⟩ (u : Y.functionField) = 0) :
    ∃ t : (Y.presheaf.stalk z)ˣ,
      algebraMap (Y.presheaf.stalk z) Y.functionField (t : Y.presheaf.stalk z)
        = (u : Y.functionField) := by
  by_cases hc : IsCodimOnePt Y z
  · exact exists_unit_of_ordPt_eq_zero hnorm ⟨z, hc⟩ u.ne_zero (hord hc)
  · exact exists_unit_of_ringKrullDim_le_zero
      (withBot_le_zero_of_le_one_of_ne hdim hc) u

/-! ## ★2. 局所方程式の取り替えに依らない -/

variable {X Y : Scheme.{u}} (g : X ⟶ Y) [IsDominant g] [IsIntegral X] [IsIntegral Y]
  [IsLocallyNoetherian X] [IsLocallyNoetherian Y]

/-- ★★★★**引き戻しの係数は局所方程式の取り替えに依らない**。 -/
theorem ordPt_ffMap_congr (hnormX : IsNormalScheme X) (hnormY : IsNormalScheme Y)
    (w : PrimeDivisorPt X)
    (hdim : ringKrullDim (Y.presheaf.stalk (g.base w.1)) ≤ 1)
    {f₁ f₂ : Y.functionField} (h₁ : f₁ ≠ 0) (h₂ : f₂ ≠ 0)
    (hord : ∀ hc : IsCodimOnePt Y (g.base w.1),
      ordPt Y hnormY ⟨g.base w.1, hc⟩ f₁ = ordPt Y hnormY ⟨g.base w.1, hc⟩ f₂) :
    ordPt X hnormX w (ffMap g f₁) = ordPt X hnormX w (ffMap g f₂) := by
  set u : (Y.functionField)ˣ := (Units.mk0 f₁ h₁) * (Units.mk0 f₂ h₂)⁻¹ with hu
  have huval : (u : Y.functionField) = f₁ * f₂⁻¹ := rfl
  obtain ⟨t, ht⟩ := exists_unit_of_ordPt_of_dim_le_one hnormY (g.base w.1) hdim u (by
    intro hc
    rw [huval, ordPt_mul hnormY _ h₁ (inv_ne_zero h₂), ordPt_inv hnormY _ h₂, hord hc]
    ring)
  have hzero : ordPt X hnormX w (ffMap g (u : Y.functionField)) = 0 := by
    rw [← ht]
    exact ordPt_ffMap_eq_zero_of_isUnit g hnormX w t
  have hff2 : ffMap g f₂ ≠ 0 := by
    intro h
    exact h₂ ((ffMap g).hom.injective (by simpa using h))
  have hsplit : ffMap g (u : Y.functionField) = ffMap g f₁ * (ffMap g f₂)⁻¹ := by
    rw [huval, map_mul, map_inv₀]
  rw [hsplit, ordPt_mul hnormX w
    (by intro h; exact h₁ ((ffMap g).hom.injective (by simpa using h))) (inv_ne_zero hff2),
    ordPt_inv hnormX w hff2] at hzero
  omega

/-! ## ★3. 引き戻しの係数 -/

/-- ★局所方程式を組にして取り出す(選択のため)。 -/
theorem cartier_local (hnormY : IsNormalScheme Y) {D : WeilDiv Y}
    (hD : IsCartierDiv hnormY D) (y : Y) :
    ∃ p : Y.Opens × Y.functionField, y ∈ p.1 ∧ p.2 ≠ 0 ∧
      ∀ v : PrimeDivisorPt Y, (v : Y) ∈ p.1 → D v = ordPt Y hnormY v p.2 := by
  obtain ⟨U, hy, f, hf, hDU⟩ := hD y
  exact ⟨⟨U, f⟩, hy, hf, hDU⟩

/-- ★★**引き戻しの係数** `ord_w(g^* f)`(`f` は `g(w)` の近くの局所方程式)。 -/
noncomputable def pullCoeff (hnormX : IsNormalScheme X) (hnormY : IsNormalScheme Y)
    {D : WeilDiv Y} (hD : IsCartierDiv hnormY D) (w : PrimeDivisorPt X) : ℤ :=
  ordPt X hnormX w (ffMap g (cartier_local hnormY hD (g.base w.1)).choose.2)

/-- ★★★**係数はどの局所方程式で測っても同じ**。 -/
theorem pullCoeff_eq (hnormX : IsNormalScheme X) (hnormY : IsNormalScheme Y)
    {D : WeilDiv Y} (hD : IsCartierDiv hnormY D) (w : PrimeDivisorPt X)
    (hdim : ringKrullDim (Y.presheaf.stalk (g.base w.1)) ≤ 1)
    {U : Y.Opens} {f : Y.functionField} (hf : f ≠ 0) (hw : g.base w.1 ∈ U)
    (hDU : ∀ v : PrimeDivisorPt Y, (v : Y) ∈ U → D v = ordPt Y hnormY v f) :
    pullCoeff g hnormX hnormY hD w = ordPt X hnormX w (ffMap g f) := by
  obtain ⟨hy0, hf0, hDU0⟩ := (cartier_local hnormY hD (g.base w.1)).choose_spec
  refine ordPt_ffMap_congr g hnormX hnormY w hdim hf0 hf (fun hc => ?_)
  rw [← hDU0 ⟨g.base w.1, hc⟩ hy0, ← hDU ⟨g.base w.1, hc⟩ hw]

/-! ## ★4. 台の有限性 -/

/-- ★★★★**引き戻しの係数の台は有限**(`X` は準コンパクト)。

★各点で「局所方程式が定義されているアフィン開」を取り、
そこでは係数が `ord_w(g^*f)` になるので `div-finite` が効く。 -/
theorem finite_pullCoeff_ne_zero (hnormX : IsNormalScheme X) (hnormY : IsNormalScheme Y)
    {D : WeilDiv Y} (hD : IsCartierDiv hnormY D) [CompactSpace X]
    (hdim : ∀ w : PrimeDivisorPt X,
      ringKrullDim (Y.presheaf.stalk (g.base w.1)) ≤ 1) :
    {w : PrimeDivisorPt X | pullCoeff g hnormX hnormY hD w ≠ 0}.Finite := by
  have key : ∀ x : X, ∃ W : X.Opens, x ∈ W ∧
      {w : PrimeDivisorPt X | (w : X) ∈ W ∧
        pullCoeff g hnormX hnormY hD w ≠ 0}.Finite := by
    intro x
    obtain ⟨hy0, hf0, hDU0⟩ := (cartier_local hnormY hD (g.base x)).choose_spec
    obtain ⟨_, ⟨W, hWaff, rfl⟩, hxW, hWU⟩ := X.isBasis_affineOpens.exists_subset_of_mem_open
      (show x ∈ (g ⁻¹ᵁ (cartier_local hnormY hD (g.base x)).choose.1 : X.Opens) from hy0)
      (g ⁻¹ᵁ (cartier_local hnormY hD (g.base x)).choose.1).2
    refine ⟨W, hxW, ?_⟩
    haveI : Nonempty W := ⟨⟨x, hxW⟩⟩
    have hffne : ffMap g (cartier_local hnormY hD (g.base x)).choose.2 ≠ 0 := by
      intro h
      exact hf0 ((ffMap g).hom.injective (by simpa using h))
    refine Set.Finite.subset (finite_ordPt_ne_zero_on_affine hnormX hWaff hffne) ?_
    rintro w ⟨hwW, hne⟩
    refine ⟨hwW, ?_⟩
    rw [← pullCoeff_eq g hnormX hnormY hD w (hdim w) hf0 (hWU hwW) hDU0]
    exact hne
  choose W hxW hWfin using key
  obtain ⟨s, hs⟩ := isCompact_univ.elim_finite_subcover (fun x : X => (W x : Set X))
    (fun x => (W x).2) (fun x _ => Set.mem_iUnion.mpr ⟨x, hxW x⟩)
  refine Set.Finite.subset (Set.Finite.biUnion s.finite_toSet (fun x _ => hWfin x)) ?_
  intro w hw
  obtain ⟨x, hx⟩ := Set.mem_iUnion.mp (hs (Set.mem_univ (w : X)))
  obtain ⟨hxs, hwW⟩ := Set.mem_iUnion.mp hx
  exact Set.mem_biUnion hxs ⟨hwW, hw⟩

open scoped Classical in
/-- ★★★★★★**Cartier 因子の引き戻し** `g^* D : WeilDiv X`。 -/
noncomputable def cartierPullback (hnormX : IsNormalScheme X) (hnormY : IsNormalScheme Y)
    {D : WeilDiv Y} (hD : IsCartierDiv hnormY D) [CompactSpace X]
    (hdim : ∀ w : PrimeDivisorPt X,
      ringKrullDim (Y.presheaf.stalk (g.base w.1)) ≤ 1) : WeilDiv X :=
  Finsupp.onFinset (finite_pullCoeff_ne_zero g hnormX hnormY hD hdim).toFinset
    (pullCoeff g hnormX hnormY hD) (fun w hw => by simpa using hw)

@[simp] theorem cartierPullback_apply (hnormX : IsNormalScheme X) (hnormY : IsNormalScheme Y)
    {D : WeilDiv Y} (hD : IsCartierDiv hnormY D) [CompactSpace X]
    (hdim : ∀ w : PrimeDivisorPt X,
      ringKrullDim (Y.presheaf.stalk (g.base w.1)) ≤ 1) (w : PrimeDivisorPt X) :
    cartierPullback g hnormX hnormY hD hdim w = pullCoeff g hnormX hnormY hD w := rfl

/-! ## ★5. 引き戻しは Cartier で、加法的 -/

/-- ★★★★★**Cartier 因子の引き戻しは Cartier**。

★局所方程式 `f` を引き戻した `g^*f` がそのまま局所方程式になる。 -/
theorem isCartierDiv_cartierPullback (hnormX : IsNormalScheme X) (hnormY : IsNormalScheme Y)
    {D : WeilDiv Y} (hD : IsCartierDiv hnormY D) [CompactSpace X]
    (hdim : ∀ w : PrimeDivisorPt X,
      ringKrullDim (Y.presheaf.stalk (g.base w.1)) ≤ 1) :
    IsCartierDiv hnormX (cartierPullback g hnormX hnormY hD hdim) := by
  intro x
  obtain ⟨hy0, hf0, hDU0⟩ := (cartier_local hnormY hD (g.base x)).choose_spec
  refine ⟨g ⁻¹ᵁ (cartier_local hnormY hD (g.base x)).choose.1, hy0,
    ffMap g (cartier_local hnormY hD (g.base x)).choose.2, ?_, fun v hv => ?_⟩
  · intro h
    exact hf0 ((ffMap g).hom.injective (by simpa using h))
  · rw [cartierPullback_apply]
    exact pullCoeff_eq g hnormX hnormY hD v (hdim v) hf0 hv hDU0

/-- ★★★★**引き戻しは加法的**。 -/
theorem pullCoeff_add (hnormX : IsNormalScheme X) (hnormY : IsNormalScheme Y)
    {D₁ D₂ : WeilDiv Y} (hD₁ : IsCartierDiv hnormY D₁) (hD₂ : IsCartierDiv hnormY D₂)
    (hDs : IsCartierDiv hnormY (D₁ + D₂)) (w : PrimeDivisorPt X)
    (hdim : ringKrullDim (Y.presheaf.stalk (g.base w.1)) ≤ 1) :
    pullCoeff g hnormX hnormY hDs w
      = pullCoeff g hnormX hnormY hD₁ w + pullCoeff g hnormX hnormY hD₂ w := by
  obtain ⟨U₁, hw₁, f₁, hf₁, hD₁U⟩ := hD₁ (g.base w.1)
  obtain ⟨U₂, hw₂, f₂, hf₂, hD₂U⟩ := hD₂ (g.base w.1)
  have hsum : ∀ v : PrimeDivisorPt Y, (v : Y) ∈ U₁ ⊓ U₂ →
      (D₁ + D₂) v = ordPt Y hnormY v (f₁ * f₂) := by
    rintro v ⟨hv₁, hv₂⟩
    rw [Finsupp.add_apply, hD₁U v hv₁, hD₂U v hv₂, ordPt_mul hnormY v hf₁ hf₂]
  rw [pullCoeff_eq g hnormX hnormY hDs w hdim (U := U₁ ⊓ U₂) (mul_ne_zero hf₁ hf₂)
      ⟨hw₁, hw₂⟩ hsum,
    pullCoeff_eq g hnormX hnormY hD₁ w hdim hf₁ hw₁ hD₁U,
    pullCoeff_eq g hnormX hnormY hD₂ w hdim hf₂ hw₂ hD₂U, map_mul,
    ordPt_mul hnormX w
      (by intro h; exact hf₁ ((ffMap g).hom.injective (by simpa using h)))
      (by intro h; exact hf₂ ((ffMap g).hom.injective (by simpa using h)))]


/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Example 6.1` の Cartier 因子の引き戻し。 -/
def cartierPullback.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Cartier 因子の引き戻し(局所方程式を引き戻す)",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
