/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ComapOnAffineOpen
import ABC3.Found.GenEll.DivIdealChart
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★段 E0 が閉じた —— `div(s₀) = ψ^*超平面`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5-7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★★★★★これは何か —— 段 E0 の到達点

    **`divisorOfSection M s₀ = (超平面).comap ψ`**

★これは「**切断 `s₀` の零点因子は、貼った射に沿った超平面の引き戻しである**」という、
原文が `Proposition 1.4, (iv)` の証明で暗黙に使っている等式そのものである。

★★★★これで **`§9-922` の `hht`（点ごとの高さの等式）が因子の等式から出る**
（`htArith_comap_eq_nsmul`）。

## ★★★機構 —— 4 段

| 段 | 内容 | 道具 |
|---|---|---|
| 1 | 自明化つきアフィン `V` で `(ψ^*超平面).ideal V = divIdeal M s₀ V` | `ext_of_iSup_eq_top` を `{V ⊓ X_{s_i}}` に |
| 2 | よって `(ψ^*超平面).ideal ≤ divIdeal`（自明化が無い開では `divIdeal = ⊤`） | `divIdeal_eq_top` |
| 3 | `le_ofIdeals_iff` で `ψ^*超平面 ≤ divisorOfSection` | mathlib |
| 4 | 逆向きは `le_of_iSup_eq_top` を `{X_{s_i} ⊓ U_j}` に | mathlib |

★★段 2 の「自明化が無い開では `divIdeal = ⊤`」——**段 E0 で長らく障害と見えていたもの**が、
ここでは**逆に効く**（`⊤` なら `≤` は自明）。

## ★★★★測定の記録

★★★★★**`divIdeal` の `⊤` は障害ではなく助けだった**。
`divisorOfSection = ofIdeals divIdeal` は「`divIdeal` 以下の**最大の**イデアル層」であり、
`ψ^*超平面` が `divIdeal` 以下であることさえ言えれば `≤` が出る。
★そして「以下」は自明化が無い開では**ただ**である。
★★残るのは自明化つき開での**等式**だけで、それは `§9-925`（幾何の側）と
`§9-923`（`divIdeal` の側）が同じ生成元を出すことから従う。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★`divIdeal` を任意の小さいアフィン開で読む -/

/-- ★★★★★★**`divIdeal M s₀ W = span {(s₀/t) の制限}`**
（`W` が自明化つき `V` の中の、`X_t` に含まれるアフィン開なら）。

★`§9-923` の一般化（そこでは `W = X_t ⊓ V` に固定していた）。
★★機構は同じ——`trivValue(s₀)` と `sectionRatio(s₀,t)` は単元倍しか違わない。 -/
theorem divIdeal_eq_span_globalRatio_of_le (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (V : X.Opens) (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s₀ t : (M.obj (op ⊤) : Type))
    {W : X.Opens} (hW : IsAffineOpen W) (hWV : W ≤ V) (hWt : W ≤ nonVanishing M t) :
    divIdeal M s₀ ⟨W, hW⟩
      = Ideal.span {(X.presheaf.map (homOfLE hWt).op).hom (globalRatio M hM s₀ t)} := by
  have hle : W ≤ nonVanishing M t ⊓ V := le_inf hWt hWV
  have h1 : (X.presheaf.map (homOfLE hWt).op).hom (globalRatio M hM s₀ t)
      = (X.presheaf.map (homOfLE hle).op).hom (sectionRatio M V e s₀ t) := by
    rw [← globalRatio_res M hM s₀ t ⟨V, e⟩, res_trans]
  rw [divIdeal_eq M s₀ ⟨W, hW⟩ (trivialOfLe M hWV e), trivValue_restrict M hWV e s₀,
    h1, sectionRatio, map_mul, res_trans]
  refine (Ideal.span_singleton_mul_right_unit ?_ _).symm
  exact ((isUnit_trivValue_res M V e t).unit⁻¹.isUnit).map
    (X.presheaf.map (homOfLE hle).op).hom

/-! ## ★★★★★★★★★★★★★★段 1 —— 自明化つきアフィン開での等式 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★**自明化つきアフィン開の上で
`(ψ^*超平面).ideal V = divIdeal M s₀ V`**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★`V` を `{V ⊓ X_{s_i}}` で覆い、`ext_of_iSup_eq_top` で貼る。
★★どの成分でも両側が `span {(s₀/s_i) の制限}` になる
——`§9-925`（幾何の側）と `§9-923`（`divIdeal` の側）である。 -/
theorem ideal_comap_globalToProj_triv (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (haff : ∀ i, IsAffineOpen (nonVanishing M (s i)))
    (V : X.Opens) (hV : IsAffineOpen V)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    ((hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov)).ideal ⟨V, hV⟩
      = divIdeal M (s 0) ⟨V, hV⟩ := by
  haveI : IsAffine V.toScheme := hV
  set K := (hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov) with hK
  set A : Ideal (Γ(V.toScheme, (⊤ : V.toScheme.Opens)) : Type) :=
    Ideal.span {V.topIso.inv.hom (trivValue M V e (s 0))} with hA
  have hVaff : ∀ i, IsAffineOpen (V.ι ⁻¹ᵁ (nonVanishing M (s i))) := by
    intro i
    rw [← Scheme.Hom.isAffineOpen_iff_of_isOpenImmersion V.ι,
      Scheme.Hom.image_preimage_eq_opensRange_inf, Scheme.Opens.opensRange_ι, inf_comm,
      nonVanishing_inf M V e (s i)]
    exact hV.basicOpen _
  have hcov' : (⨆ i, V.ι ⁻¹ᵁ (nonVanishing M (s i))) = ⊤ := by
    rw [← Scheme.Hom.preimage_iSup, hcov]; rfl
  have hmain : K.comap V.ι = Scheme.IdealSheafData.ofIdealTop A := by
    refine Scheme.IdealSheafData.ext_of_iSup_eq_top
      (fun i => ⟨V.ι ⁻¹ᵁ (nonVanishing M (s i)), hVaff i⟩) hcov' (fun i => ?_)
    have hW'V : V.ι ''ᵁ (V.ι ⁻¹ᵁ (nonVanishing M (s i))) ≤ V := Scheme.Opens.ι_image_le _ _
    have hW'i : V.ι ''ᵁ (V.ι ⁻¹ᵁ (nonVanishing M (s i))) ≤ nonVanishing M (s i) := by
      rw [Scheme.Hom.image_preimage_eq_opensRange_inf]; exact inf_le_right
    rw [ideal_comap_ι K V ⟨_, hVaff i⟩,
      ideal_comap_globalToProj_of_le M hM φ s hcov i (haff i)
        ((hVaff i).image_of_isOpenImmersion V.ι) hW'i,
      Scheme.IdealSheafData.ofIdealTop_ideal, hA, Ideal.map_span, Set.image_singleton,
      res_topIso_inv V _ le_top hW'V,
      ← divIdeal_eq_span_globalRatio_of_le M hM V e (s 0) (s i)
        ((hVaff i).image_of_isOpenImmersion V.ι) hW'V hW'i,
      divIdeal_eq M (s 0) ⟨_, (hVaff i).image_of_isOpenImmersion V.ι⟩
        (trivialOfLe M hW'V e), trivValue_restrict M hW'V e (s 0)]
    rfl
  have key : ∀ (W' : X.Opens) (hW' : IsAffineOpen W'), W' = V →
      (K.ideal ⟨W', hW'⟩ = divIdeal M (s 0) ⟨W', hW'⟩) →
      (K.ideal ⟨V, hV⟩ = divIdeal M (s 0) ⟨V, hV⟩) := by
    rintro W' hW' rfl h
    exact h
  refine key _ ((isAffineOpen_top V.toScheme).image_of_isOpenImmersion V.ι)
    V.ι_image_top ?_
  have h2 := ideal_comap_ι K V ⟨⊤, isAffineOpen_top V.toScheme⟩
  rw [hmain, Scheme.IdealSheafData.ofIdealTop_ideal, map_res_self_ideal, hA] at h2
  have hW'V : V.ι ''ᵁ (⊤ : V.toScheme.Opens) ≤ V := Scheme.Opens.ι_image_le _ _
  rw [← h2, divIdeal_eq M (s 0) ⟨_, (isAffineOpen_top V.toScheme).image_of_isOpenImmersion V.ι⟩
      (trivialOfLe M hW'V e), trivValue_restrict M hW'V e (s 0),
    ← res_topIso_inv V ⊤ le_rfl hW'V (trivValue M V e (s 0))]
  congr 1
  rw [show (homOfLE (le_rfl : (⊤ : V.toScheme.Opens) ≤ ⊤) : (⊤ : V.toScheme.Opens) ⟶ ⊤)
      = 𝟙 _ from Subsingleton.elim _ _, op_id, V.toScheme.presheaf.map_id]
  rfl

/-! ## ★★★★★★★★★★★★★★★★★★段 2-4 —— 到達 -/

/-- ★★★★★★★★★★★★★★★★★★**`div(s₀) = ψ^*超平面`** —— 段 E0 の到達点。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★「切断 `s₀` の零点因子は、貼った射に沿った超平面の引き戻しである」。
★★★★これで `§9-922` の `hht`（点ごとの高さの等式）が
`htArith_comap_eq_nsmul` を通して**因子の等式から出る**。

★★★測定: **`divIdeal` の `⊤`（自明化が無い開）は障害ではなく助けだった**
——`divisorOfSection = ofIdeals divIdeal` は「`divIdeal` 以下の最大」なので、
`ψ^*超平面 ≤ divisorOfSection` を出すには `≤ divIdeal` さえ言えばよく、
`⊤` の所ではそれが**ただ**である。 -/
theorem divisorOfSection_eq_comap_globalToProj (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    {N : ℕ} (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (haff : ∀ i, IsAffineOpen (nonVanishing M (s i)))
    {ι : Type} (U : ι → X.Opens) (hcovU : (⨆ j, U j) = ⊤)
    (hU : ∀ j, IsAffineOpen (U j))
    (eU : ∀ j, (restrictPresheafFunctor X (U j)).obj M ≅ 𝟙_ (PresheafModulesOn X (U j))) :
    divisorOfSection M (s 0)
      = (hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov) := by
  classical
  have hle1 : ((hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov)).ideal
      ≤ divIdeal M (s 0) := by
    intro V
    by_cases h : Nonempty ((restrictPresheafFunctor X V.1).obj M ≅ 𝟙_ (PresheafModulesOn X V.1))
    · obtain ⟨e⟩ := h
      exact le_of_eq (ideal_comap_globalToProj_triv M hM φ s hcov haff V.1 V.2 e)
    · rw [divIdeal_eq_top M (s 0) V (not_nonempty_iff.mp h)]
      exact le_top
  refine le_antisymm ?_ (Scheme.IdealSheafData.le_ofIdeals_iff.2 hle1)
  refine Scheme.IdealSheafData.le_of_iSup_eq_top
    (fun p : Fin (N + 1) × ι =>
      (⟨nonVanishing M (s p.1) ⊓ U p.2,
        isAffineOpen_nonVanishing_inf M (U p.2) (hU p.2) (eU p.2) (s p.1)⟩ : X.affineOpens))
    (iSup_nonVanishing_inf_eq_top M s U hcovU hcov) (fun p => ?_)
  refine le_trans (divisorOfSection_ideal_le M (s 0) _) (le_of_eq ?_)
  exact (ideal_comap_globalToProj_triv M hM φ s hcov haff _ _
    (trivialOfLe M (inf_le_right : nonVanishing M (s p.1) ⊓ U p.2 ≤ U p.2) (eU p.2))).symm

/-! ## ★出典の紐付け(`.src`) -/

def divIdeal_eq_span_globalRatio_of_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2(divIdeal は任意の小さいアフィン開で (s₀/t) の制限が生成する)",
    sectionId := "genell-def-1-2" }

def ideal_comap_globalToProj_triv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(自明化つきアフィン開の上で ψ^*超平面 = divIdeal)",
    sectionId := "genell-prop-1-4" }

def divisorOfSection_eq_comap_globalToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E0——div(s₀) = ψ^*超平面)",
    sectionId := "genell-prop-1-4" }

def divisorOfSection_eq_comap_globalToProj.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ideal_comap_globalToProj_of_le(幾何の側、§9-925)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ideal_comap_globalToProj_of_le") 3,
    .citation "[ABC3]" "divIdeal_eq(自明化があればどれで測ってもよい、§9-809)"
      (.inProject "ABC3" "ABC3.Found.GenEll.divIdeal_eq") 1,
    .citation "[mathlib]" "Scheme.IdealSheafData.ext_of_iSup_eq_top・le_of_iSup_eq_top"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.ext_of_iSup_eq_top") 1,
    .implicitStep
      ("★★★★測定: **divIdeal の ⊤(自明化が無い開)は障害ではなく助けだった**" ++
       "——divisorOfSection = ofIdeals divIdeal は「divIdeal 以下の最大」なので、" ++
       "ψ^*超平面 ≤ divisorOfSection を出すには ≤ divIdeal さえ言えばよく、" ++
       "⊤ の所ではそれが**ただ**である") 4,
    .implicitStep
      ("★★★★これで §9-922 の hht(点ごとの高さの等式)が " ++
       "htArith_comap_eq_nsmul を通して**因子の等式から出る**。" ++
       "★Proposition 1.4, (iv) に残るのは (2) 点がチャートを通ること、" ++
       "(3) 座標が点を分けること(hinj)の 2 つだけになった") 5 ]

end ABC3.Found.GenEll
