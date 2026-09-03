/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ComapOnAffineOpen
import ABC3.Found.GenEll.DivIdealChart
import ABC3.Meta.Claim
import ABC3.Found.GenEll.DivisorOfSectionComap.Definition12

/-!
# DivisorOfSectionComap —— `[GenEll] Proposition 1.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
variable {X : Scheme.{0}}

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
