/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottGlobalToProj
import ABC3.Found.GenEll.DivIdealChart
import ABC3.Found.GenEll.HyperplanePullback
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★引き戻しをチャートの上で計算する（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★これは何か —— 段 E0 の最後の計算

`§9-923` で `divIdeal M s₀` は自明化つきアフィン開の上で `span {s₀/s_i}` だと分かった。
★残っていたのは **`(ψ^*超平面).ideal` を同じ被覆の上で計算すること**である。

★★★★本ファイルはその**チャート全体の上での値**を取る:

    `((超平面).comap ψ).comap ι_i = ofIdealTop (span {s₀/s_i})`   （`X_{s_i}` の上で）

——`X_{s_i}` は（`IsAmple` から）アフィンなので、イデアル層は `⊤` での値だけで決まる。

## ★★★機構 —— アフィンなら `pullbackIdealOf` に落ちる

★`Z` がアフィンなら `Z.toSpecΓ` は同型なので、任意の `χ : Z ⟶ Y` は
`χ = Z.toSpecΓ ≫ x`（`x : Spec Γ(Z,⊤) ⟶ Y`）と書ける。そのとき

    **`(I.comap χ).ideal ⊤ = pullbackIdealOf Γ(Z,⊤) I x`**   （`ideal_comap_toSpecΓ_comp`）

★★これで `§9-865` の `pullbackIdealOf_hyperplane_point`（点に沿った超平面の引き戻し）が
**イデアル層の水準でそのまま使える**。

★★★`globalChartMorphism`（`§9-848`）はまさに `Z.toSpecΓ ≫ Spec.map β` の形なので、
`globalChartToProj i = Z.toSpecΓ ≫ (Spec.map β ≫ chartA i)` と括り直すだけでよい。

## ★配管の記録

★★★★同一視が**ちょうど打ち消し合う**のが要点である:

* `appIso (Z.toSpecΓ) ⊤` と `ΓSpecIso Γ(Z,⊤)` は
  `Scheme.toSpecΓ_appTop`（`toSpecΓ.appTop = (ΓSpecIso _).hom`）で繋がる
* `Ideal.comap (α.inv)` と `Ideal.map (α.hom)` は同型なら同じ（`comap_inv_eq_map_hom`）

★したがって `pullbackIdealOf` の定義に現れる `ΓSpecIso` の `comap` と、
`ideal_comap_of_isOpenImmersion` に現れる `appIso` の `comap` が**相殺する**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★同型に沿った `comap` は `map` である -/

/-- ★**同型に沿った `comap` は逆向きの `map` である**。 -/
theorem comap_inv_eq_map_hom {A B : CommRingCat.{0}} (α : A ≅ B) (J : Ideal A) :
    Ideal.comap α.inv.hom J = Ideal.map α.hom.hom J := by
  apply le_antisymm
  · intro b hb
    have hb' : α.hom.hom (α.inv.hom b) = b := by
      rw [← CommRingCat.comp_apply, α.inv_hom_id]; rfl
    rw [← hb']
    exact Ideal.mem_map_of_mem _ hb
  · rw [Ideal.map_le_iff_le_comap]
    intro a ha
    show α.inv.hom (α.hom.hom a) ∈ J
    rw [← CommRingCat.comp_apply, α.hom_inv_id]
    simpa using ha

/-- ★**同じ開への制限に沿った `Ideal.map` は恒等**。 -/
theorem map_res_self_ideal {Z : Scheme.{0}} {U : Z.Opens} (h : U ≤ U)
    (J : Ideal (Γ(Z, U) : Type)) :
    Ideal.map (CommRingCat.Hom.hom (Z.presheaf.map (homOfLE h).op)) J = J := by
  rw [show (homOfLE h : U ⟶ U) = 𝟙 U from Subsingleton.elim _ _, op_id, Z.presheaf.map_id,
    CommRingCat.hom_id, Ideal.map_id]

/-! ## ★★★★アフィンなら `pullbackIdealOf` に落ちる -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★**`appIso (toSpecΓ)` と `ΓSpecIso` は繋がる**。 -/
theorem res_comp_appIso_toSpecΓ (Z : Scheme.{0}) [IsAffine Z] :
    (Spec Γ(Z, (⊤ : Z.Opens))).presheaf.map
        (homOfLE (show Z.toSpecΓ ''ᵁ (⊤ : Z.Opens) ≤ ⊤ from fun ⦃_⦄ _ => trivial)).op
      ≫ (Scheme.Hom.appIso Z.toSpecΓ (⊤ : Z.Opens)).hom
      = (Scheme.ΓSpecIso Γ(Z, (⊤ : Z.Opens))).hom := by
  rw [Scheme.Hom.appIso_hom', Scheme.Hom.map_appLE, ← Scheme.toSpecΓ_appTop]
  simp [Scheme.Hom.appLE, Scheme.Hom.appTop]

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★**アフィンな始域からの引き戻しは `pullbackIdealOf` である**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `(I.comap (Z.toSpecΓ ≫ x)).ideal ⊤ = pullbackIdealOf Γ(Z,⊤) I x`

★★これで点の水準の計算（`§9-865`）が**イデアル層の水準でそのまま使える**。 -/
theorem ideal_comap_toSpecΓ_comp (Z : Scheme.{0}) [IsAffine Z] {Y : Scheme.{0}}
    (I : Y.IdealSheafData) (x : Spec Γ(Z, (⊤ : Z.Opens)) ⟶ Y) :
    (I.comap (Z.toSpecΓ ≫ x)).ideal ⟨⊤, isAffineOpen_top Z⟩
      = pullbackIdealOf Γ(Z, (⊤ : Z.Opens)) I x := by
  rw [Scheme.IdealSheafData.comap_comp,
    Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion _ Z.toSpecΓ ⟨⊤, isAffineOpen_top Z⟩,
    ← Scheme.IdealSheafData.map_ideal (I.comap x)
      (Subtype.mk_le_mk.mpr (fun ⦃_⦄ _ => trivial) :
        (⟨Z.toSpecΓ ''ᵁ (⊤ : Z.Opens),
            (isAffineOpen_top Z).image_of_isOpenImmersion Z.toSpecΓ⟩
          : (Spec Γ(Z, (⊤ : Z.Opens))).affineOpens) ≤ ⟨⊤, isAffineOpen_top _⟩),
    pullbackIdealOf, Scheme.IdealSheafData.equivOfIsAffine_apply,
    comap_inv_eq_map_hom, comap_inv_eq_map_hom, Ideal.map_map, ← CommRingCat.hom_comp]
  exact congrArg
    (fun (g : (Spec Γ(Z, (⊤ : Z.Opens))).presheaf.obj (op ⊤) ⟶ Γ(Z, (⊤ : Z.Opens))) =>
      Ideal.map g.hom ((I.comap x).ideal ⟨⊤, isAffineOpen_top _⟩))
    (res_comp_appIso_toSpecΓ Z)

/-! ## ★★★★★★★★★★★★★★チャートの上での超平面の引き戻し -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★**チャートの上で超平面の引き戻しは `(s₀/s_i)` が生成する**
（イデアル層の水準）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★`globalChartToProj i = Z.toSpecΓ ≫ (Spec.map β ≫ chartA i)` と括り直し、
`ideal_comap_toSpecΓ_comp` と `§9-865` を繋ぐだけである。 -/
theorem ideal_comap_globalChartToProj (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (haff : IsAffineOpen (nonVanishing M (s i))) :
    haveI : IsAffine (nonVanishing M (s i)).toScheme := haff
    ((hyperplaneIdeal N ℤ).comap (globalChartToProj M hM φ s i)).ideal
        ⟨⊤, isAffineOpen_top _⟩
      = Ideal.span {((nonVanishing M (s i)).topIso.inv.hom)
          (globalRatio M hM (s 0) (s i))} := by
  haveI : IsAffine (nonVanishing M (s i)).toScheme := haff
  have hfac : globalChartToProj M hM φ s i
      = (nonVanishing M (s i)).toScheme.toSpecΓ ≫
        (Spec.map (CommRingCat.ofHom
          (((nonVanishing M (s i)).topIso.inv.hom).comp (globalAwayHom M hM φ s i)))
          ≫ chartA N ℤ i) := by
    rw [globalChartToProj, globalChartMorphism, Category.assoc]
  rw [hfac, ideal_comap_toSpecΓ_comp, pullbackIdealOf_hyperplane_point]
  show Ideal.span {((nonVanishing M (s i)).topIso.inv.hom)
    (globalAwayHom M hM φ s i (projCoord N ℤ i 0))} = _
  rw [globalAwayHom_projCoord]

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★**`ψ^*超平面` はチャートの上で `(s₀/s_i)` の主イデアル層である**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `((超平面).comap ψ).comap ι_i = ofIdealTop (span {s₀/s_i})`

★`X_{s_i}` は（`IsAmple` から）アフィンなので、イデアル層は `⊤` での値だけで決まる。
★★これが `§9-923`（`divIdeal M s₀` は `span {s₀/s_i}`）と**同じ生成元**である
——段 E0 の両側が出会う点。 -/
theorem comap_globalToProj_chart (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (i : Fin (N + 1))
    (haff : IsAffineOpen (nonVanishing M (s i))) :
    haveI : IsAffine (nonVanishing M (s i)).toScheme := haff
    ((hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov)).comap
        (nonVanishing M (s i)).ι
      = Scheme.IdealSheafData.ofIdealTop
          (Ideal.span {((nonVanishing M (s i)).topIso.inv.hom)
            (globalRatio M hM (s 0) (s i))}) := by
  haveI : IsAffine (nonVanishing M (s i)).toScheme := haff
  refine Scheme.IdealSheafData.ext_of_isAffine ?_
  rw [← Scheme.IdealSheafData.comap_comp, ι_globalToProj,
    ideal_comap_globalChartToProj M hM φ s i haff,
    Scheme.IdealSheafData.ofIdealTop_ideal, map_res_self_ideal]

/-! ## ★出典の紐付け(`.src`) -/

def comap_inv_eq_map_hom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(同型に沿った comap は逆向きの map である)",
    sectionId := "genell-prop-1-4" }

def res_comp_appIso_toSpecΓ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(appIso (toSpecΓ) と ΓSpecIso は繋がる)",
    sectionId := "genell-prop-1-4" }

def ideal_comap_toSpecΓ_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(アフィンな始域からの引き戻しは pullbackIdealOf である)",
    sectionId := "genell-prop-1-4" }

def ideal_comap_globalChartToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの上で超平面の引き戻しは (s₀/s_i)——イデアル層の水準)",
    sectionId := "genell-prop-1-4" }

def comap_globalToProj_chart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(ψ^*超平面はチャートの上で (s₀/s_i) の主イデアル層である)",
    sectionId := "genell-prop-1-4" }

def comap_globalToProj_chart.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pullbackIdealOf_hyperplane_point(点に沿った超平面の引き戻し、§9-865)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pullbackIdealOf_hyperplane_point") 2,
    .citation "[mathlib]" "Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion") 1,
    .citation "[mathlib]" "Scheme.toSpecΓ_appTop(toSpecΓ.appTop = (ΓSpecIso _).hom)"
      (.inMathlib "AlgebraicGeometry.Scheme.toSpecΓ_appTop") 1,
    .implicitStep
      ("★★★★配管の要点は**同一視がちょうど打ち消し合う**ことである" ++
       "——pullbackIdealOf の定義に現れる ΓSpecIso の comap と、" ++
       "ideal_comap_of_isOpenImmersion に現れる appIso の comap が相殺する。" ++
       "繋ぐのは Scheme.toSpecΓ_appTop の 1 本だけ") 3,
    .implicitStep
      ("★★これで §9-923(divIdeal M s₀ は span {s₀/s_i})と**同じ生成元**が出た" ++
       "——段 E0 の両側が出会う点である。" ++
       "★残るのは divisorOfSection M s₀ = ψ^*超平面 を " ++
       "Scheme.IdealSheafData.ext_of_iSup_eq_top で貼ることだけである") 4 ]

end ABC3.Found.GenEll
