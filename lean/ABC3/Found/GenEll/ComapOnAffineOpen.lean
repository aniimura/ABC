/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ComapChartIdeal
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★`ψ^*超平面` の値をアフィン開ごとに読む（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★★これは何か —— 両側が同じ生成元になった

`§9-923`（段 E0 の側）:

    `divIdeal M s₀ (X_t ⊓ V) = span { (s₀/t) の制限 }`

`§9-924`（幾何の側）:

    `((超平面).comap ψ).comap ι_i = ofIdealTop (span {s₀/s_i})`

★★★★本ファイルは後者を**アフィン開 `W ≤ X_{s_i}` の上の値**に開く:

    **`(ψ^*超平面).ideal W = span { (s₀/s_i) の `W` への制限 }`**

——`§9-923` と**文字どおり同じ形**である。

## ★★★機構 —— `ι` の `appIso` は恒等である

★mathlib の `Scheme.Opens.ι_appIso : U.ι.appIso V = Iso.refl _` が効く。
これで `ideal_comap_of_isOpenImmersion` の `comap` が**消える**:

    `(K.comap U.ι).ideal V = K.ideal ⟨U.ι ''ᵁ V, _⟩`

★★あとは `ofIdealTop` の値（`Ideal.map` の制限）と `topIso.inv` の合成が
**素直な制限になる**こと（`res_topIso_inv`）だけである。
どちらも `X.presheaf.map` の合成なので `Opens` の射が `Subsingleton` であることで閉じる。

## ★配管の記録

★★★★**開集合が `U.ι ''ᵁ (U.ι ⁻¹ᵁ W)` の形で出てくる**のが摩擦の源である。
`congrArg` は `ideal` が依存型なので使えない。
★**汎化してから `rintro rfl`** すれば、証明無関係で残りは自動で合う:

```lean
have key : ∀ W' hW' hW'le, W' = W → (W' についての主張) → (W についての主張) := by
  rintro W' hW' hW'le rfl h; exact h
```

## ★残っている段（明示）

★★★残るのは `divisorOfSection M s₀ = ψ^*超平面` である。道筋は:

1. 自明化つきアフィン `V` について `(ψ^*超平面).ideal V = divIdeal M s₀ V`
   （`V` を `{V ⊓ X_{s_i}}` で覆って `ext_of_iSup_eq_top`）
2. したがって `(ψ^*超平面).ideal ≤ divIdeal`（自明化が無い開では `divIdeal = ⊤`）
3. `le_ofIdeals_iff` で `ψ^*超平面 ≤ divisorOfSection`
4. 逆向きは `le_of_iSup_eq_top` を `{X_{s_i} ⊓ U_j}` に当てる
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★`ι` に沿った `comap` は素通しである -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★**開部分スキームの `ι` に沿った `comap` は素通しである**。

★mathlib の `Scheme.Opens.ι_appIso : U.ι.appIso V = Iso.refl _` が効く。 -/
theorem ideal_comap_ι (K : X.IdealSheafData) (U : X.Opens)
    (V : (U.toScheme).affineOpens) :
    (K.comap U.ι).ideal V
      = K.ideal ⟨U.ι ''ᵁ V.1, V.2.image_of_isOpenImmersion U.ι⟩ := by
  rw [Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion K U.ι V, Scheme.Opens.ι_appIso]
  exact Submodule.toSubMulAction_inj.mp rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★★**`topIso.inv` のあとの制限は素直な制限である**。 -/
theorem res_topIso_inv (U : X.Opens) (V : (U.toScheme).Opens)
    (hle2 : V ≤ (⊤ : (U.toScheme).Opens)) (h' : U.ι ''ᵁ V ≤ U)
    (g : (Γ(X, U) : Type)) :
    ((U.toScheme).presheaf.map (homOfLE hle2).op).hom (U.topIso.inv.hom g)
      = (X.presheaf.map (homOfLE h').op).hom g := by
  rw [Scheme.Opens.topIso_inv, Scheme.Opens.toScheme_presheaf_map,
    ← CommRingCat.comp_apply, ← Functor.map_comp]
  congr 1

/-! ## ★★★★★★★★★★★★★★★アフィン開の上での値 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★**`ψ^*超平面` はアフィン開 `W ≤ X_{s_i}` の上で
`(s₀/s_i)` の制限が生成する**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★これは `§9-923`（`divIdeal M s₀ (X_t ⊓ V) = span {(s₀/t) の制限}`）と
**文字どおり同じ形**である——段 E0 の両側が同じ生成元になった。 -/
theorem ideal_comap_globalToProj_of_le (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (i : Fin (N + 1))
    (haff : IsAffineOpen (nonVanishing M (s i)))
    {W : X.Opens} (hW : IsAffineOpen W) (hWle : W ≤ nonVanishing M (s i)) :
    ((hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov)).ideal ⟨W, hW⟩
      = Ideal.span {(X.presheaf.map (homOfLE hWle).op).hom
          (globalRatio M hM (s 0) (s i))} := by
  haveI : IsAffine (nonVanishing M (s i)).toScheme := haff
  have hrange : W ≤ (nonVanishing M (s i)).ι.opensRange := by
    rw [Scheme.Opens.opensRange_ι]; exact hWle
  have hVaff : IsAffineOpen ((nonVanishing M (s i)).ι ⁻¹ᵁ W) :=
    hW.preimage_of_isOpenImmersion (nonVanishing M (s i)).ι hrange
  have himg : (nonVanishing M (s i)).ι ''ᵁ ((nonVanishing M (s i)).ι ⁻¹ᵁ W) = W := by
    rw [Scheme.Hom.image_preimage_eq_opensRange_inf, Scheme.Opens.opensRange_ι,
      inf_eq_right.mpr hWle]
  have key : ∀ (W' : X.Opens) (hW' : IsAffineOpen W') (hW'le : W' ≤ nonVanishing M (s i)),
      W' = W →
      (((hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov)).ideal ⟨W', hW'⟩
        = Ideal.span {(X.presheaf.map (homOfLE hW'le).op).hom
            (globalRatio M hM (s 0) (s i))}) →
      (((hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov)).ideal ⟨W, hW⟩
        = Ideal.span {(X.presheaf.map (homOfLE hWle).op).hom
            (globalRatio M hM (s 0) (s i))}) := by
    rintro W' hW' hW'le rfl h
    exact h
  refine key _ (hVaff.image_of_isOpenImmersion (nonVanishing M (s i)).ι)
    (Scheme.Opens.ι_image_le _ _) himg ?_
  have h1 := ideal_comap_ι ((hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov))
    (nonVanishing M (s i)) ⟨_, hVaff⟩
  rw [comap_globalToProj_chart M hM φ s hcov i haff,
    Scheme.IdealSheafData.ofIdealTop_ideal, Ideal.map_span, Set.image_singleton] at h1
  rw [← h1, res_topIso_inv]

/-! ## ★出典の紐付け(`.src`) -/

def ideal_comap_ι.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(開部分スキームの ι に沿った comap は素通しである)",
    sectionId := "genell-prop-1-4" }

def res_topIso_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(topIso.inv のあとの制限は素直な制限である)",
    sectionId := "genell-prop-1-4" }

def ideal_comap_globalToProj_of_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(ψ^*超平面はアフィン開の上で (s₀/s_i) の制限が生成する)",
    sectionId := "genell-prop-1-4" }

def ideal_comap_globalToProj_of_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "comap_globalToProj_chart(チャートの上での値、§9-924)"
      (.inProject "ABC3" "ABC3.Found.GenEll.comap_globalToProj_chart") 3,
    .citation "[mathlib]" "Scheme.Opens.ι_appIso(ι の appIso は恒等)"
      (.inMathlib "AlgebraicGeometry.Scheme.Opens.ι_appIso") 1,
    .citation "[mathlib]" "IsAffineOpen.preimage_of_isOpenImmersion"
      (.inMathlib "AlgebraicGeometry.IsAffineOpen.preimage_of_isOpenImmersion") 1,
    .implicitStep
      ("★★これは §9-923(divIdeal M s₀ (X_t ⊓ V) = span {(s₀/t) の制限})と" ++
       "**文字どおり同じ形**である——段 E0 の両側が同じ生成元になった") 4,
    .implicitStep
      ("★★★★配管: 開集合が U.ι ''ᵁ (U.ι ⁻¹ᵁ W) の形で出てくるのが摩擦の源である。" ++
       "congrArg は ideal が依存型なので使えない。" ++
       "**汎化してから rintro rfl** すれば、証明無関係で残りは自動で合う") 2,
    .implicitStep
      ("★★★残るのは divisorOfSection M s₀ = ψ^*超平面 である。道筋は " ++
       "(1) 自明化つきアフィン V について (ψ^*超平面).ideal V = divIdeal M s₀ V" ++
       "(V を {V ⊓ X_{s_i}} で覆って ext_of_iSup_eq_top)、" ++
       "(2) したがって (ψ^*超平面).ideal ≤ divIdeal(自明化が無い開では divIdeal = ⊤)、" ++
       "(3) le_ofIdeals_iff で ψ^*超平面 ≤ divisorOfSection、" ++
       "(4) 逆向きは le_of_iSup_eq_top を {X_{s_i} ⊓ U_j} に当てる") 4 ]

end ABC3.Found.GenEll
