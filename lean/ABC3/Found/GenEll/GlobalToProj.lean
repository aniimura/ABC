/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.OverlapAwayHomSymm
import ABC3.Found.GenEll.GlobalChartToProj
import ABC3.Found.GenEll.GlueOpens
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★段 E1d が閉じた —— `ψ : X ⟶ ℙᴺ_R`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★★★これは何か —— 貼り合わせが完了した

`§9-848`（`GlobalChartToProj.lean`）でチャートごとの射 `X_{s_i} ⟶ ℙᴺ_R` が入り、
`§9-836`（`GlueOpens.lean`）で貼り合わせの機構が入り、
`§9-908`〜`§9-910` で「重なりでの一致」の**代数の側**が入った。

★★★★本ファイルはそれを**幾何の側の配管で繋いで貼り合わせる**:

    `⨆_i X_{s_i} = ⊤`  ⟹  `ψ : X ⟶ ℙᴺ_R`

## ★★★機構 —— 3 段

| 段 | 内容 | 道具 |
|---|---|---|
| 1 | チャート射を `W ⊆ X_{s_i}` へ制限する | `Scheme.toSpecΓ_naturality`・`Scheme.homOfLE_appTop` |
| 2 | 両方を同じ `awayι 𝒜 (x_i·x_j)` へ落とす | `Proj.SpecMap_awayMap_awayι` ＋ `§9-910` |
| 3 | 貼る | `glueOpens`（`§9-836`） |

★段 2 が要点である。`§9-910` の `overlapAwayHom_awayMap`（`i` 側）と
`overlapAwayHom_awayMap_symm`（`j` 側）で、**どちらの合成も**

    `W ⟶ Spec A⁰_{x_i x_j} ⟶ ℙᴺ_R`

という同じ形になる。

## ★★測定の記録

★★★★**`topIso` は `eqToHom` の言い換えである**（mathlib の `Scheme.Opens.topIso` は
`X.presheaf.mapIso (eqToIso U.ι_image_top.symm).op`）。
★したがって `V.topIso.inv ≫ (X.homOfLE h).appTop = (制限) ≫ W.topIso.inv` は
**両辺とも `X.presheaf.map` の合成**であり、`Opens` の射が `Subsingleton` なので
`congr 1` だけで閉じる。

★★★★★**`instances` 透明度で 4 回はまった**——`X.presheaf` を
`(Opens X)ᵒᵖ ⥤ CommRingCat` に落とそうとして型が合わなくなる。
`set_option backward.isDefEq.respectTransparency false in` を
本ファイルの `rw` を使う定理すべてに付けた（`tools/lean-idioms.md` の既知の穴）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★段 1 —— チャート射を小さい開へ制限する -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★**`topIso` と制限は交換する**。

★両辺とも `X.presheaf.map` の合成なので、`Opens` の射が `Subsingleton` であることから
`congr 1` だけで閉じる。 -/
theorem topIso_inv_comp_appTop {W V : X.Opens} (h : W ≤ V) :
    V.topIso.inv ≫ Scheme.Hom.appTop (X.homOfLE h)
      = X.presheaf.map (homOfLE h).op ≫ W.topIso.inv := by
  rw [Scheme.homOfLE_appTop, Scheme.Opens.topIso_inv, Scheme.Opens.topIso_inv,
    ← X.presheaf.map_comp, ← X.presheaf.map_comp]
  congr 1

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★**`Spec` への射を小さい開へ制限する**。

★`Scheme.toSpecΓ_naturality` と `topIso_inv_comp_appTop` を繋いだだけである。 -/
theorem homOfLE_comp_specMap {W V : X.Opens} (h : W ≤ V) {A : Type} [CommRing A]
    (α : A →+* (Γ(X, V) : Type)) :
    X.homOfLE h ≫ (V.toScheme.toSpecΓ ≫ Spec.map (CommRingCat.ofHom (V.topIso.inv.hom.comp α)))
      = W.toScheme.toSpecΓ ≫ Spec.map (CommRingCat.ofHom
          (W.topIso.inv.hom.comp (((X.presheaf.map (homOfLE h).op).hom).comp α))) := by
  rw [← Category.assoc, Scheme.toSpecΓ_naturality, Category.assoc, ← Spec.map_comp]
  refine congrArg (fun f => W.toScheme.toSpecΓ ≫ Spec.map f) ?_
  calc CommRingCat.ofHom (V.topIso.inv.hom.comp α) ≫ Scheme.Hom.appTop (X.homOfLE h)
      = CommRingCat.ofHom α ≫ (V.topIso.inv ≫ Scheme.Hom.appTop (X.homOfLE h)) := by
        rw [← Category.assoc]; rfl
    _ = CommRingCat.ofHom α ≫ (X.presheaf.map (homOfLE h).op ≫ W.topIso.inv) := by
        rw [topIso_inv_comp_appTop]
    _ = CommRingCat.ofHom (W.topIso.inv.hom.comp
          (((X.presheaf.map (homOfLE h).op).hom).comp α)) := by
        rw [← Category.assoc]; rfl

/-- ★★★★**チャート射の制限**——`α` を `globalAwayHom` に取っただけである。 -/
theorem homOfLE_comp_globalChartMorphism (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) {W : X.Opens}
    (h : W ≤ nonVanishing M (s i)) :
    X.homOfLE h ≫ globalChartMorphism M hM φ s i
      = W.toScheme.toSpecΓ ≫ Spec.map (CommRingCat.ofHom
          (W.topIso.inv.hom.comp
            (((X.presheaf.map (homOfLE h).op).hom).comp (globalAwayHom M hM φ s i)))) :=
  homOfLE_comp_specMap h (globalAwayHom M hM φ s i)

/-! ## ★★★★★★★★★★★★★★★段 2 —— チャート射は重なりで一致する -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★**段 E1d 本体** ——
2 つのチャート射は重なりで一致する。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★機構は「両方を同じ `awayι 𝒜 (x_i·x_j)` へ落とす」ことである:
`§9-910` の `overlapAwayHom_awayMap`（`i` 側）と `overlapAwayHom_awayMap_symm`（`j` 側）で
どちらの合成も `W ⟶ Spec A⁰_{x_i x_j} ⟶ ℙᴺ_R` という同じ形になり、
`Proj.SpecMap_awayMap_awayι` で `awayι` が揃う。 -/
theorem globalChartToProj_agree (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) {W : X.Opens}
    (hi : W ≤ nonVanishing M (s i)) (hj : W ≤ nonVanishing M (s j)) :
    X.homOfLE hi ≫ globalChartToProj M hM φ s i
      = X.homOfLE hj ≫ globalChartToProj M hM φ s j := by
  have h1 : ((X.presheaf.map (homOfLE hi).op).hom).comp (globalAwayHom M hM φ s i)
      = (overlapAwayHom M hM φ s i j hi hj).comp
        (HomogeneousLocalization.awayMap
          (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
          (MvPolynomial.isHomogeneous_X R j)
          (rfl : (MvPolynomial.X i * MvPolynomial.X j : MvPolynomial (Fin (N + 1)) R)
            = MvPolynomial.X i * MvPolynomial.X j)) :=
    RingHom.ext (fun z => (overlapAwayHom_awayMap M hM φ s i j hi hj z).symm)
  have h2 : ((X.presheaf.map (homOfLE hj).op).hom).comp (globalAwayHom M hM φ s j)
      = (overlapAwayHom M hM φ s i j hi hj).comp
        (HomogeneousLocalization.awayMap
          (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
          (MvPolynomial.isHomogeneous_X R i)
          (mul_comm (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R) (MvPolynomial.X j))) :=
    RingHom.ext (fun z => (overlapAwayHom_awayMap_symm M hM φ s i j hi hj z).symm)
  rw [globalChartToProj, globalChartToProj, ← Category.assoc, ← Category.assoc,
    homOfLE_comp_globalChartMorphism, homOfLE_comp_globalChartMorphism,
    h1, h2, Category.assoc, Category.assoc]
  rw [show CommRingCat.ofHom (W.topIso.inv.hom.comp
        ((overlapAwayHom M hM φ s i j hi hj).comp
          (HomogeneousLocalization.awayMap _ (MvPolynomial.isHomogeneous_X R j) rfl)))
      = CommRingCat.ofHom (HomogeneousLocalization.awayMap
          (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
          (MvPolynomial.isHomogeneous_X R j) rfl)
        ≫ CommRingCat.ofHom (W.topIso.inv.hom.comp (overlapAwayHom M hM φ s i j hi hj)) from rfl,
    show CommRingCat.ofHom (W.topIso.inv.hom.comp
        ((overlapAwayHom M hM φ s i j hi hj).comp
          (HomogeneousLocalization.awayMap _ (MvPolynomial.isHomogeneous_X R i) (mul_comm _ _))))
      = CommRingCat.ofHom (HomogeneousLocalization.awayMap
          (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
          (MvPolynomial.isHomogeneous_X R i)
          (mul_comm (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R) (MvPolynomial.X j)))
        ≫ CommRingCat.ofHom (W.topIso.inv.hom.comp (overlapAwayHom M hM φ s i j hi hj)) from rfl]
  rw [Spec.map_comp, Spec.map_comp, Category.assoc, Category.assoc,
    Proj.SpecMap_awayMap_awayι, Proj.SpecMap_awayMap_awayι]

/-! ## ★★★★★★★★★★★★★★★★段 3 —— 貼る -/

/-- ★★★★★★★★★★★★★★★★**大域の射 `ψ : X ⟶ ℙᴺ_R`**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★★これが原文の「`L_ℚ` の或る正冪が埋め込みを与える」の**射の側**である
——切断 `s_0, …, s_N` が `X` を覆うなら（`⨆_i X_{s_i} = ⊤`）、
それらは `X` から射影空間への射を定める。

★★埋め込みであること（閉埋め込み・単射）は
`§9-848` の `isImmersion_globalChartToProj` が**チャートごとに**与えている。 -/
noncomputable def globalToProj (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ i, nonVanishing M (s i)) = ⊤) :
    X ⟶ projSpace N R :=
  glueOpens (fun i => nonVanishing M (s i)) hcov
    (fun i => globalChartToProj M hM φ s i)
    (fun i j => globalChartToProj_agree M hM φ s i j _ _)

/-- ★★★★★★**貼った射はチャートの上で元の射に戻る**。 -/
theorem ι_globalToProj (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ i, nonVanishing M (s i)) = ⊤) (i : Fin (N + 1)) :
    (nonVanishing M (s i)).ι ≫ globalToProj M hM φ s hcov
      = globalChartToProj M hM φ s i :=
  ι_glueOpens _ hcov _ _ i

/-! ## ★出典の紐付け(`.src`) -/

def topIso_inv_comp_appTop.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(topIso と制限は交換する)",
    sectionId := "genell-prop-1-4" }

def homOfLE_comp_specMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(Spec への射を小さい開へ制限する)",
    sectionId := "genell-prop-1-4" }

def homOfLE_comp_globalChartMorphism.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャート射の制限)",
    sectionId := "genell-prop-1-4" }

def globalChartToProj_agree.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E1d——チャート射は重なりで一致する)",
    sectionId := "genell-prop-1-4" }

def globalToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(大域の射 ψ : X ⟶ ℙᴺ_R)",
    sectionId := "genell-prop-1-4" }

def ι_globalToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(貼った射はチャートの上で元の射に戻る)",
    sectionId := "genell-prop-1-4" }

def globalToProj.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "globalChartToProj(チャートごとの射、§9-848)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalChartToProj") 2,
    .citation "[ABC3]" "glueOpens(開被覆に沿って射を貼る、§9-836)"
      (.inProject "ABC3" "ABC3.Found.GenEll.glueOpens") 2,
    .citation "[ABC3]" "overlapAwayHom_awayMap_symm(重なりの環準同型の対称性、§9-910)"
      (.inProject "ABC3" "ABC3.Found.GenEll.overlapAwayHom_awayMap_symm") 3,
    .citation "[mathlib]" "Proj.SpecMap_awayMap_awayι(awayι を x_i·x_j へ落とす)"
      (.inMathlib "AlgebraicGeometry.Proj.SpecMap_awayMap_awayι") 1,
    .citation "[mathlib]" "Scheme.toSpecΓ_naturality"
      (.inMathlib "AlgebraicGeometry.Scheme.toSpecΓ_naturality") 1,
    .implicitStep
      ("★★★★測定: topIso は eqToHom の言い換えである" ++
       "(mathlib の Scheme.Opens.topIso は X.presheaf.mapIso (eqToIso U.ι_image_top.symm).op)。" ++
       "したがって V.topIso.inv ≫ (X.homOfLE h).appTop = (制限) ≫ W.topIso.inv は" ++
       "両辺とも X.presheaf.map の合成であり、Opens の射が Subsingleton なので " ++
       "congr 1 だけで閉じる") 2,
    .implicitStep
      ("★★★★★測定: instances 透明度で 4 回はまった" ++
       "——X.presheaf を (Opens X)ᵒᵖ ⥤ CommRingCat に落とそうとして型が合わなくなる。" ++
       "set_option backward.isDefEq.respectTransparency false in を rw を使う定理すべてに付けた" ++
       "(tools/lean-idioms.md の既知の穴)") 2,
    .implicitStep
      ("★★残るのは (1) 被覆条件 ⨆_i X_{s_i} = ⊤(段 E2)を ample から出すことと、" ++
       "(2) 貼った射 ψ が埋め込みであること" ++
       "(チャートごとには §9-848 の isImmersion_globalChartToProj がある)") 5 ]

end ABC3.Found.GenEll
