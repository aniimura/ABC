/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop16
import ABC3.Found.FrdI.Thm51Cls

/-!
# [FrdI] Proposition 1.6, (v) の穴を isotropic 型で埋める

`Gap/FrdI/Section1.lean` の `Gap_1_6_v` は、`Proposition 1.6, (v)` の
`base-trivial` と `metrically trivial` の 2 件の `⟸` に **`Aut-ample`** が
足りない、という記録である。

★★**`Theorem 5.1, (iii)` が `base-trivial` の枝を埋める**:

原文 (FrdI p.99):
> and that all Frobenius-trivial objects of C are Aut-ample. In light of these

すなわち `𝒞` が **isotropic 型**なら、Frobenius-trivial な対象は Aut-ample である
(`isAutAmple_of_frobTrivial`)。そして `base-trivial` な対象は Frobenius-trivial
(`isFrobeniusTrivial_of_baseTrivial`)。★したがって

  isotropic 型 ＋ base-trivial ⟹ Aut-ample

が出て、`cfp_baseTrivial_mpr` の仮定 `haa` が**消える**。

★★**残る穴は `metrically trivial` の枝だけ**である。
そちらは「pre-step 自己射の底が自分の `Div` を動かさない」(`htwist`)に戻ってしまう
—— `Theorem 5.1, (iii)` の議論は `A` が Frobenius-trivial であることを使うが、
metrically trivial からそれを出す道が `htwist` を経由するためである。
★この 1 点が `Gap_1_6_v` の現在地である(`Gap/FrdI/Section1.lean` に記録)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {D' : Type u3} [Category.{v3} D']
  {C : Type u2} [Category.{v2} C] {Φ : MonoidOn.{v, u, w} D}

/-- ★★★★**isotropic 型では、base-trivial な対象は Aut-ample**。

★`Theorem 5.1, (iii)` の帰結。★これが `Gap_1_6_v` の半分を埋める。 -/
theorem isAutAmple_of_baseTrivial {P : PreFrobenioid C Φ} (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y) {A : C} (h : IsBaseTrivial P A) :
    IsAutAmple P A :=
  isAutAmple_of_frobTrivial G hiso (isFrobeniusTrivial_of_baseTrivial P G.core h)

variable (P : PreFrobenioid C Φ) (Gf : D' ⥤ D)
  (hGf : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (Gf.map α))
  (hD' : IsTotallyEpimorphic D')
  (hcC : IsConnected (CfpCat P Gf)) (hcD' : IsConnected D')

/-- ★★★★**[FrdI] Proposition 1.6, (v)** の `base-trivial` の `⟸`
—— **`Aut-ample` の仮定なしで**(isotropic 型の下で)。

原文 (FrdI p.28):
> (v) A object of C is Frobenius-trivial (respectively, quasi-Frobenius-trivial; sub-quasi-Frobenius-trivial; metrically trivial; base-trivial; perfect; grouplike; unit-trivial; Frobenius-normalized; isotropic; Frobenius-isotropic) if and only if it projects to such an object of C.
-/
theorem cfp_baseTrivial_mpr_of_isotropic (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y) (A : CfpCat P Gf)
    (h : IsBaseTrivial P A.obj.left) :
    IsBaseTrivial (cfpPreFrobenioid P Gf hGf hD' hcC hcD') A :=
  cfp_baseTrivial_mpr P Gf hGf hD' hcC hcD' A (isAutAmple_of_baseTrivial G hiso h) h

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Proposition 1.6, (v)` の `base-trivial` の条(isotropic 型の下)。 -/
def cfp_baseTrivial_mpr_of_isotropic.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 28,
    item := "Proposition 1.6, (v) — base-trivial の ⟸(isotropic 型)",
    sectionId := "frdi-prop-1-6" }

end ABC3.Found.FrdI
