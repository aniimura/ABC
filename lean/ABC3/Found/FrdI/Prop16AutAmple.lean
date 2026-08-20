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

/-! ## ★★(vi) の `Aut^sub-ample`

原文 (FrdI p.28):
> (vi) A object of C is Aut-ample (respectively, Autsub-ample; End-ample) if

★`Aut-ample` と `End-ample` は **`Prop16.lean` で閉じている**
(`cfp_autAmple_of` / `cfp_endAmple_of`)。残っていたのは `Aut^sub-ample` である。

## ★★★障壁の正体(2026-08-20 に式で特定した)

`𝒞'` の sub-automorphism の**証人**は 3 組である:

* (L) `β_C ≫ χ = χ ≫ φ₀`(`𝒞` 側の証人、`β_C` は同型)
* (R) `β_D ≫ ψ = ψ ≫ g`(`𝒟'` 側の証人、`β_D` は同型)
* (S) 両者の**四角形**。すなわち `Base B ≅ G(e)` で、
  `Base χ` が `G(ψ)` に、`Base β_C` が `G(β_D)` に対応すること

★★**`Aut^sub-ample` は (L) の証人を 1 つ与えるだけで、その底が `G` の像に入るとは
言わない**。引き戻しで作り直そうとすると、普遍性が与える一意の持ち上げ `γ₀` は
`Div γ₀ = Φ(Base χ)(Div φ₀)` を持ち、**これが 0 になる保証が無い**
(`SubAutInvariants.lean` の `div_relation_of_isSubAutomorphism`)。

## ★★★★抜け道 —— `𝒟'` が Aut-saturated なら閉じる

★`𝒟'` の sub-automorphism が**すべて同型**なら(`IsAutSaturatedObj`)、
`Aut^sub_{𝒟'}(a) = Aut_{𝒟'}(a)` なので `Aut-ample` から直に出る。
★★**これは応用では真である** —— `Example 6.1` / `Example 6.3` の底の圈
`B(G)⁰`(`Sec6GaloisCat.lean` の `(FinSub K K̄)ᵒᵖ`)では、連結対象の
自己射がすべて同型だからである。

## ★逸脱(記録)

★原文の (vi) は「`Aut^sub-ample` な対象に射影するなら」と言うが、
下の定理は代わりに **`Aut-ample`** を仮定する。
★`𝒟'` が Aut-saturated のとき、結論が要求するのは**同型の持ち上げ**なので、
この強めの仮定が自然である(`Aut^sub-ample` だけでは
持ち上げた `φ` が同型とは限らず、上の障壁に戻る)。 -/

/-- ★★★★**[FrdI] Proposition 1.6, (vi)** の `Aut^sub-ample` ——
`𝒟'` の sub-automorphism がすべて同型なら、`Aut-ample` から降りる。 -/
theorem cfp_autSubAmple_of_autSaturated (A : CfpCat P Gf)
    (hsat : ∀ (d : D') (α : End d), IsSubAutomorphism α → IsIso α)
    (h : IsAutAmple P A.obj.left) :
    IsAutSubAmple (cfpPreFrobenioid P Gf hGf hD' hcC hcD') A := by
  intro g0 hg0
  obtain ⟨φ, hiso, hb⟩ :=
    cfp_autAmple_of P Gf hGf hD' hcC hcD' A h g0 (hsat _ g0 hg0)
  exact ⟨φ, isSubAutomorphism_of_isIso φ hiso, hb⟩

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Proposition 1.6, (v)` の `base-trivial` の条(isotropic 型の下)。 -/
def cfp_baseTrivial_mpr_of_isotropic.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 28,
    item := "Proposition 1.6, (v) — base-trivial の ⟸(isotropic 型)",
    sectionId := "frdi-prop-1-6-v" }

/-- ★locator —— `Proposition 1.6, (vi)` の `Aut^sub-ample`
(★**条つき**: `𝒟'` が Aut-saturated 、かつ `Aut-ample` を仮定)。 -/
def cfp_autSubAmple_of_autSaturated.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 28,
    item := "Proposition 1.6, (vi) — Aut^sub-ample",
    sectionId := "frdi-prop-1-6-v" }

end ABC3.Found.FrdI
