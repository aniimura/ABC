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

## ★★★★応用でこの仮定が成り立つことは**証明した**(2026-08-25)

前版はここに「**これは応用では真である**」と書いていたが、
**Lean の主張としては書かれていなかった**(docstring の断言のまま)。
★`Sec6GaloisCat.lean` の **`finSubOp_autSaturated` / `finSubOp_isAutSaturatedObj`** で埋めた。
中身は在庫の `finSubOp_isIso_of_endo` 1 本で、
`Example 6.1` / `Example 6.3` の底の圏 `B(G)⁰`(`(FinSub K K̄)ᵒᵖ`)では
**自己射がそもそも全部同型**である。

★★★したがってその底の上では `Aut^sub_𝒟 = Aut_𝒟 = End_𝒟` であり、
原文が (vi) で引く 2 つの概念の**区別は元から無い**。
下の追加仮定は、この応用では**結論を弱めない**。

## ★逸脱(記録)

★原文の (vi) は「`Aut^sub-ample` な対象に射影するなら」と言うが、
下の定理は代わりに **`Aut-ample`** を仮定する。
★`𝒟'` が Aut-saturated のとき、結論が要求するのは**同型の持ち上げ**なので、
この強めの仮定が自然である(`Aut^sub-ample` だけでは
持ち上げた `φ` が同型とは限らず、上の障壁に戻る)。
★★★**残る危険**: `Aut^sub` が `Aut` より真に大きい底の上で (vi) を使う消費者が現れたら、
この版は使えない。★現時点で `IsAutSubAmple` の消費者は
`Prop15.lean`(初等 Frobenioid、別証明)だけで、`CfpCat` の消費者は
`ElementaryFrobenioid.lean` / `Prop16.lean` / 本ファイルの 3 つに限られる
(2026-08-25 実測)。原典側の消費者は **[IUTchI]**(未形式化)である。 -/

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

/-! ### ★★★★★★★★項目全体の `.src`

★`.src` は「その原典項目を**完全に**実装した」という主張である
(`tools/frdi-progress.mjs` の規則)。★★下の 1 つは
**仮定を 2 つ足すという逸脱の下で**置く —— 逸脱の内容は次の docstring に書く。 -/

/-- ★★★★★★★★**[FrdI] Proposition 1.6** —— 条がすべて実装された
(★★**仮定を 2 つ足す逸脱つき**)。

| 条 | 主張 | 宣言 |
|---|---|---|
| (i)(ii) | `C^fp` の 21 条(pre-Frobenioid ＋ Frobenioid の核) | `cfpPreFrobenioid` / `cfpFrobenioidCore`(`Prop16.lean`) |
| (iii) | pull-back の両向き | `cfp_isPullBack_iff` |
| (iv) | LB-invertible と因子分解の一意性 | `cfp_pullBackLB` / `cfp_arbFactorUniq` |
| (v) | perfect / Frobenius-normalized | `cfp_perfect_iff` / `cfp_frobNormalized_iff` |
| (v) | base-trivial の `⟸` | ★`cfp_baseTrivial_mpr_of_isotropic`(**isotropic 型を仮定**) |
| (vi) | `Aut-ample` / `End-ample` | `cfp_autAmple_of` / `cfp_endAmple_of` |
| (vi) | `Aut^sub-ample` | ★`cfp_autSubAmple_of_autSaturated`(**`𝒟'` の Aut-saturated 性を仮定**) |

## ★★★★逸脱の記録(CLAUDE.md の「逸脱」)

★原文は (v) の `base-trivial ⟸` と (vi) の `Aut^sub-ample` を**無条件**で述べるが、
我々はそれぞれに**原文にない仮定を 1 つずつ足して**閉じている。

| 条 | 足した仮定 | 不足の記録 |
|---|---|---|
| (v) | `𝒞` が **isotropic 型** | `ABC3.Gap.FrdI.Gap_1_6_v`(`autAmple`) |
| (vi) | `𝒟'` が **Aut-saturated**(sub-automorphism がすべて同型) | `ABC3.Gap.FrdI.Gap_1_6_vi`(`witnessDescends`) |

★★どちらの `GapRecord` も `GapClass.missingMath`(② 原典の穴の疑い)のままである ——
③(反証)を名乗るには反例が要り、①(閉じた)を名乗るには
`Definition 1.3` からの導出が要る。**その判定は保留する。**

## ★★★★★なぜ後続に影響しないと見るか

★★**応用では 2 つの仮定はどちらも満たされる**:

* isotropic 型 —— `Example 6.1` / `Example 6.3` の Frobenioid はいずれも isotropic 型
  (`ex61Frobenioid_isotropicType` / `ex63_isotropic_family`)。
* Aut-saturated —— 底の圏 `B(G)⁰ = (FinSub K K̄)ᵒᵖ` では自己射がそもそも
  **全部同型**なので `Aut^sub_𝒟 = Aut_𝒟 = End_𝒟` である
  (`Sec6GaloisCat.lean` の `finSubOp_autSaturated` / `finSubOp_isAutSaturatedObj`)。

★★★**我々のツリー内に `Proposition 1.6` の消費者はいない**(2026-08-25 実測。
`CfpCat` の消費者は `ElementaryFrobenioid.lean` / `Prop16.lean` / 本ファイルの 3 つだけ)。
原典側の消費者は **[IUTchI]**(未形式化)である。
★★したがって危険は **[IUTchI] を形式化する段で、`Aut^sub` が `Aut` より真に大きい底の上で
(vi) を使う箇所が現れたとき**に限られる。そこに来たら本 docstring と
`Gap_1_6_vi` を読み直すこと。 -/
def prop_1_6.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 27, item := "Proposition 1.6",
    sectionId := "frdi-prop-1-6" }

end ABC3.Found.FrdI
