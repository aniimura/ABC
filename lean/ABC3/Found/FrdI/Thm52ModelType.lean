/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Model
import ABC3.Found.FrdI.Prop33Coa

/-!
# [FrdI] Theorem 5.2, (iv) の使い方 —— model 型ならば model Frobenioid

原文 (FrdI p.101):
> that is 1-compatible with the functors C → FΦ, C → FΦ.

★`thm_5_2_iv`(`Thm52Model.lean`)は base-Frobenius pair `(𝒫, F)` を**引数に取る**形で
書いてある。★原文の仮定は「`𝒞` が model 型」なので、そこから pair を取り出す形に
まとめ直しておく(`modelType_equiv`)。

★★あわせて **`𝒞^un-tr` はつねに model 型**(`Proposition 5.5, (iii)` の一節)を
`Theorem 5.1, (iv)` から出す —— `𝒞^un-tr` は isotropic かつ unit-trivial だからである。

原文 (FrdI p.104):
> (iii) If C is of standard (respectively, rationally standard; model) type, then so is Cpf.
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} [IsConnected D]

/-- ★★★★★**`Theorem 5.2, (iv)` を原文の仮定の形で** ——
`𝒞` が **model 型**で、有理関数の単系 `B` の interface が与えられていれば、
`𝒞` は model Frobenioid と圏同値。

★`IsOfModelType` の第 1 成分(pre-model 型)から base-Frobenius pair を取り出し、
第 2 成分(birationally Frobenius-normalized 型)をそのまま `hfn` に渡すだけ。 -/
noncomputable def modelType_equiv {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y) (hm : IsOfModelType C P G) :
    C ≌ ModelData.Obj R.model :=
  thm_5_2_iv R hiso (fun Z => hm.2 Z) (Classical.choice hm.1).sec
    (Classical.choice hm.1).frob (Classical.choice hm.1).isFrobSection

/-- ★★★**[FrdI] Proposition 5.5, (iii) の一節** —— `𝒞^un-tr` は**つねに model 型**。

原文 (FrdI p.104):
> (iii) If C is of standard (respectively, rationally standard; model) type, then so is Cpf.

★`𝒞^un-tr` は isotropic(`unTr_isotropic`)かつ unit-trivial(`unTr_unitTrivial`)なので、
`Theorem 5.1, (iv)` がそのまま当たる。 -/
theorem unTr_isOfModelType (Fc : FrobenioidCore P) (G : Frobenioid P) :
    IsOfModelType (UnTr P) (unTrPre P Fc) (unTr_frobenioid P Fc G) :=
  thm_5_1_iv (unTr_frobenioid P Fc G) (fun B => unTr_isotropic P Fc B)
    (fun A => unTr_unitTrivial P Fc A)

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Proposition 5.5, (iii)` の「`𝒞^un-tr` はつねに model 型」の条。 -/
def unTr_isOfModelType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (iii) — 𝒞^un-tr はつねに model 型",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
