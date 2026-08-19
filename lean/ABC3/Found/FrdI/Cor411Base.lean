/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop113
import ABC3.Found.FrdI.Def45

/-!
# [FrdI] Corollary 4.11 —— **Div-slim による剛性**

原文 (FrdI p.93):
> Thus, since Di is Div-slim, it follows that every automorphism [of an object of

## ★★★★★測って分かった —— `Proposition 1.13, (i)` の骨がそのまま使える

在庫の `isRigidFunctor_comp_of_equivalence`(`Prop113.lean`)は
「`G` が rigid ⟹ `e.functor ⋙ G` が rigid」を、
**共役 ＋ `e.unit` に沿った運搬**で示している。

★★`Corollary 4.11` が要るのは **rigid ではなく「`Φ` が恒等へ送るなら恒等」**である
(`𝒟` は slim とは限らず **Div-slim** しか仮定しない)。
★そこで骨だけを取り出して**条件つきの形**に一般化する
(`eq_refl_of_conj_eq_refl`)——`hG` を仮定として受け取るだけの違いである。

## ★★`overPhiAut` の成分
`overPhiAut Φ A η = isoWhiskerRight (NatIso.op η) Φ.functor` の
`Z : (𝒟_A)ᵒᵖ` での成分は `Φ.functor.map ((η.hom.app Z.unop).op)` であり、
その `.hom` は **`Φ.map (η.hom.app Z.unop)`** そのものである。
★したがって「`Φ` が恒等へ送る」は `∀ Y, Φ.map (η.hom.app Y) = id` と書ける。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

/-! ## ★1. 圏論の骨(条件つき版) -/

/-- ★★★**`isRigidFunctor_comp_of_equivalence` の条件つき版** ——
`e.functor ⋙ G` の自己同型は、`G ≅ G` へ共役したものが恒等なら恒等。

★在庫の証明の骨をそのまま使う(共役 ＋ `e.unit` に沿った運搬)。 -/
theorem eq_refl_of_conj_eq_refl {A₁ : Type u} [Category.{v} A₁]
    {B₁ : Type u2} [Category.{v2} B₁] {E₁ : Type u3} [Category.{v3} E₁]
    (e : A₁ ≌ B₁) {G : B₁ ⥤ E₁} (η : e.functor ⋙ G ≅ e.functor ⋙ G)
    (hG : ((e.invFunIdAssoc G).symm ≪≫ Functor.isoWhiskerLeft e.inverse η
      ≪≫ (e.invFunIdAssoc G)) = Iso.refl G) :
    η = Iso.refl _ := by
  set X : e.inverse ⋙ e.functor ⋙ G ≅ G := e.invFunIdAssoc G with hX
  have key : ∀ b : B₁,
      η.hom.app (e.inverse.obj b) = 𝟙 ((e.functor ⋙ G).obj (e.inverse.obj b)) := by
    intro b
    have h : X.inv.app b ≫ η.hom.app (e.inverse.obj b) ≫ X.hom.app b = 𝟙 (G.obj b) :=
      congrArg (fun t : G ≅ G => t.hom.app b) hG
    have h2 : η.hom.app (e.inverse.obj b) ≫ X.hom.app b = X.hom.app b := by
      have h3 := congrArg (fun t => X.hom.app b ≫ t) h
      simp only [Iso.hom_inv_id_app_assoc, Category.comp_id] at h3
      exact h3
    have hXi : IsIso (X.hom.app b) :=
      ⟨⟨X.inv.app b, X.hom_inv_id_app b, X.inv_hom_id_app b⟩⟩
    exact eq_id_of_comp_eq_of_isIso _ (X.hom.app b) hXi h2
  apply Iso.ext
  apply NatTrans.ext
  funext a
  show η.hom.app a = 𝟙 ((e.functor ⋙ G).obj a)
  set u : (e.functor ⋙ G).obj (e.inverse.obj (e.functor.obj a)) ⟶ (e.functor ⋙ G).obj a :=
    (e.functor ⋙ G).map (e.unitInv.app a) with hu
  set w : (e.functor ⋙ G).obj a ⟶ (e.functor ⋙ G).obj (e.inverse.obj (e.functor.obj a)) :=
    (e.functor ⋙ G).map (e.unit.app a) with hw
  have hwu : w ≫ u = 𝟙 ((e.functor ⋙ G).obj a) := by
    rw [hw, hu, ← Functor.map_comp]
    simp
    rfl
  have hnat : u ≫ η.hom.app a = η.hom.app (e.inverse.obj (e.functor.obj a)) ≫ u :=
    η.hom.naturality (e.unitInv.app a)
  rw [key (e.functor.obj a), Category.id_comp] at hnat
  calc η.hom.app a = (w ≫ u) ≫ η.hom.app a := by rw [hwu, Category.id_comp]
    _ = w ≫ u ≫ η.hom.app a := by rw [Category.assoc]
    _ = w ≫ u := by rw [hnat]
    _ = 𝟙 _ := hwu

/-! ## ★2. `overPhiAut` の成分 -/

section DivSlim

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

variable (Φ) in
/-- ★`overPhiAut` は恒等を恒等へ送る。 -/
theorem overPhiAut_one (A : D) :
    overPhiAut Φ A (1 : Aut (Over.forget A)) = 1 := by
  apply Iso.ext
  apply NatTrans.ext
  funext Z
  show Φ.functor.map ((𝟙 (Z.unop.left)).op) = 𝟙 _
  rw [op_id, Φ.functor.map_id]

variable (Φ) in
/-- ★★**`Φ` がすべての成分を恒等へ送るなら `overPhiAut` の像は恒等**。

★成分は `Φ.functor.map ((η.hom.app Y).op)` であり、その `.hom` は `Φ.map (η.hom.app Y)`。 -/
theorem overPhiAut_eq_one {A : D} (η : Aut (Over.forget A))
    (h : ∀ (Y : Over A) (x : Φ.val Y.left), Φ.map (η.hom.app Y) x = x) :
    overPhiAut Φ A η = 1 := by
  apply Iso.ext
  apply NatTrans.ext
  funext Z
  show Φ.functor.map ((η.hom.app Z.unop).op) = 𝟙 _
  exact AddCommMonCat.hom_ext (AddMonoidHom.ext fun x => h Z.unop x)

/-! ## ★3. Div-slim 版の `Proposition 1.13, (i)` -/

variable (P : PreFrobenioid C Φ)

include P in
/-- ★★★★★**Div-slim ⟹ pull-back の対象の上で自明**。

★在庫の `prop_1_13_i_pullBack_trivial` は `IsSlimCat D` を使うが、
ここは **`IsDivSlim Φ` ＋「`Φ` が恒等へ送る」**だけで済ませる。

★★手筋は同じ ——`Definition 1.3, (i), (c)` の圏同値 `(𝒞^pl-bk)_A ≌ 𝒟_{A_𝒟}` で
`𝒟_{A_𝒟} → 𝒟` の自己同型へ共役し、そこで **Div-slim の単射性**を当てる。 -/
theorem divSlim_pullBack_trivial (F : FrobenioidCore P) {Φ₀ : MonoidOn.{v, u, w} D}
    (hds : IsDivSlim Φ₀) (A : C)
    (η : (Over.forget A ⋙ P.proj) ≅ (Over.forget A ⋙ P.proj))
    (hdiv : ∀ (Y : Over A) (x : Φ₀.val (P.proj.obj Y.left)),
      Φ₀.map (η.hom.app Y) x = x)
    (Y : Over A) (hpb : IsPullBack P Y.hom) :
    η.hom.app Y = 𝟙 ((Over.forget A ⋙ P.proj).obj Y) := by
  haveI := F.plBkEquiv A
  set e := (plBkOverFunctor P A).asEquivalence with he
  set θ : (plBkOverFunctor P A ⋙ Over.forget ((P.toElem.obj A).base))
      ≅ (plBkOverFunctor P A ⋙ Over.forget ((P.toElem.obj A).base)) :=
    Functor.isoWhiskerLeft (plBkToOver P A) η with hθ
  set X : e.inverse ⋙ e.functor ⋙ Over.forget ((P.toElem.obj A).base)
      ≅ Over.forget ((P.toElem.obj A).base) :=
    e.invFunIdAssoc (Over.forget ((P.toElem.obj A).base)) with hX
  -- ★共役した自己同型を `Φ` が恒等へ送ることを見る
  have hconj : overPhiAut Φ₀ ((P.toElem.obj A).base)
      (X.symm ≪≫ Functor.isoWhiskerLeft e.inverse θ ≪≫ X) = 1 := by
    refine overPhiAut_eq_one Φ₀ _ (fun Z x => ?_)
    show Φ₀.map (X.inv.app Z ≫ θ.hom.app (e.inverse.obj Z) ≫ X.hom.app Z) x = x
    have hkey := hdiv ((plBkToOver P A).obj (e.inverse.obj Z)) (Φ₀.map (X.hom.app Z) x)
    calc Φ₀.map (X.inv.app Z ≫ θ.hom.app (e.inverse.obj Z) ≫ X.hom.app Z) x
        = Φ₀.map (X.inv.app Z) (Φ₀.map (θ.hom.app (e.inverse.obj Z))
            (Φ₀.map (X.hom.app Z) x)) := by
            rw [Φ₀.map_comp, Φ₀.map_comp]
            rfl
      _ = Φ₀.map (X.inv.app Z) (Φ₀.map (X.hom.app Z) x) :=
            congrArg (fun t => Φ₀.map (X.inv.app Z) t) hkey
      _ = Φ₀.map (X.inv.app Z ≫ X.hom.app Z) x :=
            (Φ₀.map_comp (X.hom.app Z) (X.inv.app Z) x).symm
      _ = x := by rw [X.inv_hom_id_app, Φ₀.map_id]
  have h1 : (X.symm ≪≫ Functor.isoWhiskerLeft e.inverse θ ≪≫ X)
      = (1 : Aut (Over.forget ((P.toElem.obj A).base))) :=
    hds _ (hconj.trans (overPhiAut_one Φ₀ _).symm)
  have hθrefl : θ = Iso.refl _ := eq_refl_of_conj_eq_refl e θ h1
  exact congrArg (fun t => t.hom.app (Over.mk (⟨Y.hom, hpb⟩ : (⟨Y.left⟩ : PlBk P) ⟶ ⟨A⟩)))
    hθrefl

include P in
/-- ★★★★★★**Div-slim 版の `Proposition 1.13, (i)`** ——
`𝒞_A → 𝒟` の自己同型で「`Φ` が恒等へ送る」ものは**恒等**。

原文 (FrdI p.93):
> Thus, since Di is Div-slim, it follows that every automorphism [of an object of -/
theorem divSlim_over_aut_eq_id (F : FrobenioidCore P) {Φ₀ : MonoidOn.{v, u, w} D}
    (hds : IsDivSlim Φ₀) (A : C)
    (η : (Over.forget A ⋙ P.proj) ≅ (Over.forget A ⋙ P.proj))
    (hdiv : ∀ (Y : Over A) (x : Φ₀.val (P.proj.obj Y.left)),
      Φ₀.map (η.hom.app Y) x = x) :
    η = Iso.refl _ :=
  prop_1_13_i_from_pullBack P F A η
    (fun Y hpb => divSlim_pullBack_trivial P F hds A η hdiv Y hpb)

def divSlim_over_aut_eq_id.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (i) — Div-slim による剛性",
    sectionId := "frdi-cor-4-11" }

end DivSlim

end ABC3.Found.FrdI
