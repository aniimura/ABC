/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Frob

/-!
# [FrdI] `PfCat`(根 1 の部分)から `PfRootObj`(原文の `𝒞^pf`)への関手

原文 (FrdI p.57):
> Also, we obtain a natural functor C → Cpf [by mapping

★★我々の木では `𝒞^pf` を 2 段で作ってある:

* `PfCat P F` —— **根 1 の部分**(対象は `𝒞` の対象、射は `Hom^pf`)
* `PfRootObj P F` —— **原文の `𝒞^pf`**(対象は対 `(A,n)`)

★`𝒞 → PfCat`(`toPfCat`)と `𝒞 → PfRootObj`(`toPfRoot`)は在庫にあるが、
**`PfCat → PfRootObj`** が無かった。★`Proposition 5.3` の図式で
`𝒞^pf` から `𝒞^rlf` へ渡るときに要る(model Frobenioid との圏同値は
`PfRootObj` の側で立ててあるため)。

## ★中身

`rtExt A 1 : A ⟶ A^{(1)}` は**同型**(`isIso_rtExt_one`)なので、
`Hom^pf(A,B)` を `Hom^pf(A^{(1)}, B^{(1)}) = Hom_{𝒞^pf}((A,1),(B,1))` へ
**共役で移す**だけである。

★★関手性は「根 1 どうしでは `compRoot` は `compPf` そのもの」
(`compRoot_root_one`、`compRoot_eq_lift` ＋ `rtRootIso_*_eq_self`)に落ちる。

-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

section PfCatRoot

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-- ★★**根 1 の対象どうしでは `compRoot` は `compPf` そのもの**。

★`compRoot_eq_lift` を `c = PA = PB = PE = 1` で当て、
`rtRootIso_hom_eq_self` / `rtRootIso_inv_eq_self`(在庫)で 3 つの輸送を消す。 -/
theorem compRoot_root_one {A B E : C}
    (f : HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩)
    (g : HomRoot P F (⟨B, 1⟩ : PfRootObj P F) ⟨E, 1⟩) :
    compRoot P F f g = compPf P F f g := by
  have hlift := compRoot_eq_lift (P := P) (F := F) f g
    (c := 1) (PA := 1) (PB := 1) (PE := 1) (hcA := rfl) (hcB := rfl) (hcE := rfl)
    (ef := 1) (eg := 1) (er := 1) (hfA := rfl) (hfB := rfl) (hgA := rfl) (hgE := rfl)
    (hrA := rfl) (hrE := rfl)
  rw [hlift, rtRootIso_inv_eq_self, rtRootIso_inv_eq_self, rtRootIso_hom_eq_self]

/-- ★★共役が同型なら元も同型(任意の圏での代数)。

★★**抽象の対象で述べておく**のが要点 —— 型の同義語を跨ぐ具体形で
`simp` / `infer_instance` を回すと `instances` 透明度で落ちる
(`lean-idioms.md` の 14)。 -/
theorem isIso_of_conj {X : Type u3} [Category.{v3} X] {a b a' b' : X}
    (i : a' ⟶ a) [IsIso i] (φ : a ⟶ b) (j : b ⟶ b') [IsIso j]
    (h : IsIso (i ≫ φ ≫ j)) : IsIso φ := by
  haveI := h
  have hfac : φ = inv i ≫ (i ≫ φ ≫ j) ≫ inv j := by simp
  rw [hfac]
  infer_instance

/-- ★共役の同型性(3 つ組)。 -/
theorem isIso_comp₃ {X : Type u3} [Category.{v3} X] {a b c d : X}
    (f : a ⟶ b) [IsIso f] (g : b ⟶ c) [IsIso g] (k : c ⟶ d) [IsIso k] :
    IsIso (f ≫ g ≫ k) :=
  inferInstance

variable {P F} in
/-- ★`PfCat` の合成は `compPf`(`rw` のための橋)。 -/
theorem pfCat_comp_eq {A B E : PfCat P F} (f : A ⟶ B) (g : B ⟶ E) :
    f ≫ g = compPf P F f g := rfl

variable {P F} in
/-- ★`toPfCat` の射の部分は `toHomPf`(`rw` のための橋)。 -/
theorem toPfCat_map_eq {A B : C} (ψ : A ⟶ B) :
    (toPfCat P F).map ψ = toHomPf (F := F) ψ := rfl

/-- ★共役の合成(任意の圏での代数)—— `j ≫ i' = 𝟙` なら真ん中が消える。 -/
theorem conj_comp_aux {X : Type u2} [Category.{v2} X] {a b c a' b' c' : X}
    (i : a' ⟶ a) (f : a ⟶ b) (j : b ⟶ b') (i' : b' ⟶ b) (g : b ⟶ c) (k : c ⟶ c')
    (h : j ≫ i' = 𝟙 b) :
    i ≫ (f ≫ g) ≫ k = (i ≫ f ≫ j) ≫ (i' ≫ g ≫ k) := by
  simp only [Category.assoc]
  rw [← Category.assoc j, h, Category.id_comp]

/-- ★★★★**`PfCat P F ⥤ PfRootObj P F`** —— 根 1 の部分を原文の `𝒞^pf` へ入れる。

★`rtExt A 1` が同型なので、共役で移すだけである。 -/
noncomputable def pfCatToRoot : PfCat P F ⥤ PfRootObj P F where
  obj A := ⟨(pfDown P F A), 1⟩
  map {A B} f :=
    haveI := isIso_rtExt_one P F (pfDown P F A)
    (show HomRoot P F (⟨(pfDown P F A), 1⟩ : PfRootObj P F) ⟨(pfDown P F B), 1⟩ from
      (toPfCat P F).map (inv (rtExt P F (pfDown P F A) 1)) ≫ f
        ≫ (toPfCat P F).map (rtExt P F (pfDown P F B) 1))
  map_id A := by
    haveI := isIso_rtExt_one P F (pfDown P F A)
    show (toPfCat P F).map (inv (rtExt P F (pfDown P F A) 1))
        ≫ 𝟙 _ ≫ (toPfCat P F).map (rtExt P F (pfDown P F A) 1)
      = idRoot P F (⟨(pfDown P F A), 1⟩ : PfRootObj P F)
    rw [Category.id_comp, ← (toPfCat P F).map_comp, IsIso.inv_hom_id, (toPfCat P F).map_id]
    rfl
  map_comp {A B E} f g := by
    haveI := isIso_rtExt_one P F (pfDown P F A)
    haveI := isIso_rtExt_one P F (pfDown P F B)
    haveI := isIso_rtExt_one P F (pfDown P F E)
    refine Eq.trans ?_ (compRoot_root_one P F
      ((toPfCat P F).map (inv (rtExt P F (pfDown P F A) 1)) ≫ f
        ≫ (toPfCat P F).map (rtExt P F (pfDown P F B) 1))
      ((toPfCat P F).map (inv (rtExt P F (pfDown P F B) 1)) ≫ g
        ≫ (toPfCat P F).map (rtExt P F (pfDown P F E) 1))).symm
    have hcancel : (toPfCat P F).map (rtExt P F (pfDown P F B) 1)
        ≫ (toPfCat P F).map (inv (rtExt P F (pfDown P F B) 1))
      = 𝟙 B := by
      rw [← (toPfCat P F).map_comp, IsIso.hom_inv_id, (toPfCat P F).map_id]
      rfl
    exact conj_comp_aux _ f _ _ g _ hcancel

variable {P F} in
/-- ★`𝒞 → PfCat → PfRootObj` は `𝒞 → PfRootObj` に等しい(射の側)。 -/
theorem toPfCat_comp_pfCatToRoot_map {A B : C} (φ : A ⟶ B) :
    (toPfCat P F ⋙ pfCatToRoot P F).map φ = (toPfRoot P F).map φ := by
  haveI := isIso_rtExt_one P F A
  haveI := isIso_rtExt_one P F B
  show (toPfCat P F).map (inv (rtExt P F A 1)) ≫ (toPfCat P F).map φ
      ≫ (toPfCat P F).map (rtExt P F B 1)
    = toRootHom (F := F) φ
  rw [← (toPfCat P F).map_comp, ← (toPfCat P F).map_comp]
  rfl


/-- ★★★★**`pfCatToRoot` は射の上で全単射**(同型による共役だから)。

★★`Iso.homCongr`(`f ↦ α⁻¹ ≫ f ≫ β`)が**定義的にそのもの**なので 1 行で出る。
★これが `Proposition 5.5, (ii)` で `Θ.map` の全単射性を根 1 で出すのに要る。 -/
theorem pfCatToRoot_map_bijective (X Y : PfCat P F) :
    Function.Bijective ((pfCatToRoot P F).map (X := X) (Y := Y)) := by
  haveI h1 : IsIso (rtExt P F (pfDown P F X) 1) := isIso_rtExt_one P F _
  haveI h2 : IsIso (rtExt P F (pfDown P F Y) 1) := isIso_rtExt_one P F _
  exact (Iso.homCongr
    ((toPfCat P F).mapIso (@asIso _ _ _ _ (rtExt P F (pfDown P F X) 1) h1))
    ((toPfCat P F).mapIso (@asIso _ _ _ _ (rtExt P F (pfDown P F Y) 1) h2))).bijective

/-- ★★★★locator —— `pfCatToRoot` の射の全単射性。 -/
def pfCatToRoot_map_bijective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 57,
    item := "Definition 3.1, (iii) — 根 1 の埋め込みは射の上で全単射",
    sectionId := "frdi-def-3-1" }


/-! ## ★2. `𝒞^pf` の isotropic 性を根 1 の部分の言葉で


★★型の同義語 `PfCat P F := C` を跨ぐので、**在庫の `pfDown` に揃える**こと
(自前の同義語ほどき関数を書くと `pfDiv` などの暗黙引数と噛み合わない)、
**目標の側で `rw` せず自分で述べた `have` の側で `rw` する**こと、
**同型性の代数は抽象の圏で述べた補題に投げる**こと —— の 3 点が要点である
(`lean-idioms.md` の 14)。 -/

/-- ★`rtExt A 1` の逆射も Frobenius 型(同型だから)。 -/
theorem isFrobeniusType_inv_rtExt_one (A : C) :
    haveI := isIso_rtExt_one P F A
    IsFrobeniusType P (inv (rtExt P F A 1)) :=
  haveI := isIso_rtExt_one P F A
  isFrobeniusType_of_isIso P (inv (rtExt P F A 1))

set_option maxHeartbeats 1600000 in
variable {P F} in
/-- ★★★★★**`𝒞^pf`(根 1 の部分)も isotropic 型** ——
`Proposition 3.2, (iii)`(`pfRoot_isOfIsotropicType`)を `pfCatToRoot` で降ろしたもの。

原文 (FrdI p.59):
> (i), is a Frobenioid of perfect and isotropic type. Moreover, there is a natural
-/
theorem pf_isOfIsotropicType (hfi : IsOfFrobeniusIsotropicType P) (X : PfCat P F) :
    IsIsotropic (pfPre P F) X := by
  intro Dd φ hisom hstep
  haveI hiX := isIso_rtExt_one P F (pfDown P F X)
  haveI hiD := isIso_rtExt_one P F (pfDown P F Dd)
  haveI hbX : IsIso (P.Base (rtExt P F (pfDown P F X) 1)) :=
    (rtExt_frobType P F (pfDown P F X) 1).2
  haveI hbD : IsIso (P.Base (rtExt P F (pfDown P F Dd) 1)) :=
    (rtExt_frobType P F (pfDown P F Dd) 1).2
  haveI hbXi : IsIso (P.Base (inv (rtExt P F (pfDown P F X) 1))) :=
    (isFrobeniusType_inv_rtExt_one P F (pfDown P F X)).2
  haveI hbφ : IsIso (pfBase φ) := hstep.2
  have hDivD : P.Div (rtExt P F (pfDown P F Dd) 1) = 0 :=
    (rtExt_frobType P F (pfDown P F Dd) 1).1.2
  have hDivX : P.Div (inv (rtExt P F (pfDown P F X) 1)) = 0 :=
    (isFrobeniusType_inv_rtExt_one P F (pfDown P F X)).1.2
  -- ★型の同義語を跨ぐ補題は**この場で具体化**しておく
  have eDivD : pfDiv (toHomPf (F := F) (rtExt P F (pfDown P F Dd) 1)) = 0 := by
    rw [pfDiv_toHomPf, hDivD, Pf.mk_eq_zero_iff]
    exact ⟨1, by simp⟩
  have eDivX : pfDiv (toHomPf (F := F) (inv (rtExt P F (pfDown P F X) 1))) = 0 := by
    rw [pfDiv_toHomPf, hDivX, Pf.mk_eq_zero_iff]
    exact ⟨1, by simp⟩
  have eDegD : pfDeg (toHomPf (F := F) (rtExt P F (pfDown P F Dd) 1)) = 1 := by
    rw [pfDeg_toHomPf]; exact rtExt_degFr P F _ 1
  have eDegX : pfDeg (toHomPf (F := F) (inv (rtExt P F (pfDown P F X) 1))) = 1 := by
    rw [pfDeg_toHomPf]; exact degFr_inv_eq_one _ (rtExt_degFr P F _ 1)
  have eBaseD : pfBase (toHomPf (F := F) (rtExt P F (pfDown P F Dd) 1))
      = P.Base (rtExt P F (pfDown P F Dd) 1) := pfBase_toHomPf _
  have eBaseX : pfBase (toHomPf (F := F) (inv (rtExt P F (pfDown P F X) 1)))
      = P.Base (inv (rtExt P F (pfDown P F X) 1)) := pfBase_toHomPf _
  -- ★★3 つの不変量(**自分で述べた `have` の側で `rw`**)
  have hφj : pfDiv (compPf P F φ (toHomPf (F := F) (rtExt P F (pfDown P F Dd) 1))) = 0 := by
    rw [pfDiv_comp, eDivD, map_zero, show pfDiv φ = 0 from hisom, smul_zero, add_zero]
  have hdiv0 : pfDiv ((pfCatToRoot P F).map φ) = 0 := by
    have h0 : pfDiv (compPf P F (toHomPf (F := F) (inv (rtExt P F (pfDown P F X) 1)))
        (compPf P F φ (toHomPf (F := F) (rtExt P F (pfDown P F Dd) 1)))) = 0 := by
      rw [pfDiv_comp, hφj, map_zero, eDivX, smul_zero, add_zero]
    exact h0
  have hdegφj : pfDeg (compPf P F φ (toHomPf (F := F) (rtExt P F (pfDown P F Dd) 1)))
      = 1 := by
    rw [pfDeg_comp, eDegD, show pfDeg φ = 1 from hstep.1, one_mul]
  have hdeg1 : pfDeg ((pfCatToRoot P F).map φ) = 1 := by
    have h0 : pfDeg (compPf P F (toHomPf (F := F) (inv (rtExt P F (pfDown P F X) 1)))
        (compPf P F φ (toHomPf (F := F) (rtExt P F (pfDown P F Dd) 1)))) = 1 := by
      rw [pfDeg_comp, hdegφj, eDegX, one_mul]
    exact h0
  haveI hbg : IsIso (pfBase ((pfCatToRoot P F).map φ)) := by
    have h0 : pfBase (compPf P F (toHomPf (F := F) (inv (rtExt P F (pfDown P F X) 1)))
        (compPf P F φ (toHomPf (F := F) (rtExt P F (pfDown P F Dd) 1))))
        = P.Base (inv (rtExt P F (pfDown P F X) 1)) ≫ pfBase φ
          ≫ P.Base (rtExt P F (pfDown P F Dd) 1) := by
      rw [pfBase_comp, pfBase_comp, eBaseX, eBaseD]
    have hb : pfBase ((pfCatToRoot P F).map φ)
        = P.Base (inv (rtExt P F (pfDown P F X) 1)) ≫ pfBase φ
          ≫ P.Base (rtExt P F (pfDown P F Dd) 1) := h0
    rw [hb]
    exact @isIso_comp₃ _ _ _ _ _ _ _ hbXi _ hbφ _ hbD
  -- ★根つきの側で isotropic 性を使う
  have hisom' : IsIsometric (pfRootPre P F) ((pfCatToRoot P F).map φ) := by
    show rootDiv ((pfCatToRoot P F).map φ) = 0
    rw [rootDiv, hdiv0, map_zero, Pf.divBy_zero]
  have hstep' : IsPreStep (pfRootPre P F) ((pfCatToRoot P F).map φ) := by
    refine ⟨hdeg1, ?_⟩
    show IsIso (rootBase ((pfCatToRoot P F).map φ))
    haveI e1 : IsIso (P.Base (rtExt P F ((pfCatToRoot P F).obj X).obj
      ((pfCatToRoot P F).obj Dd).root)) := hbX
    haveI e2 : IsIso (P.Base (rtExt P F ((pfCatToRoot P F).obj Dd).obj
      ((pfCatToRoot P F).obj X).root)) := hbD
    rw [rootBase]
    exact isIso_comp₃ _ _ _
  haveI hgiso : IsIso ((pfCatToRoot P F).map φ) :=
    pfRoot_isOfIsotropicType (F := F) hfi _ _ _ hisom' hstep'
  -- ★根つきの同型は根 1 では `PfCat` の同型
  obtain ⟨h, h1, h2⟩ := hgiso.out
  haveI hi' : IsIso ((toPfCat P F).map (inv (rtExt P F (pfDown P F X) 1))) := inferInstance
  haveI hj' : IsIso ((toPfCat P F).map (rtExt P F (pfDown P F Dd) 1)) := by
    haveI := isIso_rtExt_one P F (pfDown P F Dd)
    exact ⟨((toPfCat P F).mapIso (asIso (rtExt P F (pfDown P F Dd) 1))).inv,
      ((toPfCat P F).mapIso (asIso (rtExt P F (pfDown P F Dd) 1))).hom_inv_id,
      ((toPfCat P F).mapIso (asIso (rtExt P F (pfDown P F Dd) 1))).inv_hom_id⟩
  haveI hgiso' : IsIso ((toPfCat P F).map (inv (rtExt P F (pfDown P F X) 1)) ≫ φ
      ≫ (toPfCat P F).map (rtExt P F (pfDown P F Dd) 1)) :=
    ⟨h, (compRoot_root_one P F ((pfCatToRoot P F).map φ) h).symm.trans h1,
      (compRoot_root_one P F h ((pfCatToRoot P F).map φ)).symm.trans h2⟩
  exact @isIso_of_conj _ _ _ _ _ _
    ((toPfCat P F).map (inv (rtExt P F (pfDown P F X) 1))) hi' φ
    ((toPfCat P F).map (rtExt P F (pfDown P F Dd) 1)) hj' hgiso'

end PfCatRoot

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Definition 3.1, (iii)` の「根 1 の部分から `𝒞^pf` へ」。 -/
def pfCatToRoot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 57,
    item := "Definition 3.1, (iii) — 根 1 の部分 PfCat から 𝒞^pf(PfRootObj)への関手",
    sectionId := "frdi-def-3-1" }

end ABC3.Found.FrdI
