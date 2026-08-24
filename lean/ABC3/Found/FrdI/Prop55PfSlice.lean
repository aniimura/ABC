/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfDiv
import ABC3.Found.FrdI.Prop44Gp
import ABC3.Found.FrdI.Prop55PfKappa
import ABC3.Found.FrdI.Prop55BiratPf

/-!
# [FrdI] `𝒞 → 𝒞^pf` は双有理射の `Div^gp` を `Pf.of` で移す

原文 (FrdI p.105):
> Finally, assertion (iv) is immediate from the definitions [cf. also assertions (i), (ii);

★★`Proposition 5.5, (iv)` の残り(`𝒞^pf` の有理関数の単系 `(Φ^pf)^birat` を
原文の `(Φ^birat)^pf = ℚ·Φ^birat` と同定する段)を、**両側から**攻める。

* `Prop55PfDiv.lean` —— `ℚ` の出どころ(`𝒞^pf` の零因子の分母)。
* **本ファイル** —— **`Φ^birat` の像が `(Φ^pf)^birat` に入る**側。

## ★中身

`Φ^birat(A)` の元は `Div^gp` の代表元の形 `sliceDivGpOf a _ φ`
(`a` は pre-step、`a⁻¹ ≫ φ` の零因子)で書ける(`Proposition 4.4, (iii)`)。
★★自然な関手 `𝒞 → 𝒞^pf`(根 1、`toRootHom`)は 3 つの不変量を

  `Base` そのまま / `deg_Fr` そのまま / `Div ↦ Pf.of (Div)`

と移す(在庫 `rootBase_toRootHom` / `rootDeg_toRootHom` / `rootDiv_toRootHom`)ので、
`sliceDivGpOf` は **`Pf.of` を後ろに掛けるだけ**で移る(`sliceDivGpOf_toRootHom`)。

★代表元の対応(`IdxBirat P G A → IdxBirat (pfRootPre P F) Gpf ⟨A,1⟩`)も通したので、
`biratDivGp` の言葉でも言える(`biratDivGp_toRootHom`)。

## ★★残り(記録)

* **像が `𝒪^×` に入る**こと —— 代表元の対応は関手ではないので、
  「単元が単元へ移る」は `𝒞^birat ⥤ (𝒞^birat)^pf ⥤ (𝒞^pf)^birat` の合成関手
  (`Proposition 5.5, (ii)` の `biratPfHom`)を通して言う必要がある。
* 逆向き(**分母を払えば `Φ^birat` に戻る**)—— こちらが
  `(Φ^pf)^birat ⊆ ℚ·Φ^birat` であり、やはり `Proposition 5.5, (ii)` を使う。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-! ## ★1. 底の同型は移る -/

/-- ★`Base` はそのままなので、底が同型なら `𝒞^pf` でも同型。 -/
theorem isIso_rootBase_toRootHom {A B : C} {a : A ⟶ B} (ha : IsIso (P.Base a)) :
    IsIso ((pfRootPre P F).Base (toRootHom (F := F) a)) := by
  rw [show (pfRootPre P F).Base (toRootHom (F := F) a) = P.Base a from rootBase_toRootHom a]
  exact ha

/-- ★逆射も同じもの。 -/
theorem inv_rootBase_toRootHom {A B : C} {a : A ⟶ B} (ha : IsIso (P.Base a)) :
    @inv _ _ _ _ ((pfRootPre P F).Base (toRootHom (F := F) a))
        (isIso_rootBase_toRootHom (F := F) ha)
      = @inv _ _ _ _ (P.Base a) ha := by
  haveI hB : IsIso ((pfRootPre P F).Base (toRootHom (F := F) a)) :=
    isIso_rootBase_toRootHom (F := F) ha
  refine IsIso.inv_eq_of_hom_inv_id ?_
  rw [show (pfRootPre P F).Base (toRootHom (F := F) a) = P.Base a from rootBase_toRootHom a]
  exact @IsIso.hom_inv_id _ _ _ _ (P.Base a) ha

/-- ★★零因子は `Pf.of` で移る(`rootDiv_toRootHom` の言い換え)。 -/
theorem div_toRootHom {A B : C} (f : A ⟶ B) :
    (pfRootPre P F).Div (toRootHom (F := F) f) = Pf.of (P.Div f) := by
  rw [show (pfRootPre P F).Div (toRootHom (F := F) f) = rootDiv (toRootHom (F := F) f) from rfl,
    rootDiv_toRootHom, Pf.of_apply]

/-! ## ★2. `Div^gp` の代表元も `Pf.of` で移る -/

/-- ★★★★★★**`𝒞 → 𝒞^pf` は双有理射の `Div^gp` を `Pf.of` で移す**。

★★`sliceDivGpOf a _ φ = (Base a)⁻¹^*(Div φ − deg_Fr(φ)·Div a)` の 3 つの材料が
どれも `toRootHom` でそのまま(`Base`・`deg_Fr`)か `Pf.of` を掛けるだけ(`Div`)なので、
**`Pf.of` を外に括り出せる**。

★★★これが `Φ^birat` の像が `(Φ^pf)^birat` に入ることの中身である。 -/
theorem sliceDivGpOf_toRootHom {A B A' : C} (a : A' ⟶ A) (ha : IsIso (P.Base a))
    (φ : A' ⟶ B) :
    sliceDivGpOf (P := pfRootPre P F) (toRootHom (F := F) a)
        (isIso_rootBase_toRootHom (F := F) ha) (toRootHom (F := F) φ)
      = gpMap _ (Pf.of (M := Φ.val (P.toElem.obj A).base)) (sliceDivGpOf a ha φ) := by
  rw [sliceDivGpOf_eq, sliceDivGpOf_eq, div_toRootHom, div_toRootHom,
    show (pfRootPre P F).degFr (toRootHom (F := F) φ) = P.degFr φ from rootDeg_toRootHom φ,
    inv_rootBase_toRootHom ha]
  show (gpMap (Pf (Φ.val (P.toElem.obj A').base))
        (Pf.map (Φ.map (@inv _ _ _ _ (P.Base a) ha))))
      (toGp (Pf (Φ.val (P.toElem.obj A').base)) (Pf.of (P.Div φ))
        - ((P.degFr φ : ℕ+) : ℕ) • toGp (Pf (Φ.val (P.toElem.obj A').base)) (Pf.of (P.Div a)))
    = gpMap _ (Pf.of (M := Φ.val (P.toElem.obj A).base))
      ((gpMap (Φ.val (P.toElem.obj A').base) (Φ.map (@inv _ _ _ _ (P.Base a) ha)))
        (toGp (Φ.val (P.toElem.obj A').base) (P.Div φ)
          - ((P.degFr φ : ℕ+) : ℕ) • toGp (Φ.val (P.toElem.obj A').base) (P.Div a)))
  have hcomp : (Pf.map (Φ.map (@inv _ _ _ _ (P.Base a) ha))).comp
        (Pf.of (M := Φ.val (P.toElem.obj A').base))
      = (Pf.of (M := Φ.val (P.toElem.obj A).base)).comp
        (Φ.map (@inv _ _ _ _ (P.Base a) ha)) := by
    ext z; rfl
  have hz : toGp (Pf (Φ.val (P.toElem.obj A').base)) (Pf.of (P.Div φ))
        - ((P.degFr φ : ℕ+) : ℕ) • toGp (Pf (Φ.val (P.toElem.obj A').base)) (Pf.of (P.Div a))
      = gpMap _ (Pf.of (M := Φ.val (P.toElem.obj A').base))
        (toGp (Φ.val (P.toElem.obj A').base) (P.Div φ)
          - ((P.degFr φ : ℕ+) : ℕ) • toGp (Φ.val (P.toElem.obj A').base) (P.Div a)) := by
    rw [map_sub, map_nsmul, gpMap_toGp, gpMap_toGp]
  rw [hz, ← AddMonoidHom.comp_apply, ← gpMap_comp, hcomp, gpMap_comp, AddMonoidHom.comp_apply]

/-! ## ★3. 代表元の対応 —— `biratDivGp` の言葉で -/

/-- ★★**pre-step は `𝒞^pf` へ移る**(`deg_Fr` はそのまま、底同型は `toRootHom_baseIso`)。 -/
theorem preStep_toRootHom {A B : C} {f : A ⟶ B} (hf : IsPreStep P f) :
    IsPreStep (pfRootPre P F) (toRootHom (F := F) f) :=
  ⟨by
    show rootDeg (toRootHom (F := F) f) = 1
    rw [rootDeg_toRootHom]
    exact hf.1,
   toRootHom_baseIso f hf.2⟩

/-- ★★**双有理化の添字対象を `𝒞^pf` へ移す**。

★co-angular 性は `𝒞^pf` では**無条件**である(`pfRoot_isCoAngular`、
`Proposition 1.4, (i)` ＋ `𝒞^pf` の全対象 isotropic 性)。 -/
noncomputable def idxBiratToRootHom (Gpf : Frobenioid (pfRootPre P F))
    (hfi : IsOfFrobeniusIsotropicType P)
    {A A' : C} (a : A' ⟶ A) (has : IsPreStep P a) :
    IdxBirat (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F) :=
  idxBiratMk (pfRootPre P F) Gpf (toRootHom (F := F) a)
    (pfRoot_isCoAngular hfi _) (preStep_toRootHom has)

/-- ★★★★★★**代表元の `Div^gp` も `Pf.of` で移る**。

★`biratDivGp` は代表元では `sliceDivGpOf` そのもの(`biratDivGp_mk`)なので、
`sliceDivGpOf_toRootHom` がそのまま効く。 -/
theorem biratDivGp_toRootHom {G : Frobenioid P} (Gpf : Frobenioid (pfRootPre P F))
    (hfi : IsOfFrobeniusIsotropicType P)
    {A B A' : C} (a : A' ⟶ A) (hac : IsCoAngular P a) (has : IsPreStep P a) (φ : A' ⟶ B) :
    biratDivGp (P := pfRootPre P F) (G := Gpf)
        (HomBirat.mk (idxBiratToRootHom Gpf hfi a has) (toRootHom (F := F) φ))
      = gpMap _ (Pf.of (M := Φ.val (P.toElem.obj A).base))
        (biratDivGp (HomBirat.mk (idxBiratMk P G a hac has) φ)) := by
  rw [biratDivGp_mk, biratDivGp_mk]
  exact sliceDivGpOf_toRootHom a has.2 φ

/-! ## ★3-b. 底と次数も移る —— `𝒪^×` に入ることの材料 -/

/-- ★★**双有理射の底も `toRootHom` でそのまま**。

★`sliceBaseOf a _ φ = (Base a)⁻¹ ≫ Base φ` の 2 つの材料がどちらもそのままだから。 -/
theorem sliceBaseOf_toRootHom {A B A' : C} (a : A' ⟶ A) (ha : IsIso (P.Base a))
    (φ : A' ⟶ B) :
    sliceBaseOf (P := pfRootPre P F) (toRootHom (F := F) a)
        (isIso_rootBase_toRootHom (F := F) ha) (toRootHom (F := F) φ)
      = sliceBaseOf (P := P) a ha φ := by
  rw [sliceBaseOf_eq, sliceBaseOf_eq, inv_rootBase_toRootHom ha,
    show (pfRootPre P F).Base (toRootHom (F := F) φ) = P.Base φ from rootBase_toRootHom φ]
  rfl

/-- ★★代表元の底は移っても同じ。 -/
theorem biratBase_toRootHom {G : Frobenioid P} (Gpf : Frobenioid (pfRootPre P F))
    (hfi : IsOfFrobeniusIsotropicType P)
    {A B A' : C} (a : A' ⟶ A) (hac : IsCoAngular P a) (has : IsPreStep P a) (φ : A' ⟶ B) :
    biratBase (P := pfRootPre P F) (G := Gpf)
        (HomBirat.mk (idxBiratToRootHom Gpf hfi a has) (toRootHom (F := F) φ))
      = biratBase (HomBirat.mk (idxBiratMk P G a hac has) φ) := by
  rw [biratBase_mk, biratBase_mk]
  exact sliceBaseOf_toRootHom a has.2 φ

/-- ★★代表元の次数は移っても同じ。 -/
theorem biratDeg_toRootHom {G : Frobenioid P} (Gpf : Frobenioid (pfRootPre P F))
    (hfi : IsOfFrobeniusIsotropicType P)
    {A B A' : C} (a : A' ⟶ A) (hac : IsCoAngular P a) (has : IsPreStep P a) (φ : A' ⟶ B) :
    biratDeg (P := pfRootPre P F) (G := Gpf)
        (HomBirat.mk (idxBiratToRootHom Gpf hfi a has) (toRootHom (F := F) φ))
      = biratDeg (HomBirat.mk (idxBiratMk P G a hac has) φ) := by
  rw [biratDeg_mk, biratDeg_mk]
  exact rootDeg_toRootHom φ

/-! ## ★4. 自然な関手 `𝒞^birat ⟶ (𝒞^pf)^birat` を代表元の上で計算する

★★`Proposition 5.5, (ii)` の `biratPfHom` と、`𝒞^birat` の完全化(根 1)への
標準射 `toHomPf` を合成したものが、原文の言う自然な関手である。
★★★**その代表元での値がちょうど `toRootHom` の像になる**ことを見る ——
これで「`Φ^birat` の像が `(Φ^pf)^birat` に入る」が閉じる。

★鍵は `W := idxOne` のとき `biratPfDeg = 1` で、
`biratPfIsoA'` / `biratPfIsoB'` が**次数 1 の `rootLift`＝恒等射**になること。 -/

section Birat

variable {G : Frobenioid P} [IsConnected D]

/-- ★★根 1 では `rootMap` は `toRootHom` そのもの(`κ_A = 𝟙`)。 -/
theorem rootMap_one (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (f : A ⟶ B) :
    rootMap (F := F) hfi f 1 = toRootHom (F := F) f := by
  refine rootMap_ext hfi _ _ ?_
  rw [rootMap_spec, pfKappa_one, pfKappa_one, Category.comp_id, Category.id_comp]

/-- ★★次数 1 の恒等射の `rootLift` は恒等射。 -/
theorem rootLift_id_one (hfi : IsOfFrobeniusIsotropicType P) {A : C}
    (h : P.degFr (𝟙 A) = 1) :
    rootLift (F := F) hfi (𝟙 A) 1 h = 𝟙 (⟨A, 1⟩ : PfRootObj P F) := by
  have hs := rootLift_spec (F := F) hfi (𝟙 A) 1 h
  rw [pfKappa_one, Category.comp_id, toRootHom_id] at hs
  exact hs

theorem biratPfIsoA'_idxOne_inv (hfi : IsOfFrobeniusIsotropicType P) (A B : C)
    (h : P.degFr (idxOne P F A B).hom.hom.1 = 1) :
    (biratPfIsoA' (F := F) hfi (idxOne P F A B) 1 h).inv = 𝟙 (⟨A, 1⟩ : PfRootObj P F) := by
  have hh : (biratPfIsoA' (F := F) hfi (idxOne P F A B) 1 h).hom
      = 𝟙 (⟨A, 1⟩ : PfRootObj P F) := rootLift_id_one hfi h
  have h2 := (biratPfIsoA' (F := F) hfi (idxOne P F A B) 1 h).hom_inv_id
  rw [hh] at h2
  exact (Category.id_comp _).symm.trans h2

theorem biratPfIsoB'_idxOne_inv (hfi : IsOfFrobeniusIsotropicType P) (A B : C)
    (h : P.degFr (idxOne P F A B).hom.hom.2 = 1) :
    (biratPfIsoB' (F := F) hfi (idxOne P F A B) 1 h).inv = 𝟙 (⟨B, 1⟩ : PfRootObj P F) := by
  have hh : (biratPfIsoB' (F := F) hfi (idxOne P F A B) 1 h).hom
      = 𝟙 (⟨B, 1⟩ : PfRootObj P F) := rootLift_id_one hfi h
  have h2 := (biratPfIsoB' (F := F) hfi (idxOne P F A B) 1 h).hom_inv_id
  rw [hh] at h2
  exact (Category.id_comp _).symm.trans h2

/-- ★★★`W := idxOne` のときの `biratPfMk'` は `toRootHom` の像そのもの。 -/
theorem biratPfMk'_idxOne (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) {A B : C}
    (hA : P.degFr (idxOne P F A B).hom.hom.1 = 1)
    (hB : P.degFr (idxOne P F A B).hom.hom.2 = 1)
    (Z : IdxBirat P G A) (ψ : Z.unop.left.obj ⟶ B) :
    biratPfMk' hfi Gpf (idxOne P F A B) 1 hA hB Z ψ
      = HomBirat.mk (idxBiratToRootHom Gpf hfi Z.unop.hom.hom Z.unop.hom.property.2)
        (toRootHom (F := F) ψ) := by
  refine homBirat_mk_congr ?_ _ _ _ _ ?_
  · rw [rootMap_one, biratPfIsoA'_idxOne_inv]
    exact Category.comp_id _
  · rw [rootMap_one, biratPfIsoB'_idxOne_inv]
    exact Category.comp_id _

/-- ★★★★★**自然な関手 `𝒞^birat ⟶ (𝒞^pf)^birat` の代表元での値**。 -/
theorem biratPfHom_toHomPf_mk (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) {A B : C}
    (Z : IdxBirat P G A) (ψ : Z.unop.left.obj ⟶ B) :
    biratPfHom hfi Gpf F' A B (toHomPf (F := F') (HomBirat.mk Z ψ))
      = HomBirat.mk (idxBiratToRootHom Gpf hfi Z.unop.hom.hom Z.unop.hom.property.2)
        (toRootHom (F := F) ψ) :=
  (biratPfHom_mk hfi Gpf F' A B (idxOne P F A B) (HomBirat.mk Z ψ)).trans
    (((biratPf_mk hfi Gpf (idxOne P F A B) Z ψ).trans
      (biratPfMk'_eq hfi Gpf (idxOne P F A B) 1 (P.degFr_id A) (P.degFr_id B) Z ψ).symm).trans
      (biratPfMk'_idxOne hfi Gpf (P.degFr_id A) (P.degFr_id B) Z ψ))

/-- ★★★★★★**自然な関手は `Div^gp` を `Pf.of` で移す**(代表元の版)。

★★★これが `Proposition 5.5, (iv)` の単系の同定の「片側」——
`Φ^birat` の元 `Div^gp(δ)` は `(Φ^pf)^birat` の元 `Div^gp(Ξ δ)` に
`Pf.of` を掛けただけで移る。 -/
theorem biratDivGp_biratPfHom_toHomPf (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) {A B : C}
    (Z : IdxBirat P G A) (ψ : Z.unop.left.obj ⟶ B) :
    biratDivGp (P := pfRootPre P F) (G := Gpf)
        (biratPfHom hfi Gpf F' A B (toHomPf (F := F') (HomBirat.mk Z ψ)))
      = gpMap _ (Pf.of (M := Φ.val (P.toElem.obj A).base))
        (biratDivGp (HomBirat.mk Z ψ)) := by
  rw [biratPfHom_toHomPf_mk]
  exact biratDivGp_toRootHom Gpf hfi Z.unop.hom.hom Z.unop.hom.property.1
    Z.unop.hom.property.2 ψ

/-- ★★★**自然な関手は底を動かさない**(代表元を取って `biratBase_toRootHom`)。 -/
theorem biratBase_biratPfHom (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A B : C)
    (ε : HomBirat P G A B) :
    biratBase (P := pfRootPre P F) (G := Gpf)
        (biratPfHom hfi Gpf F' A B (toHomPf (F := F') ε)) = biratBase ε := by
  obtain ⟨Z, ψ, hZ⟩ := HomBirat.exists_rep ε
  subst hZ
  exact (congrArg (biratBase (P := pfRootPre P F) (G := Gpf))
      (biratPfHom_toHomPf_mk hfi Gpf F' Z ψ)).trans
    (biratBase_toRootHom Gpf hfi Z.unop.hom.hom Z.unop.hom.property.1
      Z.unop.hom.property.2 ψ)

/-- ★★★**自然な関手は次数を動かさない**。 -/
theorem biratDeg_biratPfHom (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A B : C)
    (ε : HomBirat P G A B) :
    biratDeg (P := pfRootPre P F) (G := Gpf)
        (biratPfHom hfi Gpf F' A B (toHomPf (F := F') ε)) = biratDeg ε := by
  obtain ⟨Z, ψ, hZ⟩ := HomBirat.exists_rep ε
  subst hZ
  exact (congrArg (biratDeg (P := pfRootPre P F) (G := Gpf))
      (biratPfHom_toHomPf_mk hfi Gpf F' Z ψ)).trans
    (biratDeg_toRootHom Gpf hfi Z.unop.hom.hom Z.unop.hom.property.1
      Z.unop.hom.property.2 ψ)

/-- ★★★**自然な関手は `Div^gp` を `Pf.of` で移す**(一般の元について)。 -/
theorem biratDivGp_biratPfHom (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A B : C)
    (ε : HomBirat P G A B) :
    biratDivGp (P := pfRootPre P F) (G := Gpf)
        (biratPfHom hfi Gpf F' A B (toHomPf (F := F') ε))
      = gpMap _ (Pf.of (M := Φ.val (P.toElem.obj A).base)) (biratDivGp ε) := by
  obtain ⟨Z, ψ, hZ⟩ := HomBirat.exists_rep ε
  subst hZ
  exact biratDivGp_biratPfHom_toHomPf hfi Gpf F' Z ψ

end Birat

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Proposition 5.5, (iv)` の単系の同定のうち、
「`Φ^birat` の像は `(Φ^pf)^birat` に入る」側(`Div^gp` の代表元の計算)。 -/
def sliceDivGpOf_toRootHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iv) — 𝒞 → 𝒞^pf は双有理射の Div^gp を Pf.of で移す",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
