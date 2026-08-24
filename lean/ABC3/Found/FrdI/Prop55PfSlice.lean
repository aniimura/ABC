/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfDiv
import ABC3.Found.FrdI.Prop44Gp
import ABC3.Found.FrdI.Prop55PfKappa
import ABC3.Found.FrdI.Prop55BiratPf
import ABC3.Found.FrdI.Prop53Birat
import ABC3.Found.FrdI.Prop53QPhi

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

★★★**像が `𝒪^×` に入る**ことも通した(`otimes_biratPfHom`)——
代表元の対応は関手ではないので、単元性は
`Ξ := biratPfHom ∘ toHomPf : 𝒞^birat ⟶ (𝒞^pf)^birat`(`Proposition 5.5, (ii)`)が
恒等射と合成を保つこと(`biratPfHom_unit`)を通して運ぶ。
★★★★これで **`Φ^birat` の像 ⊆ `(Φ^pf)^birat`** が閉じた
(`phiBiratAt_pf_mem` / `phiBiratOn_pf_mem`)。

★★★★★**逆向きも閉じた**(`biratDivGp_nsmul_mem_pfImage`)——
`(Φ^pf)^birat` の元は**正の整数倍で `Φ^birat` の像に入る**。筋は 3 段:

1. `biratPfHom` は全単射なので、`(𝒞^pf)^birat` の単元 `δ'` は
   `x ∈ Hom^pf_{𝒞^birat}` から来る(`biratPfEndHom`、自己射の単系の準同型)。
2. `x` は `Proposition 5.5, (i)` の完全化の構造から「`α` の `k` 乗根」であり、
   `x^k` は `𝒞^birat` の元 `α` の像である(`biratPfHom_pow_eq_Xi`)。
   ★根つき/根なしの橋は `endRootOneEquiv`(`Prop55PfDiv.lean`)。
3. `Div^gp` は `𝒪^×` の上で加法的(`biratDivGp_pow_otimes`)なので
   `k · Div^gp(δ') = Div^gp(δ'^k) = Pf.of (Div^gp α)` となり、分母が払える。
   ★`𝒞^birat` の isotropic 対象では `𝒪^▷ = 𝒪^×`(`otri_mem_otimes_birat`)。

★★★★★★**仮定は `δ' ∈ 𝒪^×` だけになった**(`otri_of_otimes_biratPfHom`)——
`biratPfHom` が `Base` と `deg_Fr` を保つこと
(`biratBase_biratPfHom_eq` / `biratDeg_biratPfHom_eq`、原文が (ii) で
「compatible with the functors to the respective elementary Frobenioids」と書く条)
から、`𝒪^▷` の 2 条件がそのまま `x` へ降りる。
★その律速だった「`rootMap` / `rootLift` の `rootBase`」は
**`𝒞^pf` の標準射の底が恒等**(`rootBase_pfKappa`)から出た。

## ★★残り(記録)

両側の包含を `qPhiBiratOn`(`Prop53QPhi.lean` の飽和)の言葉で束ねて
   **部分群の等号**にすること。★`⊆` の側には
   「`Φ^birat` の元の `k` 乗根が `(Φ^pf)^birat` に在る」＋
   `Gp (Pf M)` の torsion-free 性が要る。
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

/-! ## ★5. `𝒪^×` は `𝒪^×` へ移る -/

/-- ★`𝒞^birat` の恒等射の底は `𝟙`。 -/
theorem biratBase_id (A : C) :
    biratBase (P := P) (G := G) (𝟙 (biratUp P G A)) = 𝟙 ((P.toElem.obj A).base) := by
  show biratBase (toHomBirat (P := P) (G := G) (𝟙 A)) = _
  rw [biratBase_toHomBirat, P.Base_id]

/-- ★★★★★**自然な関手は `𝒪^×` を `𝒪^×` へ移す**。

★底恒等性と次数 1 は `biratBase_biratPfHom` / `biratDeg_biratPfHom`、
単元性は `biratPfHom_unit`(恒等射と合成の保存)から。 -/
theorem otimes_biratPfHom (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C)
    {δ : biratUp P G A ⟶ biratUp P G A}
    (hδ : δ ∈ OTimes (biratPre P G) (biratUp P G A)) :
    biratPfHom hfi Gpf F' A A (toHomPf (F := F') δ)
      ∈ OTimes (biratPre (pfRootPre P F) Gpf)
        (biratUp (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F)) := by
  obtain ⟨⟨hb, hl⟩, hu⟩ := hδ
  obtain ⟨u, hu'⟩ := hu
  refine ⟨⟨?_, ?_⟩, ?_⟩
  · show biratBase (P := pfRootPre P F) (G := Gpf)
        (biratPfHom hfi Gpf F' A A (toHomPf (F := F') δ))
      = biratBase (P := pfRootPre P F) (G := Gpf)
        (𝟙 (biratUp (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F)))
    rw [biratBase_biratPfHom, biratBase_id]
    have hb' : biratBase (P := P) (G := G) δ
        = biratBase (P := P) (G := G) (𝟙 (biratUp P G A)) := hb
    exact hb'.trans (biratBase_id A)
  · show biratDeg (P := pfRootPre P F) (G := Gpf)
        (biratPfHom hfi Gpf F' A A (toHomPf (F := F') δ)) = 1
    rw [biratDeg_biratPfHom]
    exact hl
  · have hv : (u : End (biratUp P G A)) = δ := hu'
    have e1 : δ ≫ (u.inv : biratUp P G A ⟶ biratUp P G A) = 𝟙 _ := by
      have := u.inv_val
      rw [hv] at this
      exact this
    have e2 : (u.inv : biratUp P G A ⟶ biratUp P G A) ≫ δ = 𝟙 _ := by
      have := u.val_inv
      rw [hv] at this
      exact this
    have h1 : compPf (biratPre P G) F' (toHomPf (F := F') δ)
        (toHomPf (F := F') (u.inv : biratUp P G A ⟶ biratUp P G A))
        = toHomPf (F := F') (𝟙 (biratUp P G A)) := by
      rw [← toHomPf_comp, e1]
    have h2 : compPf (biratPre P G) F'
        (toHomPf (F := F') (u.inv : biratUp P G A ⟶ biratUp P G A)) (toHomPf (F := F') δ)
        = toHomPf (F := F') (𝟙 (biratUp P G A)) := by
      rw [← toHomPf_comp, e2]
    have k1 := biratPfHom_unit hfi Gpf F' A _ _ h1
    have k2 := biratPfHom_unit hfi Gpf F' A _ _ h2
    exact ⟨⟨biratPfHom hfi Gpf F' A A (toHomPf (F := F') δ),
      biratPfHom hfi Gpf F' A A
        (toHomPf (F := F') (u.inv : biratUp P G A ⟶ biratUp P G A)), k2, k1⟩, rfl⟩

/-- ★★★★★★**`Φ^birat` の像は `(Φ^pf)^birat` に入る**。

★★★これが `Proposition 5.5, (iv)` の単系の同定の**片側**である ——
`Φ^birat(A)` の元 `Div^gp(δ)` は、自然な関手で移した `Div^gp(Ξ δ)` に
`Pf.of` を掛けただけのものであり、`Ξ δ` は `𝒪^×` に入る。 -/
theorem phiBiratAt_pf_mem (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C)
    {y : Gp (Φ.val (P.toElem.obj A).base)} (hy : y ∈ phiBiratAt P G (biratUp P G A)) :
    gpMap _ (Pf.of (M := Φ.val (P.toElem.obj A).base)) y
      ∈ phiBiratAt (pfRootPre P F) Gpf
        (biratUp (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F)) := by
  obtain ⟨δ, hδ, rfl⟩ := hy
  exact ⟨biratPfHom hfi Gpf F' A A (toHomPf (F := F') (δ : biratUp P G A ⟶ biratUp P G A)),
    otimes_biratPfHom hfi Gpf F' A hδ,
    biratDivGp_biratPfHom hfi Gpf F' A A (δ : biratUp P G A ⟶ biratUp P G A)⟩

/-- ★★★★★★同上を **`𝒟` で添字づけた形**で(`Φ^birat` は `𝒟` 上の部分関手なので、
model Frobenioid が使うのはこちらの形である)。 -/
theorem phiBiratOn_pf_mem (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C)
    {y : Gp (Φ.val (P.toElem.obj A).base)}
    (hy : y ∈ phiBiratOn G (P.toElem.obj A).base) :
    gpMap _ (Pf.of (M := Φ.val (P.toElem.obj A).base)) y
      ∈ phiBiratOn Gpf ((pfRootPre P F).toElem.obj (⟨A, 1⟩ : PfRootObj P F)).base := by
  rw [phiBiratOn_base G hiso A] at hy
  rw [phiBiratOn_base Gpf (pfRoot_isOfIsotropicType (F := F) hfi) (⟨A, 1⟩ : PfRootObj P F)]
  exact phiBiratAt_pf_mem hfi Gpf F' A hy

/-! ## ★5-b. `𝒞^pf` の標準射の**底**は恒等

★★★**測って分かった**(2026-08-24): `κ_{A,k} = rtExt(A,1)⁻¹ ≫ rtExt(A,k)` の底は
`rootBase` の共役(`Base(rtExt X.obj Y.root) ≫ − ≫ Base(rtExt Y.obj X.root)⁻¹`)と
**ちょうど打ち消し合って恒等になる**。
★★これで `rootMap` と `rootLift` の底も明示式になり、
`𝒞^pf` の中の `Base` の計算がすべて `𝒞` の `Base` に落ちる。

★これらは本来 `Prop55PfKappa.lean` の在庫だが、下流の再ビルドを避けて本ファイルに置く。 -/

/-- ★★★★**標準射 `κ_{A,k}` の底は恒等**。 -/
theorem rootBase_pfKappa (A : C) (k : ℕ+) :
    rootBase (pfKappa (F := F) A k) = 𝟙 ((P.toElem.obj A).base) := by
  haveI h1 : IsIso (P.Base (rtExt P F A 1)) := (rtExt_frobType P F A 1).2
  haveI hk : IsIso (P.Base (rtExt P F A k)) := (rtExt_frobType P F A k).2
  haveI := isIso_rtExt_one P F A
  show P.Base (rtExt P F A 1) ≫ pfBase (toHomPf (F := F) (kappaRep (P := P) (F := F) A k))
      ≫ inv (P.Base (rtExt P F A k)) = _
  rw [pfBase_toHomPf]
  show P.Base (rtExt P F A 1) ≫ P.Base (inv (rtExt P F A 1) ≫ rtExt P F A k)
      ≫ inv (P.Base (rtExt P F A k)) = _
  rw [P.Base_comp, ← Category.assoc, ← Category.assoc, ← P.Base_comp, IsIso.hom_inv_id,
    P.Base_id, Category.id_comp, IsIso.hom_inv_id]

/-- ★★★**根を上げた射の底はもとの射の底**。 -/
theorem rootBase_rootMap (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (f : A ⟶ B) (k : ℕ+) :
    rootBase (rootMap (F := F) hfi f k) = P.Base f := by
  have h := congrArg (pfRootPre P F).Base (rootMap_spec (F := F) hfi f k)
  rw [(pfRootPre P F).Base_comp, (pfRootPre P F).Base_comp] at h
  have h2 : rootBase (rootMap (F := F) hfi f k) ≫ rootBase (pfKappa (F := F) B k)
      = rootBase (pfKappa (F := F) A k) ≫ rootBase (toRootHom (F := F) f) := h
  rw [rootBase_pfKappa, rootBase_pfKappa, Category.comp_id, Category.id_comp,
    rootBase_toRootHom] at h2
  exact h2

/-- ★★★**持ち上げた射の底ももとの射の底**。 -/
theorem rootBase_rootLift (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    (k : ℕ+) (hk : P.degFr α = k) :
    rootBase (rootLift (F := F) hfi α k hk) = P.Base α := by
  have h := congrArg (pfRootPre P F).Base (rootLift_spec (F := F) hfi α k hk)
  rw [(pfRootPre P F).Base_comp] at h
  have h2 : rootBase (rootLift (F := F) hfi α k hk) ≫ rootBase (pfKappa (F := F) A₁ k)
      = rootBase (toRootHom (F := F) α) := h
  rw [rootBase_pfKappa, Category.comp_id, rootBase_toRootHom] at h2
  exact h2

/-! ## ★5-c. `biratPfHom` は `Base` と `deg_Fr` を保つ

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the func-

★★原文が「compatible with the functors to the respective elementary Frobenioids」と
書く条である。★代表元(`homPf_birat_exists_rep` で `idxToBirat` の像に取る)の上で
両辺を計算すると、**どちらも同じ 3 つの `Base` の合成**になる。 -/

section Compat

variable {G : Frobenioid P} [IsConnected D]

/-- ★共通の終域を持つ 2 つの同型についての恒等式(`Base` の計算で 2 度使う)。 -/
theorem comp_inv_comp_inv_id {X Y Z : D} (f : X ⟶ Z) (g : Y ⟶ Z)
    (hf : IsIso f) (hg : IsIso g) :
    (f ≫ @inv _ _ _ _ g hg) ≫ (g ≫ @inv _ _ _ _ f hf) = 𝟙 X := by
  haveI := hf
  haveI := hg
  exact (Category.assoc _ _ _).trans
    ((congrArg (fun t : Z ⟶ X => f ≫ t)
      (IsIso.inv_hom_id_assoc g (@inv _ _ _ _ f hf))).trans (IsIso.hom_inv_id f))

theorem rootBase_biratPfIsoA_inv (hfi : IsOfFrobeniusIsotropicType P) {A B : C}
    (W : IdxPf P F A B) (h1 : IsIso (P.Base W.hom.hom.1)) :
    rootBase ((biratPfIsoA (F := F) hfi W).inv) = @inv _ _ _ _ (P.Base W.hom.hom.1) h1 := by
  have hh : (pfRootPre P F).Base ((biratPfIsoA (F := F) hfi W).hom) = P.Base W.hom.hom.1 :=
    rootBase_rootLift hfi W.hom.hom.1 (biratPfDeg W) rfl
  have h2 : (pfRootPre P F).Base ((biratPfIsoA (F := F) hfi W).hom)
      ≫ (pfRootPre P F).Base ((biratPfIsoA (F := F) hfi W).inv) = 𝟙 _ := by
    rw [← (pfRootPre P F).Base_comp, (biratPfIsoA (F := F) hfi W).hom_inv_id,
      (pfRootPre P F).Base_id]
  rw [hh] at h2
  exact (@IsIso.inv_eq_of_hom_inv_id _ _ _ _ (P.Base W.hom.hom.1) h1 _ h2).symm

theorem rootBase_biratPfIsoB_inv (hfi : IsOfFrobeniusIsotropicType P) {A B : C}
    (W : IdxPf P F A B) (h2i : IsIso (P.Base W.hom.hom.2)) :
    rootBase ((biratPfIsoB (F := F) hfi W).inv) = @inv _ _ _ _ (P.Base W.hom.hom.2) h2i := by
  have hh : (pfRootPre P F).Base ((biratPfIsoB (F := F) hfi W).hom) = P.Base W.hom.hom.2 :=
    rootBase_rootLift hfi W.hom.hom.2 (biratPfDeg W) W.hom.property.2.2.symm
  have h2 : (pfRootPre P F).Base ((biratPfIsoB (F := F) hfi W).hom)
      ≫ (pfRootPre P F).Base ((biratPfIsoB (F := F) hfi W).inv) = 𝟙 _ := by
    rw [← (pfRootPre P F).Base_comp, (biratPfIsoB (F := F) hfi W).hom_inv_id,
      (pfRootPre P F).Base_id]
  rw [hh] at h2
  exact (@IsIso.inv_eq_of_hom_inv_id _ _ _ _ (P.Base W.hom.hom.2) h2i _ h2).symm

/-- ★★★代表元の側の `deg_Fr`。 -/
theorem biratDeg_biratPfMk (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) {A B : C} (W : IdxPf P F A B)
    (Z : IdxBirat P G W.right.obj.1) (ψ : Z.unop.left.obj ⟶ W.right.obj.2) :
    biratDeg (P := pfRootPre P F) (G := Gpf) (biratPfMk hfi Gpf W Z ψ) = P.degFr ψ := by
  show biratDeg (P := pfRootPre P F) (G := Gpf)
      (HomBirat.mk (biratPfIdx hfi Gpf W Z)
        (rootMap (F := F) hfi ψ (biratPfDeg W) ≫ (biratPfIsoB hfi W).inv)) = _
  rw [biratDeg_mk]
  show (pfRootPre P F).degFr
      (rootMap (F := F) hfi ψ (biratPfDeg W) ≫ (biratPfIsoB hfi W).inv) = _
  rw [(pfRootPre P F).degFr_comp, degFr_of_isIso (pfRootPre P F) (biratPfIsoB hfi W).inv,
    one_mul, rootMap_degFr]

/-- ★★★代表元の側の `Base`(`𝒞^pf` の中で計算したもの)。 -/
theorem biratBase_biratPfMk (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) {A B : C} (W : IdxPf P F A B)
    (Z : IdxBirat P G W.right.obj.1) (ψ : Z.unop.left.obj ⟶ W.right.obj.2)
    (h1 : IsIso (P.Base W.hom.hom.1)) (h2 : IsIso (P.Base W.hom.hom.2))
    (hZi : IsIso (P.Base Z.unop.hom.hom)) :
    biratBase (P := pfRootPre P F) (G := Gpf) (biratPfMk hfi Gpf W Z ψ)
      = P.Base W.hom.hom.1 ≫ @inv _ _ _ _ (P.Base Z.unop.hom.hom) hZi
        ≫ P.Base ψ ≫ @inv _ _ _ _ (P.Base W.hom.hom.2) h2 := by
  have hbi : IsIso ((pfRootPre P F).Base ((biratPfIdx hfi Gpf W Z).unop.hom.hom)) :=
    (biratPfIdx hfi Gpf W Z).unop.hom.property.2.2
  have ha : (pfRootPre P F).Base ((biratPfIdx hfi Gpf W Z).unop.hom.hom)
      = P.Base Z.unop.hom.hom ≫ @inv _ _ _ _ (P.Base W.hom.hom.1) h1 := by
    show (pfRootPre P F).Base (rootMap (F := F) hfi Z.unop.hom.hom (biratPfDeg W)
        ≫ (biratPfIsoA (F := F) hfi W).inv) = _
    rw [(pfRootPre P F).Base_comp]
    exact congrArg₂ (· ≫ ·) (rootBase_rootMap hfi _ _) (rootBase_biratPfIsoA_inv hfi W h1)
  have hphi : (pfRootPre P F).Base (rootMap (F := F) hfi ψ (biratPfDeg W)
        ≫ (biratPfIsoB (F := F) hfi W).inv)
      = P.Base ψ ≫ @inv _ _ _ _ (P.Base W.hom.hom.2) h2 := by
    rw [(pfRootPre P F).Base_comp]
    exact congrArg₂ (· ≫ ·) (rootBase_rootMap hfi _ _) (rootBase_biratPfIsoB_inv hfi W h2)
  have hinv : @inv _ _ _ _ ((pfRootPre P F).Base ((biratPfIdx hfi Gpf W Z).unop.hom.hom)) hbi
      = P.Base W.hom.hom.1 ≫ @inv _ _ _ _ (P.Base Z.unop.hom.hom) hZi := by
    refine @IsIso.inv_eq_of_hom_inv_id _ _ _ _ _ hbi _ ?_
    exact ha ▸ comp_inv_comp_inv_id (P.Base Z.unop.hom.hom) (P.Base W.hom.hom.1) hZi h1
  rw [show biratPfMk hfi Gpf W Z ψ
      = HomBirat.mk (biratPfIdx hfi Gpf W Z)
        (rootMap (F := F) hfi ψ (biratPfDeg W) ≫ (biratPfIsoB (F := F) hfi W).inv) from rfl,
    biratBase_mk, sliceBaseOf_eq]
  exact (congrArg₂ (· ≫ ·) hinv hphi).trans (Category.assoc _ _ _)

/-- ★★★代表元の側の `Base`(`(𝒞^birat)^pf` の中で計算したもの)—— **同じ式になる**。 -/
theorem pfBase_mk_idxToBirat (F' : FrobenioidCore (biratPre P G)) {A B : C} (W : IdxPf P F A B)
    (Z : IdxBirat P G W.right.obj.1) (ψ : Z.unop.left.obj ⟶ W.right.obj.2)
    (h2 : IsIso (P.Base W.hom.hom.2)) (hZi : IsIso (P.Base Z.unop.hom.hom)) :
    pfBase (HomPf.mk ((idxToBirat P F G F' A B).obj W) (HomBirat.mk Z ψ))
      = P.Base W.hom.hom.1 ≫ @inv _ _ _ _ (P.Base Z.unop.hom.hom) hZi
        ≫ P.Base ψ ≫ @inv _ _ _ _ (P.Base W.hom.hom.2) h2 := by
  have hbi2 : IsIso ((biratPre P G).Base
      (((idxToBirat P F G F' A B).obj W).hom.hom.2)) :=
    ((idxToBirat P F G F' A B).obj W).hom.property.2.1.2
  have hb1 : (biratPre P G).Base (((idxToBirat P F G F' A B).obj W).hom.hom.1)
      = P.Base W.hom.hom.1 := biratBase_toHomBirat _
  have hb2 : (biratPre P G).Base (((idxToBirat P F G F' A B).obj W).hom.hom.2)
      = P.Base W.hom.hom.2 := biratBase_toHomBirat _
  have hbz : (biratPre P G).Base (HomBirat.mk Z ψ)
      = @inv _ _ _ _ (P.Base Z.unop.hom.hom) hZi ≫ P.Base ψ := by
    show biratBase (HomBirat.mk Z ψ) = _
    rw [biratBase_mk, sliceBaseOf_eq]
  have hinv2 : @inv _ _ _ _ ((biratPre P G).Base
        (((idxToBirat P F G F' A B).obj W).hom.hom.2)) hbi2
      = @inv _ _ _ _ (P.Base W.hom.hom.2) h2 := by
    refine @IsIso.inv_eq_of_hom_inv_id _ _ _ _ _ hbi2 _ ?_
    rw [hb2]
    exact @IsIso.hom_inv_id _ _ _ _ (P.Base W.hom.hom.2) h2
  rw [pfBase_mk]
  show (biratPre P G).Base (((idxToBirat P F G F' A B).obj W).hom.hom.1)
      ≫ (biratPre P G).Base (HomBirat.mk Z ψ)
      ≫ @inv _ _ _ _ ((biratPre P G).Base
        (((idxToBirat P F G F' A B).obj W).hom.hom.2)) hbi2 = _
  rw [hb1]
  exact congrArg (P.Base W.hom.hom.1 ≫ ·)
    ((congrArg₂ (· ≫ ·) hbz hinv2).trans (Category.assoc _ _ _))

/-- ★★★★★**`biratPfHom` は `deg_Fr` を保つ**。 -/
theorem biratDeg_biratPfHom_eq (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A B : C)
    (x : HomPf (biratPre P G) F' (biratUp P G A) (biratUp P G B)) :
    biratDeg (P := pfRootPre P F) (G := Gpf) (biratPfHom hfi Gpf F' A B x) = pfDeg x := by
  obtain ⟨W, z, hW⟩ := homPf_birat_exists_rep (F := F) F' A B x
  subst hW
  obtain ⟨Z, ψ, hZ⟩ := HomBirat.exists_rep z
  subst hZ
  rw [biratPfHom_mk, biratPf_mk, biratDeg_biratPfMk, pfDeg_mk]
  show _ = biratDeg (HomBirat.mk Z ψ)
  rw [biratDeg_mk]

/-- ★★★★★**`biratPfHom` は `Base` を保つ**。

★★★これが原文の「compatible with the functors to the respective elementary
Frobenioids」の中身であり、`𝒪^▷` / `𝒪^×` の所属を両側で行き来させる鍵である。 -/
theorem biratBase_biratPfHom_eq (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A B : C)
    (x : HomPf (biratPre P G) F' (biratUp P G A) (biratUp P G B)) :
    biratBase (P := pfRootPre P F) (G := Gpf) (biratPfHom hfi Gpf F' A B x) = pfBase x := by
  obtain ⟨W, z, hW⟩ := homPf_birat_exists_rep (F := F) F' A B x
  subst hW
  obtain ⟨Z, ψ, hZ⟩ := HomBirat.exists_rep z
  subst hZ
  have h1 : IsIso (P.Base W.hom.hom.1) := W.hom.property.1.2
  have h2 : IsIso (P.Base W.hom.hom.2) := W.hom.property.2.1.2
  have hZi : IsIso (P.Base Z.unop.hom.hom) := Z.unop.hom.property.2.2
  rw [biratPfHom_mk, biratPf_mk, biratBase_biratPfMk hfi Gpf W Z ψ h1 h2 hZi]
  exact (pfBase_mk_idxToBirat F' W Z ψ h2 hZi).symm

end Compat

/-! ## ★5-d. 仮定 `hx` を `hδ` から導く -/

/-- ★★**同型による共役は `𝒪^▷` を保つ**(底恒等性も次数 1 も共役で不変)。 -/
theorem mem_otri_conj {Cc : Type u2} [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} D}
    (P₀ : PreFrobenioid Cc Φ₀) {X Y : Cc} (u : X ≅ Y) (f : End X)
    (hf : f ∈ OTri P₀ X) : (u.conj f) ∈ OTri P₀ Y := by
  obtain ⟨hb, hl⟩ := hf
  haveI : IsIso u.hom := u.isIso_hom
  haveI : IsIso u.inv := u.isIso_inv
  refine ⟨?_, ?_⟩
  · show P₀.Base (u.inv ≫ ((f : End X) : X ⟶ X) ≫ u.hom) = P₀.Base (𝟙 Y)
    rw [P₀.Base_comp, P₀.Base_comp, show P₀.Base ((f : End X) : X ⟶ X) = P₀.Base (𝟙 X) from hb,
      P₀.Base_id, Category.id_comp, ← P₀.Base_comp, u.inv_hom_id]
  · show P₀.degFr (u.inv ≫ ((f : End X) : X ⟶ X) ≫ u.hom) = 1
    rw [P₀.degFr_comp, P₀.degFr_comp, degFr_of_isIso P₀ u.hom, degFr_of_isIso P₀ u.inv,
      show P₀.degFr ((f : End X) : X ⟶ X) = 1 from hl, one_mul, one_mul]

/-- ★★★**根なしの `𝒪^▷` は根つき化身の `𝒪^▷` に移る**(橋 `endRootOneEquiv` で)。 -/
theorem mem_otri_endRootOneEquiv_symm (A : C) (x : End (pfObjOf P F A))
    (hx : x ∈ OTri (pfPre P F) (pfObjOf P F A)) :
    (endRootOneEquiv (F := F) A).symm x
      ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F) := by
  refine (mem_otri_rootSelfIso_inv A 1 _).mpr ?_
  have h : endRootMulEquiv (F := F) A 1 ((endRootOneEquiv (F := F) A).symm x)
      = endPfCatRtOne (F := F) A x :=
    (endRootMulEquiv (F := F) A 1).apply_symm_apply _
  show endRootMulEquiv (F := F) A 1 ((endRootOneEquiv (F := F) A).symm x)
      ∈ OTri (pfPre P F) (rtObjPf P F A 1)
  rw [h]
  exact mem_otri_conj (pfPre P F)
    ((toPfCat P F).mapIso (@asIso _ _ _ _ (rtExt P F A 1) (isIso_rtExt_one P F A))) x hx

/-! ## ★6. 逆向き —— 分母を払えば `Φ^birat` に戻る -/

/-- ★★★★**`biratPfHom` は自己射の単系の準同型**(`Proposition 5.5, (ii)` の
「恒等射と合成を保つ」を `MonoidHom` に束ねたもの)。

★`End` の積は `x * y = y ≫ x` なので、`biratPfHom_comp` を引数を入れ替えて渡す。 -/
noncomputable def biratPfEndHom (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C) :
    End (pfObjOf (biratPre P G) F' (biratUp P G A))
      →* End (biratUp (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F)) where
  toFun x := biratPfHom hfi Gpf F' A A x
  map_one' := biratPfHom_id hfi Gpf F' A
  map_mul' x y := biratPfHom_comp hfi Gpf F' A A A y x

@[simp] theorem biratPfEndHom_apply (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C)
    (x : End (pfObjOf (biratPre P G) F' (biratUp P G A))) :
    biratPfEndHom hfi Gpf F' A x = biratPfHom hfi Gpf F' A A x := rfl

/-- ★★★★★**`δ' ∈ 𝒪^×` から仮定 `hx` が出る**。

★★`biratPfHom` が `Base` と `deg_Fr` を保つこと(`biratBase_biratPfHom_eq` /
`biratDeg_biratPfHom_eq`)から、`𝒪^▷` の 2 条件がそのまま `x` へ降りる。
★★★これで逆向きの包含が**単元性の仮定だけ**から出るようになった。 -/
theorem otri_of_otimes_biratPfHom (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C)
    {x : End (pfObjOf (biratPre P G) F' (biratUp P G A))}
    (hδ : biratPfEndHom hfi Gpf F' A x
      ∈ OTimes (biratPre (pfRootPre P F) Gpf)
        (biratUp (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F))) :
    (endRootOneEquiv (F := F') (biratUp P G A)).symm x
      ∈ OTri (pfRootPre (biratPre P G) F')
        (⟨biratUp P G A, 1⟩ : PfRootObj (biratPre P G) F') := by
  refine mem_otri_endRootOneEquiv_symm (F := F') (biratUp P G A) x ?_
  obtain ⟨⟨hb, hl⟩, -⟩ := hδ
  refine ⟨?_, ?_⟩
  · refine (isBaseIdentity_pf_iff (P := biratPre P G) (F := F') x).mpr ?_
    have hb' : biratBase (P := pfRootPre P F) (G := Gpf)
        (biratPfHom hfi Gpf F' A A x)
        = biratBase (P := pfRootPre P F) (G := Gpf)
          (𝟙 (biratUp (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F))) := hb
    rw [biratBase_id] at hb'
    exact (biratBase_biratPfHom_eq hfi Gpf F' A A x).symm.trans hb'
  · show pfDeg x = 1
    have hl' : biratDeg (P := pfRootPre P F) (G := Gpf)
        (biratPfHom hfi Gpf F' A A x) = 1 := hl
    exact (biratDeg_biratPfHom_eq hfi Gpf F' A A x).symm.trans hl'

/-- ★★★★★★**`(𝒞^pf)^birat` の自己射は、正の冪で `𝒞^birat` から来る**。

★★★これが `(Φ^pf)^birat ⊆ ℚ·Φ^birat`(逆向きの包含)の本体である ——
`Proposition 5.5, (ii)` の全単射で `𝒞^birat` の完全化へ移し、
そこで `Proposition 5.5, (i)` の「正の冪で戻る」(`hom_pf_pow_toRootHom`)を当てる。 -/
theorem biratPfHom_pow_eq_Xi (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C)
    (hisoB : ∀ X : BiratCat P G, IsIsotropic (biratPre P G) X)
    (hAB : IsFrobeniusTrivial (biratPre P G) (biratUp P G A))
    (hfnB : IsFrobeniusNormalized (biratPre P G) (biratUp P G A))
    (hfnB' : IsFrobeniusNormalized (biratPre P G)
      (rtObj (biratPre P G) F' (biratUp P G A) 1))
    (ζ : ℕ+ →* End (biratUp P G A))
    (hdeg : ∀ n : ℕ+, (biratPre P G).degFr ((ζ n : End (biratUp P G A)) :
      biratUp P G A ⟶ biratUp P G A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity (biratPre P G) (ζ n)
      ∧ IsFrobeniusType (biratPre P G) ((ζ n : End (biratUp P G A)) : _ ⟶ _))
    (x : End (pfObjOf (biratPre P G) F' (biratUp P G A)))
    (hx : (endRootOneEquiv (F := F') (biratUp P G A)).symm x
      ∈ OTri (pfRootPre (biratPre P G) F')
        (⟨biratUp P G A, 1⟩ : PfRootObj (biratPre P G) F')) :
    ∃ (k : ℕ+) (α : OTri (biratPre P G) (biratUp P G A)),
      (biratPfEndHom hfi Gpf F' A x) ^ ((k : ℕ+) : ℕ)
        = biratPfEndHom hfi Gpf F' A
          (toHomPf (F := F') ((α : End (biratUp P G A)) : _ ⟶ _)) := by
  obtain ⟨k, α, h⟩ := hom_pf_pow_toRootHom (F := F') hisoB hAB hfnB hfnB' ζ hdeg hprop x hx
  refine ⟨k, α, ?_⟩
  rw [← map_pow, h, endRootOneEquiv_toRootHom]
  rfl

/-- ★★**`𝒪^×` の上では `Div^gp` は冪を倍にする**(`biratDivGp_mul_otimes` の帰納)。 -/
theorem biratDivGp_pow_otimes {X : BiratCat P G} {δ : End X}
    (hδ : δ ∈ OTimes (biratPre P G) X) (n : ℕ) :
    biratDivGp ((δ ^ n : End X) : X ⟶ X) = n • biratDivGp ((δ : End X) : X ⟶ X) := by
  induction n with
  | zero =>
      show biratDivGp (((1 : End X)) : X ⟶ X) = 0
      exact biratDivGp_id X
  | succ n ih =>
      rw [pow_succ, biratDivGp_mul_otimes ((OTimes (biratPre P G) X).pow_mem hδ n) hδ, ih,
        succ_nsmul]

/-- ★★★**`𝒞^birat` の isotropic 対象では `𝒪^▷ = 𝒪^×`**。

★pre-step は `𝒞^birat` で同型になる(`birat_isIso_of_preStep_of_isotropic`)ので、
`𝒪^▷` の元は自動的に単元である。 -/
theorem otri_mem_otimes_birat
    (hfnBirat : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {X : BiratCat P G} (hX : IsIsotropic (biratPre P G) X) {α : End X}
    (hα : α ∈ OTri (biratPre P G) X) : α ∈ OTimes (biratPre P G) X := by
  have hb : (biratPre P G).Base ((α : End X) : X ⟶ X) = 𝟙 _ := by
    have h := hα.1
    rw [show (biratPre P G).Base ((α : End X) : X ⟶ X)
      = (biratPre P G).Base (𝟙 X) from h, (biratPre P G).Base_id]
  have hps : IsPreStep (biratPre P G) ((α : End X) : X ⟶ X) := by
    refine ⟨hα.2, ?_⟩
    show IsIso ((biratPre P G).Base ((α : End X) : X ⟶ X))
    rw [hb]
    infer_instance
  haveI := birat_isIso_of_preStep_of_isotropic P G hfnBirat hX hps
  refine ⟨hα, ⟨⟨α, @inv _ _ _ _ ((α : End X) : X ⟶ X) _, ?_, ?_⟩, rfl⟩⟩
  · exact IsIso.inv_hom_id _
  · exact IsIso.hom_inv_id _

/-- ★★★★★★**逆向きの包含** —— `(Φ^pf)^birat` の元は**分母を払えば `Φ^birat` に戻る**。

★★★これで `Proposition 5.5, (iv)` の単系の同定が両側から挟まれた:

* 順向き(`phiBiratOn_pf_mem`)—— `Φ^birat` の像は `(Φ^pf)^birat` に入る。
* 逆向き(**本定理**)—— `(Φ^pf)^birat` の元は正の整数倍で `Φ^birat` の像に入る。

★筋は 3 段: `Proposition 5.5, (ii)` の全単射で `𝒞^birat` の完全化へ移し(`x`)、
そこで `Proposition 5.5, (i)` の「正の冪で戻る」を当て(`biratPfHom_pow_eq_Xi`)、
`Div^gp` が `𝒪^×` の上で加法的であること(`biratDivGp_pow_otimes`)で分母を払う。 -/
theorem biratDivGp_nsmul_mem_pfImage (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C)
    (hisoB : ∀ X : BiratCat P G, IsIsotropic (biratPre P G) X)
    (hfnBirat : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hAB : IsFrobeniusTrivial (biratPre P G) (biratUp P G A))
    (hfnB : IsFrobeniusNormalized (biratPre P G) (biratUp P G A))
    (hfnB' : IsFrobeniusNormalized (biratPre P G)
      (rtObj (biratPre P G) F' (biratUp P G A) 1))
    (ζ : ℕ+ →* End (biratUp P G A))
    (hdeg : ∀ n : ℕ+, (biratPre P G).degFr ((ζ n : End (biratUp P G A)) :
      biratUp P G A ⟶ biratUp P G A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity (biratPre P G) (ζ n)
      ∧ IsFrobeniusType (biratPre P G) ((ζ n : End (biratUp P G A)) : _ ⟶ _))
    (x : End (pfObjOf (biratPre P G) F' (biratUp P G A)))
    (hδ : biratPfEndHom hfi Gpf F' A x
      ∈ OTimes (biratPre (pfRootPre P F) Gpf)
        (biratUp (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F))) :
    ∃ k : ℕ+, ((k : ℕ+) : ℕ) • biratDivGp (P := pfRootPre P F) (G := Gpf)
        ((biratPfEndHom hfi Gpf F' A x :
          End (biratUp (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F))) : _ ⟶ _)
      ∈ AddSubgroup.map (gpMap _ (Pf.of (M := Φ.val (P.toElem.obj A).base)))
        (phiBiratAt P G (biratUp P G A)) := by
  have hx := otri_of_otimes_biratPfHom hfi Gpf F' A hδ
  obtain ⟨k, α, hk⟩ := biratPfHom_pow_eq_Xi hfi Gpf F' A hisoB hAB hfnB hfnB' ζ hdeg hprop x hx
  refine ⟨k, ?_⟩
  rw [← biratDivGp_pow_otimes hδ ((k : ℕ+) : ℕ), hk]
  refine ⟨biratDivGp ((α : End (biratUp P G A)) : biratUp P G A ⟶ biratUp P G A),
    ⟨(α : End (biratUp P G A)), otri_mem_otimes_birat hfnBirat (hisoB _) α.2, rfl⟩, ?_⟩
  exact (biratDivGp_biratPfHom hfi Gpf F' A A
    ((α : End (biratUp P G A)) : biratUp P G A ⟶ biratUp P G A)).symm

/-- ★★★★★★**部分群の包含として** —— `Φ^birat` の像は `(Φ^pf)^birat` の部分群に入る。

★`Prop53QPhi.lean` の `phiBiratPfImage`(`Φ^birat` の `Pf` への像)の言葉で述べたもの。
★★これで `Proposition 5.3` の `ℚ·Φ^birat = qPhiBiratOn`(飽和)との比較が
**部分群の層で**できるようになる。 -/
theorem phiBiratPfImage_le_phiBiratOn (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C) :
    phiBiratPfImage P G ((P.toElem.obj A).base)
      ≤ phiBiratOn Gpf (((pfRootPre P F).toElem.obj (⟨A, 1⟩ : PfRootObj P F)).base) := by
  rintro _ ⟨y, hy, rfl⟩
  exact phiBiratOn_pf_mem hfi hiso Gpf F' A hy

end Birat

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Proposition 5.5, (iv)` の単系の同定のうち、
「`Φ^birat` の像は `(Φ^pf)^birat` に入る」側(`Div^gp` の代表元の計算)。 -/
def sliceDivGpOf_toRootHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iv) — 𝒞 → 𝒞^pf は双有理射の Div^gp を Pf.of で移す",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★★locator —— `Proposition 5.5, (iv)` の単系の同定の**順向き**
(`Φ^birat` の像は `(Φ^pf)^birat` に入る)。 -/
def phiBiratOn_pf_mem.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iv) — Φ^birat の像は (Φ^pf)^birat に入る",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★★locator —— `Proposition 5.5, (iv)` の単系の同定の**逆向き**
(`(Φ^pf)^birat` の元は分母を払えば `Φ^birat` に戻る)。
★仮定は `δ' ∈ 𝒪^×` だけになった(根つき化身の `𝒪^▷` は `otri_of_otimes_biratPfHom` で出る)。 -/
def biratDivGp_nsmul_mem_pfImage.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iv) — (Φ^pf)^birat の元は分母を払えば Φ^birat に戻る",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★locator —— `Proposition 5.5, (ii)` の全単射を自己射の単系の準同型に束ねたもの。 -/
def biratPfEndHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 射の全単射は自己射の単系の準同型",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
