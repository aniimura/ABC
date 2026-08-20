/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm51Slice

/-!
# [FrdI] Theorem 5.1, (i) の中心 —— pre-step の span と `Φ^birat`

`Theorem 5.1, (i)` は `Pic_Φ(A) ≅ Pic_𝒞(A)` という**全単射**を主張する。
その中身は、突き詰めると次の 1 個の**判定条件**である:

> `X` から出る 2 本の pre-step `φ : X → A'`、`ψ : X → A` に対し、
> `A'` と `A` が(底で `α = Base(φ)⁻¹ ∘ Base(ψ)` と整合して)同型であることと、
> `Div(ψ) − Div(φ) ∈ Φ^birat(X)` であることは同値。

原文 (FrdI p.96):
> (φ : B → A, ψ : B → C)

原文 (FrdI p.96):
> to the object

★★**我々の `Φ^birat` は `Proposition 4.4, (iii)` で `𝒪^×(A^birat)` の像として
定義してある**(`phiBiratAt`)ので、原文が p.97–98 で長く展開する
「`Φ^birat` の定義に戻る」議論は、**`𝒞^birat` の同型 1 本**に置き換わる。
これは原文の意図(「`Φ^birat` は `𝒞^birat` の単元が測るずれ」)そのままである。

★本ファイルは `Φ^birat` の**2 つの向きの辞書**を作るところから始める:

| 向き | 補題 |
|---|---|
| base-equivalent な pre-step の対 から `Φ^birat` の元 | `mem_phiBiratAt_of_preStepPair` |
| `Φ^birat` の元 から そのような対 | `exists_preStepPair_of_mem_phiBiratAt` |
| `Φ^birat` は co-angular pre-step で輸送される | `mem_phiBiratAt_transport` |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★0. `Gp` と `Φ.map` の小道具 -/

section GpTools

variable {D : Type u} [Category.{v} D] {Φ : MonoidOn.{v, u, w} D}

theorem gpMap_phi_id {X : D} (x : Gp (Φ.val X)) : gpMap _ (Φ.map (𝟙 X)) x = x := by
  have hid : Φ.map (𝟙 X) = AddMonoidHom.id _ := by ext z; exact Φ.map_id _ z
  rw [hid]; simp [gpMap_id]

theorem gpMap_phi_inv_left {X Y : D} (g : X ⟶ Y) [IsIso g] (x : Gp (Φ.val X)) :
    gpMap _ (Φ.map g) (gpMap _ (Φ.map (inv g)) x) = x := by
  rw [gpMap_phi_comp, IsIso.hom_inv_id]; exact gpMap_phi_id x

theorem gpMap_phi_inv_right {X Y : D} (g : X ⟶ Y) [IsIso g] (x : Gp (Φ.val Y)) :
    gpMap _ (Φ.map (inv g)) (gpMap _ (Φ.map g) x) = x := by
  rw [gpMap_phi_comp, IsIso.inv_hom_id]; exact gpMap_phi_id x

end GpTools

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. `Div` と `Base` の小道具 -/

/-- pre-step を後ろに合成したときの `Div`。 -/
theorem Div_comp_preStep {W X A : C} (δ : W ⟶ X) (ψ : X ⟶ A) (h : P.degFr ψ = 1) :
    P.Div (δ ≫ ψ) = Φ.map (P.Base δ) (P.Div ψ) + P.Div δ := by
  rw [P.Div_comp, h]
  show _ + ((1 : ℕ+) : ℕ) • P.Div δ = _
  rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]

/-- 同型の底は同型。 -/
theorem base_isIso_of_iso {A B : C} (θ : A ≅ B) : IsIso (P.Base θ.hom) :=
  ⟨P.Base θ.inv, by rw [← P.Base_comp, θ.hom_inv_id]; exact P.Base_id A,
    by rw [← P.Base_comp, θ.inv_hom_id]; exact P.Base_id B⟩

/-- 同型の逆射の底は、底の逆射。 -/
theorem base_inv_of_iso {A B : C} (θ : A ≅ B) [inst : IsIso (P.Base θ.hom)] :
    P.Base θ.inv = @inv _ _ _ _ (P.Base θ.hom) inst := by
  refine (IsIso.inv_eq_of_hom_inv_id ?_).symm
  rw [← P.Base_comp, θ.hom_inv_id]; exact P.Base_id A

/-! ## ★2. `Φ^birat` の辞書(その 1)—— 対から元へ -/

/-- ★`[δ₁]⁻¹ ≫ [δ₂]` は `𝒞^birat` の同型。 -/
theorem birat_mk_isIso (G : Frobenioid P) {A W : C}
    (δ₁ δ₂ : W ⟶ A) (hc₁ : IsCoAngular P δ₁) (hs₁ : IsPreStep P δ₁)
    (hc₂ : IsCoAngular P δ₂) (hs₂ : IsPreStep P δ₂) :
    IsIso (show (toBiratCat P G).obj A ⟶ (toBiratCat P G).obj A from
      HomBirat.mk (idxBiratMk P G δ₁ hc₁ hs₁) δ₂) := by
  haveI h1 : IsIso ((toBiratCat P G).map δ₁) := birat_isIso_of_coaPre δ₁ hc₁ hs₁
  haveI h2 : IsIso ((toBiratCat P G).map δ₂) := birat_isIso_of_coaPre δ₂ hc₂ hs₂
  haveI : IsIso ((toBiratCat P G).map δ₁ ≫ (show (toBiratCat P G).obj A ⟶ (toBiratCat P G).obj A
      from HomBirat.mk (idxBiratMk P G δ₁ hc₁ hs₁) δ₂)) := by
    have hcomp : ((toBiratCat P G).map δ₁ ≫ (show (toBiratCat P G).obj A ⟶ (toBiratCat P G).obj A
        from HomBirat.mk (idxBiratMk P G δ₁ hc₁ hs₁) δ₂)) = (toBiratCat P G).map δ₂ :=
      birat_toHom_comp_mk (G := G) δ₁ hc₁ hs₁ δ₂
    rw [hcomp]; exact h2
  exact IsIso.of_isIso_comp_left ((toBiratCat P G).map δ₁) _

/-- ★★★★**base-equivalent な co-angular pre-step の対は `Φ^birat` の元を定める**。

原文 (FrdI p.97):
> there exists a pair of base-equivalent pre-steps δ1, δ2 : D → A such that
-/
theorem mem_phiBiratAt_of_preStepPair (G : Frobenioid P) {A W : C}
    (δ₁ δ₂ : W ⟶ A) (hc₁ : IsCoAngular P δ₁) (hs₁ : IsPreStep P δ₁)
    (hc₂ : IsCoAngular P δ₂) (hs₂ : IsPreStep P δ₂)
    (hbase : P.Base δ₁ = P.Base δ₂) :
    sliceDivGpOf (P := P) δ₁ hs₁.2 δ₂ ∈ phiBiratAt P G A := by
  haveI : IsIso (P.Base δ₁) := hs₁.2
  have key : inv (P.Base δ₁) ≫ P.Base δ₂ = 𝟙 ((P.toElem.obj A).base) := by
    rw [← hbase]; exact IsIso.inv_hom_id _
  refine ⟨HomBirat.mk (idxBiratMk P G δ₁ hc₁ hs₁) δ₂, ⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
  · refine Eq.trans ?_ ((biratPre P G).Base_id A).symm
    refine (biratBase_mk _ _).trans ?_
    refine (sliceBaseOf_eq _ _ _).trans ?_
    exact key
  · show biratDeg (HomBirat.mk (idxBiratMk P G δ₁ hc₁ hs₁) δ₂) = 1
    rw [biratDeg_mk]
    exact hs₂.1
  · exact (isUnit_iff_isIso _).mpr (birat_mk_isIso G δ₁ δ₂ hc₁ hs₁ hc₂ hs₂)
  · rw [biratDivGp_mk]
    exact sliceDivGpOf_congr rfl _ _ δ₂

/-! ## ★3. `Φ^birat` の辞書(その 2)—— 元から対へ -/

/-- ★★★★**`Φ^birat` の元は base-equivalent な pre-step の対から来る**。

★`𝒪^×(A^birat)` の元の代表 `(a, φ)` について、
base-identity から `Base φ = Base a`、linear から `degFr φ = 1` が出るので、
`φ` は `a` に base-equivalent な pre-step になる。 -/
theorem exists_preStepPair_of_mem_phiBiratAt (G : Frobenioid P) {A : C}
    {x : Gp (Φ.val (P.toElem.obj A).base)} (hx : x ∈ phiBiratAt P G A) :
    ∃ (W : C) (δ₁ δ₂ : W ⟶ A) (hs₁ : IsPreStep P δ₁),
      IsCoAngular P δ₁ ∧ IsPreStep P δ₂ ∧ P.Base δ₁ = P.Base δ₂ ∧
      sliceDivGpOf (P := P) δ₁ hs₁.2 δ₂ = x := by
  obtain ⟨δ, hδ, rfl⟩ := hx
  obtain ⟨Z, φ, hrep⟩ := HomBirat.exists_rep (↑δ : HomBirat P G A A)
  have hb := hδ.1.1
  have hd := hδ.1.2
  rw [← hrep] at hb hd
  haveI ha : IsIso (P.Base Z.unop.hom.hom) := Z.unop.hom.property.2.2
  have h0 : biratBase (HomBirat.mk Z φ) = 𝟙 ((P.toElem.obj A).base) :=
    hb.trans ((biratPre P G).Base_id A)
  rw [biratBase_mk, sliceBaseOf_eq] at h0
  have hbeq : P.Base Z.unop.hom.hom = P.Base φ :=
    (((IsIso.inv_comp_eq _).mp h0).trans (Category.comp_id _)).symm
  have hlin : P.degFr φ = 1 := (biratDeg_mk Z φ) ▸ hd
  refine ⟨Z.unop.left.obj, Z.unop.hom.hom, φ, Z.unop.hom.property.2, Z.unop.hom.property.1,
    ⟨hlin, show IsIso (P.Base φ) from hbeq ▸ ha⟩, hbeq, ?_⟩
  rw [← hrep, biratDivGp_mk]
  exact (sliceDivGpOf_congr rfl _ _ φ).symm

/-! ## ★4. `Φ^birat` の輸送 -/

/-- ★`𝒞^birat` の同型による輸送(両向き)。 -/
theorem mem_phiBiratAt_iso (G : Frobenioid P) {X A : BiratCat P G} (e : X ≅ A)
    (y : Gp (Φ.val (P.toElem.obj (biratDown P G A)).base)) :
    gpMap _ (Φ.map (biratBase e.hom)) y ∈ phiBiratAt P G X ↔ y ∈ phiBiratAt P G A := by
  have hinv : biratBase e.inv ≫ biratBase e.hom
      = 𝟙 ((P.toElem.obj (biratDown P G A)).base) := by
    show (biratPre P G).Base e.inv ≫ (biratPre P G).Base e.hom = _
    rw [← (biratPre P G).Base_comp, e.inv_hom_id]
    exact (biratPre P G).Base_id A
  constructor
  · intro h
    have h2 := phiBiratAt_map_le (P := P) (G := G) e.symm ⟨_, h, rfl⟩
    simp only [Iso.symm_hom] at h2
    have heq : gpMap _ (Φ.map (biratBase e.inv)) (gpMap _ (Φ.map (biratBase e.hom)) y) = y := by
      rw [gpMap_phi_comp, hinv]
      exact gpMap_phi_id y
    rwa [heq] at h2
  · intro h
    exact phiBiratAt_map_le (P := P) (G := G) e ⟨y, h, rfl⟩

/-- ★★★**`Φ^birat` は co-angular pre-step に沿って輸送される**。 -/
theorem mem_phiBiratAt_transport (G : Frobenioid P) {X A : C} (ν : X ⟶ A)
    (hc : IsCoAngular P ν) (hs : IsPreStep P ν)
    (y : Gp (Φ.val (P.toElem.obj A).base)) :
    gpMap _ (Φ.map (P.Base ν)) y ∈ phiBiratAt P G X ↔ y ∈ phiBiratAt P G A := by
  haveI hiso : IsIso ((toBiratCat P G).map ν) := birat_isIso_of_coaPre ν hc hs
  have h := mem_phiBiratAt_iso (P := P) (G := G) (X := X) (A := A)
    (@asIso (BiratCat P G) _ X A ((toBiratCat P G).map ν) hiso) y
  have hbh : biratBase ((toBiratCat P G).map ν) = P.Base ν := biratBase_toHomBirat ν
  exact hbh ▸ h

/-! ## ★5. ★★★★span の判定条件 —— `Theorem 5.1, (i)` の中心 -/

/-- ★span `(φ : X → A', ψ : X → A)` の**類** —— `Φ^gp(A')` の元。

原文 (FrdI p.96):
> Φ(φ)−1(Div(ψ) − Div(φ)) ∈ Φgp(A)
-/
noncomputable def spanCls {X A' A : C} (φ : X ⟶ A') (hφ : IsIso (P.Base φ)) (ψ : X ⟶ A) :
    Gp (Φ.val (P.toElem.obj A').base) :=
  haveI := hφ
  gpMap _ (Φ.map (inv (P.Base φ))) (toGp _ (P.Div ψ) - toGp _ (P.Div φ))

theorem spanCls_eq {X A' A : C} (φ : X ⟶ A') (hφ : IsIso (P.Base φ)) (ψ : X ⟶ A) :
    spanCls φ hφ ψ = haveI := hφ
      gpMap _ (Φ.map (inv (P.Base φ))) (toGp _ (P.Div ψ) - toGp _ (P.Div φ)) := rfl

/-- ★`spanCls` を `Base φ` で押し出すと差そのものに戻る。 -/
theorem gpMap_base_spanCls {X A' A : C} (φ : X ⟶ A') (hφ : IsIso (P.Base φ)) (ψ : X ⟶ A) :
    gpMap _ (Φ.map (P.Base φ)) (spanCls φ hφ ψ)
      = toGp _ (P.Div ψ) - toGp _ (P.Div φ) := by
  haveI := hφ
  exact gpMap_phi_inv_left (P.Base φ) _

/-- ★`degFr ψ = 1` なら `spanCls` は `sliceDivGpOf` と同じ。 -/
theorem spanCls_eq_sliceDivGpOf {X A' A : C} (φ : X ⟶ A') (hφ : IsIso (P.Base φ))
    (ψ : X ⟶ A) (hψ : P.degFr ψ = 1) :
    spanCls φ hφ ψ = sliceDivGpOf (P := P) φ hφ ψ := by
  rw [spanCls_eq, sliceDivGpOf_eq, hψ]
  show _ = gpMap _ (Φ.map (inv (P.Base φ)))
    (toGp _ (P.Div ψ) - ((1 : ℕ+) : ℕ) • toGp _ (P.Div φ))
  rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]

/-- ★★**底で整合する同型があれば、差は `Φ^birat(X)` に入る**。 -/
theorem span_mem_phiBirat_of_iso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {X A' A : C} (φ : X ⟶ A') (ψ : X ⟶ A)
    (hsφ : IsPreStep P φ) (hsψ : IsPreStep P ψ)
    (θ : A' ≅ A) (hθ : P.Base θ.hom = @inv _ _ _ _ (P.Base φ) hsφ.2 ≫ P.Base ψ) :
    toGp _ (P.Div ψ) - toGp _ (P.Div φ) ∈ phiBiratAt P G X := by
  haveI hbφ : IsIso (P.Base φ) := hsφ.2
  haveI hbψ : IsIso (P.Base ψ) := hsψ.2
  haveI hbθ : IsIso (P.Base θ.hom) := base_isIso_of_iso θ
  have hstθ : IsPreStep P θ.hom := ⟨degFr_of_isIso P θ.hom, hbθ⟩
  have hs₁ : IsPreStep P (φ ≫ θ.hom) := IsPreStep.comp P hsφ hstθ
  haveI hb₁ : IsIso (P.Base (φ ≫ θ.hom)) := hs₁.2
  have hb : P.Base (φ ≫ θ.hom) = P.Base ψ := by
    rw [P.Base_comp, hθ, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  have hd : P.Div (φ ≫ θ.hom) = P.Div φ := by
    rw [Div_comp_preStep _ _ (degFr_of_isIso P θ.hom),
      show P.Div θ.hom = 0 from isIsometric_of_isIso P θ.hom, map_zero, zero_add]
  have hmem := mem_phiBiratAt_of_preStepPair G (φ ≫ θ.hom) ψ
    (prop_1_4_i P _ (fun Y _ => hiso Y)) hs₁
    (prop_1_4_i P _ (fun Y _ => hiso Y)) hsψ hb
  have hinvb : (@inv _ _ _ _ (P.Base (φ ≫ θ.hom)) hs₁.2) = @inv _ _ _ _ (P.Base ψ) hbψ := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [hb, IsIso.hom_inv_id]
  have hval : sliceDivGpOf (P := P) (φ ≫ θ.hom) hs₁.2 ψ
      = gpMap _ (Φ.map (inv (P.Base ψ))) (toGp _ (P.Div ψ) - toGp _ (P.Div φ)) := by
    rw [sliceDivGpOf_eq, hinvb, hd, show P.degFr ψ = 1 from hsψ.1]
    show gpMap _ (Φ.map (inv (P.Base ψ)))
      (toGp _ (P.Div ψ) - ((1 : ℕ+) : ℕ) • toGp _ (P.Div φ)) = _
    rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  rw [hval] at hmem
  have h2 := (mem_phiBiratAt_transport G ψ (prop_1_4_i P _ (fun Y _ => hiso Y)) hsψ _).mpr hmem
  rwa [gpMap_phi_inv_left] at h2

/-- ★対から作った 2 本の合成の `Div` が一致する、という計算。 -/
theorem span_div_eq_of_pair {X A' A : C} (φ : X ⟶ A') (ψ : X ⟶ A)
    (hsφ : IsPreStep P φ) (hsψ : IsPreStep P ψ)
    {W : C} (δ₁ δ₂ : W ⟶ X)
    (hbeq : P.Base δ₁ = P.Base δ₂)
    (hstep2 : toGp _ (P.Div δ₂) - toGp _ (P.Div δ₁)
      = gpMap _ (Φ.map (P.Base δ₁)) (toGp _ (P.Div ψ) - toGp _ (P.Div φ))) :
    P.Div (δ₁ ≫ ψ) = P.Div (δ₂ ≫ φ) := by
  refine (P.divisorial (P.toElem.obj W).base).1.1 ?_
  rw [Div_comp_preStep _ _ hsψ.1, Div_comp_preStep _ _ hsφ.1, ← hbeq]
  rw [toGp_add, toGp_add, ← gpMap_toGp, ← gpMap_toGp]
  rw [map_sub] at hstep2
  have he : toGp _ (P.Div δ₂)
      = toGp _ (P.Div δ₁)
        + (gpMap _ (Φ.map (P.Base δ₁)) (toGp _ (P.Div ψ))
          - gpMap _ (Φ.map (P.Base δ₁)) (toGp _ (P.Div φ))) := by
    rw [← hstep2]; abel
  rw [he]
  abel

/-- ★★**差が `Φ^birat(X)` に入れば、底で整合する同型がある**。 -/
theorem span_iso_of_mem_phiBirat (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {X A' A : C} (φ : X ⟶ A') (ψ : X ⟶ A)
    (hsφ : IsPreStep P φ) (hsψ : IsPreStep P ψ)
    (hmem : toGp _ (P.Div ψ) - toGp _ (P.Div φ) ∈ phiBiratAt P G X) :
    ∃ θ : A' ≅ A, P.Base θ.hom = @inv _ _ _ _ (P.Base φ) hsφ.2 ≫ P.Base ψ := by
  haveI hbφ : IsIso (P.Base φ) := hsφ.2
  haveI hbψ : IsIso (P.Base ψ) := hsψ.2
  obtain ⟨W, δ₁, δ₂, hs₁, hc₁, hs₂, hbeq, hval⟩ := exists_preStepPair_of_mem_phiBiratAt G hmem
  haveI hb₁ : IsIso (P.Base δ₁) := hs₁.2
  have hstep2 : toGp _ (P.Div δ₂) - toGp _ (P.Div δ₁)
      = gpMap _ (Φ.map (P.Base δ₁)) (toGp _ (P.Div ψ) - toGp _ (P.Div φ)) := by
    rw [← hval, sliceDivGpOf_eq, show P.degFr δ₂ = 1 from hs₂.1]
    rw [show (((1 : ℕ+) : ℕ) • toGp (Φ.val (P.toElem.obj W).base) (P.Div δ₁))
      = toGp _ (P.Div δ₁) from by rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]]
    exact (gpMap_phi_inv_left (P.Base δ₁) _).symm
  have hdiv := span_div_eq_of_pair φ ψ hsφ hsψ δ₁ δ₂ hbeq hstep2
  obtain ⟨θ', heq⟩ := coaPreStep_iso_of_div_eq G (δ₁ ≫ ψ) (δ₂ ≫ φ)
    ⟨prop_1_4_i P _ (fun Y _ => hiso Y), IsPreStep.comp P hs₁ hsψ⟩
    ⟨prop_1_4_i P _ (fun Y _ => hiso Y), IsPreStep.comp P hs₂ hsφ⟩ hdiv
  haveI hbθ : IsIso (P.Base θ'.hom) := base_isIso_of_iso θ'
  have hbase : P.Base ψ ≫ P.Base θ'.hom = P.Base φ := by
    have h := congrArg P.Base heq
    rw [P.Base_comp, P.Base_comp, P.Base_comp, ← hbeq, Category.assoc] at h
    exact (cancel_epi (P.Base δ₁)).mp h
  refine ⟨θ'.symm, ?_⟩
  show P.Base θ'.inv = _
  rw [base_inv_of_iso θ']
  refine IsIso.inv_eq_of_hom_inv_id ?_
  have hθ' : P.Base θ'.hom = inv (P.Base ψ) ≫ P.Base φ := by
    rw [← hbase, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  rw [hθ', Category.assoc, ← Category.assoc (P.Base φ), IsIso.hom_inv_id, Category.id_comp,
    IsIso.inv_hom_id]

/-- ★★★★**span の判定条件**(`Theorem 5.1, (i)` の中心)。

原文 (FrdI p.96):
> on the one hand and to the element
-/
theorem span_iso_iff (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {X A' A : C} (φ : X ⟶ A') (ψ : X ⟶ A)
    (hsφ : IsPreStep P φ) (hsψ : IsPreStep P ψ) :
    (∃ θ : A' ≅ A, P.Base θ.hom = @inv _ _ _ _ (P.Base φ) hsφ.2 ≫ P.Base ψ)
      ↔ spanCls φ hsφ.2 ψ ∈ phiBiratAt P G A' := by
  have htr := mem_phiBiratAt_transport G φ (prop_1_4_i P _ (fun Y _ => hiso Y)) hsφ
    (spanCls φ hsφ.2 ψ)
  rw [gpMap_base_spanCls] at htr
  constructor
  · rintro ⟨θ, hθ⟩
    exact htr.mp (span_mem_phiBirat_of_iso G hiso φ ψ hsφ hsψ θ hθ)
  · intro h
    exact span_iso_of_mem_phiBirat G hiso φ ψ hsφ hsψ (htr.mpr h)

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.1, (i)` の中心(判定条件)。 -/
def span_iso_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 96, item := "Theorem 5.1, (i) — span の判定条件",
    sectionId := "frdi-thm-5-1" }

end ABC3.Found.FrdI
