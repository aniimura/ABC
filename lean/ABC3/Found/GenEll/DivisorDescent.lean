/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SectionDescent

/-!
# [GenEll] Remark 1.5.1 —— **因子の降下（閉部分スキームの言葉で）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★★★mathlib の欠落を迂回する

`SectionDescent.lean` で測ったとおり、**mathlib には一般の射に対する
`ideal_comap`（`(I.comap f).ideal U = (I.ideal V).map (f.appLE V U)`）が無い**。
★イデアル層の言葉で因子を降ろすには、まずその補題を作らねばならない。

★★★**しかし因子を「閉部分スキーム」で書けば、その補題は要らない。**

因子 `D ⊆ X` は閉埋め込み `iZ : Z ⟶ X` であり、
「`D` が一緒に降りる」とは**可換正方形が有限段で立つ**ことである:

    Z_n ─ψ→ Z′_n
     │       │
    iZ_n    iZ′_n
     ↓       ↓
    X_n ─φ→ X′_n

★★`ψ` と `φ` は `IsoDescent.lean` の `bcDescHom`、
`iZ_n` は本ファイルの `bcBC`（`ℤ` 上の射の底変換）である。
★★★正方形が立つことは `descent_unique_baseChangeRatTower`（Stacks 01ZC 単射側）で出る
——**`ideal_comap` を経由しない**。

## ★★★★★★機構

`Spec ℤ[1/n!]` 成分は引き戻しの普遍性（`bcBC_fst` / `bcDescHom_fst`）だけで自動。
★残る `X′` 成分だけが降下の対象であり、そこは前ファイルまでの道具でそのまま流れる。

## ★残っている段

★★同型（`IsoDescent.lean`）と正方形（本ファイル）を**同じ段に揃える**組み立て。
★★★`Σ` の外での conductor の一致。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

/-! ## ★★★★★★`ℤ` 上の射の底変換 -/

/-- ★★★★**`ℤ` 上の射 `iZ : Z ⟶ X` を段 `n` へ底変換したもの**。

★因子 `D ⊆ X` の閉埋め込みを `X ×_ℤ ℤ[1/n!]` へ運ぶのがこれである。 -/
noncomputable def bcBC {Z X : Scheme.{0}}
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (f : X ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (n : ℕᵒᵖ) : bcObj g n ⟶ bcObj f n :=
  liftDescent g f (pullback.snd (overRatTowerDiagram.obj n).hom g ≫ iZ)

@[simp] theorem bcBC_snd {Z X : Scheme.{0}}
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (f : X ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (n : ℕᵒᵖ) :
    bcBC g f iZ n ≫ pullback.snd (overRatTowerDiagram.obj n).hom f
      = pullback.snd (overRatTowerDiagram.obj n).hom g ≫ iZ :=
  liftDescent_snd _ _ _

@[simp] theorem bcBC_fst {Z X : Scheme.{0}}
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (f : X ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (n : ℕᵒᵖ) :
    bcBC g f iZ n ≫ pullback.fst (overRatTowerDiagram.obj n).hom f
      = pullback.fst (overRatTowerDiagram.obj n).hom g :=
  liftDescent_fst _ _ _

/-- ★★段を下げても底変換は同じもの。 -/
theorem bcBC_map {Z X : Scheme.{0}}
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (f : X ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) {m n : ℕᵒᵖ} (h : m ⟶ n) :
    bcMap g h ≫ bcBC g f iZ n = bcBC g f iZ m ≫ bcMap f h := by
  apply pullback.hom_ext
  · rw [Category.assoc, bcBC_fst, bcMap_fst, Category.assoc, bcMap_fst, ← Category.assoc, bcBC_fst]
  · rw [Category.assoc, bcBC_snd, ← Category.assoc, bcMap_snd, Category.assoc,
      bcMap_snd, bcBC_snd]

/-- ★★★★**`ℤ` 上の射 `iZ : Z ⟶ X` を `X_ℚ` の水準へ底変換したもの**。 -/
noncomputable def bcBCpt {Z X : Scheme.{0}}
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (f : X ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) : bcPt g ⟶ bcPt f :=
  pullback.lift (pullback.fst (overRatTowerCone.pt).hom g)
    (pullback.snd (overRatTowerCone.pt).hom g ≫ iZ) (specZIsTerminal.hom_ext _ _)

@[simp] theorem bcBCpt_snd {Z X : Scheme.{0}}
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (f : X ⟶ Spec (CommRingCat.of ℤ)) (iZ : Z ⟶ X) :
    bcBCpt g f iZ ≫ pullback.snd (overRatTowerCone.pt).hom f
      = pullback.snd (overRatTowerCone.pt).hom g ≫ iZ := pullback.lift_snd _ _ _

@[simp] theorem bcBCpt_fst {Z X : Scheme.{0}}
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (f : X ⟶ Spec (CommRingCat.of ℤ)) (iZ : Z ⟶ X) :
    bcBCpt g f iZ ≫ pullback.fst (overRatTowerCone.pt).hom f
      = pullback.fst (overRatTowerCone.pt).hom g := pullback.lift_fst _ _ _

/-- ★★★★★★底変換は錐の脚と両立する。 -/
theorem bcLeg_bcBC {Z X : Scheme.{0}}
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (f : X ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (n : ℕᵒᵖ) :
    bcLeg g n ≫ bcBC g f iZ n = bcBCpt g f iZ ≫ bcLeg f n := by
  apply pullback.hom_ext
  · rw [Category.assoc, bcBC_fst, bcLeg_fst, Category.assoc, bcLeg_fst, ← Category.assoc,
      bcBCpt_fst]
  · rw [Category.assoc, bcBC_snd, ← Category.assoc, bcLeg_snd, Category.assoc,
      bcLeg_snd, bcBCpt_snd]

/-! ## ★★★★★★★★★★可換正方形の降下 -/

/-- ★★★★★★★★★★**可換正方形の降下** ——
`X_ℚ` の水準で立っている正方形は、**有限段で既に立っている**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

    Z_ℚ ─u→ Z′_ℚ            Z_n ─ψ→ Z′_n
     │       │       ⟹       │       │
     ↓       ↓                ↓       ↓
    X_ℚ ─v→ X′_ℚ            X_n ─φ→ X′_n

★★★**これが原文の「因子 `D` も一緒に `ℤ[Σ^{-1}]` へ延びる」段である**
——因子を閉部分スキームで書いたので、mathlib に無い `ideal_comap` を**経由しない**。

★機構: `Spec ℤ[1/n!]` 成分は引き戻しの普遍性だけで自動、
`X′` 成分だけを `descent_unique_baseChangeRatTower`（Stacks 01ZC 単射側）で降ろす。 -/
theorem exists_square_bcDescHom {Z Z' X X' : Scheme.{0}}
    [CompactSpace Z] [QuasiSeparatedSpace Z]
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (g' : Z' ⟶ Spec (CommRingCat.of ℤ))
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    [LocallyOfFiniteType f']
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X')
    (u : bcPt g ⟶ bcPt g') (v : bcPt f ⟶ bcPt f')
    (hsq : u ≫ bcBCpt g' f' iZ' = bcBCpt g f iZ ≫ v)
    {i : ℕᵒᵖ} (aZ : bcObj g i ⟶ Z') (a : bcObj f i ⟶ X')
    (haZ : bcLeg g i ≫ aZ = u ≫ pullback.snd (overRatTowerCone.pt).hom g')
    (ha : bcLeg f i ≫ a = v ≫ pullback.snd (overRatTowerCone.pt).hom f') :
    ∃ (n : ℕᵒᵖ) (h : n ⟶ i),
      bcDescHom g g' aZ h ≫ bcBC g' f' iZ' n = bcBC g f iZ n ≫ bcDescHom f f' a h := by
  have hL : bcLeg g i ≫ (aZ ≫ iZ')
      = bcBCpt g f iZ ≫ v ≫ pullback.snd (overRatTowerCone.pt).hom f' := by
    rw [← Category.assoc, haZ, Category.assoc, ← bcBCpt_snd g' f' iZ', ← Category.assoc, hsq,
      Category.assoc]
  have hR : bcLeg g i ≫ (bcBC g f iZ i ≫ a)
      = bcBCpt g f iZ ≫ v ≫ pullback.snd (overRatTowerCone.pt).hom f' := by
    rw [← Category.assoc, bcLeg_bcBC, Category.assoc, ha]
  obtain ⟨n, h1, h2, hn⟩ := descent_unique_baseChangeRatTower g f'
    (aZ ≫ iZ') (specZIsTerminal.hom_ext _ _)
    (bcBC g f iZ i ≫ a) (specZIsTerminal.hom_ext _ _) (hL.trans hR.symm)
  -- ★`ℕᵒᵖ` は前順序なので 2 本の道は等しい
  have h12 : h1 = h2 := Subsingleton.elim _ _
  subst h12
  have hn' : bcMap g h1 ≫ (aZ ≫ iZ') = bcMap g h1 ≫ (bcBC g f iZ i ≫ a) := hn
  refine ⟨n, h1, ?_⟩
  apply pullback.hom_ext
  · rw [Category.assoc, bcBC_fst, bcDescHom_fst, Category.assoc, bcDescHom_fst, bcBC_fst]
  · rw [Category.assoc, bcBC_snd, ← Category.assoc, bcDescHom_snd, Category.assoc,
      Category.assoc, bcDescHom_snd, hn', ← Category.assoc, bcBC_map, Category.assoc]

/-! ## ★★★★★★★★★★対 `(X, D)` の spreading out -/

/-- ★★正方形の等式を下の段へ運ぶ。 -/
theorem square_bcDescHom_of_le {Z Z' X X' : Scheme.{0}}
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (g' : Z' ⟶ Spec (CommRingCat.of ℤ))
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X')
    {iz ix m n : ℕᵒᵖ} (aZ : bcObj g iz ⟶ Z') (a : bcObj f ix ⟶ X')
    (hz : m ⟶ iz) (hx : m ⟶ ix) (h : n ⟶ m)
    (H : bcDescHom g g' aZ hz ≫ bcBC g' f' iZ' m = bcBC g f iZ m ≫ bcDescHom f f' a hx) :
    bcDescHom g g' aZ (h ≫ hz) ≫ bcBC g' f' iZ' n
      = bcBC g f iZ n ≫ bcDescHom f f' a (h ≫ hx) := by
  have H2 : bcMap g hz ≫ aZ ≫ iZ' = bcBC g f iZ m ≫ bcMap f hx ≫ a := by
    calc bcMap g hz ≫ aZ ≫ iZ'
        = (bcDescHom g g' aZ hz ≫ pullback.snd (overRatTowerDiagram.obj m).hom g') ≫ iZ' := by
          rw [bcDescHom_snd, Category.assoc]
      _ = bcDescHom g g' aZ hz ≫
            (bcBC g' f' iZ' m ≫ pullback.snd (overRatTowerDiagram.obj m).hom f') := by
          rw [bcBC_snd, Category.assoc]
      _ = (bcDescHom g g' aZ hz ≫ bcBC g' f' iZ' m) ≫
            pullback.snd (overRatTowerDiagram.obj m).hom f' := by rw [Category.assoc]
      _ = (bcBC g f iZ m ≫ bcDescHom f f' a hx) ≫
            pullback.snd (overRatTowerDiagram.obj m).hom f' := by rw [H]
      _ = bcBC g f iZ m ≫
            (bcDescHom f f' a hx ≫ pullback.snd (overRatTowerDiagram.obj m).hom f') := by
          rw [Category.assoc]
      _ = bcBC g f iZ m ≫ bcMap f hx ≫ a := by rw [bcDescHom_snd]
  apply pullback.hom_ext
  · rw [Category.assoc, bcBC_fst, bcDescHom_fst, Category.assoc, bcDescHom_fst, bcBC_fst]
  · rw [Category.assoc, bcBC_snd, ← Category.assoc, bcDescHom_snd, Category.assoc,
      Category.assoc, bcDescHom_snd, bcMap_comp, bcMap_comp, Category.assoc,
      Category.assoc, H2, ← Category.assoc, bcBC_map, Category.assoc]

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— 対 `(X, D)` の spreading out**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

`(X_ℚ, D_ℚ) ≅ (X′_ℚ, D′_ℚ)`（`Spec ℚ` 上、閉埋め込みと両立）なら、
**有限段 `n` で両方が同型になり、その同型は閉埋め込みと両立する**。

    Z_ℚ ─u→ Z′_ℚ            Z_n ─ψ→ Z′_n     (ψ は同型)
     │       │       ⟹       │       │
     ↓       ↓                ↓       ↓
    X_ℚ ─v→ X′_ℚ            X_n ─φ→ X′_n     (φ は同型)

★★★**これが原文の「同型がある有限素数集合 `Σ` の上で `ℤ[Σ⁻¹]` へ延びる」の全体である**
——因子ごと延びることを含む。

★★★★★組み立ての鍵は、同型の側（`exists_comp_id_of_desc`）と
正方形の側（`exists_square_bcDescHom`）を**同じ材料 `a`, `aZ` から作る**ことである。
そうすれば `ℕᵒᵖ` が前順序であることから**添字の道が自動的に一致する**
（`Subsingleton.elim`）。 -/
theorem exists_pair_descent {Z Z' X X' : Scheme.{0}}
    [CompactSpace Z] [QuasiSeparatedSpace Z] [CompactSpace Z'] [QuasiSeparatedSpace Z']
    [CompactSpace X] [QuasiSeparatedSpace X] [CompactSpace X'] [QuasiSeparatedSpace X']
    (g : Z ⟶ Spec (CommRingCat.of ℤ)) (g' : Z' ⟶ Spec (CommRingCat.of ℤ))
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    [LocallyOfFinitePresentation g] [LocallyOfFinitePresentation g']
    [LocallyOfFinitePresentation f] [LocallyOfFinitePresentation f']
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X')
    (e : bcPt f ≅ bcPt f')
    (he : e.hom ≫ pullback.fst (overRatTowerCone.pt).hom f'
      = pullback.fst (overRatTowerCone.pt).hom f)
    (eZ : bcPt g ≅ bcPt g')
    (heZ : eZ.hom ≫ pullback.fst (overRatTowerCone.pt).hom g'
      = pullback.fst (overRatTowerCone.pt).hom g)
    (hsq : eZ.hom ≫ bcBCpt g' f' iZ' = bcBCpt g f iZ ≫ e.hom) :
    ∃ (n : ℕᵒᵖ) (φ : bcObj f n ⟶ bcObj f' n) (φ' : bcObj f' n ⟶ bcObj f n)
      (ψ : bcObj g n ⟶ bcObj g' n) (ψ' : bcObj g' n ⟶ bcObj g n),
      φ ≫ φ' = 𝟙 (bcObj f n) ∧ φ' ≫ φ = 𝟙 (bcObj f' n) ∧
      ψ ≫ ψ' = 𝟙 (bcObj g n) ∧ ψ' ≫ ψ = 𝟙 (bcObj g' n) ∧
      ψ ≫ bcBC g' f' iZ' n = bcBC g f iZ n ≫ φ := by
  classical
  obtain ⟨ix, a, hai, -⟩ := exists_factor_baseChangeRatTower f f'
    (e.hom ≫ pullback.snd (overRatTowerCone.pt).hom f')
  obtain ⟨jx, b, hbj, -⟩ := exists_factor_baseChangeRatTower f' f
    (e.inv ≫ pullback.snd (overRatTowerCone.pt).hom f)
  obtain ⟨iz, aZ, haZ, -⟩ := exists_factor_baseChangeRatTower g g'
    (eZ.hom ≫ pullback.snd (overRatTowerCone.pt).hom g')
  obtain ⟨jz, bZ, hbZ, -⟩ := exists_factor_baseChangeRatTower g' g
    (eZ.inv ≫ pullback.snd (overRatTowerCone.pt).hom g)
  have hai2 : bcLeg f ix ≫ a = e.hom ≫ pullback.snd (overRatTowerCone.pt).hom f' := hai
  have hbj2 : bcLeg f' jx ≫ b = e.inv ≫ pullback.snd (overRatTowerCone.pt).hom f := hbj
  have haZ2 : bcLeg g iz ≫ aZ = eZ.hom ≫ pullback.snd (overRatTowerCone.pt).hom g' := haZ
  have hbZ2 : bcLeg g' jz ≫ bZ = eZ.inv ≫ pullback.snd (overRatTowerCone.pt).hom g := hbZ
  obtain ⟨nX, pX, pX', hXAB, hXBA⟩ := exists_comp_id_of_desc f f' e he a b hai2 hbj2
  obtain ⟨nZ, pZ, pZ', hZAB, hZBA⟩ := exists_comp_id_of_desc g g' eZ heZ aZ bZ haZ2 hbZ2
  -- ★正方形は共通の段で受ける
  have qz : IsCofiltered.min iz ix ⟶ iz := IsCofiltered.minToLeft iz ix
  have qx : IsCofiltered.min iz ix ⟶ ix := IsCofiltered.minToRight iz ix
  obtain ⟨nS, pS, hS⟩ := exists_square_bcDescHom g g' f f' iZ iZ' eZ.hom e.hom hsq
    (bcMap g qz ≫ aZ) (bcMap f qx ≫ a)
    (by rw [← Category.assoc, bcLeg_map, haZ2])
    (by rw [← Category.assoc, bcLeg_map, hai2])
  rw [bcDescHom_bcMap, bcDescHom_bcMap] at hS
  -- ★★共通の段へ落とす
  have r12 : IsCofiltered.min (IsCofiltered.min nX nZ) nS ⟶ IsCofiltered.min nX nZ :=
    IsCofiltered.minToLeft _ _
  have r3 : IsCofiltered.min (IsCofiltered.min nX nZ) nS ⟶ nS := IsCofiltered.minToRight _ _
  have r1 : IsCofiltered.min (IsCofiltered.min nX nZ) nS ⟶ nX :=
    r12 ≫ IsCofiltered.minToLeft nX nZ
  have r2 : IsCofiltered.min (IsCofiltered.min nX nZ) nS ⟶ nZ :=
    r12 ≫ IsCofiltered.minToRight nX nZ
  set n : ℕᵒᵖ := IsCofiltered.min (IsCofiltered.min nX nZ) nS
  have H1 := bcDescHom_comp_id_of_le f f' a b pX pX' r1 hXAB
  have H2 := bcDescHom_comp_id_of_le f' f b a pX' pX r1 hXBA
  have H3 := bcDescHom_comp_id_of_le g g' aZ bZ pZ pZ' r2 hZAB
  have H4 := bcDescHom_comp_id_of_le g' g bZ aZ pZ' pZ r2 hZBA
  have H5 := square_bcDescHom_of_le g g' f f' iZ iZ' aZ a (pS ≫ qz) (pS ≫ qx) r3 hS
  -- ★★★`ℕᵒᵖ` は前順序なので添字の道は一意
  have e1 : r3 ≫ pS ≫ qx = r1 ≫ pX := Subsingleton.elim _ _
  have e2 : r3 ≫ pS ≫ qz = r2 ≫ pZ := Subsingleton.elim _ _
  rw [e1, e2] at H5
  exact ⟨n, _, _, _, _, H1, H2, H3, H4, H5⟩

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Remark 1.5.1` には
★同型と正方形を同じ段に揃える組み立て、★★`Σ` の外での conductor の一致、が残っている。 -/

def exists_pair_descent.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(対 (X, D) の spreading out——conductor の一致は含まない)",
    sectionId := "genell-rem-1-5-1" }

def exists_square_bcDescHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(因子の降下——可換正方形が有限段で立つ。同型との段揃えは含まない)",
    sectionId := "genell-rem-1-5-1" }

def exists_square_bcDescHom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "descent_unique_baseChangeRatTower(降下の一意性、Stacks 01ZC 単射側)"
      (.inProject "ABC3" "ABC3.Found.GenEll.descent_unique_baseChangeRatTower") 9,
    .citation "[ABC3]" "liftDescent(引き戻しの普遍性で段を保って持ち上げる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.liftDescent") 9,
    .citation "[mathlib]" "Scheme.IdealSheafData.comapIso((I.comap f).subscheme ≅ pullback f I.subschemeι)"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.comapIso") 9,
    .implicitStep
      ("★★★★★因子を**閉部分スキーム**で書いたので、mathlib に無い ideal_comap を" ++
       "経由しない。★イデアル層の言葉で書くと " ++
       "(I.comap f).ideal U = (I.ideal V).map (f.appLE V U) が要るが、" ++
       "mathlib にあるのは ideal_comap_of_isOpenImmersion だけである" ++
       "(2026-08-27 実測、ResearchPaper/mathlib-gap.json の ideal-comap-on-affine-opens)") 9,
    .implicitStep
      ("★Spec ℤ[1/n!] 成分は引き戻しの普遍性(bcBC_fst / bcDescHom_fst)だけで自動。" ++
       "降下が要るのは X′ 成分だけである") 9,
    .implicitStep
      ("★★残る段: 同型(IsoDescent.lean)と正方形(本ファイル)を同じ段に揃える組み立て。" ++
       "★機構は min を 2 回取って Subsingleton.elim で道を同一視するという、" ++
       "exists_iso_baseChangeRatTower と同じ形である") 9,
    .implicitStep
      ("★★★残る段: Σ の外での conductor の一致。" ++
       "★これが揃えば Skeleton/GenEll/Section1.lean の remark_1_5_1 が受けている" ++
       "仮定 hagree を証明で置き換えられる") 9 ]

end ABC3.Found.GenEll
