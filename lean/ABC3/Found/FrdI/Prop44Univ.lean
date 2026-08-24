/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44

/-!
# [FrdI] Proposition 4.4 —— `Hom^birat` の**普遍性**

原文 (FrdI p.82):
> (Birationalization of a Frobenioid I) For A, B

★★在庫の `HomBirat` には `mk` / `exists_rep` / `sound` / `eq_iff`
(帰納極限**への**射と等号判定)しかなく、**極限からの射**(普遍性)が無かった。
本ファイルはそれを置く。

## ★なぜ要るか

`Corollary 4.10` の `Ψ^birat`(`psiBiratCor`)は**手で組んである**
(`biratPsiMap` / `biratPsiMap_id` / `biratPsiMap_comp`)。
★同じことを毎回手で組むのは高くつく ——
`Proposition 5.5, (ii)` の梱包でも同じ形が要る
(`Prop55ScaleRootCoa.lean` の「逃げ道 B」)。

★★普遍性があれば「**co-angular pre-step を同型に送る写像は `Hom^birat` を経由する**」
と 1 行で書ける。`HomBirat = HomColim (homFunctorBirat …)` で
`HomColim` は mathlib の `colimit` なので、`HomColim.desc` がそのまま乗る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 uE

section BiratUniv

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {G : Frobenioid P}

/-- ★遷移射との両立を `HomColim` の形に持ち上げたもの。 -/
theorem HomBirat.descAux {A B : C} {T : Type (max v u2 v2)}
    (f : ∀ Z : IdxBirat P G A, (Z.unop.left.obj ⟶ B) → T)
    (hf : ∀ {Z W : IdxBirat P G A} (u : Z ⟶ W) (φ : Z.unop.left.obj ⟶ B),
      f W (u.unop.left.hom ≫ φ) = f Z φ)
    {Z W : IdxBirat P G A} (u : Z ⟶ W) (x : (homFunctorBirat P G A B).obj Z) :
    f W ((homFunctorBirat P G A B).map u x).down = f Z x.down := by
  rw [homFunctorBirat_map]
  exact hf u x.down

/-- ★★★★★**`Hom^birat` の普遍性** ——
各添字での写像が遷移射(＝前合成)と両立すれば、`Hom^birat` から降りる。 -/
noncomputable def HomBirat.desc {A B : C} {T : Type (max v u2 v2)}
    (f : ∀ Z : IdxBirat P G A, (Z.unop.left.obj ⟶ B) → T)
    (hf : ∀ {Z W : IdxBirat P G A} (u : Z ⟶ W) (φ : Z.unop.left.obj ⟶ B),
      f W (u.unop.left.hom ≫ φ) = f Z φ) :
    HomBirat P G A B → T :=
  HomColim.desc (homFunctorBirat P G A B) (fun Z x => f Z x.down)
    (fun {_ _} u x => HomBirat.descAux f hf u x)

@[simp] theorem HomBirat.desc_mk {A B : C} {T : Type (max v u2 v2)}
    (f : ∀ Z : IdxBirat P G A, (Z.unop.left.obj ⟶ B) → T)
    (hf : ∀ {Z W : IdxBirat P G A} (u : Z ⟶ W) (φ : Z.unop.left.obj ⟶ B),
      f W (u.unop.left.hom ≫ φ) = f Z φ)
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) :
    HomBirat.desc f hf (HomBirat.mk Z φ) = f Z φ :=
  HomColim.desc_mk (homFunctorBirat P G A B) (fun Z x => f Z x.down)
    (fun {_ _} u x => HomBirat.descAux f hf u x) Z (ULift.up φ)

/-- ★★**`Hom^birat` からの 2 つの写像は代表元で決まる**。 -/
theorem HomBirat.desc_ext {A B : C} {T : Type (max v u2 v2)}
    (g h : HomBirat P G A B → T)
    (hgh : ∀ (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B),
      g (HomBirat.mk Z φ) = h (HomBirat.mk Z φ)) : g = h :=
  HomColim.desc_ext (homFunctorBirat P G A B) g h (fun Z x => hgh Z x.down)

/-! ## ★2. `co-angular pre-step` を同型に送る関手は `Hom^birat` を経由する -/

variable {E : Type uE} [Category.{max v u2 v2} E]

/-! ★分数 `φ ∘ (δ)⁻¹`(`δ` は co-angular pre-step)を
`Ω φ ∘ (Ω δ)⁻¹` に送るだけ。★well-defined 性は遷移射(前合成)との両立で、
`Over.w` と「`Ω` が同型に送る」ことから出る。 -/

/-- ★代表元での値 —— 分数 `φ ∘ δ⁻¹` を `Ω φ ∘ (Ω δ)⁻¹` に送る。 -/
noncomputable def biratDescFun (Ω : C ⥤ E)
    (hΩ : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso (Ω.map f))
    {A B : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) : Ω.obj A ⟶ Ω.obj B :=
  @inv _ _ _ _ (Ω.map Z.unop.hom.hom) (hΩ _ Z.unop.hom.property) ≫ Ω.map φ

/-- ★遷移射(前合成)との両立。 -/
theorem biratDescFun_compat (Ω : C ⥤ E)
    (hΩ : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso (Ω.map f))
    {A B : C} {Z W : IdxBirat P G A} (u : Z ⟶ W) (φ : Z.unop.left.obj ⟶ B) :
    biratDescFun (G := G) Ω hΩ W (u.unop.left.hom ≫ φ)
      = biratDescFun (G := G) Ω hΩ Z φ := by
  haveI hZ : IsIso (Ω.map Z.unop.hom.hom) := hΩ _ Z.unop.hom.property
  haveI hW : IsIso (Ω.map W.unop.hom.hom) := hΩ _ W.unop.hom.property
  -- ★`Over.w` —— `u.left ≫ Z.hom = W.hom`
  have hw : u.unop.left.hom ≫ Z.unop.hom.hom = W.unop.hom.hom :=
    congrArg (fun t : W.unop.left ⟶ coaPreObj P G A => t.hom) (Over.w u.unop)
  have hΩw : Ω.map u.unop.left.hom ≫ Ω.map Z.unop.hom.hom = Ω.map W.unop.hom.hom := by
    rw [← Ω.map_comp, hw]
  show @inv _ _ _ _ (Ω.map W.unop.hom.hom) hW ≫ Ω.map (u.unop.left.hom ≫ φ)
    = @inv _ _ _ _ (Ω.map Z.unop.hom.hom) hZ ≫ Ω.map φ
  refine (cancel_epi (Ω.map W.unop.hom.hom)).mp ?_
  have hL : Ω.map W.unop.hom.hom ≫ (@inv _ _ _ _ (Ω.map W.unop.hom.hom) hW
        ≫ Ω.map (u.unop.left.hom ≫ φ))
      = Ω.map (u.unop.left.hom ≫ φ) := by
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  have hR : Ω.map W.unop.hom.hom ≫ (@inv _ _ _ _ (Ω.map Z.unop.hom.hom) hZ ≫ Ω.map φ)
      = Ω.map u.unop.left.hom ≫ Ω.map φ := by
    rw [← hΩw, Category.assoc, ← Category.assoc (Ω.map Z.unop.hom.hom),
      IsIso.hom_inv_id, Category.id_comp]
  rw [hL, hR, Ω.map_comp]

/-- ★★★★★★**`Hom^birat` の普遍性(本体)** ——
`Ω : 𝒞 ⥤ E` が co-angular pre-step を**同型に送る**なら、
`Hom^birat_𝒞(A,B) → Hom_E(Ω A, Ω B)` が定まる。 -/
noncomputable def biratDescHom (Ω : C ⥤ E)
    (hΩ : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso (Ω.map f))
    {A B : C} : HomBirat P G A B → (Ω.obj A ⟶ Ω.obj B) :=
  HomBirat.desc (biratDescFun (G := G) Ω hΩ)
    (fun {_ _} u φ => biratDescFun_compat Ω hΩ u φ)

@[simp] theorem biratDescHom_mk (Ω : C ⥤ E)
    (hΩ : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso (Ω.map f))
    {A B : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) :
    biratDescHom (G := G) Ω hΩ (HomBirat.mk Z φ) = biratDescFun (G := G) Ω hΩ Z φ :=
  HomBirat.desc_mk (biratDescFun (G := G) Ω hΩ)
    (fun {_ _} u φ => biratDescFun_compat Ω hΩ u φ) Z φ

/-! ## ★3. 恒等と合成 —— 関手にする

★★配管メモ(idiom 14)—— `inv` を含む式を `rw` で動かすと
「パターンが見つからない」「motive が型付かない」に必ず当たる。
`IsIso` の引数がメタ変数だとインスタンス探索が走らないため。
★★逃げ方は「**圏 `E` の中だけの抽象補題**にして、`IsIso` を
インスタンス束縛にしてしまう」こと。下の 2 つがそれ。 -/

section FracAux

universe uv

variable {E : Type uE} [Category.{uv} E]

/-- ★分数の要 —— 四角形 `g ≫ p = a ≫ w` から `g⁻¹ ≫ a = p ≫ w⁻¹`。 -/
theorem frac_key_aux {X Y T U : E} (g : X ⟶ Y) [IsIso g] (a : X ⟶ T) (p : Y ⟶ U)
    (w : T ⟶ U) [IsIso w] (hsq : g ≫ p = a ≫ w) :
    inv g ≫ a = p ≫ inv w := by
  rw [IsIso.inv_comp_eq, ← Category.assoc, hsq, Category.assoc, IsIso.hom_inv_id,
    Category.comp_id]

/-- ★分数の合成 —— `(g ≫ z)⁻¹ ≫ (a ≫ q) = (z⁻¹ ≫ p) ≫ (w⁻¹ ≫ q)`。 -/
theorem frac_comp_aux {X Y Z T U V : E} (g : X ⟶ Y) [IsIso g] (z : Y ⟶ Z) [IsIso z]
    (a : X ⟶ T) (p : Y ⟶ U) (w : T ⟶ U) [IsIso w] (q : T ⟶ V)
    (gz : X ⟶ Z) [IsIso gz] (aq : X ⟶ V)
    (hgz : gz = g ≫ z) (haq : aq = a ≫ q)
    (key : inv g ≫ a = p ≫ inv w) :
    inv gz ≫ aq = (inv z ≫ p) ≫ (inv w ≫ q) := by
  subst hgz
  subst haq
  rw [IsIso.inv_comp]
  simp only [Category.assoc]
  rw [← Category.assoc (inv g) a q, key, Category.assoc]

/-- ★`e = 𝟙` なら `e⁻¹ ≫ k = k`。 -/
theorem inv_id_comp_aux {X Y : E} (e : X ⟶ X) [IsIso e] (he : e = 𝟙 X) (k : X ⟶ Y) :
    inv e ≫ k = k := by
  subst he
  rw [IsIso.inv_id, Category.id_comp]

end FracAux

theorem biratDescHom_id (Ω : C ⥤ E)
    (hΩ : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso (Ω.map f)) (A : C) :
    biratDescHom (G := G) Ω hΩ (toHomBirat (P := P) (G := G) (𝟙 A)) = 𝟙 (Ω.obj A) := by
  show biratDescHom (G := G) Ω hΩ (HomBirat.mk (idxBiratOne P G A) (𝟙 A)) = _
  rw [biratDescHom_mk]
  show @inv _ _ _ _ (Ω.map (𝟙 A)) (hΩ _ (idxBiratOne P G A).unop.hom.property)
      ≫ Ω.map (𝟙 A) = 𝟙 (Ω.obj A)
  exact IsIso.inv_hom_id _

theorem biratDescHom_comp (Ω : C ⥤ E)
    (hΩ : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso (Ω.map f))
    {A B E' : C} (f : HomBirat P G A B) (g : HomBirat P G B E') :
    biratDescHom (G := G) Ω hΩ (compBirat P G G.core f g)
      = biratDescHom (G := G) Ω hΩ f ≫ biratDescHom (G := G) Ω hΩ g := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  rw [compBirat_mk, biratDescHom_mk, biratDescHom_mk, biratDescHom_mk]
  haveI hZ : IsIso (Ω.map Z.unop.hom.hom) := hΩ _ Z.unop.hom.property
  haveI hW : IsIso (Ω.map W.unop.hom.hom) := hΩ _ W.unop.hom.property
  haveI hγ : IsIso (Ω.map (biratPullGamma G.core Z φ W)) :=
    hΩ _ ⟨biratPullGamma_coAngular G.core Z φ W, biratPullGamma_preStep G.core Z φ W⟩
  haveI hP : IsIso (Ω.map (biratPullGamma G.core Z φ W ≫ Z.unop.hom.hom)) :=
    hΩ _ (biratPullIdx G.core Z φ W).unop.hom.property
  -- ★引き戻しの四角形を `E` へ送る
  have hsq : Ω.map (biratPullGamma G.core Z φ W) ≫ Ω.map φ
      = Ω.map (biratPullAlpha G.core Z φ W) ≫ Ω.map W.unop.hom.hom := by
    rw [← Ω.map_comp, ← Ω.map_comp]
    exact congrArg Ω.map (biratPull_sq G.core Z φ W)
  -- ★★インスタンス引数は `@` で明示する ——
  -- 暗黙の対象がメタ変数のうちに探索が走って落ちるため(idiom 14 の続き)。
  have key : @inv _ _ _ _ (Ω.map (biratPullGamma G.core Z φ W)) hγ
        ≫ Ω.map (biratPullAlpha G.core Z φ W)
      = Ω.map φ ≫ @inv _ _ _ _ (Ω.map W.unop.hom.hom) hW :=
    @frac_key_aux _ _ _ _ _ _ (Ω.map (biratPullGamma G.core Z φ W)) hγ
      (Ω.map (biratPullAlpha G.core Z φ W)) (Ω.map φ) (Ω.map W.unop.hom.hom) hW hsq
  exact @frac_comp_aux _ _ _ _ _ _ _ _
    (Ω.map (biratPullGamma G.core Z φ W)) hγ (Ω.map Z.unop.hom.hom) hZ
    (Ω.map (biratPullAlpha G.core Z φ W)) (Ω.map φ) (Ω.map W.unop.hom.hom) hW (Ω.map ψ)
    (Ω.map (biratPullGamma G.core Z φ W ≫ Z.unop.hom.hom)) hP
    (Ω.map (biratPullAlpha G.core Z φ W ≫ ψ))
    (Ω.map_comp _ _) (Ω.map_comp _ _) key

/-- ★★★★★★★**`𝒞^birat` の普遍性** ——
co-angular pre-step を同型に送る関手は `𝒞^birat` を経由する。

★★これで `Corollary 4.10` の `Ψ^birat` や `Proposition 5.5, (ii)` の梱包を
**手で組まずに済む**。 -/
noncomputable def biratDescFunctor (Ω : C ⥤ E)
    (hΩ : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso (Ω.map f)) :
    BiratCat P G ⥤ E where
  obj A := Ω.obj (biratDown P G A)
  map f := biratDescHom (G := G) Ω hΩ f
  map_id A := biratDescHom_id Ω hΩ (biratDown P G A)
  map_comp f g := biratDescHom_comp Ω hΩ f g

/-- ★★**`𝒞 → 𝒞^birat` と 1-可換**。 -/
noncomputable def toBiratCatCompBiratDescFunctor (Ω : C ⥤ E)
    (hΩ : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso (Ω.map f)) :
    toBiratCat P G ⋙ biratDescFunctor (G := G) Ω hΩ ≅ Ω :=
  NatIso.ofComponents (fun _ => Iso.refl _) (by
    intro X Y f
    have hmap : biratDescHom (G := G) Ω hΩ (toHomBirat (P := P) (G := G) f) = Ω.map f := by
      show biratDescHom (G := G) Ω hΩ (HomBirat.mk (idxBiratOne P G X) f) = _
      rw [biratDescHom_mk]
      show @inv _ _ _ _ (Ω.map (𝟙 X)) (hΩ _ (idxBiratOne P G X).unop.hom.property)
          ≫ Ω.map f = Ω.map f
      exact @inv_id_comp_aux _ _ _ _ (Ω.map (𝟙 X))
        (hΩ _ (idxBiratOne P G X).unop.hom.property) (Ω.map_id X) (Ω.map f)
    show biratDescHom (G := G) Ω hΩ (toHomBirat (P := P) (G := G) f) ≫ 𝟙 (Ω.obj Y)
      = 𝟙 (Ω.obj X) ≫ Ω.map f
    rw [Category.comp_id, Category.id_comp, hmap])

end BiratUniv
