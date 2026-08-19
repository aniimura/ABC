/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.AppA
import ABC3.Found.FrdI.Prop19

/-!
# [FrdI] Theorem 3.4, (v) —— 底の射の 3 分解

原文 (FrdI p.68):
> φD = Base(ψ) ◦ Base(γ) ◦ Base(α)−1

★★`Ψ_Base : 𝒟₁ → 𝒟₂` の構成には、`𝒟` の任意の射を
**pre-step 2 本と pull-back 1 本**に分解する必要がある。

## ★★★在庫 2 本で出る(2026-08-19 に測った)

| 段 | 在庫 |
|---|---|
| pull-back を取る | `plBkEquiv`(`Definition 1.3, (i), (c)`)の**本質的全射性** |
| 残りの同型を span にする | `preStepSpan`(`Definition 1.3, (i), (b)`) |

★`plBkOverFunctor A' : Over (⟨A'⟩ : PlBk P) ⥤ Over (Base A')` が圏同値なので、
`φ_𝒟 : Base A ⟶ Base A'` を `Over (Base A')` の対象と見て本質的全射性を使えば、
pull-back 射 `ψ : B ⟶ A'` と同型 `θ : Base B ≅ Base A` で
`Base ψ = θ ≫ φ_𝒟` なるものが取れる。
★あとは `preStepSpan` を `θ^{-1}` に当てるだけ。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★1. pull-back で `φ_𝒟` を実現する -/

include F in
/-- ★★★**`𝒟` の任意の射は pull-back の底に同型を継いだもの**。

★`plBkEquiv` の**本質的全射性**そのもの。 -/
theorem exists_pullBack_realizing {A A' : C}
    (φ𝒟 : (P.toElem.obj A).base ⟶ (P.toElem.obj A').base) :
    ∃ (B : C) (ψ : B ⟶ A') (_ : IsPullBack P ψ)
      (θ : (P.toElem.obj B).base ≅ (P.toElem.obj A).base),
      P.Base ψ = θ.hom ≫ φ𝒟 := by
  haveI := F.plBkEquiv A'
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := plBkOverFunctor P A')
    (Over.mk φ𝒟)
  refine ⟨Z.left.obj, Z.hom.hom, Z.hom.property,
    (Over.forget ((P.toElem.obj A').base)).mapIso e, ?_⟩
  exact (Over.w e.hom).symm

/-! ## ★2. 3 分解 -/

include F in
/-- ★★★★★**[FrdI] Theorem 3.4, (v) の 3 分解**。

原文 (FrdI p.68):
> φD = Base(ψ) ◦ Base(γ) ◦ Base(α)−1

★`α : X ⟶ A`、`γ : X ⟶ B` は pre-step、`ψ : B ⟶ A'` は pull-back。 -/
theorem base_three_factor {A A' : C}
    (φ𝒟 : (P.toElem.obj A).base ⟶ (P.toElem.obj A').base) :
    ∃ (X B : C) (α : X ⟶ A) (γ : X ⟶ B) (ψ : B ⟶ A')
      (hα : IsPreStep P α) (_ : IsPreStep P γ) (_ : IsPullBack P ψ),
      φ𝒟 = @inv _ _ _ _ (P.Base α) hα.2 ≫ P.Base γ ≫ P.Base ψ := by
  obtain ⟨B, ψ, hψ, θ, hθ⟩ := exists_pullBack_realizing P F φ𝒟
  haveI : IsIso θ.inv := inferInstance
  obtain ⟨X, α, γ, hα, hγ, hspan⟩ := F.preStepSpan A B θ.inv (by infer_instance)
  refine ⟨X, B, α, γ, ψ, hα, hγ, hψ, ?_⟩
  have h1 : θ.inv ≫ P.Base ψ = φ𝒟 := by
    rw [hθ, ← Category.assoc, θ.inv_hom_id, Category.id_comp]
  rw [← h1, hspan, Category.assoc]

def base_three_factor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — 底の射の 3 分解",
    sectionId := "frdi-thm-3-4" }

/-! ## ★3. 自然性は `Base` の像から `𝒟` 全体へ延びる

★★★3 分解があるので、**`Base` の像で自然なら `𝒟` のすべての射で自然**になる ——
`Base(α)^{-1}` の側は自然性の四角形を逆にするだけ。

★★これが `Ψ_Base` の **1-一意性**の芯である
(`Prop311.projPrecompIso` は `Φ = 0` の仮定を使うが、ここは一般の `Φ` で通る)。 -/

section Naturality

universe v3 u3

variable {D₂ : Type u3} [Category.{v3} D₂]

/-- ★同型の像での自然性は逆向きにも延びる。 -/
theorem naturality_inv {G G' : D ⥤ D₂} {Y Y' : D} (e : Y ⟶ Y') [IsIso e]
    {a : G.obj Y ⟶ G'.obj Y} {a' : G.obj Y' ⟶ G'.obj Y'}
    (h : G.map e ≫ a' = a ≫ G'.map e) :
    G.map (inv e) ≫ a = a' ≫ G'.map (inv e) := by
  haveI : IsIso (G.map e) := G.map_isIso e
  haveI : IsIso (G'.map e) := G'.map_isIso e
  have h1 : G.map (inv e) = inv (G.map e) := by
    rw [← CategoryTheory.Functor.map_inv]
  have h2 : G'.map (inv e) = inv (G'.map e) := by
    rw [← CategoryTheory.Functor.map_inv]
  rw [h1, h2]
  rw [IsIso.inv_comp_eq, ← Category.assoc, IsIso.eq_comp_inv]
  exact h.symm

include F in
/-- ★★★★★**`Base` の像で自然なら `𝒟` のすべての射で自然**。 -/
theorem naturality_of_base {G G' : D ⥤ D₂}
    (app : ∀ A : C, G.obj ((P.toElem.obj A).base) ⟶ G'.obj ((P.toElem.obj A).base))
    (hnat : ∀ {A A' : C} (f : A ⟶ A'),
      G.map (P.Base f) ≫ app A' = app A ≫ G'.map (P.Base f))
    {A A' : C} (φ𝒟 : (P.toElem.obj A).base ⟶ (P.toElem.obj A').base) :
    G.map φ𝒟 ≫ app A' = app A ≫ G'.map φ𝒟 := by
  obtain ⟨X, B, α, γ, ψ, hα, hγ, hψ, hdec⟩ := base_three_factor P F φ𝒟
  haveI : IsIso (P.Base α) := hα.2
  have hαinv : G.map (@inv _ _ _ _ (P.Base α) hα.2) ≫ app X
      = app A ≫ G'.map (@inv _ _ _ _ (P.Base α) hα.2) :=
    naturality_inv (P.Base α) (hnat α)
  subst hdec
  rw [CategoryTheory.Functor.map_comp, CategoryTheory.Functor.map_comp,
    CategoryTheory.Functor.map_comp, CategoryTheory.Functor.map_comp]
  simp only [Category.assoc]
  calc G.map (@inv _ _ _ _ (P.Base α) hα.2) ≫ G.map (P.Base γ) ≫ G.map (P.Base ψ) ≫ app A'
      = G.map (@inv _ _ _ _ (P.Base α) hα.2) ≫ G.map (P.Base γ) ≫ app B ≫ G'.map (P.Base ψ) := by
        rw [hnat ψ]
    _ = G.map (@inv _ _ _ _ (P.Base α) hα.2) ≫ app X ≫ G'.map (P.Base γ) ≫ G'.map (P.Base ψ) := by
        rw [← Category.assoc (G.map (P.Base γ)), hnat γ, Category.assoc]
    _ = app A ≫ G'.map (@inv _ _ _ _ (P.Base α) hα.2) ≫ G'.map (P.Base γ)
          ≫ G'.map (P.Base ψ) := by
        rw [← Category.assoc, hαinv, Category.assoc]

def naturality_of_base.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — 自然性は Base の像から 𝒟 全体へ延びる",
    sectionId := "frdi-thm-3-4" }

end Naturality

/-! ## ★4. ★★★★★一般の `projPrecompIso`

★★`Prop311.projPrecompIso` は `𝒞 ≌ 𝒟 × 𝒩`(`Φ = 0` の場合)を使うが、
**3 分解があれば一般の `Φ` で通る**。
★これが `Theorem 3.4, (v)` の `Ψ_Base` の **1-一意性**である。 -/

section Uniq

universe v3 u3

variable {D₂ : Type u3} [Category.{v3} D₂]

/-- ★`baseSurj` が選ぶ対象。 -/
noncomputable def chosenObj (Y : D) : C := (F.baseSurj Y).choose

/-- ★`baseSurj` が選ぶ同型。 -/
noncomputable def chosenIso (Y : D) :
    (P.toElem.obj (chosenObj P F Y)).base ≅ Y :=
  ((F.baseSurj Y).choose_spec.2).some

include F in
/-- ★★★★成分の自然性 —— `naturality_of_base` を同型で挟むだけ。 -/
theorem projPrecompIso_naturality {G G' : D ⥤ D₂} (h : P.proj ⋙ G ≅ P.proj ⋙ G')
    {A A' : C} {Y Y' : D}
    (eY : (P.toElem.obj A).base ≅ Y) (eY' : (P.toElem.obj A').base ≅ Y') (φ : Y ⟶ Y') :
    G.map φ ≫ (G.map eY'.inv ≫ h.hom.app A' ≫ G'.map eY'.hom)
      = (G.map eY.inv ≫ h.hom.app A ≫ G'.map eY.hom) ≫ G'.map φ := by
  have key := naturality_of_base P F (fun Z => h.hom.app Z)
    (fun {_ _} f => h.hom.naturality f) (A := A) (A' := A') (eY.hom ≫ φ ≫ eY'.inv)
  haveI : IsIso (G.map eY.hom) := G.map_isIso _
  refine (cancel_epi (G.map eY.hom)).mp ?_
  have hL : G.map eY.hom ≫ (G.map φ ≫ (G.map eY'.inv ≫ h.hom.app A' ≫ G'.map eY'.hom))
      = h.hom.app A ≫ G'.map eY.hom ≫ G'.map φ := by
    have s1 : G.map eY.hom ≫ (G.map φ ≫ (G.map eY'.inv ≫ h.hom.app A' ≫ G'.map eY'.hom))
        = (G.map (eY.hom ≫ φ ≫ eY'.inv) ≫ h.hom.app A') ≫ G'.map eY'.hom := by
      rw [CategoryTheory.Functor.map_comp, CategoryTheory.Functor.map_comp]
      simp only [Category.assoc]
    have hc : G'.map eY'.inv ≫ G'.map eY'.hom = 𝟙 _ := by
      rw [← CategoryTheory.Functor.map_comp, eY'.inv_hom_id, CategoryTheory.Functor.map_id]
    have s2inner : (G'.map eY.hom ≫ G'.map φ ≫ G'.map eY'.inv) ≫ G'.map eY'.hom
        = G'.map eY.hom ≫ G'.map φ := by
      rw [Category.assoc, Category.assoc, hc, Category.comp_id]
    have s2 : (h.hom.app A ≫ G'.map (eY.hom ≫ φ ≫ eY'.inv)) ≫ G'.map eY'.hom
        = h.hom.app A ≫ G'.map eY.hom ≫ G'.map φ := by
      rw [CategoryTheory.Functor.map_comp, CategoryTheory.Functor.map_comp, Category.assoc]
      exact congrArg (fun t => h.hom.app A ≫ t) s2inner
    exact s1.trans ((congrArg (fun t => t ≫ G'.map eY'.hom) key).trans s2)
  have hc2 : G.map eY.hom ≫ G.map eY.inv = 𝟙 _ := by
    rw [← CategoryTheory.Functor.map_comp, eY.hom_inv_id, CategoryTheory.Functor.map_id]
  have hR : G.map eY.hom ≫ ((G.map eY.inv ≫ h.hom.app A ≫ G'.map eY.hom) ≫ G'.map φ)
      = h.hom.app A ≫ G'.map eY.hom ≫ G'.map φ := by
    rw [Category.assoc, ← Category.assoc (G.map eY.hom) (G.map eY.inv), hc2,
      Category.id_comp, Category.assoc]
    rfl
  rw [hL, hR]

include F in
/-- ★★★★★**一般の 1-一意性** —— `P.proj` との前合成は関手の同型を反射する。 -/
noncomputable def projPrecompIsoGen {G G' : D ⥤ D₂}
    (h : P.proj ⋙ G ≅ P.proj ⋙ G') : G ≅ G' :=
  NatIso.ofComponents
    (fun Y => (G.mapIso (chosenIso P F Y)).symm
      ≪≫ h.app (chosenObj P F Y)
      ≪≫ G'.mapIso (chosenIso P F Y))
    (fun {Y Y'} φ =>
      projPrecompIso_naturality P F h (chosenIso P F Y) (chosenIso P F Y') φ)

def projPrecompIsoGen.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — Ψ_Base の 1-一意性(一般の Φ)",
    sectionId := "frdi-thm-3-4" }

end Uniq

/-! ## ★5. `Ψ` は pull-back のスライスに関手を誘導する

原文 (FrdI p.68):
> [i.e., which maps A → (Pi)A, (A → A) → {(Pi)A → (Pi)A}]. Similarly, the

★★★原典の `(P_i)_A` は **`A` へ入る pull-back 射のなす圏**であり、
`Definition 1.3, (i), (c)`(`plBkEquiv`)で `𝒟_{A_𝒟}` と同値になる。
★`Ψ` が pull-back 射を保てば、この圏の間の関手が立つ。 -/

section PlBkPsi

universe v4 u4 w4 u5 v5

variable {D₂ : Type u4} [Category.{v4} D₂] {C₂ : Type u5} [Category.{v5} C₂]
  {Φ₂ : MonoidOn.{v4, u4, w4} D₂} (P₂ : PreFrobenioid C₂ Φ₂)
  (Ψ : C ⥤ C₂)
  (hPB : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P f → IsPullBack P₂ (Ψ.map f))

include hPB in
/-- ★★**`Ψ` が誘導する `𝒞^pl-bk` の間の関手**。 -/
def plBkPsi : PlBk P ⥤ PlBk P₂ where
  obj X := ⟨Ψ.obj X.obj⟩
  map {X Y} f := ⟨Ψ.map f.hom, hPB f.hom f.property⟩
  map_id X := by
    refine WideSubcategory.hom_ext _ ?_
    show Ψ.map (𝟙 X.obj) = 𝟙 _
    exact CategoryTheory.Functor.map_id _ _
  map_comp {X Y Z} f g := by
    refine WideSubcategory.hom_ext _ ?_
    show Ψ.map (f.hom ≫ g.hom) = Ψ.map f.hom ≫ Ψ.map g.hom
    exact CategoryTheory.Functor.map_comp _ _ _

include hPB in
/-- ★★★**`(P₁)_A → (P₂)_{ΨA}`** —— 原典の「`(A → A) → {(P_i)_A → (P_i)_A}`」。 -/
def plBkSlicePsi (A : C) :
    Over (⟨A⟩ : PlBk P) ⥤ Over (⟨Ψ.obj A⟩ : PlBk P₂) :=
  Over.post (plBkPsi P P₂ Ψ hPB)

def plBkPsi.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — Ψ が pull-back のスライスに誘導する関手",
    sectionId := "frdi-thm-3-4" }

/-! ## ★6. 共役関手 `𝒟₁_{Base A} ⥤ 𝒟₂_{Base ΨA}`

原文 (FrdI p.68):
> [i.e., by inverting the equivalence of categories induced by α and then composing

★★★`plBkEquiv`(`Definition 1.3, (i), (c)`)で `𝒟_{A_𝒟}` を `(𝒞^pl-bk)_A` に置き換え、
`Ψ` で移して戻す。★**選択を含まない**(`plBkOverFunctor` の逆関手は
`IsEquivalence` インスタンスから来る canonical なもの)ので、
これ自体は well-defined である。 -/

include F hPB in
/-- ★★★★★**共役関手** —— `𝒟₁` のスライスを `𝒟₂` のスライスへ移す。

★これが原典の `Q_i → R_i` の「新しい関手」の実体である。 -/
noncomputable def basePsiSlice (A : C) :
    Over ((P.toElem.obj A).base) ⥤ Over ((P₂.toElem.obj (Ψ.obj A)).base) :=
  haveI := F.plBkEquiv A
  (plBkOverFunctor P A).inv ⋙ plBkSlicePsi P P₂ Ψ hPB A
    ⋙ plBkOverFunctor P₂ (Ψ.obj A)

def basePsiSlice.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — 𝒟 のスライスの間の共役関手",
    sectionId := "frdi-thm-3-4" }

/-! ## ★7. ★★★★★同型の同定は**終対象**から出る —— 選択が要らない

★★★`Over X` の対象 `f : Y ⟶ X` は **`f` が同型 ⟺ 終対象**である。
★`basePsiSlice` は圏同値 2 つと `Over.post` の合成なので**終対象を保つ**。
★★★したがって `θ : Base B ≅ Base A` を `Over (Base A)` の対象と見ると、
その像の構造射は**自動的に同型**になり、
`Base (ΨB) ≅ Base (ΨA)` の同定が**選択なしで**得られる。

★★これが `preStepSpan` の取り方への依存を消す鍵である。 -/

/-- ★`Over X` の対象は構造射が同型なら `Over.mk (𝟙 X)` と同型。 -/
noncomputable def overIsoMkId {X : D} (Z : Over X) (h : IsIso Z.hom) : Z ≅ Over.mk (𝟙 X) :=
  Over.isoMk (@asIso _ _ _ _ Z.hom h) (Category.comp_id Z.hom)

/-- ★★逆に `Over.mk (𝟙 X)` と同型なら構造射は同型。 -/
theorem over_isIso_hom_of_iso {X : D} {Z : Over X} (e : Z ≅ Over.mk (𝟙 X)) :
    IsIso Z.hom := by
  have h2 : Z.hom = e.hom.left := (Over.w e.hom).symm.trans (Category.comp_id _)
  rw [h2]
  exact (Over.forget X).map_isIso e.hom

def overIsoMkId.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — Over X で構造射が同型 ⟺ 終対象と同型",
    sectionId := "frdi-thm-3-4" }

/-! ## ★8. `plBkPsi` と `basePsiSlice` は圏同値

★★`Cor410Birat.lean` の `coaPrePsi_isEquivalence` と**同じ形**である ——
`Ψ` の忠実・充満・本質的全射を `WideSubcategory` に包み直すだけ。 -/

include hPB in
/-- ★忠実。 -/
theorem plBkPsi_faithful [Ψ.Faithful] : (plBkPsi P P₂ Ψ hPB).Faithful where
  map_injective {_ _ f g} h :=
    WideSubcategory.hom_ext _ (Ψ.map_injective
      (congrArg (fun t : (plBkPsi P P₂ Ψ hPB).obj _ ⟶ _ => t.hom) h))

include hPB in
/-- ★★充満 —— `Ψ` の充満性 ＋ pull-back の**反射**。 -/
theorem plBkPsi_full [Ψ.Full]
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f) :
    (plBkPsi P P₂ Ψ hPB).Full where
  map_surjective {X Y} g := by
    obtain ⟨f₀, hf₀⟩ := Ψ.map_surjective g.hom
    refine ⟨⟨f₀, hPB' f₀ ?_⟩, WideSubcategory.hom_ext _ hf₀⟩
    rw [hf₀]
    exact g.property

include hPB in
/-- ★本質的全射。 -/
theorem plBkPsi_essSurj [Ψ.EssSurj] : (plBkPsi P P₂ Ψ hPB).EssSurj where
  mem_essImage Y := by
    obtain ⟨A, ⟨ε⟩⟩ := Functor.EssSurj.mem_essImage (F := Ψ) Y.obj
    refine ⟨(⟨A⟩ : PlBk P), ⟨?_⟩⟩
    exact { hom := ⟨ε.hom, isPullBack_of_isIso P₂ ε.hom⟩,
            inv := ⟨ε.inv, isPullBack_of_isIso P₂ ε.inv⟩,
            hom_inv_id := WideSubcategory.hom_ext _ ε.hom_inv_id,
            inv_hom_id := WideSubcategory.hom_ext _ ε.inv_hom_id }

include hPB in
/-- ★★★★`plBkPsi` は圏同値。 -/
theorem plBkPsi_isEquivalence [Ψ.IsEquivalence]
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f) :
    (plBkPsi P P₂ Ψ hPB).IsEquivalence where
  faithful := plBkPsi_faithful P P₂ Ψ hPB
  full := plBkPsi_full P P₂ Ψ hPB hPB'
  essSurj := plBkPsi_essSurj P P₂ Ψ hPB

include F hPB in
/-- ★★★★★**共役関手は圏同値** —— 圏同値 3 つの合成。 -/
theorem basePsiSlice_isEquivalence [Ψ.IsEquivalence] (F₂ : FrobenioidCore P₂)
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f) (A : C) :
    (basePsiSlice P F P₂ Ψ hPB A).IsEquivalence := by
  haveI := F.plBkEquiv A
  haveI := F₂.plBkEquiv (Ψ.obj A)
  haveI := plBkPsi_isEquivalence P P₂ Ψ hPB hPB'
  haveI : (plBkSlicePsi P P₂ Ψ hPB A).IsEquivalence :=
    inferInstanceAs (Functor.IsEquivalence (Over.post (X := (⟨A⟩ : PlBk P))
      (plBkPsi P P₂ Ψ hPB)))
  exact inferInstanceAs (Functor.IsEquivalence
    ((plBkOverFunctor P A).inv ⋙ plBkSlicePsi P P₂ Ψ hPB A
      ⋙ plBkOverFunctor P₂ (Ψ.obj A)))

def basePsiSlice_isEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — 共役関手は圏同値",
    sectionId := "frdi-thm-3-4" }

/-! ## ★9. `plBkOverFunctor` は `Over.map` と可換(pull-back の場合)

★★★`Ψ_Base` の組み上げでは
`conj(Base ψ) ≅ Over.map (Base (Ψψ))`(`ψ` は pull-back)を示す必要がある。
★その芯が本補題である —— `plBkOverFunctor` は構造射に `Base` を当てるだけなので、
`Base` の合成則がそのまま自然性になる。 -/

include hPB in
/-- ★★★★**`plBkOverFunctor` は pull-back による `Over.map` と可換**。

★どちらも `Z ↦ Over.mk (Base Z.hom.hom ≫ Base ψ)` になる(`Base_comp`)。 -/
theorem plBkOverFunctor_over_map_obj {A A' : C} (ψ : A ⟶ A')
    (Z : Over (⟨A⟩ : PlBk P)) :
    (Over.map (P.Base ψ)).obj ((plBkOverFunctor P A).obj Z)
      = Over.mk (P.Base (Z.hom.hom ≫ ψ)) :=
  congrArg Over.mk (P.Base_comp Z.hom.hom ψ).symm

/-! ### ★関手としての可換性は `eqToHom` の処理が残る(2026-08-19 実測)

★対象の側(`plBkOverFunctor_over_map_obj`)は `Base_comp` そのもので出るが、
`Functor.ext` の射の側は `Over.Hom.left (eqToHom _)` を潰す補題が要る。
★`Comma.eqToHom_left` では `simp` が閉じなかった。
★★**object 版があれば `Ψ_Base` の組み上げには足りる見込み**なので、
関手版は必要になった時点で書く。 -/


end PlBkPsi

end ABC3.Found.FrdI
