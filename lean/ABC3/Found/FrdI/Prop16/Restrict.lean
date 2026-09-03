/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop19
import ABC3.Found.FrdI.PlBkShuffle

/-!
# Prop16 —— `Φ` の `𝒟'` への制限・(i) の `𝔽_{Φ'}`

☆もとの 1 枚を**入れ子の切れ目**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {D' : Type u3} [Category.{v3} D']



/-! ## ★`Φ` の `𝒟'` への制限

原文 (FrdI p.27):
> monoid on a connected, totally epimorphic category D; C →FΦ a Frobenioid.

★**`FSM ↦ FSM` の仮定はここでだけ使う** —— `Definition 1.1, (ii), (b)` を
`Φ'` について確かめるのに要る。`(ii), (a)` は `Φ` のものがそのまま降りる。 -/

/-- `Φ` を関手 `G : 𝒟' ⥤ 𝒟` に沿って制限した monoid。 -/
def MonoidOn.restrict (Φ : MonoidOn.{v, u, w} D) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α)) :
    MonoidOn.{v3, u3, w} D' where
  functor := G.op ⋙ Φ.functor
  charInj α := Φ.charInj (G.map α)
  fsmIso α hα := Φ.fsmIso (G.map α) (hG α hα)

@[simp] theorem MonoidOn.restrict_val (Φ : MonoidOn.{v, u, w} D) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α)) (A : D') :
    (Φ.restrict G hG).val A = Φ.val (G.obj A) := rfl

theorem MonoidOn.restrict_map (Φ : MonoidOn.{v, u, w} D) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α))
    {A B : D'} (α : B ⟶ A) (x : Φ.val (G.obj A)) :
    (Φ.restrict G hG).map α x = Φ.map (G.map α) x := rfl

/-! ## ★(i) —— `𝔽_{Φ'} ≅ 𝔽_Φ ×_𝒟 𝒟'`

原文 (FrdI p.27):
> (i) There is a natural equivalence of categories

★原文は「follows formally from the definitions」と言う。**測る。** -/

/-- **(i)** の関手 `𝔽_{Φ'} ⥤ 𝔽_Φ ×_𝒟 𝒟'`。

★対象 `W` を `(⟨G W⟩, W, 𝟙)` に送る。`𝒟` 成分の同型は **`𝟙` に取れる**。 -/
def cfpElemFunctor (Φ : MonoidOn.{v, u, w} D) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α)) :
    ElemFrobCat (Φ.restrict G hG) ⥤ CFP (ElemFrobCat.proj (Φ := Φ)) G where
  obj W := ⟨⟨⟨G.obj W.base⟩, W.base, 𝟙 _⟩, inferInstanceAs (IsIso (𝟙 _))⟩
  map {W₁ W₂} f := InducedCategory.homMk ⟨⟨G.map f.base, f.div, f.deg⟩, f.base, by
    show G.map f.base ≫ 𝟙 (G.obj W₂.base) = 𝟙 (G.obj W₁.base) ≫ G.map f.base
    rw [Category.comp_id, Category.id_comp]⟩
  map_id W := by
    refine InducedCategory.hom_ext (CommaMorphism.ext ?_ rfl)
    refine ElemFrobCat.Hom.ext ?_ rfl rfl
    show G.map (𝟙 W.base) = 𝟙 _
    rw [G.map_id]
  map_comp {W₁ W₂ W₃} f g := by
    refine InducedCategory.hom_ext (CommaMorphism.ext ?_ rfl)
    refine ElemFrobCat.Hom.ext ?_ rfl rfl
    show G.map (f.base ≫ g.base) = G.map f.base ≫ G.map g.base
    rw [G.map_comp]

/-- **(i)** —— ★**`𝔽_{Φ'} ⥤ 𝔽_Φ ×_𝒟 𝒟'` は圏同値**。

★原文の「follows formally from the definitions」は**正しい**。測ると:
* **忠実性**: 像の 2 成分から `base` / `div` / `deg` がそのまま読める
* **充満性**: 四角形が `h₁.base = G h₂` を強制するので `h₂` が原像を与える
* **本質的全射性**: `(Y, W, α)` に対し `W` を取り、同型を **`⟨α⁻¹, 0, 1⟩`** で作る
★**3 成分とも「四角形が第2成分を第1成分から決める」ことに帰着する。** -/
theorem cfpElemFunctor_isEquivalence (Φ : MonoidOn.{v, u, w} D) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α)) :
    (cfpElemFunctor Φ G hG).IsEquivalence := by
  haveI hfaith : (cfpElemFunctor Φ G hG).Faithful := by
    constructor
    intro W₁ W₂ f g hfg
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · exact congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) hfg
    · exact congrArg
        (fun t => ElemFrobCat.Hom.div (CommaMorphism.left (InducedCategory.Hom.hom t))) hfg
    · exact congrArg
        (fun t => ElemFrobCat.Hom.deg (CommaMorphism.left (InducedCategory.Hom.hom t))) hfg
  haveI hfull : (cfpElemFunctor Φ G hG).Full := by
    constructor
    intro W₁ W₂ h
    have hw := (InducedCategory.Hom.hom h).w
    have hb := (Category.comp_id (InducedCategory.Hom.hom h).left.base).symm.trans
      (hw.trans (Category.id_comp _))
    refine ⟨⟨(InducedCategory.Hom.hom h).right, (InducedCategory.Hom.hom h).left.div,
      (InducedCategory.Hom.hom h).left.deg⟩, ?_⟩
    refine InducedCategory.hom_ext (CommaMorphism.ext ?_ rfl)
    exact ElemFrobCat.Hom.ext hb.symm rfl rfl
  haveI hess : (cfpElemFunctor Φ G hG).EssSurj := by
    constructor
    intro X
    haveI hXi : IsIso X.obj.hom := X.property
    obtain ⟨v, hv1, hv2⟩ := hXi.out
    refine ⟨⟨X.obj.right⟩, ⟨?_⟩⟩
    refine ⟨InducedCategory.homMk ⟨⟨v, 0, 1⟩, 𝟙 _, ?_⟩,
      InducedCategory.homMk ⟨⟨X.obj.hom, 0, 1⟩, 𝟙 _, ?_⟩, ?_, ?_⟩
    · show v ≫ X.obj.hom = 𝟙 _ ≫ G.map (𝟙 X.obj.right)
      rw [hv2, G.map_id, Category.comp_id]
    · show X.obj.hom ≫ 𝟙 _ = X.obj.hom ≫ G.map (𝟙 X.obj.right)
      rw [G.map_id]
    · refine InducedCategory.hom_ext (CommaMorphism.ext ?_ (Category.comp_id _))
      refine ElemFrobCat.Hom.ext ?_ ?_ ?_
      · show v ≫ X.obj.hom = 𝟙 _
        exact hv2
      · simp [ElemFrobCat.comp_div]
      · simp [ElemFrobCat.comp_deg]
    · refine InducedCategory.hom_ext (CommaMorphism.ext ?_ (Category.comp_id _))
      refine ElemFrobCat.Hom.ext ?_ ?_ ?_
      · show X.obj.hom ≫ v = 𝟙 _
        exact hv1
      · simp [ElemFrobCat.comp_div]
      · simp [ElemFrobCat.comp_deg]
  exact ⟨hfaith, hfull, hess⟩

/-! ### ★★教訓 —— **「補題が無い」ではなく「名前を間違えていた」**

(i) の圏同値は一度「`ElemFrobCat` の合成の `div` 成分が `simp` で閉じない、原因未特定」
として切った。★**原因は単純で、`simp [ElemFrobCat.Hom.comp]` と
存在しない名前を指定していた**だけだった。正しくは
**`ElemFrobCat.comp_div` / `comp_deg`(どちらも既に `@[simp]`)**である。

★**「無い」と言う前に S1–S4(ファイル名の列挙を含む)** —— この規律は
**補題の名前にも当てはまる**。★**推測した名前が通らないことを「原因未特定」と書いてはいけない。**
まず `grep` する。 -/

/-! ## ★`𝒞' = 𝒞 ×_𝒟 𝒟'` -/

variable {C : Type u2} [Category.{v2} C] {Φ : MonoidOn.{v, u, w} D}

/-- **`𝒞' = 𝒞 ×_𝒟 𝒟'`**。 -/
abbrev CfpCat (P : PreFrobenioid C Φ) (G : D' ⥤ D) : Type _ := CFP P.proj G

/-- ★`𝒞'` の射の `𝒞` 成分。 -/
abbrev CfpCat.fst {P : PreFrobenioid C Φ} {G : D' ⥤ D} {X Y : CfpCat P G} (f : X ⟶ Y) :
    X.obj.left ⟶ Y.obj.left := f.hom.left

/-- ★`𝒞'` の射の `𝒟'` 成分。 -/
abbrev CfpCat.snd {P : PreFrobenioid C Φ} {G : D' ⥤ D} {X Y : CfpCat P G} (f : X ⟶ Y) :
    X.obj.right ⟶ Y.obj.right := f.hom.right

/-- ★**`𝒞'` は totally epimorphic**。

原文 (FrdI p.28):
> well; similarly, [in light of the various properties of the natural projection functor

★原文は `𝒟'` の側しか挙げないが、**`𝒞` の側も要る**(射は対なので、
両成分がそれぞれ epi でなければならない)。`𝒞` は Frobenioid なので既に totally epimorphic。 -/
theorem cfp_totEpi (P : PreFrobenioid C Φ) (G : D' ⥤ D) (hD' : IsTotallyEpimorphic D') :
    IsTotallyEpimorphic (CfpCat P G) := by
  intro X Y f
  constructor
  intro Z g h hgh
  haveI h1 : Epi (CfpCat.fst f) := P.totEpiC _ _ _
  haveI h2 : Epi (CfpCat.snd f) := hD' _ _ _
  have e1 : CfpCat.fst f ≫ CfpCat.fst g = CfpCat.fst f ≫ CfpCat.fst h :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hgh
  have e2 : CfpCat.snd f ≫ CfpCat.snd g = CfpCat.snd f ≫ CfpCat.snd h :=
    congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) hgh
  exact InducedCategory.hom_ext
    (CommaMorphism.ext ((cancel_epi (CfpCat.fst f)).mp e1) ((cancel_epi (CfpCat.snd f)).mp e2))

/-! ## ★`𝒞' → 𝔽_{Φ'}`

★**書き方の注意(手順4として記録)**: `P.Base f` の型は
`(P.toElem.obj A).base ⟶ (P.toElem.obj B).base`、`Comma` の `hom` の型は
`P.proj.obj A ⟶ G.obj A'` である。この2つは**定義的には等しいが字面が違う**ので、
`rw` が通らない。★**`P.proj.map` の綴りに統一する**のが正しい対処である
(`Istr` の `rw` 問題と同じ形の、しかし別の原因の問題)。 -/

/-- ★`𝒞'` の射の `𝒟'` 成分と `𝒞` 成分を結ぶ四角形(`Comma` の `w`)。 -/
theorem cfp_square {P : PreFrobenioid C Φ} {G : D' ⥤ D} {X Y : CfpCat P G} (f : X ⟶ Y) :
    P.proj.map (CfpCat.fst f) ≫ Y.obj.hom = X.obj.hom ≫ G.map (CfpCat.snd f) :=
  f.hom.w

/-- ★上の四角形を「`α⁻¹` を通す」形に直したもの。 -/
theorem cfp_square_inv {P : PreFrobenioid C Φ} {G : D' ⥤ D} {X Y : CfpCat P G} (f : X ⟶ Y)
    [IsIso X.obj.hom] [IsIso Y.obj.hom] :
    inv X.obj.hom ≫ P.proj.map (CfpCat.fst f) = G.map (CfpCat.snd f) ≫ inv Y.obj.hom := by
  rw [IsIso.eq_comp_inv, Category.assoc, cfp_square f, ← Category.assoc, IsIso.inv_hom_id,
    Category.id_comp]

/-- ★同じものを `P.Base` の綴りで言い換えたもの(型は定義的に等しい)。 -/
theorem cfp_square_inv' {P : PreFrobenioid C Φ} {G : D' ⥤ D} {X Y : CfpCat P G} (f : X ⟶ Y)
    [IsIso X.obj.hom] [IsIso Y.obj.hom] :
    inv X.obj.hom ≫ P.Base (CfpCat.fst f) = G.map (CfpCat.snd f) ≫ inv Y.obj.hom :=
  cfp_square_inv f

/-- ★**`𝒞' → 𝔽_{Φ'}`** —— `𝒞'` の pre-Frobenioid 構造の本体。

対象 `(A, A', α : Base A ≅ G A')` を `A' ∈ 𝒟'` に、
射 `(γ, γ')` を `⟨γ', Φ(α⁻¹)(Div γ), deg_Fr γ⟩` に送る。

★`Div` の付け替え `Φ(α⁻¹)` が要るのは、`Div γ ∈ Φ(Base A)` であって
`Φ'(A') = Φ(G A')` ではないからである。**この付け替えが関手性の中身**であり、
それを支えるのが `cfp_square_inv` である。 -/
noncomputable def cfpToElem (P : PreFrobenioid C Φ) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α)) :
    CfpCat P G ⥤ ElemFrobCat (Φ.restrict G hG) where
  obj X := ⟨X.obj.right⟩
  map {X Y} f :=
    haveI : IsIso X.obj.hom := X.property
    { base := CfpCat.snd f
      div := Φ.map (inv X.obj.hom) (P.Div (CfpCat.fst f))
      deg := P.degFr (CfpCat.fst f) }
  map_id X := by
    haveI hX : IsIso X.obj.hom := X.property
    have hid := P.toElem.map_id X.obj.left
    have hdiv0 : P.Div (𝟙 X.obj.left) = 0 := congrArg ElemFrobCat.Hom.div hid
    have hdeg0 : P.degFr (𝟙 X.obj.left) = 1 := congrArg ElemFrobCat.Hom.deg hid
    apply ElemFrobCat.Hom.ext
    · rfl
    · show Φ.map (inv X.obj.hom) (P.Div (𝟙 X.obj.left)) = 0
      rw [hdiv0]
      exact map_zero _
    · exact hdeg0
  map_comp {X Y Z} f g := by
    haveI hX : IsIso X.obj.hom := X.property
    haveI hY : IsIso Y.obj.hom := Y.property
    have hcomp := P.toElem.map_comp (CfpCat.fst f) (CfpCat.fst g)
    have hdiv : P.Div (CfpCat.fst f ≫ CfpCat.fst g)
        = Φ.map (P.Base (CfpCat.fst f)) (P.Div (CfpCat.fst g))
          + (P.degFr (CfpCat.fst g) : ℕ) • P.Div (CfpCat.fst f) :=
      congrArg ElemFrobCat.Hom.div hcomp
    have hdeg : P.degFr (CfpCat.fst f ≫ CfpCat.fst g)
        = P.degFr (CfpCat.fst g) * P.degFr (CfpCat.fst f) :=
      congrArg ElemFrobCat.Hom.deg hcomp
    apply ElemFrobCat.Hom.ext
    · rfl
    · show Φ.map (inv X.obj.hom) (P.Div (CfpCat.fst f ≫ CfpCat.fst g))
        = Φ.map (G.map (CfpCat.snd f)) (Φ.map (inv Y.obj.hom) (P.Div (CfpCat.fst g)))
          + (P.degFr (CfpCat.fst g) : ℕ) • Φ.map (inv X.obj.hom) (P.Div (CfpCat.fst f))
      -- ★`rw` は使わず**項で繋ぐ**。`P.proj.obj A` と `(P.toElem.obj A).base` は
      -- 定義的に等しいだけで字面が違うので、`rw` の照合が通らない。
      refine Eq.trans (congrArg (Φ.map (inv X.obj.hom)) hdiv) ?_
      refine Eq.trans ((Φ.map (inv X.obj.hom)).map_add _ _) ?_
      refine congrArg₂ (· + ·) ?_ ?_
      · exact (((Φ.map_comp (P.Base (CfpCat.fst f)) (@inv _ _ _ _ X.obj.hom hX)
            (P.Div (CfpCat.fst g))).symm.trans
          (congrArg (fun t => Φ.map t (P.Div (CfpCat.fst g))) (cfp_square_inv' f))).trans
          (Φ.map_comp (@inv _ _ _ _ Y.obj.hom hY) (G.map (CfpCat.snd f))
            (P.Div (CfpCat.fst g))))
      · exact (Φ.map (inv X.obj.hom)).map_nsmul _ _
    · exact hdeg

/-- ★★**`𝒞' = 𝒞 ×_𝒟 𝒟'` は pre-Frobenioid**。

原文 (FrdI p.28):
> Frobenioid. Now assertion (vi) follows immediately from the definitions; one checks
-/
noncomputable def cfpPreFrobenioid (P : PreFrobenioid C Φ) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α))
    (hD' : IsTotallyEpimorphic D') (hcC : IsConnected (CfpCat P G))
    (hcD' : IsConnected D') : PreFrobenioid (CfpCat P G) (Φ.restrict G hG) where
  toElem := cfpToElem P G hG
  divisorial A := P.divisorial (G.obj A)
  totEpiC := cfp_totEpi P G hD'
  totEpiD := hD'
  connectedC := hcC
  connectedD := hcD'

/-! ## ★辞書 —— `𝒞'` の `Base` / `Div` / `deg_Fr`

★`Istr` のときの `istr_compat_*` に当たるもの。ただし **`Div` だけは `rfl` ではない**
(`Φ(α⁻¹)` の付け替えが挟まる)。そこが `Istr`(充満部分圏)との違いである。 -/

end ABC3.Found.FrdI
