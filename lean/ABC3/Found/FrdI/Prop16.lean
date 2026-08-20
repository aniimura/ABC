import ABC3.Found.FrdI.Prop19
import ABC3.Found.FrdI.PlBkShuffle

/-!
# [FrdI] Proposition 1.6 —— Categorical Fiber Products

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.27–p.28(**400 dpi 目視確認 2026-08-15、実装のための再読**)。
§0 の `CFP` は物理 p.17(同じく目視)。

原文 (FrdI p.27):
> Proposition 1.6.

原文 (FrdI p.27):
> (Categorical Fiber Products) Let Φ be a divisorial

★原文の「`D′` を取る」行は **`′`(prime) を含むため逐語照合できない**
(pdftotext の layout 抽出で `′` が拾えない)。**書き換えず、照合できない事実として記す。**

## ★規模の測定(目視)

| 条 | 内容 | 主張の数 |
|---|---|---|
| (i) | `𝔽_Φ' ≅ 𝔽_Φ ×_𝒟 𝒟'` | 1 |
| (ii) | `𝒞' = 𝒞 ×_𝒟 𝒟'` が Frobenioid | 1(＋`Definition 1.3` の 21 条) |
| (iii) | 5 クラスが射影で判定できる | 5 |
| (iv) | base-iso について 3 クラス ＋ `𝒪^▷` の全単射 | 4 |
| (v) | 11 個の対象タイプ | 11 |
| (vi) | 3 タイプ(★**片向きのみ**) | 3 |

★**証明は 1 段落**(11 行)。21 条は「routine verification」の一言で片づけられている。

## ★★仕分け —— `Istr` とは**保証の出どころが違う**

`Proposition 1.9, (v)` では `Definition 1.3, (vii), (b)`(isotropic からの射の終域は
isotropic)が**閉性**を与え、前向きの構成が自動で `𝒞^istr` に留まった。
**CFP に閉性は無い。** 代わりに効くのは:

★**新しく作る対象・射は必ず base-isomorphism に沿って現れ、
`𝒞'` の base-isomorphism とは「`𝒟'` 成分が同型」という意味**なので、
`𝒟'` 成分を `𝟙` に取れる(新規対象)か、入力から与えられる(`plBkEquiv` など)。

★したがって **`G : 𝒟' ⥤ 𝒟` の充満性は要らない**。
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

section Dict

variable (P : PreFrobenioid C Φ) (G : D' ⥤ D)
  (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α))
  (hD' : IsTotallyEpimorphic D')
  (hcC : IsConnected (CfpCat P G)) (hcD' : IsConnected D')

theorem cfp_compat_Base {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD' hcC hcD').Base f = CfpCat.snd f := rfl

theorem cfp_compat_degFr {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD' hcC hcD').degFr f = P.degFr (CfpCat.fst f) := rfl

theorem cfp_compat_Div {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD' hcC hcD').Div f
      = Φ.map (@inv _ _ _ _ X.obj.hom X.property) (P.Div (CfpCat.fst f)) := rfl

/-- **(iii)** —— **Frobenius 次数**は射影で決まる。 -/
theorem cfp_degFr_eq {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD' hcC hcD').degFr f = P.degFr (CfpCat.fst f) := rfl

/-- **(iv)** —— **base-isomorphism** は `𝒟'` 成分が同型であること。

★★**これが CFP の移送を支える一点である** —— `𝒟'` 成分の情報は
base-isomorphism の定義そのものに入っている。 -/
theorem cfp_baseIso_iff {X Y : CfpCat P G} (f : X ⟶ Y) :
    IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsIso (CfpCat.snd f) := Iff.rfl

/-- **(iii)** —— **linear** は射影で決まる。 -/
theorem cfp_linear_iff {X Y : CfpCat P G} (f : X ⟶ Y) :
    IsLinear (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsLinear P (CfpCat.fst f) := Iff.rfl

/-- **(iii)** —— **isometry** は射影で決まる。

★中身は「`Φ(α⁻¹)` が単射」の一点(`Definition 1.1, (ii), (a)`)。 -/
theorem cfp_isometric_iff {X Y : CfpCat P G} (f : X ⟶ Y) :
    IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsIsometric P (CfpCat.fst f) := by
  constructor
  · intro h
    exact Φ.map_injective (@inv _ _ _ _ X.obj.hom X.property) (h.trans (map_zero _).symm)
  · intro h
    show Φ.map (@inv _ _ _ _ X.obj.hom X.property) (P.Div (CfpCat.fst f)) = 0
    rw [show P.Div (CfpCat.fst f) = 0 from h]
    exact map_zero _

/-! ### ★**手順5**(2つの綴り問題への対処、確定版)

★`rw` を項で置き換えるだけでは足りない。**部分型の成分や `congrArg` の射影**は
綴りが食い違ったままなので、
**「綴りの決まった変数を `obtain ⟨h', rfl⟩ : ∃ h' : 正しい型, h' = h` で先に導入する」**
のが正しい対処である(`Prop19` の `Over.Hom.left` に使った定型と同じもの)。 -/

/-- ★`𝒞'` の射の `𝒞` 成分の底射は、`𝒟'` 成分と両端の同型で**完全に決まる**。

★★**これが CFP の移送を支える一点**である —— `𝒞` 成分の底射は自由ではない。 -/
theorem cfp_base_fst {X Y : CfpCat P G} (f : X ⟶ Y) [IsIso Y.obj.hom] :
    P.proj.map (CfpCat.fst f) = X.obj.hom ≫ G.map (CfpCat.snd f) ≫ inv Y.obj.hom := by
  rw [← Category.assoc, ← cfp_square f, Category.assoc, IsIso.hom_inv_id, Category.comp_id]

/-- **(iii)** —— `𝒞` の pull-back は `𝒞'` の pull-back(★**構成の向き**)。

★中身は「`𝒟'` 成分 `h` を与えると、`𝒞` 側で使うべき底射
`u = α_Z ≫ G(h) ≫ α_X⁻¹` が**一意に決まる**」の一点。 -/
theorem cfp_isPullBack_of {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsPullBack P (CfpCat.fst φ)) : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') φ := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  intro Z
  haveI hZ : IsIso Z.obj.hom := Z.property
  constructor
  · intro f₁ f₂ hf
    have hp := Subtype.ext_iff.mp hf
    have hs : CfpCat.snd f₁ = CfpCat.snd f₂ := congrArg Prod.snd hp
    have hcomp : (f₁ ≫ φ : Z ⟶ Y) = f₂ ≫ φ := congrArg Prod.fst hp
    have hc : CfpCat.fst f₁ ≫ CfpCat.fst φ = CfpCat.fst f₂ ≫ CfpCat.fst φ :=
      congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hcomp
    have hb : P.Base (CfpCat.fst f₁) = P.Base (CfpCat.fst f₂) := by
      show P.proj.map (CfpCat.fst f₁) = P.proj.map (CfpCat.fst f₂)
      rw [cfp_base_fst P G f₁, cfp_base_fst P G f₂, hs]
    exact InducedCategory.hom_ext
      (CommaMorphism.ext ((h Z.obj.left).1 (Subtype.ext (Prod.ext hc hb))) hs)
  · rintro ⟨⟨g, hh0⟩, hcond⟩
    -- ★手順5: 綴りの決まった変数を先に導入する
    obtain ⟨hh, rfl⟩ : ∃ hh : Z.obj.right ⟶ X.obj.right, hh = hh0 := ⟨hh0, rfl⟩
    have hcond' : CfpCat.snd g = hh ≫ CfpCat.snd φ := hcond
    obtain ⟨u, hu⟩ : ∃ u : P.proj.obj Z.obj.left ⟶ P.proj.obj X.obj.left,
        u = Z.obj.hom ≫ G.map hh ≫ inv X.obj.hom := ⟨_, rfl⟩
    have hbase : P.proj.map (CfpCat.fst g) = u ≫ P.proj.map (CfpCat.fst φ) := by
      rw [hu, cfp_base_fst P G g, cfp_base_fst P G φ, hcond', G.map_comp]
      simp only [Category.assoc]
      rw [← Category.assoc (inv X.obj.hom) X.obj.hom, IsIso.inv_hom_id, Category.id_comp]
    obtain ⟨f₁, hf₁⟩ := (h Z.obj.left).2 ⟨(CfpCat.fst g, u), hbase⟩
    have hp := Subtype.ext_iff.mp hf₁
    have h1 : (f₁ ≫ CfpCat.fst φ) = CfpCat.fst g := congrArg Prod.fst hp
    have h2 : P.proj.map f₁ = u := congrArg Prod.snd hp
    have hw : P.proj.map f₁ ≫ X.obj.hom = Z.obj.hom ≫ G.map hh := by
      rw [h2, hu]
      simp only [Category.assoc]
      rw [IsIso.inv_hom_id, Category.comp_id]
    refine ⟨InducedCategory.homMk ⟨f₁, hh, hw⟩, Subtype.ext (Prod.ext ?_ rfl)⟩
    exact InducedCategory.hom_ext (CommaMorphism.ext h1 hcond'.symm)

/-- ★`𝒞'` の射は、両成分が同型なら同型。 -/
theorem cfp_isIso_of {X Y : CfpCat P G} (f : X ⟶ Y) (h1 : IsIso (CfpCat.fst f))
    (h2 : IsIso (CfpCat.snd f)) : IsIso f := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  haveI := h1
  haveI := h2
  have hw : P.proj.map (inv (CfpCat.fst f)) ≫ X.obj.hom
      = Y.obj.hom ≫ G.map (inv (CfpCat.snd f)) := by
    rw [P.proj.map_inv, G.map_inv, IsIso.inv_comp_eq, ← Category.assoc, IsIso.eq_comp_inv]
    exact (cfp_square f).symm
  refine ⟨InducedCategory.homMk ⟨inv (CfpCat.fst f), inv (CfpCat.snd f), hw⟩, ?_, ?_⟩
  · exact InducedCategory.hom_ext (CommaMorphism.ext (IsIso.hom_inv_id _) (IsIso.hom_inv_id _))
  · exact InducedCategory.hom_ext (CommaMorphism.ext (IsIso.inv_hom_id _) (IsIso.inv_hom_id _))

/-- ★`𝒞'` の base-isomorphism は `𝒞` の base-isomorphism。

★`cfp_base_fst` により `𝒞` 成分の底射は同型3つの合成になる。 -/
theorem cfp_baseIso_fst {X Y : CfpCat P G} (f : X ⟶ Y)
    (h : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f) :
    IsBaseIsomorphism P (CfpCat.fst f) := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  haveI : IsIso (CfpCat.snd f) := h
  show IsIso (P.proj.map (CfpCat.fst f))
  rw [cfp_base_fst P G f]
  infer_instance

/-- **(iii)** —— **co-angular** は `𝒞` から `𝒞'` へ降りる(★構成の向き)。 -/
theorem cfp_coAngular_of {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsCoAngular P (CfpCat.fst φ)) :
    IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ := by
  intro Z W γ β α hfac hαl hβi hβs hdisj
  have hfac' : CfpCat.fst φ = CfpCat.fst γ ≫ CfpCat.fst β ≫ CfpCat.fst α :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hfac
  have hβsC : IsPreStep P (CfpCat.fst β) :=
    ⟨hβs.1, cfp_baseIso_fst P G hG hD' hcC hcD' β hβs.2⟩
  have hdisjC : IsBaseIsomorphism P (CfpCat.fst α) ∨ IsBaseIsomorphism P (CfpCat.fst γ) :=
    hdisj.imp (cfp_baseIso_fst P G hG hD' hcC hcD' α) (cfp_baseIso_fst P G hG hD' hcC hcD' γ)
  have := h _ _ (CfpCat.fst γ) (CfpCat.fst β) (CfpCat.fst α) hfac' hαl
    ((cfp_isometric_iff P G hG hD' hcC hcD' β).mp hβi) hβsC hdisjC
  exact cfp_isIso_of P G β this hβs.2

/-- ★`𝒞` の対象を、底が `G` の像と同型であるときに `𝒞'` へ持ち上げる。 -/
def cfpMk (A : C) (W : D') (e : P.proj.obj A ⟶ G.obj W) [IsIso e] : CfpCat P G :=
  ⟨⟨A, W, e⟩, inferInstanceAs (IsIso e)⟩

/-- ★`𝒞'` の射を両成分と四角形から作る。 -/
def cfpHom {X Y : CfpCat P G} (u : X.obj.left ⟶ Y.obj.left) (v : X.obj.right ⟶ Y.obj.right)
    (w : P.proj.map u ≫ Y.obj.hom = X.obj.hom ≫ G.map v) : X ⟶ Y :=
  InducedCategory.homMk ⟨u, v, w⟩

/-- ★`𝒞'` の同型の `𝒞` 成分は同型。 -/
theorem cfp_isIso_fst {X Y : CfpCat P G} (f : X ⟶ Y) (h : IsIso f) : IsIso (CfpCat.fst f) := by
  obtain ⟨g, h1, h2⟩ := h.out
  refine ⟨CfpCat.fst g, ?_, ?_⟩
  · exact congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) h1
  · exact congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) h2

/-! ### ★**手順6**(2つの綴り問題への対処、3つ目)

★**射影が簡約されない `def` を作らず、構造リテラルを直接書く**。
`cfpMk` のような補助定義を挟むと `(cfpMk …).obj.right` が `Y.obj.right` に
**簡約されず**、`rw` も インスタンス合成も落ちる。
★手順5(部分型の成分)では届かない、3つ目の出方である。 -/

/-- **(iii)** —— **co-angular** は `𝒞'` から `𝒞` へ上がる(★**難しい向き**)。

★★**要点**: `co-angular` の定義に入っている**選言**
「`α` か `γ` が base-isomorphism」が、**中間対象を `𝒞'` へ持ち上げる橋**になる。
`β` は pre-step なので必ず base-isomorphism であり、
選言の側からもう一方の端まで **base-isomorphism の鎖**が繋がるので、
`Base Z₀` も `Base W₀` も `G` の像と同型になる。

★**`G` の本質的全射性は要らない** —— 使うのは「与えられた分解の中の base-isomorphism」だけ。
★★**これが CFP 版の仕分け基準**である: **定義の中に base-isomorphism の鎖があるか**。 -/
theorem cfp_coAngular_to {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ) : IsCoAngular P (CfpCat.fst φ) := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  intro Z₀ W₀ γ₀ β₀ α₀ hfac hα₀l hβ₀i hβ₀s hdisj
  haveI hβb : IsIso (P.proj.map β₀) := hβ₀s.2
  have hbf : P.proj.map (CfpCat.fst φ) ≫ Y.obj.hom = X.obj.hom ≫ G.map (CfpCat.snd φ) :=
    cfp_square φ
  have hfacp : P.proj.map (CfpCat.fst φ)
      = P.proj.map γ₀ ≫ P.proj.map β₀ ≫ P.proj.map α₀ := by
    rw [hfac, P.proj.map_comp, P.proj.map_comp]
  rcases hdisj with hα | hγ
  · -- ★`α₀` が base-isomorphism: `Y` の側から鎖を辿る
    haveI hαb : IsIso (P.proj.map α₀) := hα
    have hz : IsIso (P.proj.map β₀ ≫ P.proj.map α₀ ≫ Y.obj.hom) := inferInstance
    have hw : IsIso (P.proj.map α₀ ≫ Y.obj.hom) := inferInstance
    refine cfp_isIso_fst P G
      (X := ⟨⟨Z₀, Y.obj.right, P.proj.map β₀ ≫ P.proj.map α₀ ≫ Y.obj.hom⟩, hz⟩)
      (Y := ⟨⟨W₀, Y.obj.right, P.proj.map α₀ ≫ Y.obj.hom⟩, hw⟩)
      (InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩) ?_
    refine h ⟨⟨Z₀, Y.obj.right, P.proj.map β₀ ≫ P.proj.map α₀ ≫ Y.obj.hom⟩, hz⟩
      ⟨⟨W₀, Y.obj.right, P.proj.map α₀ ≫ Y.obj.hom⟩, hw⟩
      (InducedCategory.homMk ⟨γ₀, (CfpCat.snd φ : X.obj.right ⟶ Y.obj.right), ?_⟩)
      (InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩)
      (InducedCategory.homMk ⟨α₀, 𝟙 _, by simp⟩) ?_
      hα₀l ((cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hβ₀i) ⟨hβ₀s.1, by show IsIso (𝟙 _); infer_instance⟩
      (Or.inl (by show IsIso (𝟙 _); infer_instance))
    · show P.proj.map γ₀ ≫ P.proj.map β₀ ≫ P.proj.map α₀ ≫ Y.obj.hom
        = X.obj.hom ≫ G.map (CfpCat.snd φ)
      rw [← hbf, hfacp]
      simp only [Category.assoc]
    · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
      show CfpCat.snd φ = CfpCat.snd φ ≫ 𝟙 _ ≫ 𝟙 _
      simp
  · -- ★`γ₀` が base-isomorphism: `X` の側から鎖を辿る
    haveI hγb : IsIso (P.proj.map γ₀) := hγ
    have hz : IsIso (inv (P.proj.map γ₀) ≫ X.obj.hom) := inferInstance
    have hw : IsIso (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ X.obj.hom) := inferInstance
    refine cfp_isIso_fst P G
      (X := ⟨⟨Z₀, X.obj.right, inv (P.proj.map γ₀) ≫ X.obj.hom⟩, hz⟩)
      (Y := ⟨⟨W₀, X.obj.right,
        inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ X.obj.hom⟩, hw⟩)
      (InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩) ?_
    refine h ⟨⟨Z₀, X.obj.right, inv (P.proj.map γ₀) ≫ X.obj.hom⟩, hz⟩
      ⟨⟨W₀, X.obj.right,
        inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ X.obj.hom⟩, hw⟩
      (InducedCategory.homMk ⟨γ₀, 𝟙 _, by simp⟩)
      (InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩)
      (InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : X.obj.right ⟶ Y.obj.right), ?_⟩) ?_
      hα₀l ((cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hβ₀i) ⟨hβ₀s.1, by show IsIso (𝟙 _); infer_instance⟩
      (Or.inr (by show IsIso (𝟙 _); infer_instance))
    · show P.proj.map α₀ ≫ Y.obj.hom
        = (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ X.obj.hom) ≫ G.map (CfpCat.snd φ)
      rw [Category.assoc, Category.assoc, ← hbf, hfacp]
      simp only [Category.assoc]
      rw [← Category.assoc (inv (P.proj.map γ₀)), IsIso.inv_hom_id, Category.id_comp,
        ← Category.assoc (inv (P.proj.map β₀)), IsIso.inv_hom_id, Category.id_comp]
    · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
      show CfpCat.snd φ = 𝟙 _ ≫ 𝟙 _ ≫ CfpCat.snd φ
      simp

/-- **(iii)** —— co-angular は射影で決まる(両向き)。 -/
theorem cfp_coAngular_iff {X Y : CfpCat P G} (φ : X ⟶ Y) :
    IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ ↔ IsCoAngular P (CfpCat.fst φ) :=
  ⟨cfp_coAngular_to P G hG hD' hcC hcD' φ, cfp_coAngular_of P G hG hD' hcC hcD' φ⟩

/-- **(iii)** —— **LB-invertible** は射影で決まる。

★`LB-invertible = co-angular ∧ isometric` なので、上の2つの合成。 -/
theorem cfp_lbInvertible_iff {X Y : CfpCat P G} (φ : X ⟶ Y) :
    IsLBInvertible (cfpPreFrobenioid P G hG hD' hcC hcD') φ ↔ IsLBInvertible P (CfpCat.fst φ) :=
  and_congr (cfp_coAngular_iff P G hG hD' hcC hcD' φ) (cfp_isometric_iff P G hG hD' hcC hcD' φ)

/-! ### ★★**分類表** —— 「型は等しいが字面が違う」問題の症状・原因・対処

★**個数ではなく「症状が分岐したこと」が表を作る合図**だった。
同じ `def` 由来でも対処が違うと分かった時点で分ける。
★**切り替え点(新しいファイル/section/新しい種類の補題)でこの表を引く。**

| # | 症状(出るエラー) | 原因 | 対処 |
|---|---|---|---|
| 1 | `rw` が「目標に見えているのに見つからない」 | 同じ型に2つの綴り(`def` が展開されない) | **`rw` を使わず `Eq.trans`/`congrArg` で項として繋ぐ** |
| 2 | `inv` を書くと `motive is not type correct` | `InducedCategory.Hom` に包まれて型が簡約されない | **`inv` を書かず `IsIso.out` から逆射を取る**(例外: 主張自体に `inv` が現れるときは `IsIso.eq_inv_comp` にインスタンスを明示) |
| 3 | 部分型 `{p : _ × _ // _}` の成分の型が食い違う | 成分の型が別綴りで書かれている | **`obtain ⟨h', rfl⟩ : ∃ h' : 正しい型, h' = h` で綴りの決まった変数を先に導入** |
| 4 | `(f …).obj.right` が簡約されず `rw`/合成が落ちる | 補助 `def` の射影は簡約されない | **補助 `def` を作らず構造リテラルを直接書く** |
| 5 | `failed to synthesize instance` (文脈に `haveI` があるのに) | インスタンス型が `def`(例: `cfpProp`)で書かれ展開されない | **`have h : IsIso … := inferInstance` を先に置き、その項を明示的に渡す** |
| 6 | `refine h _ _ …` の穴が埋まらない | 中間対象がメタ変数のまま残る | **対象を明示的に書く**(`_` に頼らない) |
| 7 | `include` し忘れで `Unknown identifier` | 主張に現れない変数は自動包含されない | **最初に `include F in`** |
-/

/-! ## ★(iv) —— base-isomorphism についての3クラスと `𝒪^▷`

原文 (FrdI p.28):
> tively, pre-step; step) if and only if its projection to C is. Moreover, the projection

★★**基準の4例目**: (iv) が「**base-isomorphism について**」と限定して述べられているのは、
**その仮定が鎖そのものを与える**からである。仮定を外すと
「`fst f` が base-iso でも `snd f` が同型とは限らない」(`G` は同型を反映しない)ので、
**⟸ 向きが壊れる**。★**原文の限定は必要であり、我々の基準がその理由を説明する。** -/

/-- **(iv)** —— base-isomorphism については **pre-step** が射影で決まる。 -/
theorem cfp_preStep_iff {X Y : CfpCat P G} (f : X ⟶ Y)
    (hb : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f) :
    IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsPreStep P (CfpCat.fst f) :=
  ⟨fun h => ⟨h.1, cfp_baseIso_fst P G hG hD' hcC hcD' f hb⟩, fun h => ⟨h.1, hb⟩⟩

/-- **(iv)** —— base-isomorphism については **Frobenius 型**が射影で決まる。 -/
theorem cfp_frobType_iff {X Y : CfpCat P G} (f : X ⟶ Y)
    (hb : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f) :
    IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsFrobeniusType P (CfpCat.fst f) :=
  ⟨fun h => ⟨(cfp_lbInvertible_iff P G hG hD' hcC hcD' f).mp h.1, cfp_baseIso_fst P G hG hD' hcC hcD' f hb⟩,
   fun h => ⟨(cfp_lbInvertible_iff P G hG hD' hcC hcD' f).mpr h.1, hb⟩⟩

/-- **(iv)** —— base-isomorphism については **step** が射影で決まる。

★同型性の両向きは `cfp_isIso_of` / `cfp_isIso_fst`。 -/
theorem cfp_step_iff {X Y : CfpCat P G} (f : X ⟶ Y)
    (hb : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f) :
    IsStep (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsStep P (CfpCat.fst f) := by
  constructor
  · rintro ⟨hs, hni⟩
    exact ⟨(cfp_preStep_iff P G hG hD' hcC hcD' f hb).mp hs,
      fun hi => hni (cfp_isIso_of P G f hi hb)⟩
  · rintro ⟨hs, hni⟩
    exact ⟨(cfp_preStep_iff P G hG hD' hcC hcD' f hb).mpr hs,
      fun hi => hni (cfp_isIso_fst P G f hi)⟩

/-- ★`𝒞'` の base-identity 自己射は、`𝒟'` 成分が `𝟙` であること。 -/
theorem cfp_baseIdentity_iff {A : CfpCat P G} (e : A ⟶ A) :
    IsBaseIdentity (cfpPreFrobenioid P G hG hD' hcC hcD') e ↔ CfpCat.snd e = 𝟙 A.obj.right :=
  Iff.rfl

/-- ★その `𝒞` 成分は `𝒞` の base-identity。 -/
theorem cfp_baseIdentity_fst {A : CfpCat P G} (e : A ⟶ A)
    (h : IsBaseIdentity (cfpPreFrobenioid P G hG hD' hcC hcD') e) :
    IsBaseIdentity P (CfpCat.fst e) := by
  haveI hA : IsIso A.obj.hom := A.property
  have hsq := cfp_square e
  rw [(cfp_baseIdentity_iff P G hG hD' hcC hcD' e).mp h, G.map_id, Category.comp_id] at hsq
  show P.proj.map (CfpCat.fst e) = P.proj.map (𝟙 A.obj.left)
  rw [P.proj.map_id]
  exact (cancel_mono A.obj.hom).mp (hsq.trans (Category.id_comp _).symm)

/-- **(iv)** —— 射影は **`𝒪^▷` のモノイド同型**を誘導する。

★原文の当該行(`functor C′ →C determines a bijection of monoids O▷(A′)`)は
**`′`(prime) と `▷` を含むため逐語照合できない**。書き換えず、照合できない事実として記す。

★★**基準の5例目**: `base-identity` の定義が **`𝒟'` 成分を `𝟙` に固定する**ので、
`𝒪^▷` の元は `𝒞` 成分だけで決まる。**鎖どころか、`𝒟'` 成分が一意に定まる。** -/
def cfpOTriEquiv (A : CfpCat P G) :
    OTri (cfpPreFrobenioid P G hG hD' hcC hcD') A ≃* OTri P A.obj.left where
  toFun e := ⟨CfpCat.fst (e : End A), cfp_baseIdentity_fst P G hG hD' hcC hcD' _ e.2.1, e.2.2⟩
  invFun e :=
    ⟨InducedCategory.homMk ⟨(e : End A.obj.left), 𝟙 _, by
      have h1 : P.proj.map (e : End A.obj.left) = 𝟙 _ := e.2.1.trans (P.Base_id _)
      rw [h1, G.map_id, Category.comp_id, Category.id_comp]⟩,
     by show CfpCat.snd _ = 𝟙 _; rfl, e.2.2⟩
  left_inv e := by
    refine Subtype.ext (InducedCategory.hom_ext (CommaMorphism.ext rfl ?_))
    exact ((cfp_baseIdentity_iff P G hG hD' hcC hcD' (e : End A)).mp e.2.1).symm
  right_inv e := Subtype.ext rfl
  map_mul' x y := rfl

/-! ## ★(v) —— 対象タイプ

原文 (FrdI p.28):
> and only if it projects to such an object of C.

★★**仕分け**(基準を対象タイプ向けに変形したもの):

* **(a) 鎖型** —— 定義に現れる射が base-isomorphism なので中間対象が持ち上がる。
  `isotropic`(isometric pre-step) / `sub-quasi-Frobenius-trivial`(pre-step) /
  `Frobenius-isotropic`(Frobenius 型) / `perfect`(仮定の `base-isomorphic` が鎖)
* **(b) `𝒟'` 成分固定型** —— `base-identity` が `𝒟'` 成分を `𝟙` に固定する。
  `Frobenius-trivial` / `quasi-Frobenius-trivial` / `Frobenius-normalized` /
  `unit-trivial`(`𝒪^×` 経由) / `group-like`(`Φ` の同型経由)
* ★**(c) 未解決** —— `metrically trivial` / `base-trivial`。
  **結論が「`Nonempty (X ≅ A)`」**で、`𝒞'` の同型を作るには四角形が
  `𝒞` 成分の `Base` を1つに指定してしまう。**`Aut-ample` 相当**が要り、仮定にない。
-/

/-- **(v)** —— **isotropic** は射影で決まる(★(a) 鎖型)。 -/
theorem cfp_isotropic_iff (A : CfpCat P G) :
    IsIsotropic (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔ IsIsotropic P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · intro h Dd₀ φ₀ hi hs
    haveI hφb : IsIso (P.proj.map φ₀) := hs.2
    have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
    refine cfp_isIso_fst P G
      (X := A) (Y := ⟨⟨Dd₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩)
      (InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩) ?_
    refine h ⟨⟨Dd₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩
      (InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩)
      ((cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hi) ⟨hs.1, by show IsIso (𝟙 _); infer_instance⟩
  · intro h Dd' φ hi hs
    exact cfp_isIso_of P G φ
      (h Dd'.obj.left (CfpCat.fst φ) ((cfp_isometric_iff P G hG hD' hcC hcD' φ).mp hi)
        ((cfp_preStep_iff P G hG hD' hcC hcD' φ hs.2).mp hs)) hs.2

/-! ### ★★(c) の 2 件 —— **片向きだけ実装する**(2026-08-16)

★`metrically trivial` / `base-trivial` は結論が `Nonempty (X ≅ A)` である。

★★**`⟸`(`𝒞` で成り立つ ⟹ `𝒞'` で成り立つ)は出ない** ——
`𝒞'` の同型を作るには四角形
`A.hom ∘ Base f = G(g) ∘ Dd'.hom` を満たす `(f, g)` が要り、
★**`f` の `Base` が 1 つに指定されてしまう。**
`base-trivial` が与えるのは**ある**同型であって、底を指定した同型ではない。
★**直すには `Aut-ample` 相当(または `G` の充満性)が要るが、`Proposition 1.6` の
仮定は「FSM 射を FSM 射に送る関手」だけである。**

★★**`⟹` は出る。** `𝒞` の対象 `Dd₀` を `A.obj.right` と組んで `𝒞'` の対象にし、
`𝒞'` の主張を当てて `fst` を取ればよい ——
★**この向きは「底を指定する」必要がない。**

★★**原文が (vi) を「if」(片向き)としか書いていないことに注意** ——
★**著者は (v) と (vi) で向きを書き分けている。**
★**(v) の 2 件については、我々は片向きしか出せていない。**
-/

/-- ★**(v) の `metrically trivial`(★片向き)**。 -/
theorem cfp_metricallyTrivial_mp (A : CfpCat P G)
    (h : IsMetricallyTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A) :
    IsMetricallyTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  intro Dd₀ φ₀ hc hs
  haveI hφb : IsIso (P.proj.map φ₀) := hs.2
  have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
  obtain ⟨w⟩ := h ⟨⟨Dd₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩
    (InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' _).mpr hc)
    ⟨hs.1, by show IsIso (𝟙 _); infer_instance⟩
  haveI := cfp_isIso_fst P G w.hom inferInstance
  exact ⟨asIso (CfpCat.fst w.hom)⟩

/-- ★**(v) の `base-trivial`(★片向き)**。 -/
theorem cfp_baseTrivial_mp (A : CfpCat P G)
    (h : IsBaseTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A) :
    IsBaseTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  intro Dd₀ hbi
  obtain ⟨e⟩ := hbi
  have hh : (P.toElem.obj Dd₀).base ≅ G.obj A.obj.right :=
    e.symm ≪≫ (@asIso _ _ _ _ A.obj.hom hA)
  obtain ⟨w⟩ := h ⟨⟨Dd₀, A.obj.right, hh.hom⟩, inferInstanceAs (IsIso hh.hom)⟩
    ⟨Iso.refl _⟩
  haveI := cfp_isIso_fst P G w.hom inferInstance
  exact ⟨asIso (CfpCat.fst w.hom)⟩

/-! ### ★★(c) の `⟸` —— **`Aut-ample` を足せば通る**(2026-08-16)

★**両向きが 1 点に帰着することを先に確定させた**:
`𝒞'` の同型は対 `(θ, g)` で四角形
`P.proj.map θ ≫ A.obj.hom = X.obj.hom ≫ G.map g` を満たすものなので、
★**`θ` の底が 1 つに指定される。**

* `base-trivial`: 指定される底は `X.obj.hom ≫ G.map g ≫ inv A.obj.hom`
* `metrically trivial`: `cfp_square` を解くと `inv (P.Base (fst φ))`

★★**したがって要るのはただ 1 つ**:
> **(†) 指定された同型 `v : Base X₀ ≅ Base A₀` に対し、`Base θ = v` なる同型 `θ : X₀ ≅ A₀`**

★`X₀ := A₀` と取れるので **(†) ⟺ `A₀` が `Aut-ample`** である。
★**`Definition 1.3` の 21 条からは (†) は出ない**(測定の記録は上の docstring)。
★★**ここでは原文の仮定に `Aut-ample` を足して閉じる。** -/

/-- ★★**底を指定した同型**(`Aut-ample` の言い換え) ——
「**ある**同型」を「**指定の底を持つ**同型」に直す。

★`base-trivial` / `metrically trivial` はどちらも「ある同型」しか与えないので、
★**この 1 本が (v) の残り 2 件の共通の心臓部**である。 -/
theorem cfp_iso_of_isAutAmple {A₀ X₀ : C} (haa : IsAutAmple P A₀)
    (θ₀ : X₀ ≅ A₀) (v : P.proj.obj X₀ ⟶ P.proj.obj A₀) [IsIso v] :
    ∃ θ : X₀ ⟶ A₀, IsIso θ ∧ P.proj.map θ = v := by
  haveI hb0 : IsIso (P.proj.map θ₀.hom) := by
    refine ⟨P.proj.map θ₀.inv, ?_, ?_⟩
    · rw [← P.proj.map_comp, θ₀.hom_inv_id, P.proj.map_id]
    · rw [← P.proj.map_comp, θ₀.inv_hom_id, P.proj.map_id]
  have hviso : IsIso (P.proj.map θ₀.inv ≫ v) := by
    haveI : IsIso (P.proj.map θ₀.inv) := by
      refine ⟨P.proj.map θ₀.hom, ?_, ?_⟩
      · rw [← P.proj.map_comp, θ₀.inv_hom_id, P.proj.map_id]
      · rw [← P.proj.map_comp, θ₀.hom_inv_id, P.proj.map_id]
    infer_instance
  obtain ⟨c, hci, hcb⟩ := haa (P.proj.map θ₀.inv ≫ v) hviso
  haveI := hci
  refine ⟨θ₀.hom ≫ (c : A₀ ⟶ A₀), inferInstance, ?_⟩
  rw [P.proj.map_comp,
    show P.proj.map (c : A₀ ⟶ A₀) = P.proj.map θ₀.inv ≫ v from hcb,
    ← Category.assoc, ← P.proj.map_comp, θ₀.hom_inv_id, P.proj.map_id, Category.id_comp]

/-- ★★**(v) の `metrically trivial`(`⟸`、`Aut-ample` 仮定)**。

★`𝒞'` の co-angular pre-step `φ` の `𝒞` 成分に metric triviality を当てて同型 `κ` を得、
★**`cfp_iso_of_isAutAmple` で底を `inv (Base (fst φ))` に直す。**
`𝒟'` 成分は `inv (snd φ)` を取れば四角形が `cfp_square` から出る。 -/
theorem cfp_metricallyTrivial_mpr (A : CfpCat P G) (haa : IsAutAmple P A.obj.left)
    (h : IsMetricallyTrivial P A.obj.left) :
    IsMetricallyTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
  haveI hA : IsIso A.obj.hom := A.property
  intro Dd' φ hc hs
  haveI hD : IsIso Dd'.obj.hom := Dd'.property
  have hψc : IsCoAngular P (CfpCat.fst φ) := (cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hc
  have hψs : IsPreStep P (CfpCat.fst φ) := (cfp_preStep_iff P G hG hD' hcC hcD' φ hs.2).mp hs
  haveI hψb : IsIso (P.proj.map (CfpCat.fst φ)) := hψs.2
  haveI hsb : IsIso (CfpCat.snd φ) := hs.2
  obtain ⟨κ⟩ := h Dd'.obj.left (CfpCat.fst φ) hψc hψs
  obtain ⟨θ, hθi, hθv⟩ :=
    cfp_iso_of_isAutAmple P haa κ (inv (P.proj.map (CfpCat.fst φ)))
  haveI := hθi
  have hsq : P.proj.map θ ≫ A.obj.hom
      = Dd'.obj.hom ≫ G.map (inv (CfpCat.snd φ)) := by
    rw [hθv, G.map_inv, IsIso.eq_comp_inv, Category.assoc, ← cfp_square φ,
      ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  haveI : IsIso (cfpHom P G (X := Dd') (Y := A) θ (inv (CfpCat.snd φ)) hsq) :=
    cfp_isIso_of P G _ hθi (inferInstanceAs (IsIso (inv (CfpCat.snd φ))))
  exact ⟨asIso (cfpHom P G (X := Dd') (Y := A) θ (inv (CfpCat.snd φ)) hsq)⟩

/-- ★★**(v) の `base-trivial`(`⟸`、`Aut-ample` 仮定)**。

★`𝒟'` 成分の同型 `e'` を `G` で送って `𝒞` 成分の底同型 `v` を作り、
base-triviality で**ある**同型を取ってから
★**`cfp_iso_of_isAutAmple` で底を `v` に直す。** -/
theorem cfp_baseTrivial_mpr (A : CfpCat P G) (haa : IsAutAmple P A.obj.left)
    (h : IsBaseTrivial P A.obj.left) :
    IsBaseTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
  haveI hA : IsIso A.obj.hom := A.property
  intro X hbi
  haveI hX : IsIso X.obj.hom := X.property
  obtain ⟨e'⟩ := hbi
  have he' : A.obj.right ≅ X.obj.right := e'
  haveI hv : IsIso (X.obj.hom ≫ G.map he'.inv ≫ inv A.obj.hom) := inferInstance
  obtain ⟨θ₀⟩ := h X.obj.left
    ⟨(@asIso _ _ _ _ (X.obj.hom ≫ G.map he'.inv ≫ inv A.obj.hom) hv).symm⟩
  obtain ⟨θ, hθi, hθv⟩ :=
    cfp_iso_of_isAutAmple P haa θ₀ (X.obj.hom ≫ G.map he'.inv ≫ inv A.obj.hom)
  haveI := hθi
  have hsq : P.proj.map θ ≫ A.obj.hom = X.obj.hom ≫ G.map he'.inv := by
    rw [hθv, Category.assoc, Category.assoc, IsIso.inv_hom_id, Category.comp_id]
  haveI : IsIso (cfpHom P G (X := X) (Y := A) θ he'.inv hsq) :=
    cfp_isIso_of P G _ hθi (inferInstanceAs (IsIso he'.inv))
  exact ⟨asIso (cfpHom P G (X := X) (Y := A) θ he'.inv hsq)⟩

/-- ★★★**上の 2 本の仮定は空虚ではない**。

★★**「仮定を足せば何でも証明できる」への歯止め**である ——
足した `Aut-ample` が `base-trivial` と両立しないなら上の 2 本は**空虚に真**になるが、
★**両方を満たす対象は実在する**: `𝔽_Φ` の対象はすべてそうである
(`Proposition 1.5, (i)`、`elemFrob_autAmple` / `elemFrob_baseTrivial`)。

★あわせて記録: 上の 2 本は `haa` と `h` を**どちらも使う**(片方だけでは通らない)。
★そして `.src` は**付けない** —— `Proposition 1.6` の (v) を完全に実装したという
主張ではないからである。穴は `Gap/FrdI/Section1.lean` の `Gap_1_6_v` のまま。 -/
theorem autAmple_and_baseTrivial_nonvacuous {Φ₀ : MonoidOn.{v, u, w} D}
    (hD₀ : IsTotallyEpimorphic D) [IsConnected D]
    (hpd₀ : ∀ A : D, IsPreDivisorial (Φ₀.val A)) (A : ElemFrobCat Φ₀) :
    IsAutAmple (elemPreFrobenioid Φ₀ hD₀ hpd₀) A ∧
      IsBaseTrivial (elemPreFrobenioid Φ₀ hD₀ hpd₀) A :=
  ⟨elemFrob_autAmple Φ₀ hD₀ hpd₀ A, elemFrob_baseTrivial Φ₀ hD₀ hpd₀ A⟩

/-- **(v)** —— **Frobenius-isotropic** は射影で決まる(★(a) 鎖型)。 -/
theorem cfp_frobIsotropic_iff (A : CfpCat P G) :
    IsFrobeniusIsotropic (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔
      IsFrobeniusIsotropic P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · rintro ⟨Dd', φ, hft, hiso⟩
    exact ⟨Dd'.obj.left, CfpCat.fst φ, (cfp_frobType_iff P G hG hD' hcC hcD' φ hft.2).mp hft,
      (cfp_isotropic_iff P G hG hD' hcC hcD' Dd').mp hiso⟩
  · rintro ⟨Dd₀, φ₀, hft, hiso⟩
    haveI hφb : IsIso (P.proj.map φ₀) := hft.2
    have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
    refine ⟨⟨⟨Dd₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩,
      InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩, ?_, ?_⟩
    · exact (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hft
    · exact (cfp_isotropic_iff P G hG hD' hcC hcD' _).mpr hiso

/-- **(v)** —— **quasi-Frobenius-trivial** は射影で決まる(★(b) `𝒟'` 成分固定型)。 -/
theorem cfp_quasiFrobTrivial_iff (A : CfpCat P G) :
    IsQuasiFrobeniusTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔
      IsQuasiFrobeniusTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · intro h n
    obtain ⟨φ, hbi, hdeg⟩ := h n
    exact ⟨CfpCat.fst φ, cfp_baseIdentity_fst P G hG hD' hcC hcD' φ hbi, hdeg⟩
  · intro h n
    obtain ⟨φ₀, hbi, hdeg⟩ := h n
    have hid : P.proj.map φ₀ = 𝟙 _ := hbi.trans (P.Base_id _)
    exact ⟨InducedCategory.homMk ⟨φ₀, 𝟙 _, by rw [hid, G.map_id, Category.comp_id,
      Category.id_comp]⟩, by show CfpCat.snd _ = 𝟙 _; rfl, hdeg⟩

/-- **(v)** —— **Frobenius-trivial** は射影で決まる(★(b) 型)。

★`ζ : ℕ≥1 →* End A` を運ぶとき、**`𝒟'` 成分をすべて `𝟙` に取る**ので
モノイド準同型であることが自動になる。 -/
theorem cfp_frobTrivial_iff (A : CfpCat P G) :
    IsFrobeniusTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔ IsFrobeniusTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · rintro ⟨ζ, hdeg, hprop⟩
    refine ⟨⟨⟨fun n => CfpCat.fst (ζ n), ?_⟩, ?_⟩, hdeg, fun n =>
      ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ (hprop n).1,
       (cfp_frobType_iff P G hG hD' hcC hcD' _ (hprop n).2.2).mp (hprop n).2⟩⟩
    · show CfpCat.fst (ζ 1) = 𝟙 _
      rw [ζ.map_one]; rfl
    · intro m n
      show CfpCat.fst (ζ (m * n)) = CfpCat.fst (ζ n) ≫ CfpCat.fst (ζ m)
      rw [ζ.map_mul]; rfl
  · rintro ⟨ζ₀, hdeg, hprop⟩
    have hsq : ∀ n : ℕ+, P.proj.map (ζ₀ n) ≫ A.obj.hom
        = A.obj.hom ≫ G.map (𝟙 A.obj.right) := fun n => by
      rw [show P.proj.map (ζ₀ n) = 𝟙 _ from (hprop n).1.trans (P.Base_id _),
        G.map_id, Category.comp_id, Category.id_comp]
    refine ⟨⟨⟨fun n => InducedCategory.homMk ⟨ζ₀ n, 𝟙 _, hsq n⟩, ?_⟩, ?_⟩, hdeg, fun n =>
      ⟨by show CfpCat.snd _ = 𝟙 _; rfl,
       (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr (hprop n).2⟩⟩
    · refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
      · show (ζ₀ 1 : A.obj.left ⟶ A.obj.left) = 𝟙 _
        rw [ζ₀.map_one]; rfl
      · show (𝟙 A.obj.right) = 𝟙 _; rfl
    · intro m n
      refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
      · show (ζ₀ (m * n) : A.obj.left ⟶ A.obj.left) = ζ₀ n ≫ ζ₀ m
        rw [ζ₀.map_mul]; rfl
      · show (𝟙 A.obj.right) = 𝟙 _ ≫ 𝟙 _
        simp

/-- **(v)** —— **sub-quasi-Frobenius-trivial** は射影で決まる(★(a) 鎖型)。 -/
theorem cfp_subQuasiFrobTrivial_iff (A : CfpCat P G) :
    IsSubQuasiFrobeniusTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔
      IsSubQuasiFrobeniusTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · rintro ⟨Dd', α, hca, hps, hq⟩
    exact ⟨Dd'.obj.left, CfpCat.fst α, (cfp_coAngular_iff P G hG hD' hcC hcD' α).mp hca,
      (cfp_preStep_iff P G hG hD' hcC hcD' α hps.2).mp hps,
      (cfp_quasiFrobTrivial_iff P G hG hD' hcC hcD' Dd').mp hq⟩
  · rintro ⟨Dd₀, α₀, hca, hps, hq⟩
    haveI hαb : IsIso (P.proj.map α₀) := hps.2
    have hzi : IsIso (P.proj.map α₀ ≫ A.obj.hom) := inferInstance
    refine ⟨⟨⟨Dd₀, A.obj.right, P.proj.map α₀ ≫ A.obj.hom⟩, hzi⟩,
      InducedCategory.homMk ⟨α₀, 𝟙 _, by simp⟩, ?_, ⟨hps.1, by show IsIso (𝟙 _); infer_instance⟩,
      ?_⟩
    · exact cfp_coAngular_of P G hG hD' hcC hcD' _ hca
    · exact (cfp_quasiFrobTrivial_iff P G hG hD' hcC hcD' _).mpr hq

/-- **(v)** —— **unit-trivial** は射影で決まる(★(b) `𝒟'` 成分固定型)。 -/
theorem cfp_unitTrivial_iff (A : CfpCat P G) :
    IsUnitTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔ IsUnitTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · intro h
    refine Submonoid.eq_bot_iff_forall _ |>.mpr ?_
    intro x₀ hx₀
    have hid : P.proj.map (x₀ : A.obj.left ⟶ A.obj.left) = 𝟙 _ :=
      hx₀.1.1.trans (P.Base_id _)
    have hsq : P.proj.map (x₀ : A.obj.left ⟶ A.obj.left) ≫ A.obj.hom
        = A.obj.hom ≫ G.map (𝟙 A.obj.right) := by
      rw [hid, G.map_id, Category.comp_id, Category.id_comp]
    haveI hxi : IsIso (x₀ : A.obj.left ⟶ A.obj.left) := (isUnit_iff_isIso _).mp hx₀.2
    have hmem : (InducedCategory.homMk ⟨(x₀ : A.obj.left ⟶ A.obj.left), 𝟙 _, hsq⟩ : End A)
        ∈ OTimes (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
      refine ⟨⟨by show CfpCat.snd _ = 𝟙 _; rfl, hx₀.1.2⟩, (isUnit_iff_isIso _).mpr ?_⟩
      exact cfp_isIso_of P G _ hxi (by show IsIso (𝟙 _); infer_instance)
    have := (Submonoid.eq_bot_iff_forall _).mp h _ hmem
    exact congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) this
  · intro h
    refine Submonoid.eq_bot_iff_forall _ |>.mpr ?_
    intro x hx
    have hfst : CfpCat.fst (x : End A) ∈ OTimes P A.obj.left := by
      refine ⟨⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hx.1.1, hx.1.2⟩, (isUnit_iff_isIso _).mpr ?_⟩
      exact cfp_isIso_fst P G _ ((isUnit_iff_isIso _).mp hx.2)
    have h1 : CfpCat.fst (x : End A) = 𝟙 _ := (Submonoid.eq_bot_iff_forall _).mp h _ hfst
    exact InducedCategory.hom_ext (CommaMorphism.ext h1 hx.1.1)

/-! ### ★★**第3の基準** —— モノイドの性質として定義された型

★(a) 鎖 / (b) `𝒟'` 成分固定 は「**射が絡む定義**」用の基準である。
`group-like` のように **`Φ(A_𝒟)` というモノイドの性質**として定義された型には効かない。
そこに効くのは:

★**底の同型 `α : Base A ≅ G A'` が誘導する加法モノイド同型 `Φ(α)` に沿って移す**。

★**基準が3つに増えた理由**は「定義の形が3種類ある」ことにある:
射の存在/全称(→鎖)、自己射の条件(→`𝒟'` 成分固定)、底のモノイドの条件(→`Φ(α)`)。 -/

/-- **(v)** —— **group-like** は射影で決まる(★第3の基準)。 -/
theorem cfp_groupLike_iff (A : CfpCat P G) :
    IsGroupLikeObj (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔ IsGroupLikeObj P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  show IsGroupLike (Φ.val (G.obj A.obj.right)) ↔ IsGroupLike (Φ.val (P.proj.obj A.obj.left))
  rw [isGroupLike_iff, isGroupLike_iff]
  constructor
  · intro h a
    have := (h (Φ.map (inv A.obj.hom) a)).map (Φ.map A.obj.hom)
    rwa [← Φ.map_comp (inv A.obj.hom) A.obj.hom, IsIso.hom_inv_id, Φ.map_id] at this
  · intro h a
    have := (h (Φ.map A.obj.hom a)).map (Φ.map (inv A.obj.hom))
    rwa [← Φ.map_comp A.obj.hom (inv A.obj.hom), IsIso.inv_hom_id, Φ.map_id] at this

/-- ★`𝒞'` の自己射モノイドから `𝒞` のそれへのモノイド準同型(★`End` の積は `x * y = y ≫ x`)。 -/
def cfpEndHom (A : CfpCat P G) : End A →* End A.obj.left where
  toFun f := CfpCat.fst f
  map_one' := rfl
  map_mul' _ _ := rfl

theorem cfpEndHom_pow {A : CfpCat P G} (x : End A) (k : ℕ) :
    cfpEndHom P G A (x ^ k) = (cfpEndHom P G A x) ^ k :=
  map_pow (cfpEndHom P G A) x k

/-! ★**測定**: `CfpCat.fst (x ^ k) = (CfpCat.fst x) ^ k` を
`CfpCat.fst` の綴りで述べようとすると、型注釈を付けても
`HPow (A.obj.left ⟶ A.obj.left) ℕ` が合成できない
(`End X` の `Monoid` インスタンスは `End X` の綴りにしか付かない)。
★**したがって `cfpEndHom_pow`(返り値の型が `End A.obj.left` である形)を
`show` で当てるのが正しい**。`Frobenius-normalized` はこの形で書き直せるはずだが、未着手。 -/

/-- ★`𝒟'` 成分の方のモノイド準同型。 -/
def cfpEndHomSnd (A : CfpCat P G) : End A →* End A.obj.right where
  toFun f := CfpCat.snd f
  map_one' := rfl
  map_mul' _ _ := rfl

theorem cfpEndHomSnd_pow {A : CfpCat P G} (x : End A) (k : ℕ) :
    cfpEndHomSnd P G A (x ^ k) = (cfpEndHomSnd P G A x) ^ k :=
  map_pow (cfpEndHomSnd P G A) x k

/-! ### ★(参考) `Frobenius-normalized` の以前の記録

★**数学は自明に近い**: `𝒪^▷` の元も base-identity 自己射も `𝒟'` 成分が `𝟙` に固定されるので、
`𝒞` の等式がそのまま両成分に降りる。

★**Lean で止まる理由(特定済み)**: 原文の条件は `α ^ deg_Fr(φ)` を含み、
**`End A` と `A ⟶ A` は同じ型の2つの綴りで、`Monoid`(したがって `^`)は
`End A` の綴りにしか付かない**。`CfpCat.fst` の引数型は `X ⟶ Y` なので、
`CfpCat.fst (α ^ k)` と書くと `α ^ k` が `A ⟶ A` として推論され
`HPow (A ⟶ A) ℕ` が合成できない。

★試して駄目だったもの: 型注釈 `(… : End A)` / `cfpEndHom`(返り値型 `End A.obj.left`)経由 /
`show` で目標を書き換え / `map_mul`(`≫` は `*` ではないので不適) / `one_pow`(`𝟙` と `1` の綴り)。
★**`cfpEndHom` / `cfpEndHom_pow` / `cfpEndHomSnd` / `cfpEndHomSnd_pow` は残してある。**
次に当たる人が同じ道を繰り返さずに済む形で記録する。
-/


/-- ★`𝒞'` で base-isomorphic なら `𝒞` でも base-isomorphic(★片向き)。

★★**逆は言えない**: `𝒞` の同型 `proj A.left ≅ proj B.left` から
`A.right ≅ B.right`(`𝒟'` の同型)を作るには **`G` が同型を反映する**必要があり、
仮定にない。★これも「未解決 4 件」と同じ根である。 -/
theorem cfp_baseIsomorphic_of (A B : CfpCat P G)
    (h : BaseIsomorphic (cfpPreFrobenioid P G hG hD' hcC hcD') A B) :
    BaseIsomorphic P A.obj.left B.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  obtain ⟨w⟩ := h
  exact ⟨(@asIso _ _ _ _ A.obj.hom hA) ≪≫ G.mapIso w ≪≫ (@asIso _ _ _ _ B.obj.hom hB).symm⟩

/-! ## ★(vi) —— **片向きだけ**の3タイプ

原文 (FrdI p.28):
> (vi) A object of C is Aut-ample (respectively, Autsub-ample; End-ample) if

★★**なぜ片向きなのか**が基準で説明できる:
`Aut-ample` 等は「**`𝒟` の自己射が `𝒞` から来る**」という**全射性**の条件である。
`𝒞'` について言うには `𝒟'` の自己射 `g` を `𝒟` へ送って(`G g`)`𝒞` の全射性を使えばよい ——
**`G` は関手なので送れる**。
★逆向きは「`𝒟` の自己射 `g₀` に対応する `𝒟'` の自己射」を要求するので
**`G` の充満性**が要り、仮定にない。★**原文が片向きでしか述べていないのは正しい。** -/

/-! ### ★★**#5 の原因の特定** —— 「インスタンス合成の失敗」は独立の症状ではなかった

★`IsAutAmple P' A` の `g` の型は `End ((cfpToElem …).obj A).base` であり、
**`A.obj.right` に簡約されるのは `cfpToElem` を展開したあと**である。
したがって `G.map g` の型が「もう一つの綴り」になり、
`A.obj.hom ≫ G.map g ≫ w` の `IsIso` を探すときに
文脈の `hA : IsIso A.obj.hom` と**綴りが合わず**合成が失敗する。

★★**つまり #5 は #1(2つの綴り)が「インスタンス探索」を通して現れたもの**であり、
対処も #3 と同じ ——**綴りの決まった変数を先に導入する**。
★**表は 7 行のままでよい。#5 の「原因」欄を #1 と同じに直すのが正しい。** -/

/-- **(vi)** —— **End-ample** は射影から降りる(★片向き)。 -/
theorem cfp_endAmple_of (A : CfpCat P G) (h : IsEndAmple P A.obj.left) :
    IsEndAmple (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨w, hw1, hw2⟩ := hA.out
  intro g0
  -- ★#3: 綴りの決まった変数を先に導入する
  obtain ⟨g, rfl⟩ : ∃ g : End A.obj.right, g = g0 := ⟨g0, rfl⟩
  obtain ⟨φ₀, hφ₀⟩ := h (A.obj.hom ≫ G.map g ≫ w)
  refine ⟨InducedCategory.homMk ⟨φ₀, g, ?_⟩, rfl⟩
  show P.proj.map φ₀ ≫ A.obj.hom = A.obj.hom ≫ G.map g
  rw [show P.proj.map φ₀ = A.obj.hom ≫ G.map g ≫ w from hφ₀, Category.assoc,
    Category.assoc, hw2, Category.comp_id]

/-- **(vi)** —— **Aut-ample** は射影から降りる(★片向き)。 -/
theorem cfp_autAmple_of (A : CfpCat P G) (h : IsAutAmple P A.obj.left) :
    IsAutAmple (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨w, hw1, hw2⟩ := hA.out
  haveI hwi : IsIso w := ⟨A.obj.hom, hw2, hw1⟩
  intro g0 hg0
  obtain ⟨g, rfl⟩ : ∃ g : End A.obj.right, g = g0 := ⟨g0, rfl⟩
  haveI hgi : IsIso g := hg0
  haveI hGg : IsIso (G.map g) := inferInstance
  haveI hcomp : IsIso (A.obj.hom ≫ G.map g ≫ w) := inferInstance
  obtain ⟨φ₀, hiso, hφ₀⟩ := h (A.obj.hom ≫ G.map g ≫ w) hcomp
  refine ⟨InducedCategory.homMk ⟨φ₀, g, ?_⟩, ?_, rfl⟩
  · show P.proj.map φ₀ ≫ A.obj.hom = A.obj.hom ≫ G.map g
    rw [show P.proj.map φ₀ = A.obj.hom ≫ G.map g ≫ w from hφ₀, Category.assoc,
      Category.assoc, hw2, Category.comp_id]
  · exact cfp_isIso_of P G _ hiso hg0

/-! ## ★(ii) —— `Definition 1.3` の 21 条の移送

原文 (FrdI p.28):
> equivalences, the conditions of Definition 1.3 follow via a routine verification. Thus,

★**机上の仕分けを実装で検証する段**である。まず辞書から直に出る条から。 -/

section Core

variable (F : FrobenioidCore P)

include F in
/-- **(iii)(a)** の移送 —— co-angular は合成で閉じる。 -/
theorem cfp_coAngularComp {X Y Z : CfpCat P G} (ψ : X ⟶ Y) (φ : Y ⟶ Z) :
    IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') ψ →
      IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ →
      IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') (ψ ≫ φ) := by
  intro hψ hφ
  refine cfp_coAngular_of P G hG hD' hcC hcD' _ ?_
  exact F.coAngularComp (CfpCat.fst ψ) (CfpCat.fst φ)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' ψ).mp hψ) ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hφ)

include F in
/-- **(iii)(b)** の移送。 -/
theorem cfp_coAngularOfPreStep {X Y : CfpCat P G} (α : X ⟶ Y)
    (hca : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hps : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (φ : X ⟶ Y) : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ :=
  cfp_coAngular_of P G hG hD' hcC hcD' φ
    (F.coAngularOfPreStep (CfpCat.fst α) ((cfp_coAngular_iff P G hG hD' hcC hcD' α).mp hca)
      ((cfp_preStep_iff P G hG hD' hcC hcD' α hps.2).mp hps) (CfpCat.fst φ))

include F in
/-- **(v)(a)** の移送 —— pre-step は mono。

★**両成分がそれぞれ mono** であればよい: `𝒞` 側は `F.preStepMono`、
`𝒟'` 側は **pre-step の定義から `snd` が同型**。 -/
theorem cfp_preStepMono {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hφ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ) : Mono φ := by
  haveI hm : Mono (CfpCat.fst φ) :=
    F.preStepMono (CfpCat.fst φ) ((cfp_preStep_iff P G hG hD' hcC hcD' φ hφ.2).mp hφ)
  haveI hi : IsIso (CfpCat.snd φ) := hφ.2
  constructor
  intro Z g h hgh
  have e1 : CfpCat.fst g ≫ CfpCat.fst φ = CfpCat.fst h ≫ CfpCat.fst φ :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hgh
  have e2 : CfpCat.snd g ≫ CfpCat.snd φ = CfpCat.snd h ≫ CfpCat.snd φ :=
    congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) hgh
  exact InducedCategory.hom_ext
    (CommaMorphism.ext ((cancel_mono (CfpCat.fst φ)).mp e1)
      ((cancel_mono (CfpCat.snd φ)).mp e2))

include F in
/-- **(vii)(b)** の移送 —— isotropic な対象から出る射の終域は isotropic。 -/
theorem cfp_isotropicClosed {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsIsotropic (cfpPreFrobenioid P G hG hD' hcC hcD') X) :
    IsIsotropic (cfpPreFrobenioid P G hG hD' hcC hcD') Y :=
  (cfp_isotropic_iff P G hG hD' hcC hcD' Y).mpr
    (F.isotropicClosed (CfpCat.fst φ) ((cfp_isotropic_iff P G hG hD' hcC hcD' X).mp h))

include F in
/-- **(ii)** の移送 —— 各次数の Frobenius 型射が存在する。

★Frobenius 型は base-isomorphism なので**鎖**があり、
新しい対象 `B` の `𝒟'` 成分は `A` のものを流用して `snd φ = 𝟙` に取れる。 -/
theorem cfp_frobDegSurj (A : CfpCat P G) (n : ℕ+) :
    ∃ (B : CfpCat P G) (φ : A ⟶ B),
      IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') φ ∧
        (cfpPreFrobenioid P G hG hD' hcC hcD').degFr φ = n := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨B₀, φ₀, hft, hdeg⟩ := F.frobDegSurj A.obj.left n
  haveI hφb : IsIso (P.proj.map φ₀) := hft.2
  have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
  refine ⟨⟨⟨B₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩,
    InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩, ?_, hdeg⟩
  exact (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hft

include F in
/-- **(v)(b)** の移送 —— pre-step は「co-angular pre-step ≫ isometric pre-step」に分解する。

★中間対象は**両側とも pre-step に挟まれる**ので鎖がある。 -/
theorem cfp_preStepFactor {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hφ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ) :
    ∃ (Z : CfpCat P G) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') α ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  haveI hX : IsIso X.obj.hom := X.property
  obtain ⟨Z₀, β₀, α₀, hfac, hβc, hβs, hαi, hαs⟩ :=
    F.preStepFactor (CfpCat.fst φ) ((cfp_preStep_iff P G hG hD' hcC hcD' φ hφ.2).mp hφ)
  haveI hβb : IsIso (P.proj.map β₀) := hβs.2
  have hzi : IsIso (inv (P.proj.map β₀) ≫ X.obj.hom) := inferInstance
  refine ⟨⟨⟨Z₀, X.obj.right, inv (P.proj.map β₀) ≫ X.obj.hom⟩, hzi⟩,
    InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : X.obj.right ⟶ Y.obj.right), ?_⟩, ?_, ?_,
    ⟨hβs.1, by show IsIso (𝟙 _); infer_instance⟩, ?_, ?_⟩
  · show P.proj.map α₀ ≫ Y.obj.hom
      = (inv (P.proj.map β₀) ≫ X.obj.hom) ≫ G.map (CfpCat.snd φ)
    rw [Category.assoc, ← cfp_square φ,
      show P.proj.map (CfpCat.fst φ) = P.proj.map β₀ ≫ P.proj.map α₀ from by
        rw [hfac, P.proj.map_comp],
      ← Category.assoc, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
    show CfpCat.snd φ = 𝟙 _ ≫ CfpCat.snd φ
    simp
  · exact cfp_coAngular_of P G hG hD' hcC hcD' _ hβc
  · exact (cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hαi
  · exact ⟨hαs.1, hφ.2⟩

include F in
/-- **(v)(c)** の移送 —— pre-step は「isometric pre-step ≫ co-angular pre-step」に分解する。 -/
theorem cfp_preStepFactor' {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hφ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ) :
    ∃ (Z : CfpCat P G) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') α ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  haveI hX : IsIso X.obj.hom := X.property
  obtain ⟨Z₀, β₀, α₀, hfac, hβi, hβs, hαc, hαs⟩ :=
    F.preStepFactor' (CfpCat.fst φ) ((cfp_preStep_iff P G hG hD' hcC hcD' φ hφ.2).mp hφ)
  haveI hβb : IsIso (P.proj.map β₀) := hβs.2
  have hzi : IsIso (inv (P.proj.map β₀) ≫ X.obj.hom) := inferInstance
  refine ⟨⟨⟨Z₀, X.obj.right, inv (P.proj.map β₀) ≫ X.obj.hom⟩, hzi⟩,
    InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : X.obj.right ⟶ Y.obj.right), ?_⟩, ?_, ?_,
    ⟨hβs.1, by show IsIso (𝟙 _); infer_instance⟩, ?_, ?_⟩
  · show P.proj.map α₀ ≫ Y.obj.hom
      = (inv (P.proj.map β₀) ≫ X.obj.hom) ≫ G.map (CfpCat.snd φ)
    rw [Category.assoc, ← cfp_square φ,
      show P.proj.map (CfpCat.fst φ) = P.proj.map β₀ ≫ P.proj.map α₀ from by
        rw [hfac, P.proj.map_comp],
      ← Category.assoc, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
    show CfpCat.snd φ = 𝟙 _ ≫ CfpCat.snd φ
    simp
  · exact (cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hβi
  · exact cfp_coAngular_of P G hG hD' hcC hcD' _ hαc
  · exact ⟨hαs.1, hφ.2⟩

include F in
/-- **(i)(a)** の移送 —— `𝒟'` のどの対象の上にも Frobenius-trivial な対象がある。

★`𝒞` の `baseSurj` を `G Y` に当てて得た同型を、そのまま CFP の三つ組の第3成分にする。
★**新しい対象を作るのに鎖は要らない** —— 同型が入力として与えられるから。 -/
theorem cfp_baseSurj (Y : D') :
    ∃ A : CfpCat P G, IsFrobeniusTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ∧
      Nonempty (((cfpPreFrobenioid P G hG hD' hcC hcD').toElem.obj A).base ≅ Y) := by
  obtain ⟨A₀, hft, ⟨e⟩⟩ := F.baseSurj (G.obj Y)
  haveI : IsIso e.hom := e.isIso_hom
  refine ⟨⟨⟨A₀, Y, e.hom⟩, inferInstanceAs (IsIso e.hom)⟩, ?_, ⟨Iso.refl _⟩⟩
  exact (cfp_frobTrivial_iff P G hG hD' hcC hcD' _).mpr hft

include F in
/-- **(ii)** の移送 —— 同じ次数の Frobenius 型射の本質的一意性。

★★**`𝒟'` 成分は `(snd φ)⁻¹ ≫ snd ψ` に取れる** —— Frobenius 型は base-isomorphism なので
両方の `𝒟'` 成分が同型であり、**`G` の充満性は要らない**。 -/
theorem cfp_frobDegUniq (A B E : CfpCat P G) (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hψ : IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') ψ)
    (hd : (cfpPreFrobenioid P G hG hD' hcC hcD').degFr φ = (cfpPreFrobenioid P G hG hD' hcC hcD').degFr ψ) :
    ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  haveI hE : IsIso E.obj.hom := E.property
  haveI hsφ : IsIso (CfpCat.snd φ) := hφ.2
  haveI hsψ : IsIso (CfpCat.snd ψ) := hψ.2
  haveI hpφ : IsIso (P.proj.map (CfpCat.fst φ)) := cfp_baseIso_fst P G hG hD' hcC hcD' φ hφ.2
  obtain ⟨β₀, hβiso, hβ⟩ := F.frobDegUniq A.obj.left B.obj.left E.obj.left
    (CfpCat.fst φ) (CfpCat.fst ψ)
    ((cfp_frobType_iff P G hG hD' hcC hcD' φ hφ.2).mp hφ)
    ((cfp_frobType_iff P G hG hD' hcC hcD' ψ hψ.2).mp hψ) hd
  have hsq : P.proj.map β₀ ≫ E.obj.hom
      = B.obj.hom ≫ G.map (inv (CfpCat.snd φ) ≫ CfpCat.snd ψ) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst φ))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, hβ, cfp_square ψ, ← Category.assoc,
      cfp_square φ, Category.assoc, ← G.map_comp, ← Category.assoc,
      IsIso.hom_inv_id, Category.id_comp]
  have hsnd : IsIso (inv (CfpCat.snd φ) ≫ CfpCat.snd ψ) := inferInstance
  refine ⟨InducedCategory.homMk ⟨β₀, inv (CfpCat.snd φ) ≫ CfpCat.snd ψ, hsq⟩,
    cfp_isIso_of P G _ hβiso hsnd, ?_⟩
  refine InducedCategory.hom_ext (CommaMorphism.ext hβ ?_)
  show CfpCat.snd φ ≫ inv (CfpCat.snd φ) ≫ CfpCat.snd ψ = CfpCat.snd ψ
  rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

include F in
/-- **(vii)(a)** の移送 —— isotropic hull の存在(普遍性つき)。

★hull は isometric pre-step なので**鎖**があり、新しい対象が持ち上がる。
★普遍性の `∃!` は、**`𝒟'` 成分が `snd γ` に一意に決まる**ので `𝒞` の一意性から出る。 -/
theorem cfp_isotropicHullExists (A : CfpCat P G) :
    ∃ (B : CfpCat P G) (α : A ⟶ B), IsIsotropicHull (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨B₀, α₀, hαi, hαs, hBiso, huniv⟩ := F.isotropicHullExists A.obj.left
  haveI hαb : IsIso (P.proj.map α₀) := hαs.2
  have hzi : IsIso (inv (P.proj.map α₀) ≫ A.obj.hom) := inferInstance
  refine ⟨⟨⟨B₀, A.obj.right, inv (P.proj.map α₀) ≫ A.obj.hom⟩, hzi⟩,
    InducedCategory.homMk ⟨α₀, 𝟙 _, by simp⟩,
    (cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hαi,
    ⟨hαs.1, by show IsIso (𝟙 _); infer_instance⟩,
    (cfp_isotropic_iff P G hG hD' hcC hcD' _).mpr hBiso, ?_⟩
  intro Cc hCc γ
  haveI hC : IsIso Cc.obj.hom := Cc.property
  obtain ⟨β₀, hβ₀, hβ₀u⟩ := huniv Cc.obj.left ((cfp_isotropic_iff P G hG hD' hcC hcD' Cc).mp hCc)
    (CfpCat.fst γ)
  have hsq : P.proj.map β₀ ≫ Cc.obj.hom
      = (inv (P.proj.map α₀) ≫ A.obj.hom) ≫ G.map (CfpCat.snd γ) := by
    rw [Category.assoc, ← cfp_square γ,
      show P.proj.map (CfpCat.fst γ) = P.proj.map α₀ ≫ P.proj.map β₀ from by
        rw [hβ₀, P.proj.map_comp],
      ← Category.assoc, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  refine ⟨InducedCategory.homMk ⟨β₀, (CfpCat.snd γ : A.obj.right ⟶ Cc.obj.right), hsq⟩, ?_, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext hβ₀ ?_)
    show CfpCat.snd γ = 𝟙 _ ≫ CfpCat.snd γ
    simp
  · intro β hβ
    have hf : CfpCat.fst γ = α₀ ≫ CfpCat.fst β :=
      congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hβ
    have hs : CfpCat.snd γ = 𝟙 _ ≫ CfpCat.snd β :=
      congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) hβ
    exact InducedCategory.hom_ext
      (CommaMorphism.ext (hβ₀u _ hf) (hs.trans (Category.id_comp _)).symm)

include F in
/-- ★`𝒪^▷` の元は `𝒞'` へそのまま持ち上がる(★(b) `𝒟'` 成分固定型)。 -/
theorem cfp_otri_lift (A : CfpCat P G) (x : End A.obj.left) (hx : x ∈ OTri P A.obj.left) :
    ∃ y : End A, y ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') A ∧ CfpCat.fst y = x := by
  haveI hA : IsIso A.obj.hom := A.property
  have hsq : P.proj.map (x : A.obj.left ⟶ A.obj.left) ≫ A.obj.hom
      = A.obj.hom ≫ G.map (𝟙 A.obj.right) := by
    rw [show P.proj.map (x : A.obj.left ⟶ A.obj.left) = 𝟙 _ from hx.1.trans (P.Base_id _),
      G.map_id, Category.comp_id, Category.id_comp]
  exact ⟨InducedCategory.homMk ⟨(x : A.obj.left ⟶ A.obj.left), 𝟙 _, hsq⟩,
    ⟨by show CfpCat.snd _ = 𝟙 _; rfl, hx.2⟩, rfl⟩

include F in
/-- **(iii)(c)** 順方向の移送。

★`𝒪^▷` の元は base-identity なので **`𝒟'` 成分が `𝟙` に固定**され、
`𝒞` の `∃!` がそのまま上がる。 -/
theorem cfp_otriFwd {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hca : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hps : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (α : End X) (hα : α ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') X) :
    ∃! β : End Y, β ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') Y ∧
      (φ ≫ β : X ⟶ Y) = (α : X ⟶ X) ≫ φ := by
  obtain ⟨β₀, ⟨hβ₀m, hβ₀e⟩, hβ₀u⟩ := F.otriFwd (CfpCat.fst φ)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hca)
    ((cfp_preStep_iff P G hG hD' hcC hcD' φ hps.2).mp hps)
    (CfpCat.fst (α : End X))
    ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hα.1, hα.2⟩
  obtain ⟨β, hβm, hβf⟩ := cfp_otri_lift P G hG hD' hcC hcD' F Y β₀ hβ₀m
  have hαs : CfpCat.snd (α : End X) = 𝟙 _ := hα.1
  have hβs : CfpCat.snd (β : End Y) = 𝟙 _ := hβm.1
  refine ⟨β, ⟨hβm, ?_⟩, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
    · show CfpCat.fst φ ≫ CfpCat.fst (β : End Y) = CfpCat.fst (α : End X) ≫ CfpCat.fst φ
      rw [hβf]
      exact hβ₀e
    · show CfpCat.snd φ ≫ CfpCat.snd (β : End Y) = CfpCat.snd (α : End X) ≫ CfpCat.snd φ
      rw [hβs, hαs, Category.comp_id, Category.id_comp]
  · rintro γ ⟨hγm, hγe⟩
    have hγf : CfpCat.fst φ ≫ CfpCat.fst (γ : End Y)
        = CfpCat.fst (α : End X) ≫ CfpCat.fst φ :=
      congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hγe
    refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
    · show CfpCat.fst (γ : End Y) = CfpCat.fst (β : End Y)
      rw [hβf]
      exact hβ₀u _ ⟨⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hγm.1, hγm.2⟩, hγf⟩
    · show CfpCat.snd (γ : End Y) = CfpCat.snd (β : End Y)
      rw [hβs]
      exact hγm.1

include F in
/-- **(iii)(c)** 逆方向の移送。 -/
theorem cfp_otriBwd {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hca : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hps : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (β : End Y) (hβ : β ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') Y) :
    ∃! α : End X, α ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') X ∧
      (φ ≫ β : X ⟶ Y) = (α : X ⟶ X) ≫ φ := by
  obtain ⟨α₀, ⟨hα₀m, hα₀e⟩, hα₀u⟩ := F.otriBwd (CfpCat.fst φ)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hca)
    ((cfp_preStep_iff P G hG hD' hcC hcD' φ hps.2).mp hps)
    (CfpCat.fst (β : End Y))
    ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hβ.1, hβ.2⟩
  obtain ⟨α, hαm, hαf⟩ := cfp_otri_lift P G hG hD' hcC hcD' F X α₀ hα₀m
  have hβs : CfpCat.snd (β : End Y) = 𝟙 _ := hβ.1
  have hαs : CfpCat.snd (α : End X) = 𝟙 _ := hαm.1
  refine ⟨α, ⟨hαm, ?_⟩, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
    · show CfpCat.fst φ ≫ CfpCat.fst (β : End Y) = CfpCat.fst (α : End X) ≫ CfpCat.fst φ
      rw [hαf]
      exact hα₀e
    · show CfpCat.snd φ ≫ CfpCat.snd (β : End Y) = CfpCat.snd (α : End X) ≫ CfpCat.snd φ
      rw [hβs, hαs, Category.comp_id, Category.id_comp]
  · rintro γ ⟨hγm, hγe⟩
    have hγf : CfpCat.fst φ ≫ CfpCat.fst (β : End Y)
        = CfpCat.fst (γ : End X) ≫ CfpCat.fst φ :=
      congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hγe
    refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
    · show CfpCat.fst (γ : End X) = CfpCat.fst (α : End X)
      rw [hαf]
      exact hα₀u _ ⟨⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hγm.1, hγm.2⟩, hγf⟩
    · show CfpCat.snd (γ : End X) = CfpCat.snd (α : End X)
      rw [hαs]
      exact hγm.1

include F in
/-- **(vi)** の移送 —— `𝒪^×` を除いた忠実性。

★`Div` の移送には **`Φ(α⁻¹)` の単射性**（`Definition 1.1, (ii), (a)`）を使う。 -/
theorem cfp_faithfulUpToUnits {X Y : CfpCat P G} (φ ψ : X ⟶ Y)
    (hb : BaseEquivalent (cfpPreFrobenioid P G hG hD' hcC hcD') φ ψ)
    (hm : MetricallyEquivalent (cfpPreFrobenioid P G hG hD' hcC hcD') φ ψ)
    (hφc : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hφs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hψc : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') ψ)
    (hψs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') ψ) :
    ∃ α : End Y, α ∈ OTimes (cfpPreFrobenioid P G hG hD' hcC hcD') Y ∧
      φ = ψ ≫ (α : Y ⟶ Y) := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  have hb' : CfpCat.snd φ = CfpCat.snd ψ := hb
  have hmC : MetricallyEquivalent P (CfpCat.fst φ) (CfpCat.fst ψ) :=
    Φ.map_injective (@inv _ _ _ _ X.obj.hom hX) hm
  have hbC : BaseEquivalent P (CfpCat.fst φ) (CfpCat.fst ψ) := by
    show P.proj.map (CfpCat.fst φ) = P.proj.map (CfpCat.fst ψ)
    rw [cfp_base_fst P G φ, cfp_base_fst P G ψ, hb']
  obtain ⟨α₀, hα₀, hφψ⟩ := F.faithfulUpToUnits (CfpCat.fst φ) (CfpCat.fst ψ) hbC hmC
    ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hφc)
    ((cfp_preStep_iff P G hG hD' hcC hcD' φ hφs.2).mp hφs)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' ψ).mp hψc)
    ((cfp_preStep_iff P G hG hD' hcC hcD' ψ hψs.2).mp hψs)
  obtain ⟨α, hαm, hαf⟩ := cfp_otri_lift P G hG hD' hcC hcD' F Y α₀ hα₀.1
  haveI hα₀i : IsIso (α₀ : Y.obj.left ⟶ Y.obj.left) := (isUnit_iff_isIso _).mp hα₀.2
  have hsnd1 : CfpCat.snd (α : End Y) = 𝟙 _ := hαm.1
  have hsndi : IsIso (CfpCat.snd (α : End Y)) := by rw [hsnd1]; infer_instance
  have hfsti : IsIso (CfpCat.fst (α : End Y)) := by rw [hαf]; infer_instance
  refine ⟨α, ⟨hαm, (isUnit_iff_isIso _).mpr (cfp_isIso_of P G _ hfsti hsndi)⟩, ?_⟩
  refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
  · show CfpCat.fst φ = CfpCat.fst ψ ≫ CfpCat.fst (α : End Y)
    rw [hαf]
    exact hφψ
  · show CfpCat.snd φ = CfpCat.snd ψ ≫ CfpCat.snd (α : End Y)
    rw [hsnd1, Category.comp_id]
    exact hb'

include F in
/-- **(iii)(c)** 全単射が `Base` にしか依らないことの移送。 -/
theorem cfp_otriBase {X Y : CfpCat P G} (φ φ' : X ⟶ Y)
    (hca : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hps : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hca' : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ')
    (hps' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ')
    (hbase : (cfpPreFrobenioid P G hG hD' hcC hcD').Base φ = (cfpPreFrobenioid P G hG hD' hcC hcD').Base φ')
    (α : End X) (hα : α ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') X)
    (β : End Y) (hβ : β ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') Y)
    (h : (φ ≫ β : X ⟶ Y) = (α : X ⟶ X) ≫ φ) :
    (φ' ≫ β : X ⟶ Y) = (α : X ⟶ X) ≫ φ' := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  have hb' : CfpCat.snd φ = CfpCat.snd φ' := hbase
  have hαs : CfpCat.snd (α : End X) = 𝟙 _ := hα.1
  have hβs : CfpCat.snd (β : End Y) = 𝟙 _ := hβ.1
  have hf : CfpCat.fst φ ≫ CfpCat.fst (β : End Y)
      = CfpCat.fst (α : End X) ≫ CfpCat.fst φ :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) h
  refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
  · refine F.otriBase (CfpCat.fst φ) (CfpCat.fst φ')
      ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hca)
      ((cfp_preStep_iff P G hG hD' hcC hcD' φ hps.2).mp hps)
      ((cfp_coAngular_iff P G hG hD' hcC hcD' φ').mp hca')
      ((cfp_preStep_iff P G hG hD' hcC hcD' φ' hps'.2).mp hps') ?_
      (CfpCat.fst (α : End X)) ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hα.1, hα.2⟩
      (CfpCat.fst (β : End Y)) ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hβ.1, hβ.2⟩ hf
    show P.proj.map (CfpCat.fst φ) = P.proj.map (CfpCat.fst φ')
    rw [cfp_base_fst P G φ, cfp_base_fst P G φ', hb']
  · show CfpCat.snd φ' ≫ CfpCat.snd (β : End Y) = CfpCat.snd (α : End X) ≫ CfpCat.snd φ'
    rw [hβs, hαs, Category.comp_id, Category.id_comp]

include F in
/-- **(iv)(a)** の移送 —— 任意の射の3分解。

★`γ`(Frobenius 型)と `β`(pre-step)は base-isomorphism なので**鎖**があり、
中間対象 2 つが持ち上がる。`α`(pull-back)は `𝒟'` 成分として `snd φ` を受け取る。 -/
theorem cfp_arbFactor {A B : CfpCat P G} (φ : A ⟶ B) :
    ∃ (X Y : CfpCat P G) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') γ ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  obtain ⟨X₀, Y₀, γ₀, β₀, α₀, hfac, hγ, hβ, hα⟩ := F.arbFactor (CfpCat.fst φ)
  haveI hγb : IsIso (P.proj.map γ₀) := hγ.2
  haveI hβb : IsIso (P.proj.map β₀) := hβ.2
  have hxi : IsIso (inv (P.proj.map γ₀) ≫ A.obj.hom) := inferInstance
  have hyi : IsIso (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom) := inferInstance
  have hfacp : P.proj.map (CfpCat.fst φ)
      = P.proj.map γ₀ ≫ P.proj.map β₀ ≫ P.proj.map α₀ := by
    rw [hfac, P.proj.map_comp, P.proj.map_comp]
  refine ⟨⟨⟨X₀, A.obj.right, inv (P.proj.map γ₀) ≫ A.obj.hom⟩, hxi⟩,
    ⟨⟨Y₀, A.obj.right, inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom⟩, hyi⟩,
    InducedCategory.homMk ⟨γ₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : A.obj.right ⟶ B.obj.right), ?_⟩, ?_, ?_, ?_, ?_⟩
  · show P.proj.map α₀ ≫ B.obj.hom
      = (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom) ≫ G.map (CfpCat.snd φ)
    rw [Category.assoc, Category.assoc, ← cfp_square φ, hfacp]
    simp only [Category.assoc]
    rw [← Category.assoc (inv (P.proj.map γ₀)), IsIso.inv_hom_id, Category.id_comp,
      ← Category.assoc (inv (P.proj.map β₀)), IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
    show CfpCat.snd φ = 𝟙 _ ≫ 𝟙 _ ≫ CfpCat.snd φ
    simp
  · exact (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hγ
  · exact ⟨hβ.1, by show IsIso (𝟙 _); infer_instance⟩
  · exact cfp_isPullBack_of P G hG hD' hcC hcD' _ hα

include F in
/-- **(v)(b)** の一意性の移送。

★**同型を作るのに `𝒟'` 成分は `(snd β)⁻¹ ≫ snd β'` に取る** ——
`β`, `β'` はどちらも pre-step なので `𝒟'` 成分が同型であり、
★**逆射も成分ごとに書き下せる**ので `inv` を CFP の射に対して使わずに済む(表 #2)。 -/
theorem cfp_preStepFactorUniq {A B : CfpCat P G} (X X' : CfpCat P G)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβc : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hβs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hαi : IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hαs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hβc' : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hβs' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hαi' : IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') α')
    (hαs' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hX : IsIso X.obj.hom := X.property
  haveI hX' : IsIso X'.obj.hom := X'.property
  haveI hsβ : IsIso (CfpCat.snd β) := hβs.2
  haveI hsβ' : IsIso (CfpCat.snd β') := hβs'.2
  haveI hpβ : IsIso (P.proj.map (CfpCat.fst β)) := cfp_baseIso_fst P G hG hD' hcC hcD' β hβs.2
  haveI hpβ' : IsIso (P.proj.map (CfpCat.fst β')) := cfp_baseIso_fst P G hG hD' hcC hcD' β' hβs'.2
  obtain ⟨γ₀, hγ1, hγ2⟩ := F.preStepFactorUniq X.obj.left X'.obj.left
    (CfpCat.fst β) (CfpCat.fst α) (CfpCat.fst β') (CfpCat.fst α')
    (congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) heq)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' β).mp hβc)
    ((cfp_preStep_iff P G hG hD' hcC hcD' β hβs.2).mp hβs)
    ((cfp_isometric_iff P G hG hD' hcC hcD' α).mp hαi)
    ((cfp_preStep_iff P G hG hD' hcC hcD' α hαs.2).mp hαs)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' β').mp hβc')
    ((cfp_preStep_iff P G hG hD' hcC hcD' β' hβs'.2).mp hβs')
    ((cfp_isometric_iff P G hG hD' hcC hcD' α').mp hαi')
    ((cfp_preStep_iff P G hG hD' hcC hcD' α' hαs'.2).mp hαs')
  -- ★四角形: `𝒞` 成分の底射は `𝒟'` 成分から決まる
  have hsq : P.proj.map γ₀.hom ≫ X'.obj.hom
      = X.obj.hom ≫ G.map (inv (CfpCat.snd β) ≫ CfpCat.snd β') := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, ← hγ2, cfp_square β', ← Category.assoc,
      cfp_square β, Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  have hγ3 : CfpCat.fst β' ≫ γ₀.inv = CfpCat.fst β := by
    rw [hγ2, Category.assoc, γ₀.hom_inv_id, Category.comp_id]
  have hsq' : P.proj.map γ₀.inv ≫ X.obj.hom
      = X'.obj.hom ≫ G.map (inv (CfpCat.snd β') ≫ CfpCat.snd β) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β'))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, hγ3, cfp_square β, ← Category.assoc,
      cfp_square β', Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  refine ⟨⟨InducedCategory.homMk ⟨γ₀.hom, inv (CfpCat.snd β) ≫ CfpCat.snd β', hsq⟩,
    InducedCategory.homMk ⟨γ₀.inv, inv (CfpCat.snd β') ≫ CfpCat.snd β, hsq'⟩, ?_, ?_⟩, ?_, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext γ₀.hom_inv_id ?_)
    show (inv (CfpCat.snd β) ≫ CfpCat.snd β') ≫ inv (CfpCat.snd β') ≫ CfpCat.snd β = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext γ₀.inv_hom_id ?_)
    show (inv (CfpCat.snd β') ≫ CfpCat.snd β) ≫ inv (CfpCat.snd β) ≫ CfpCat.snd β' = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext hγ1 ?_)
    show CfpCat.snd α' = (inv (CfpCat.snd β') ≫ CfpCat.snd β) ≫ CfpCat.snd α
    have hs : CfpCat.snd β ≫ CfpCat.snd α = CfpCat.snd β' ≫ CfpCat.snd α' :=
      congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) heq
    rw [Category.assoc, hs, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hγ2 ?_)
    show CfpCat.snd β' = CfpCat.snd β ≫ inv (CfpCat.snd β) ≫ CfpCat.snd β'
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

include F in
/-- **(v)(c)** の一意性の移送(★`preStepFactorUniq` と同じ形)。 -/
theorem cfp_preStepFactorUniq' {A B : CfpCat P G} (X X' : CfpCat P G)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβi : IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hβs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hαc : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hαs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hβi' : IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hβs' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hαc' : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') α')
    (hαs' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hX : IsIso X.obj.hom := X.property
  haveI hX' : IsIso X'.obj.hom := X'.property
  haveI hsβ : IsIso (CfpCat.snd β) := hβs.2
  haveI hsβ' : IsIso (CfpCat.snd β') := hβs'.2
  haveI hpβ : IsIso (P.proj.map (CfpCat.fst β)) := cfp_baseIso_fst P G hG hD' hcC hcD' β hβs.2
  haveI hpβ' : IsIso (P.proj.map (CfpCat.fst β')) := cfp_baseIso_fst P G hG hD' hcC hcD' β' hβs'.2
  obtain ⟨γ₀, hγ1, hγ2⟩ := F.preStepFactorUniq' X.obj.left X'.obj.left
    (CfpCat.fst β) (CfpCat.fst α) (CfpCat.fst β') (CfpCat.fst α')
    (congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) heq)
    ((cfp_isometric_iff P G hG hD' hcC hcD' β).mp hβi)
    ((cfp_preStep_iff P G hG hD' hcC hcD' β hβs.2).mp hβs)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' α).mp hαc)
    ((cfp_preStep_iff P G hG hD' hcC hcD' α hαs.2).mp hαs)
    ((cfp_isometric_iff P G hG hD' hcC hcD' β').mp hβi')
    ((cfp_preStep_iff P G hG hD' hcC hcD' β' hβs'.2).mp hβs')
    ((cfp_coAngular_iff P G hG hD' hcC hcD' α').mp hαc')
    ((cfp_preStep_iff P G hG hD' hcC hcD' α' hαs'.2).mp hαs')
  have hsq : P.proj.map γ₀.hom ≫ X'.obj.hom
      = X.obj.hom ≫ G.map (inv (CfpCat.snd β) ≫ CfpCat.snd β') := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, ← hγ2, cfp_square β', ← Category.assoc,
      cfp_square β, Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  have hγ3 : CfpCat.fst β' ≫ γ₀.inv = CfpCat.fst β := by
    rw [hγ2, Category.assoc, γ₀.hom_inv_id, Category.comp_id]
  have hsq' : P.proj.map γ₀.inv ≫ X.obj.hom
      = X'.obj.hom ≫ G.map (inv (CfpCat.snd β') ≫ CfpCat.snd β) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β'))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, hγ3, cfp_square β, ← Category.assoc,
      cfp_square β', Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  refine ⟨⟨InducedCategory.homMk ⟨γ₀.hom, inv (CfpCat.snd β) ≫ CfpCat.snd β', hsq⟩,
    InducedCategory.homMk ⟨γ₀.inv, inv (CfpCat.snd β') ≫ CfpCat.snd β, hsq'⟩, ?_, ?_⟩, ?_, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext γ₀.hom_inv_id ?_)
    show (inv (CfpCat.snd β) ≫ CfpCat.snd β') ≫ inv (CfpCat.snd β') ≫ CfpCat.snd β = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext γ₀.inv_hom_id ?_)
    show (inv (CfpCat.snd β') ≫ CfpCat.snd β) ≫ inv (CfpCat.snd β) ≫ CfpCat.snd β' = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext hγ1 ?_)
    show CfpCat.snd α' = (inv (CfpCat.snd β') ≫ CfpCat.snd β) ≫ CfpCat.snd α
    have hs : CfpCat.snd β ≫ CfpCat.snd α = CfpCat.snd β' ≫ CfpCat.snd α' :=
      congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) heq
    rw [Category.assoc, hs, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hγ2 ?_)
    show CfpCat.snd β' = CfpCat.snd β ≫ inv (CfpCat.snd β) ≫ CfpCat.snd β'
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

include F in
/-- **(i)(b)** の移送 —— `𝒟'` の同型を pre-step の span で実現する。

★`𝒟'` の同型 `α` を `G` で送り、両端の同型で共役して `𝒞` の `preStepSpan` に渡す。
★中間対象の `𝒟'` 成分は `A` のものを流用し、`φ` の `𝒟'` 成分は `𝟙`、
`ψ` の `𝒟'` 成分は `α` そのものに取れる。 -/
theorem cfp_preStepSpan (A B : CfpCat P G)
    (α : ((cfpPreFrobenioid P G hG hD' hcC hcD').toElem.obj A).base ⟶
      ((cfpPreFrobenioid P G hG hD' hcC hcD').toElem.obj B).base) (hα : IsIso α) :
    ∃ (X : CfpCat P G) (φ : X ⟶ A) (ψ : X ⟶ B)
      (hφ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ),
      IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') ψ ∧
        α = @inv _ _ _ _ ((cfpPreFrobenioid P G hG hD' hcC hcD').Base φ) hφ.2 ≫
          (cfpPreFrobenioid P G hG hD' hcC hcD').Base ψ := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  -- ★#3: 綴りの決まった変数を先に導入する
  obtain ⟨a, rfl⟩ : ∃ a : A.obj.right ⟶ B.obj.right, a = α := ⟨α, rfl⟩
  haveI hai : IsIso a := hα
  haveI hGa : IsIso (G.map a) := inferInstance
  obtain ⟨wB, hwB1, hwB2⟩ := hB.out
  haveI hwBi : IsIso wB := ⟨B.obj.hom, hwB2, hwB1⟩
  have hui : IsIso (A.obj.hom ≫ G.map a ≫ wB) := inferInstance
  obtain ⟨X₀, φ₀, ψ₀, hφ₀, hψ₀, heq⟩ :=
    F.preStepSpan A.obj.left B.obj.left (A.obj.hom ≫ G.map a ≫ wB) hui
  haveI hφb : IsIso (P.proj.map φ₀) := hφ₀.2
  have hxi : IsIso (P.proj.map φ₀ ≫ A.obj.hom) := inferInstance
  have h1 : A.obj.hom ≫ G.map a ≫ wB
      = @inv _ _ _ _ (P.proj.map φ₀) hφ₀.2 ≫ P.proj.map ψ₀ := heq
  have h2 : A.obj.hom ≫ G.map a
      = @inv _ _ _ _ (P.proj.map φ₀) hφ₀.2 ≫ P.proj.map ψ₀ ≫ B.obj.hom := by
    have h3 := congrArg (fun t => t ≫ B.obj.hom) h1
    simp only [Category.assoc] at h3
    rw [hwB2, Category.comp_id] at h3
    exact h3
  have hkey : P.proj.map ψ₀ ≫ B.obj.hom
      = (P.proj.map φ₀ ≫ A.obj.hom) ≫ G.map a := by
    rw [Category.assoc, h2, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  refine ⟨⟨⟨X₀, A.obj.right, P.proj.map φ₀ ≫ A.obj.hom⟩, hxi⟩,
    InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨ψ₀, a, hkey⟩,
    ⟨hφ₀.1, by show IsIso (𝟙 _); infer_instance⟩,
    ⟨hψ₀.1, hai⟩, ?_⟩
  show a = @inv _ _ _ _ (𝟙 A.obj.right) _ ≫ a
  rw [IsIso.inv_id, Category.id_comp]

/-- ★★**pull-back は左から簡約できる** —— `f ≫ w` と `w` が pull-back なら `f` も。

★★**`Definition 1.2, (ii)` の全単射条件だけから出る**(Frobenioid の公理は要らない)。
`Proposition 1.7, (v)` の pull-back の段は `FrobenioidCore` を仮定するが、
**この向きだけなら仮定なしで示せる**。
★これが無いと `plBkEquiv` の充満性が循環する(`𝒞'` が Frobenioid であることを要してしまう)。 -/
theorem isPullBack_of_comp_left {Cc : Type u2} [Category.{v2} Cc] {Ψ : MonoidOn.{v, u, w} D}
    (Q : PreFrobenioid Cc Ψ) {A B E : Cc} (f : A ⟶ B) (wm : B ⟶ E)
    (hw : IsPullBack Q wm) (hq : IsPullBack Q (f ≫ wm)) : IsPullBack Q f := by
  intro T
  constructor
  · intro f₁ f₂ hf
    have hp := Subtype.ext_iff.mp hf
    have e1 : (f₁ ≫ f) = f₂ ≫ f := congrArg Prod.fst hp
    have e2 : Q.Base f₁ = Q.Base f₂ := congrArg Prod.snd hp
    refine (hq T).1 (Subtype.ext (Prod.ext ?_ e2))
    show (f₁ ≫ f ≫ wm) = f₂ ≫ f ≫ wm
    rw [← Category.assoc, e1, Category.assoc]
  · rintro ⟨⟨g, u⟩, hcond⟩
    have hcond' : Q.Base (g ≫ wm) = u ≫ Q.Base (f ≫ wm) := by
      rw [Q.Base_comp, Q.Base_comp, hcond, Category.assoc]
    obtain ⟨h, hh⟩ := (hq T).2 ⟨(g ≫ wm, u), hcond'⟩
    have hp := Subtype.ext_iff.mp hh
    have h1 : (h ≫ f ≫ wm) = g ≫ wm := congrArg Prod.fst hp
    have h2 : Q.Base h = u := congrArg Prod.snd hp
    refine ⟨h, Subtype.ext (Prod.ext ?_ h2)⟩
    refine (hw T).1 (Subtype.ext (Prod.ext ?_ ?_))
    · show ((h ≫ f) ≫ wm) = g ≫ wm
      rw [Category.assoc]; exact h1
    · show Q.Base (h ≫ f) = Q.Base g
      rw [Q.Base_comp, h2, hcond]

include F in
/-- ★**(iv)(a) の移送に「pull-back 因子の射影がもとの pull-back である」ことを足した形**。

★これが (iii) の残っていた向き(`cfp_isPullBack_to`)の入口になる。 -/
theorem cfp_arbFactor' {A B : CfpCat P G} (φ : A ⟶ B) :
    ∃ (X Y : CfpCat P G) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') γ ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α ∧
        IsPullBack P (CfpCat.fst α) := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  obtain ⟨X₀, Y₀, γ₀, β₀, α₀, hfac, hγ, hβ, hα⟩ := F.arbFactor (CfpCat.fst φ)
  haveI hγb : IsIso (P.proj.map γ₀) := hγ.2
  haveI hβb : IsIso (P.proj.map β₀) := hβ.2
  have hxi : IsIso (inv (P.proj.map γ₀) ≫ A.obj.hom) := inferInstance
  have hyi : IsIso (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom) := inferInstance
  have hfacp : P.proj.map (CfpCat.fst φ)
      = P.proj.map γ₀ ≫ P.proj.map β₀ ≫ P.proj.map α₀ := by
    rw [hfac, P.proj.map_comp, P.proj.map_comp]
  refine ⟨⟨⟨X₀, A.obj.right, inv (P.proj.map γ₀) ≫ A.obj.hom⟩, hxi⟩,
    ⟨⟨Y₀, A.obj.right, inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom⟩, hyi⟩,
    InducedCategory.homMk ⟨γ₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : A.obj.right ⟶ B.obj.right), ?_⟩,
    ?_, ?_, ?_, ?_, hα⟩
  · show P.proj.map α₀ ≫ B.obj.hom
      = (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom) ≫ G.map (CfpCat.snd φ)
    rw [Category.assoc, Category.assoc, ← cfp_square φ, hfacp]
    simp only [Category.assoc]
    rw [← Category.assoc (inv (P.proj.map γ₀)), IsIso.inv_hom_id, Category.id_comp,
      ← Category.assoc (inv (P.proj.map β₀)), IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
    show CfpCat.snd φ = 𝟙 _ ≫ 𝟙 _ ≫ CfpCat.snd φ
    simp
  · exact (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hγ
  · exact ⟨hβ.1, by show IsIso (𝟙 _); infer_instance⟩
  · exact cfp_isPullBack_of P G hG hD' hcC hcD' _ hα

include F in
/-- ★★★★**(iii) の残っていた向き** —— **`𝒞'` の pull-back は `𝒞` の pull-back**。

原文 (FrdI p.28):
> (iii) A morphism of C is a(n) isometry (respectively, morphism of a given

★★**筋**: 任意射の 3 分解 `φ = γ ≫ β ≫ α`(`cfp_arbFactor'`)で
`γ ≫ β` は base-isomorphism、`α` は pull-back。
`φ` が pull-back なら **pull-back は左から簡約できる**(`isPullBack_of_comp_left`)ので
`γ ≫ β` も pull-back になり、**底が同型な pull-back は同型**
(`isIso_of_isPullBack_of_baseIso`)だから `γ ≫ β` は同型。
★あとは `fst φ = fst (γ ≫ β) ≫ fst α`(前者は同型、後者は `𝒞` の pull-back)である。

★★この向きが無かったために `pullBackLB` と `arbFactorUniq` の移送が止まっていた。 -/
theorem cfp_isPullBack_to {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') φ) :
    IsPullBack P (CfpCat.fst φ) := by
  obtain ⟨Z, W, γ, β, α, hfac, hγ, hβ, hα, hα₀⟩ :=
    cfp_arbFactor' P G hG hD' hcC hcD' F φ
  haveI hsγ : IsIso (CfpCat.snd γ) := hγ.2
  haveI hsβ : IsIso (CfpCat.snd β) := hβ.2
  have hbi : IsIso ((cfpPreFrobenioid P G hG hD' hcC hcD').Base (γ ≫ β)) := by
    show IsIso (CfpCat.snd γ ≫ CfpCat.snd β)
    infer_instance
  have hpb : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') ((γ ≫ β) ≫ α) := by
    rw [Category.assoc, ← hfac]; exact h
  have h1 : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') (γ ≫ β) :=
    isPullBack_of_comp_left _ (γ ≫ β) α hα hpb
  haveI : IsIso (γ ≫ β) :=
    isIso_of_isPullBack_of_baseIso (cfpPreFrobenioid P G hG hD' hcC hcD') h1 hbi
  haveI : IsIso (CfpCat.fst (γ ≫ β)) := cfp_isIso_fst P G _ inferInstance
  have hfst : CfpCat.fst φ = CfpCat.fst (γ ≫ β) ≫ CfpCat.fst α := by
    rw [hfac, ← Category.assoc]; rfl
  rw [hfst]
  exact IsPullBack.comp P (isPullBack_of_isIso P _) hα₀

/-- **(iii)** —— pull-back は射影で決まる(両向き)。 -/
theorem cfp_isPullBack_iff (F : FrobenioidCore P) {X Y : CfpCat P G} (φ : X ⟶ Y) :
    IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') φ ↔ IsPullBack P (CfpCat.fst φ) :=
  ⟨cfp_isPullBack_to P G hG hD' hcC hcD' F φ, cfp_isPullBack_of P G hG hD' hcC hcD' φ⟩

include F in
/-- **(iv)(b)** の移送 —— pull-back は LB-invertible かつ linear。 -/
theorem cfp_pullBackLB {A B : CfpCat P G} (α : A ⟶ B)
    (h : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α) :
    IsLBInvertible (cfpPreFrobenioid P G hG hD' hcC hcD') α ∧
      IsLinear (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  obtain ⟨h1, h2⟩ := F.pullBackLB (CfpCat.fst α)
    (cfp_isPullBack_to P G hG hD' hcC hcD' F α h)
  exact ⟨(cfp_lbInvertible_iff P G hG hD' hcC hcD' α).mpr h1,
    (cfp_linear_iff P G hG hD' hcC hcD' α).mpr h2⟩

include F in
/-- **(iv)(a)** の一意性の移送。

★`𝒟'` 成分は `γ`(Frobenius 型)と `β`(pre-step)の `𝒟'` 成分がどちらも同型なので、
**逆射を成分ごとに書き下せる**(分類表 #2 の回避)。 -/
theorem cfp_arbFactorUniq {A B : CfpCat P G} (X Y X' Y' : CfpCat P G)
    (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B)
    (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B)
    (heq : γ ≫ β ≫ α = γ' ≫ β' ≫ α')
    (hγ : IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') γ)
    (hβ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hα : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hγ' : IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') γ')
    (hβ' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hα' : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α') :
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hX : IsIso X.obj.hom := X.property
  haveI hX' : IsIso X'.obj.hom := X'.property
  haveI hY : IsIso Y.obj.hom := Y.property
  haveI hY' : IsIso Y'.obj.hom := Y'.property
  haveI hsγ : IsIso (CfpCat.snd γ) := hγ.2
  haveI hsγ' : IsIso (CfpCat.snd γ') := hγ'.2
  haveI hsβ : IsIso (CfpCat.snd β) := hβ.2
  haveI hsβ' : IsIso (CfpCat.snd β') := hβ'.2
  have hfst : CfpCat.fst γ ≫ CfpCat.fst β ≫ CfpCat.fst α
      = CfpCat.fst γ' ≫ CfpCat.fst β' ≫ CfpCat.fst α' :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) heq
  have hsnd : CfpCat.snd γ ≫ CfpCat.snd β ≫ CfpCat.snd α
      = CfpCat.snd γ' ≫ CfpCat.snd β' ≫ CfpCat.snd α' :=
    congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) heq
  obtain ⟨δ₀, ε₀, e1, e2, e3⟩ := F.arbFactorUniq X.obj.left Y.obj.left X'.obj.left Y'.obj.left
    (CfpCat.fst γ) (CfpCat.fst β) (CfpCat.fst α)
    (CfpCat.fst γ') (CfpCat.fst β') (CfpCat.fst α') hfst
    ((cfp_frobType_iff P G hG hD' hcC hcD' γ hγ.2).mp hγ)
    ((cfp_preStep_iff P G hG hD' hcC hcD' β hβ.2).mp hβ)
    (cfp_isPullBack_to P G hG hD' hcC hcD' F α hα)
    ((cfp_frobType_iff P G hG hD' hcC hcD' γ' hγ'.2).mp hγ')
    ((cfp_preStep_iff P G hG hD' hcC hcD' β' hβ'.2).mp hβ')
    (cfp_isPullBack_to P G hG hD' hcC hcD' F α' hα')
  haveI hpγ : IsIso (P.proj.map (CfpCat.fst γ)) :=
    ((cfp_frobType_iff P G hG hD' hcC hcD' γ hγ.2).mp hγ).2
  haveI hpγ' : IsIso (P.proj.map (CfpCat.fst γ')) :=
    ((cfp_frobType_iff P G hG hD' hcC hcD' γ' hγ'.2).mp hγ').2
  haveI hpβ : IsIso (P.proj.map (CfpCat.fst β)) :=
    ((cfp_preStep_iff P G hG hD' hcC hcD' β hβ.2).mp hβ).2
  haveI hpβ' : IsIso (P.proj.map (CfpCat.fst β')) :=
    ((cfp_preStep_iff P G hG hD' hcC hcD' β' hβ'.2).mp hβ').2
  have e3' : CfpCat.fst γ = CfpCat.fst γ' ≫ ε₀.inv := by
    rw [e3, Category.assoc, ε₀.hom_inv_id, Category.comp_id]
  have e2' : CfpCat.fst β ≫ δ₀.hom = ε₀.hom ≫ CfpCat.fst β' := by
    rw [e2, ← Category.assoc, ε₀.hom_inv_id, Category.id_comp]
  have e2'' : CfpCat.fst β' ≫ δ₀.inv = ε₀.inv ≫ CfpCat.fst β := by
    rw [e2, Category.assoc, Category.assoc, δ₀.hom_inv_id, Category.comp_id]
  have hsqε : P.proj.map ε₀.hom ≫ X'.obj.hom
      = X.obj.hom ≫ G.map (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') := by
    refine (cancel_epi (P.proj.map (CfpCat.fst γ))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, ← e3, cfp_square γ', ← Category.assoc,
      cfp_square γ, Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  have hsqε' : P.proj.map ε₀.inv ≫ X.obj.hom
      = X'.obj.hom ≫ G.map (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst γ'))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, ← e3', cfp_square γ, ← Category.assoc,
      cfp_square γ', Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  have hsqδ : P.proj.map δ₀.hom ≫ Y'.obj.hom
      = Y.obj.hom ≫ G.map (inv (CfpCat.snd β) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ')
          ≫ CfpCat.snd β') := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, e2', P.proj.map_comp, Category.assoc,
      cfp_square β', ← Category.assoc, hsqε, ← Category.assoc, cfp_square β]
    simp only [Category.assoc, ← G.map_comp]
    rw [← Category.assoc (CfpCat.snd β), IsIso.hom_inv_id, Category.id_comp]
  have hsqδ' : P.proj.map δ₀.inv ≫ Y.obj.hom
      = Y'.obj.hom ≫ G.map (inv (CfpCat.snd β') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ)
          ≫ CfpCat.snd β) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β'))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, e2'', P.proj.map_comp, Category.assoc,
      cfp_square β, ← Category.assoc, hsqε', ← Category.assoc, cfp_square β']
    simp only [Category.assoc, ← G.map_comp]
    rw [← Category.assoc (CfpCat.snd β'), IsIso.hom_inv_id, Category.id_comp]
  have hsα : CfpCat.snd α'
      = (inv (CfpCat.snd β') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β)
        ≫ CfpCat.snd α := by
    simp only [Category.assoc]
    rw [hsnd]
    simp
  refine ⟨⟨InducedCategory.homMk ⟨δ₀.hom, inv (CfpCat.snd β) ≫
      (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ CfpCat.snd β', hsqδ⟩,
    InducedCategory.homMk ⟨δ₀.inv, inv (CfpCat.snd β') ≫
      (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β, hsqδ'⟩, ?_, ?_⟩,
    ⟨InducedCategory.homMk ⟨ε₀.hom, inv (CfpCat.snd γ) ≫ CfpCat.snd γ', hsqε⟩,
    InducedCategory.homMk ⟨ε₀.inv, inv (CfpCat.snd γ') ≫ CfpCat.snd γ, hsqε'⟩, ?_, ?_⟩,
    ?_, ?_, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext δ₀.hom_inv_id ?_)
    show (inv (CfpCat.snd β) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ CfpCat.snd β')
      ≫ (inv (CfpCat.snd β') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β) = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext δ₀.inv_hom_id ?_)
    show (inv (CfpCat.snd β') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β)
      ≫ (inv (CfpCat.snd β) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ CfpCat.snd β') = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext ε₀.hom_inv_id ?_)
    show (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext ε₀.inv_hom_id ?_)
    show (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') = 𝟙 _
    simp
  · exact InducedCategory.hom_ext (CommaMorphism.ext e1 hsα)
  · refine InducedCategory.hom_ext (CommaMorphism.ext e2 ?_)
    show CfpCat.snd β' = (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β
      ≫ (inv (CfpCat.snd β) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ CfpCat.snd β')
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext e3 ?_)
    show CfpCat.snd γ' = CfpCat.snd γ ≫ inv (CfpCat.snd γ) ≫ CfpCat.snd γ'
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

/-! ### ★`plBkEquiv` —— 数学は完了、Lean は 1 箇所で止まっている(再現手順つき)

★★**数学**: 忠実性・充満性は **`𝒞'` の pull-back 性(対象の定義そのもの)を直接使う**。
充満性で要る「pull-back の左簡約」は **`isPullBack_of_comp_left`(仮定なし)** で証明済みなので
**循環しない**。本質的全射性は `𝒞` の `plBkEquiv` の本質的全射性 ＋ `cfp_isPullBack_of`
(構成の向き)だけで足りる。★**pull-back の「`𝒞' ⟹ 𝒞`」向きは要らない。**

★**Lean で止まる唯一の箇所**(本質的全射性の中):
```
hw  : e ≫ (G.map Y.hom ≫ inv A.obj.hom) = P.proj.map Z'.hom.hom
⊢    P.proj.map Z'.hom.hom ≫ A.obj.hom = e ≫ G.map Y.hom
```
`rw [← hw, Category.assoc, Category.assoc]` の後、目標は表示上
`e ≫ G.map Y.hom ≫ inv A.obj.hom ≫ A.obj.hom = e ≫ G.map Y.hom` になるが、
★**`rw [IsIso.inv_hom_id]` が「`inv ?f ≫ ?f` が見つからない」と言う。**

★試して駄目だったもの(6 通り): `simp` / `simp only [Category.assoc, hwA2]` /
`Category.assoc` を引数明示で 2 回 / `calc` で括弧を明示 /
`wA`(`hA.out` から取った逆射)版 / `@inv _ _ _ _ A.obj.hom hA` 版。
★**原因は特定できていない。** 分類表 #1 の一種と見ているが、
**「症状ではなく原因を特定する」を実行できていない**ので、そのまま記録する。

★★**親による原因の候補(2026-08-15。★これは答えではなく「判別法」である)**

第12段で我々は「もっともらしい原因を思いついても、それが原因であることは
別に確かめる必要がある」を学んだ(C1 が偽だった)。だから候補として書く。

**候補A —— `inv` のインスタンス引数違い**
  `inv f` は `[IsIso f]` を暗黙に取る。目標中の `inv A.obj.hom` が担いでいる
  インスタンスと、`rw [IsIso.inv_hom_id]` が単一化で合成するインスタンスが
  **別物**なら、項として異なるので構文的な `rw` は照合に失敗する。
  `A.obj.hom` の `IsIso` には `A.property`(`FullSubcategory` の `Prop` フィールド)由来と
  `haveI hA` 由来の 2 つがありうる。
  ★**第12段で偽と判定した C1 が、別の場所で真になっている**可能性である。

**候補B —— `rw [Category.assoc]` の適用位置**
  `rw` は**最初の出現だけ**を書き換える。2 回当てても狙った位置に行くとは限らない。

★**判別法**: `simp only [Category.assoc]` **だけ**を当てて右結合に正規化してから
`rw [IsIso.inv_hom_id]` を試す。
  - **通れば B**(位置の問題だった)
  - **なお落ちれば A**(インスタンスの問題)。そのとき `set_option pp.all true` で
    `inv` の第4引数を目視すれば確定する。

★**未試行の手が 2 つある**(子の 6 通りに含まれない):
  1. `simp only [Category.assoc, IsIso.inv_hom_id, Category.comp_id]`
     —— 子が試したのは `hwA2` 版で、`IsIso.inv_hom_id` を simp 補題として
     入れた版は試していない。
  2. ★**`slice_lhs 3 4 => rw [IsIso.inv_hom_id]`**
     —— `Mathlib/Tactic/CategoryTheory/Slice.lean` にある。
     ★**位置で指定するので、結合にも `rw` の適用順にも依存しない。**
     候補 A・B のどちらであっても B なら確実に抜ける。
-/

/-! ### ★(参考) `plBkEquiv` の構造

★★**数学は片付きました**: 忠実性と充満性は **`𝒞'` の pull-back 性を直接使う**だけでよく、
充満性で要る「pull-back の左簡約」は **`isPullBack_of_comp_left`(上)で証明済み** ——
`Definition 1.2, (ii)` の全単射条件だけから出るので**循環しません**。
本質的全射性も `cfp_isPullBack_of`(構成の向き)だけで足ります。

★Lean では本質的全射性の中の `wA ≫ A.obj.hom` の書き換えが通らず(分類表 #1 の一種)、
**規模を超えたのでここで切りました**。★**pull-back の「`𝒞'` ⟹ `𝒞`」向きは要らない**
という結論は変わりません。
-/

end Core

end Dict

end ABC3.Found.FrdI
