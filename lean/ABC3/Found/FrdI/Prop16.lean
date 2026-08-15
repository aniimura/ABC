import ABC3.Found.FrdI.Prop19

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
    (hD' : IsTotallyEpimorphic D') : PreFrobenioid (CfpCat P G) (Φ.restrict G hG) where
  toElem := cfpToElem P G hG
  divisorial A := P.divisorial (G.obj A)
  totEpiC := cfp_totEpi P G hD'
  totEpiD := hD'

/-! ## ★辞書 —— `𝒞'` の `Base` / `Div` / `deg_Fr`

★`Istr` のときの `istr_compat_*` に当たるもの。ただし **`Div` だけは `rfl` ではない**
(`Φ(α⁻¹)` の付け替えが挟まる)。そこが `Istr`(充満部分圏)との違いである。 -/

section Dict

variable (P : PreFrobenioid C Φ) (G : D' ⥤ D)
  (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α))
  (hD' : IsTotallyEpimorphic D')

theorem cfp_compat_Base {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD').Base f = CfpCat.snd f := rfl

theorem cfp_compat_degFr {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD').degFr f = P.degFr (CfpCat.fst f) := rfl

theorem cfp_compat_Div {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD').Div f
      = Φ.map (@inv _ _ _ _ X.obj.hom X.property) (P.Div (CfpCat.fst f)) := rfl

/-- **(iii)** —— **Frobenius 次数**は射影で決まる。 -/
theorem cfp_degFr_eq {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD').degFr f = P.degFr (CfpCat.fst f) := rfl

/-- **(iv)** —— **base-isomorphism** は `𝒟'` 成分が同型であること。

★★**これが CFP の移送を支える一点である** —— `𝒟'` 成分の情報は
base-isomorphism の定義そのものに入っている。 -/
theorem cfp_baseIso_iff {X Y : CfpCat P G} (f : X ⟶ Y) :
    IsBaseIsomorphism (cfpPreFrobenioid P G hG hD') f ↔ IsIso (CfpCat.snd f) := Iff.rfl

/-- **(iii)** —— **linear** は射影で決まる。 -/
theorem cfp_linear_iff {X Y : CfpCat P G} (f : X ⟶ Y) :
    IsLinear (cfpPreFrobenioid P G hG hD') f ↔ IsLinear P (CfpCat.fst f) := Iff.rfl

/-- **(iii)** —— **isometry** は射影で決まる。

★中身は「`Φ(α⁻¹)` が単射」の一点(`Definition 1.1, (ii), (a)`)。 -/
theorem cfp_isometric_iff {X Y : CfpCat P G} (f : X ⟶ Y) :
    IsIsometric (cfpPreFrobenioid P G hG hD') f ↔ IsIsometric P (CfpCat.fst f) := by
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
    (h : IsPullBack P (CfpCat.fst φ)) : IsPullBack (cfpPreFrobenioid P G hG hD') φ := by
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
    (h : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD') f) :
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
    IsCoAngular (cfpPreFrobenioid P G hG hD') φ := by
  intro Z W γ β α hfac hαl hβi hβs hdisj
  have hfac' : CfpCat.fst φ = CfpCat.fst γ ≫ CfpCat.fst β ≫ CfpCat.fst α :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hfac
  have hβsC : IsPreStep P (CfpCat.fst β) :=
    ⟨hβs.1, cfp_baseIso_fst P G hG hD' β hβs.2⟩
  have hdisjC : IsBaseIsomorphism P (CfpCat.fst α) ∨ IsBaseIsomorphism P (CfpCat.fst γ) :=
    hdisj.imp (cfp_baseIso_fst P G hG hD' α) (cfp_baseIso_fst P G hG hD' γ)
  have := h _ _ (CfpCat.fst γ) (CfpCat.fst β) (CfpCat.fst α) hfac' hαl
    ((cfp_isometric_iff P G hG hD' β).mp hβi) hβsC hdisjC
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

end Dict

end ABC3.Found.FrdI
