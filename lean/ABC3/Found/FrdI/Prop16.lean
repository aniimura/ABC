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
    (h : IsCoAngular (cfpPreFrobenioid P G hG hD') φ) : IsCoAngular P (CfpCat.fst φ) := by
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
      hα₀l ((cfp_isometric_iff P G hG hD' _).mpr hβ₀i) ⟨hβ₀s.1, by show IsIso (𝟙 _); infer_instance⟩
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
      hα₀l ((cfp_isometric_iff P G hG hD' _).mpr hβ₀i) ⟨hβ₀s.1, by show IsIso (𝟙 _); infer_instance⟩
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
    IsCoAngular (cfpPreFrobenioid P G hG hD') φ ↔ IsCoAngular P (CfpCat.fst φ) :=
  ⟨cfp_coAngular_to P G hG hD' φ, cfp_coAngular_of P G hG hD' φ⟩

/-- **(iii)** —— **LB-invertible** は射影で決まる。

★`LB-invertible = co-angular ∧ isometric` なので、上の2つの合成。 -/
theorem cfp_lbInvertible_iff {X Y : CfpCat P G} (φ : X ⟶ Y) :
    IsLBInvertible (cfpPreFrobenioid P G hG hD') φ ↔ IsLBInvertible P (CfpCat.fst φ) :=
  and_congr (cfp_coAngular_iff P G hG hD' φ) (cfp_isometric_iff P G hG hD' φ)

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
    (hb : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD') f) :
    IsPreStep (cfpPreFrobenioid P G hG hD') f ↔ IsPreStep P (CfpCat.fst f) :=
  ⟨fun h => ⟨h.1, cfp_baseIso_fst P G hG hD' f hb⟩, fun h => ⟨h.1, hb⟩⟩

/-- **(iv)** —— base-isomorphism については **Frobenius 型**が射影で決まる。 -/
theorem cfp_frobType_iff {X Y : CfpCat P G} (f : X ⟶ Y)
    (hb : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD') f) :
    IsFrobeniusType (cfpPreFrobenioid P G hG hD') f ↔ IsFrobeniusType P (CfpCat.fst f) :=
  ⟨fun h => ⟨(cfp_lbInvertible_iff P G hG hD' f).mp h.1, cfp_baseIso_fst P G hG hD' f hb⟩,
   fun h => ⟨(cfp_lbInvertible_iff P G hG hD' f).mpr h.1, hb⟩⟩

/-- **(iv)** —— base-isomorphism については **step** が射影で決まる。

★同型性の両向きは `cfp_isIso_of` / `cfp_isIso_fst`。 -/
theorem cfp_step_iff {X Y : CfpCat P G} (f : X ⟶ Y)
    (hb : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD') f) :
    IsStep (cfpPreFrobenioid P G hG hD') f ↔ IsStep P (CfpCat.fst f) := by
  constructor
  · rintro ⟨hs, hni⟩
    exact ⟨(cfp_preStep_iff P G hG hD' f hb).mp hs,
      fun hi => hni (cfp_isIso_of P G f hi hb)⟩
  · rintro ⟨hs, hni⟩
    exact ⟨(cfp_preStep_iff P G hG hD' f hb).mpr hs,
      fun hi => hni (cfp_isIso_fst P G f hi)⟩

/-- ★`𝒞'` の base-identity 自己射は、`𝒟'` 成分が `𝟙` であること。 -/
theorem cfp_baseIdentity_iff {A : CfpCat P G} (e : A ⟶ A) :
    IsBaseIdentity (cfpPreFrobenioid P G hG hD') e ↔ CfpCat.snd e = 𝟙 A.obj.right :=
  Iff.rfl

/-- ★その `𝒞` 成分は `𝒞` の base-identity。 -/
theorem cfp_baseIdentity_fst {A : CfpCat P G} (e : A ⟶ A)
    (h : IsBaseIdentity (cfpPreFrobenioid P G hG hD') e) :
    IsBaseIdentity P (CfpCat.fst e) := by
  haveI hA : IsIso A.obj.hom := A.property
  have hsq := cfp_square e
  rw [(cfp_baseIdentity_iff P G hG hD' e).mp h, G.map_id, Category.comp_id] at hsq
  show P.proj.map (CfpCat.fst e) = P.proj.map (𝟙 A.obj.left)
  rw [P.proj.map_id]
  exact (cancel_mono A.obj.hom).mp (hsq.trans (Category.id_comp _).symm)

/-- **(iv)** —— 射影は **`𝒪^▷` のモノイド同型**を誘導する。

★原文の当該行(`functor C′ →C determines a bijection of monoids O▷(A′)`)は
**`′`(prime) と `▷` を含むため逐語照合できない**。書き換えず、照合できない事実として記す。

★★**基準の5例目**: `base-identity` の定義が **`𝒟'` 成分を `𝟙` に固定する**ので、
`𝒪^▷` の元は `𝒞` 成分だけで決まる。**鎖どころか、`𝒟'` 成分が一意に定まる。** -/
def cfpOTriEquiv (A : CfpCat P G) :
    OTri (cfpPreFrobenioid P G hG hD') A ≃* OTri P A.obj.left where
  toFun e := ⟨CfpCat.fst (e : End A), cfp_baseIdentity_fst P G hG hD' _ e.2.1, e.2.2⟩
  invFun e :=
    ⟨InducedCategory.homMk ⟨(e : End A.obj.left), 𝟙 _, by
      have h1 : P.proj.map (e : End A.obj.left) = 𝟙 _ := e.2.1.trans (P.Base_id _)
      rw [h1, G.map_id, Category.comp_id, Category.id_comp]⟩,
     by show CfpCat.snd _ = 𝟙 _; rfl, e.2.2⟩
  left_inv e := by
    refine Subtype.ext (InducedCategory.hom_ext (CommaMorphism.ext rfl ?_))
    exact ((cfp_baseIdentity_iff P G hG hD' (e : End A)).mp e.2.1).symm
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
    IsIsotropic (cfpPreFrobenioid P G hG hD') A ↔ IsIsotropic P A.obj.left := by
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
      ((cfp_isometric_iff P G hG hD' _).mpr hi) ⟨hs.1, by show IsIso (𝟙 _); infer_instance⟩
  · intro h Dd' φ hi hs
    exact cfp_isIso_of P G φ
      (h Dd'.obj.left (CfpCat.fst φ) ((cfp_isometric_iff P G hG hD' φ).mp hi)
        ((cfp_preStep_iff P G hG hD' φ hs.2).mp hs)) hs.2

/-- **(v)** —— **Frobenius-isotropic** は射影で決まる(★(a) 鎖型)。 -/
theorem cfp_frobIsotropic_iff (A : CfpCat P G) :
    IsFrobeniusIsotropic (cfpPreFrobenioid P G hG hD') A ↔
      IsFrobeniusIsotropic P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · rintro ⟨Dd', φ, hft, hiso⟩
    exact ⟨Dd'.obj.left, CfpCat.fst φ, (cfp_frobType_iff P G hG hD' φ hft.2).mp hft,
      (cfp_isotropic_iff P G hG hD' Dd').mp hiso⟩
  · rintro ⟨Dd₀, φ₀, hft, hiso⟩
    haveI hφb : IsIso (P.proj.map φ₀) := hft.2
    have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
    refine ⟨⟨⟨Dd₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩,
      InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩, ?_, ?_⟩
    · exact (cfp_frobType_iff P G hG hD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hft
    · exact (cfp_isotropic_iff P G hG hD' _).mpr hiso

/-- **(v)** —— **quasi-Frobenius-trivial** は射影で決まる(★(b) `𝒟'` 成分固定型)。 -/
theorem cfp_quasiFrobTrivial_iff (A : CfpCat P G) :
    IsQuasiFrobeniusTrivial (cfpPreFrobenioid P G hG hD') A ↔
      IsQuasiFrobeniusTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · intro h n
    obtain ⟨φ, hbi, hdeg⟩ := h n
    exact ⟨CfpCat.fst φ, cfp_baseIdentity_fst P G hG hD' φ hbi, hdeg⟩
  · intro h n
    obtain ⟨φ₀, hbi, hdeg⟩ := h n
    have hid : P.proj.map φ₀ = 𝟙 _ := hbi.trans (P.Base_id _)
    exact ⟨InducedCategory.homMk ⟨φ₀, 𝟙 _, by rw [hid, G.map_id, Category.comp_id,
      Category.id_comp]⟩, by show CfpCat.snd _ = 𝟙 _; rfl, hdeg⟩

/-- **(v)** —— **Frobenius-trivial** は射影で決まる(★(b) 型)。

★`ζ : ℕ≥1 →* End A` を運ぶとき、**`𝒟'` 成分をすべて `𝟙` に取る**ので
モノイド準同型であることが自動になる。 -/
theorem cfp_frobTrivial_iff (A : CfpCat P G) :
    IsFrobeniusTrivial (cfpPreFrobenioid P G hG hD') A ↔ IsFrobeniusTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · rintro ⟨ζ, hdeg, hprop⟩
    refine ⟨⟨⟨fun n => CfpCat.fst (ζ n), ?_⟩, ?_⟩, hdeg, fun n =>
      ⟨cfp_baseIdentity_fst P G hG hD' _ (hprop n).1,
       (cfp_frobType_iff P G hG hD' _ (hprop n).2.2).mp (hprop n).2⟩⟩
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
       (cfp_frobType_iff P G hG hD' _ (by show IsIso (𝟙 _); infer_instance)).mpr (hprop n).2⟩⟩
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
    IsSubQuasiFrobeniusTrivial (cfpPreFrobenioid P G hG hD') A ↔
      IsSubQuasiFrobeniusTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · rintro ⟨Dd', α, hca, hps, hq⟩
    exact ⟨Dd'.obj.left, CfpCat.fst α, (cfp_coAngular_iff P G hG hD' α).mp hca,
      (cfp_preStep_iff P G hG hD' α hps.2).mp hps,
      (cfp_quasiFrobTrivial_iff P G hG hD' Dd').mp hq⟩
  · rintro ⟨Dd₀, α₀, hca, hps, hq⟩
    haveI hαb : IsIso (P.proj.map α₀) := hps.2
    have hzi : IsIso (P.proj.map α₀ ≫ A.obj.hom) := inferInstance
    refine ⟨⟨⟨Dd₀, A.obj.right, P.proj.map α₀ ≫ A.obj.hom⟩, hzi⟩,
      InducedCategory.homMk ⟨α₀, 𝟙 _, by simp⟩, ?_, ⟨hps.1, by show IsIso (𝟙 _); infer_instance⟩,
      ?_⟩
    · exact cfp_coAngular_of P G hG hD' _ hca
    · exact (cfp_quasiFrobTrivial_iff P G hG hD' _).mpr hq

/-- **(v)** —— **unit-trivial** は射影で決まる(★(b) `𝒟'` 成分固定型)。 -/
theorem cfp_unitTrivial_iff (A : CfpCat P G) :
    IsUnitTrivial (cfpPreFrobenioid P G hG hD') A ↔ IsUnitTrivial P A.obj.left := by
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
        ∈ OTimes (cfpPreFrobenioid P G hG hD') A := by
      refine ⟨⟨by show CfpCat.snd _ = 𝟙 _; rfl, hx₀.1.2⟩, (isUnit_iff_isIso _).mpr ?_⟩
      exact cfp_isIso_of P G _ hxi (by show IsIso (𝟙 _); infer_instance)
    have := (Submonoid.eq_bot_iff_forall _).mp h _ hmem
    exact congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) this
  · intro h
    refine Submonoid.eq_bot_iff_forall _ |>.mpr ?_
    intro x hx
    have hfst : CfpCat.fst (x : End A) ∈ OTimes P A.obj.left := by
      refine ⟨⟨cfp_baseIdentity_fst P G hG hD' _ hx.1.1, hx.1.2⟩, (isUnit_iff_isIso _).mpr ?_⟩
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
    IsGroupLikeObj (cfpPreFrobenioid P G hG hD') A ↔ IsGroupLikeObj P A.obj.left := by
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
    IsEndAmple (cfpPreFrobenioid P G hG hD') A := by
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
    IsAutAmple (cfpPreFrobenioid P G hG hD') A := by
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
    IsCoAngular (cfpPreFrobenioid P G hG hD') ψ →
      IsCoAngular (cfpPreFrobenioid P G hG hD') φ →
      IsCoAngular (cfpPreFrobenioid P G hG hD') (ψ ≫ φ) := by
  intro hψ hφ
  refine cfp_coAngular_of P G hG hD' _ ?_
  exact F.coAngularComp (CfpCat.fst ψ) (CfpCat.fst φ)
    ((cfp_coAngular_iff P G hG hD' ψ).mp hψ) ((cfp_coAngular_iff P G hG hD' φ).mp hφ)

include F in
/-- **(iii)(b)** の移送。 -/
theorem cfp_coAngularOfPreStep {X Y : CfpCat P G} (α : X ⟶ Y)
    (hca : IsCoAngular (cfpPreFrobenioid P G hG hD') α)
    (hps : IsPreStep (cfpPreFrobenioid P G hG hD') α)
    (φ : X ⟶ Y) : IsCoAngular (cfpPreFrobenioid P G hG hD') φ :=
  cfp_coAngular_of P G hG hD' φ
    (F.coAngularOfPreStep (CfpCat.fst α) ((cfp_coAngular_iff P G hG hD' α).mp hca)
      ((cfp_preStep_iff P G hG hD' α hps.2).mp hps) (CfpCat.fst φ))

include F in
/-- **(v)(a)** の移送 —— pre-step は mono。

★**両成分がそれぞれ mono** であればよい: `𝒞` 側は `F.preStepMono`、
`𝒟'` 側は **pre-step の定義から `snd` が同型**。 -/
theorem cfp_preStepMono {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hφ : IsPreStep (cfpPreFrobenioid P G hG hD') φ) : Mono φ := by
  haveI hm : Mono (CfpCat.fst φ) :=
    F.preStepMono (CfpCat.fst φ) ((cfp_preStep_iff P G hG hD' φ hφ.2).mp hφ)
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
    (h : IsIsotropic (cfpPreFrobenioid P G hG hD') X) :
    IsIsotropic (cfpPreFrobenioid P G hG hD') Y :=
  (cfp_isotropic_iff P G hG hD' Y).mpr
    (F.isotropicClosed (CfpCat.fst φ) ((cfp_isotropic_iff P G hG hD' X).mp h))

include F in
/-- **(ii)** の移送 —— 各次数の Frobenius 型射が存在する。

★Frobenius 型は base-isomorphism なので**鎖**があり、
新しい対象 `B` の `𝒟'` 成分は `A` のものを流用して `snd φ = 𝟙` に取れる。 -/
theorem cfp_frobDegSurj (A : CfpCat P G) (n : ℕ+) :
    ∃ (B : CfpCat P G) (φ : A ⟶ B),
      IsFrobeniusType (cfpPreFrobenioid P G hG hD') φ ∧
        (cfpPreFrobenioid P G hG hD').degFr φ = n := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨B₀, φ₀, hft, hdeg⟩ := F.frobDegSurj A.obj.left n
  haveI hφb : IsIso (P.proj.map φ₀) := hft.2
  have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
  refine ⟨⟨⟨B₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩,
    InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩, ?_, hdeg⟩
  exact (cfp_frobType_iff P G hG hD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hft

include F in
/-- **(v)(b)** の移送 —— pre-step は「co-angular pre-step ≫ isometric pre-step」に分解する。

★中間対象は**両側とも pre-step に挟まれる**ので鎖がある。 -/
theorem cfp_preStepFactor {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hφ : IsPreStep (cfpPreFrobenioid P G hG hD') φ) :
    ∃ (Z : CfpCat P G) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsCoAngular (cfpPreFrobenioid P G hG hD') β ∧
        IsPreStep (cfpPreFrobenioid P G hG hD') β ∧
        IsIsometric (cfpPreFrobenioid P G hG hD') α ∧
        IsPreStep (cfpPreFrobenioid P G hG hD') α := by
  haveI hX : IsIso X.obj.hom := X.property
  obtain ⟨Z₀, β₀, α₀, hfac, hβc, hβs, hαi, hαs⟩ :=
    F.preStepFactor (CfpCat.fst φ) ((cfp_preStep_iff P G hG hD' φ hφ.2).mp hφ)
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
  · exact cfp_coAngular_of P G hG hD' _ hβc
  · exact (cfp_isometric_iff P G hG hD' _).mpr hαi
  · exact ⟨hαs.1, hφ.2⟩

include F in
/-- **(v)(c)** の移送 —— pre-step は「isometric pre-step ≫ co-angular pre-step」に分解する。 -/
theorem cfp_preStepFactor' {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hφ : IsPreStep (cfpPreFrobenioid P G hG hD') φ) :
    ∃ (Z : CfpCat P G) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsIsometric (cfpPreFrobenioid P G hG hD') β ∧
        IsPreStep (cfpPreFrobenioid P G hG hD') β ∧
        IsCoAngular (cfpPreFrobenioid P G hG hD') α ∧
        IsPreStep (cfpPreFrobenioid P G hG hD') α := by
  haveI hX : IsIso X.obj.hom := X.property
  obtain ⟨Z₀, β₀, α₀, hfac, hβi, hβs, hαc, hαs⟩ :=
    F.preStepFactor' (CfpCat.fst φ) ((cfp_preStep_iff P G hG hD' φ hφ.2).mp hφ)
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
  · exact (cfp_isometric_iff P G hG hD' _).mpr hβi
  · exact cfp_coAngular_of P G hG hD' _ hαc
  · exact ⟨hαs.1, hφ.2⟩

include F in
/-- **(i)(a)** の移送 —— `𝒟'` のどの対象の上にも Frobenius-trivial な対象がある。

★`𝒞` の `baseSurj` を `G Y` に当てて得た同型を、そのまま CFP の三つ組の第3成分にする。
★**新しい対象を作るのに鎖は要らない** —— 同型が入力として与えられるから。 -/
theorem cfp_baseSurj (Y : D') :
    ∃ A : CfpCat P G, IsFrobeniusTrivial (cfpPreFrobenioid P G hG hD') A ∧
      Nonempty (((cfpPreFrobenioid P G hG hD').toElem.obj A).base ≅ Y) := by
  obtain ⟨A₀, hft, ⟨e⟩⟩ := F.baseSurj (G.obj Y)
  haveI : IsIso e.hom := e.isIso_hom
  refine ⟨⟨⟨A₀, Y, e.hom⟩, inferInstanceAs (IsIso e.hom)⟩, ?_, ⟨Iso.refl _⟩⟩
  exact (cfp_frobTrivial_iff P G hG hD' _).mpr hft

include F in
/-- **(ii)** の移送 —— 同じ次数の Frobenius 型射の本質的一意性。

★★**`𝒟'` 成分は `(snd φ)⁻¹ ≫ snd ψ` に取れる** —— Frobenius 型は base-isomorphism なので
両方の `𝒟'` 成分が同型であり、**`G` の充満性は要らない**。 -/
theorem cfp_frobDegUniq (A B E : CfpCat P G) (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : IsFrobeniusType (cfpPreFrobenioid P G hG hD') φ)
    (hψ : IsFrobeniusType (cfpPreFrobenioid P G hG hD') ψ)
    (hd : (cfpPreFrobenioid P G hG hD').degFr φ = (cfpPreFrobenioid P G hG hD').degFr ψ) :
    ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  haveI hE : IsIso E.obj.hom := E.property
  haveI hsφ : IsIso (CfpCat.snd φ) := hφ.2
  haveI hsψ : IsIso (CfpCat.snd ψ) := hψ.2
  haveI hpφ : IsIso (P.proj.map (CfpCat.fst φ)) := cfp_baseIso_fst P G hG hD' φ hφ.2
  obtain ⟨β₀, hβiso, hβ⟩ := F.frobDegUniq A.obj.left B.obj.left E.obj.left
    (CfpCat.fst φ) (CfpCat.fst ψ)
    ((cfp_frobType_iff P G hG hD' φ hφ.2).mp hφ)
    ((cfp_frobType_iff P G hG hD' ψ hψ.2).mp hψ) hd
  have hsq : P.proj.map β₀ ≫ E.obj.hom
      = B.obj.hom ≫ G.map (inv (CfpCat.snd φ) ≫ CfpCat.snd ψ) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst φ))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, hβ, cfp_square ψ, ← Category.assoc,
      cfp_square φ, Category.assoc, ← G.map_comp, ← Category.assoc,
      IsIso.hom_inv_id, Category.id_comp]
  refine ⟨InducedCategory.homMk ⟨β₀, inv (CfpCat.snd φ) ≫ CfpCat.snd ψ, hsq⟩,
    cfp_isIso_of P G _ hβiso inferInstance, ?_⟩
  refine InducedCategory.hom_ext (CommaMorphism.ext hβ ?_)
  show CfpCat.snd φ ≫ inv (CfpCat.snd φ) ≫ CfpCat.snd ψ = CfpCat.snd ψ
  rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

end Core

end Dict

end ABC3.Found.FrdI
