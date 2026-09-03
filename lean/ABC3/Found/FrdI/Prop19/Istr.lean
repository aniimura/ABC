/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop18
import ABC3.Found.FrdI.Prop19.CoaPreIso

/-!
# Prop19 —— (v) の機械・保存される 11 クラス・21 条の移送

☆もとの 1 枚を**入れ子の切れ目**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)


section Istr

/-- **`𝒞^istr`** —— isotropic な対象の充満部分圏の述語。 -/
def isotropicProp : ObjectProperty C := fun A => IsIsotropic P A

/-- **`𝒞^istr`**。 -/
abbrev Istr : Type u2 := (isotropicProp P).FullSubcategory

variable (F : FrobenioidCore P)

/-- `A` の isotropic hull の終域(選択)。★`isotropicHullExists` は `∃` なので選択が要る。 -/
noncomputable def hullObj (A : C) : C := (F.isotropicHullExists A).choose

/-- `A` から選んだ isotropic hull への射。 -/
noncomputable def hullMap (A : C) : A ⟶ hullObj P F A :=
  (F.isotropicHullExists A).choose_spec.choose

theorem hullMap_spec (A : C) : IsIsotropicHull P (hullMap P F A) :=
  (F.isotropicHullExists A).choose_spec.choose_spec

/-- `𝒞^istr` の対象としての `A^istr`。 -/
noncomputable def hullIstr (A : C) : Istr P :=
  ⟨hullObj P F A, (hullMap_spec P F A).2.2.1⟩

/-- ★★**isotropic hull の普遍性が、そのまま随伴の hom 同値になる**。

`Hom_{𝒞^istr}(A^istr, Y) ≃ Hom_𝒞(A, Y)`(`Y` は isotropic)。
★`∃!` の存在部分が `right_inv`、一意性部分が `left_inv` になる。 -/
noncomputable def hullHomEquiv (A : C) (Y : Istr P) :
    (hullIstr P F A ⟶ Y) ≃ (A ⟶ (isotropicProp P).ι.obj Y) :=
  InducedCategory.homEquiv.trans
    { toFun := fun g => hullMap P F A ≫ g
      invFun := fun h => ((hullMap_spec P F A).2.2.2 Y.obj Y.property h).choose
      left_inv := fun g =>
        (((hullMap_spec P F A).2.2.2 Y.obj Y.property
          (hullMap P F A ≫ g)).choose_spec.2 g rfl).symm
      right_inv := fun h =>
        (((hullMap_spec P F A).2.2.2 Y.obj Y.property h).choose_spec.1).symm }

/-- ★★**isotropification 関手** —— `leftAdjointOfEquiv` が
**hom 同値から関手そのものを作る**。手で組み立てる必要がない。 -/
noncomputable def isotropification : C ⥤ Istr P :=
  Adjunction.leftAdjointOfEquiv (F_obj := hullIstr P F) (G := (isotropicProp P).ι)
    (e := hullHomEquiv P F)
    (fun X _ _ g h => (Category.assoc (hullMap P F X) h.hom g.hom).symm)

/-- ★★**isotropification は包含関手の左随伴**。

原文 (FrdI p.32):
> Bistr forms a left adjoint to the inclusion functor Cistr →C, through which

★**記録(2026-08-15)**: この行の包含記号は、**PDF の描画では `↪`** だが
**`pdftotext` の抽出では `→`** になる。私は一度 PDF 画像で見たとおり `↪` と写して
ゲートに落とされた。★**「PDF 目視」と「抽出テキスト」が食い違う文字がある**という
測定であり、`▷` / `′` / `≠` が**抽出側で拾えない**のとは別の種類である
(あちらは照合不能、こちらは**別の文字に化ける**)。
★照合できる側(抽出)に合わせ、食い違いを事実として書き残す。 -/
noncomputable def isotropificationAdj : isotropification P F ⊣ (isotropicProp P).ι :=
  Adjunction.adjunctionOfEquivLeft _ _

/-- ★**`𝒞^istr` は totally epimorphic** —— 充満部分圏なので `𝒞` からそのまま移る。

★★これが「移送」の**最初の実例**である: `𝒞` の性質を1行で運ぶ。 -/
theorem istr_totEpi : IsTotallyEpimorphic (Istr P) := by
  intro A B f
  refine ⟨fun {Z} g h hgh => ?_⟩
  haveI : Epi f.hom := P.totEpiC _ _ f.hom
  refine InducedCategory.hom_ext ?_
  exact (cancel_epi f.hom).mp (congrArg InducedCategory.Hom.hom hgh)

include F in
/-- ★**`𝒞` が connected なら `𝒞^istr` も connected**（2026-08-16 に追加）。

★`isotropification` が `𝒞` の zigzag を `𝒞^istr` に送り、
どの isotropic 対象 `Y` も `hullMap` によって `Y^istr` とつながる。
★★**同型である必要はない** —— 辺が 1 本あれば連結性には十分である。 -/
theorem isConnected_istr : IsConnected (Istr P) := by
  haveI := P.connectedC
  obtain ⟨A₀⟩ := (inferInstance : Nonempty C)
  refine IsConnected.of_induct (j₀ := hullIstr P F A₀) ?_
  intro p hp0 hstep Y
  have key : ∀ A : C, hullIstr P F A ∈ p :=
    induct_on_objects (J := C) {A | hullIstr P F A ∈ p} hp0
      (fun {A B} f => hstep ((isotropification P F).map f))
  exact (hstep (⟨hullMap P F Y.obj⟩ : Y ⟶ hullIstr P F Y.obj)).mpr (key Y.obj)

/-- ★★**`𝒞^istr` の pre-Frobenioid 構造** —— 原文の
「equipped with the **restriction** to `𝒞` of the given functor `𝒞 → 𝔽_Φ`」。

★★**4フィールドのうち3つが `P` のものそのまま**である。
`totEpiC` だけが1行の議論を要した。**これが「移送」の意味である。** -/
def istrPre : PreFrobenioid (Istr P) Φ where
  toElem := (isotropicProp P).ι ⋙ P.toElem
  divisorial := P.divisorial
  totEpiC := istr_totEpi P
  totEpiD := P.totEpiD
  connectedC := isConnected_istr P F
  connectedD := P.connectedD

/-- ★`isotropification` の射を、`𝒞` の素の射として取り出したもの。

★`FullSubcategory` の射は `InducedCategory.Hom` に包まれており、
その型が簡約されないので、**素の型に落としてから使う**(第2段と同じ定型)。 -/
noncomputable def istrMap {A B : C} (f : A ⟶ B) : hullObj P F A ⟶ hullObj P F B :=
  ((isotropification P F).map f).hom

/-- ★★**`𝒞^istr` の射の性質は `𝒞` のそれと一致する**(充満部分圏だから)。

原文 (FrdI p.32):
> Cistr satisfies one of these properties with respect to Cistr if and only if it does with

★原文が「compatible with the inclusion functor」と言うのはこれで、
**`rfl` である**(`istrPre` の `toElem` が `ι ⋙ P.toElem` だから)。 -/
theorem istr_compat_degFr {X Y : Istr P} (g : X ⟶ Y) :
    (istrPre P F).degFr g = P.degFr g.hom := rfl

theorem istr_compat_Base {X Y : Istr P} (g : X ⟶ Y) :
    (istrPre P F).Base g = P.Base g.hom := rfl

theorem istr_compat_Div {X Y : Istr P} (g : X ⟶ Y) :
    (istrPre P F).Div g = P.Div g.hom := rfl

/-- ★★**isotropification の定義四角形**。

`leftAdjointOfEquiv` の作る `map f` は、**`hullMap A ≫ istrMap f = f ≫ hullMap B`
を満たす唯一の射**である。原文の「the induced [i.e., by the definition of an
"isotropic hull"!] morphism `A^istr → B^istr`」がこれ。 -/
theorem isotropification_square {A B : C} (f : A ⟶ B) :
    hullMap P F A ≫ istrMap P F f = f ≫ hullMap P F B := by
  have h := (hullHomEquiv P F A (hullIstr P F B)).apply_symm_apply
    (f ≫ hullHomEquiv P F B (hullIstr P F B) (𝟙 _))
  show hullHomEquiv P F A (hullIstr P F B) ((isotropification P F).map f) = _
  rw [show (isotropification P F).map f
      = (hullHomEquiv P F A (hullIstr P F B)).symm
        (f ≫ hullHomEquiv P F B (hullIstr P F B) (𝟙 _)) from rfl, h]
  show f ≫ (hullMap P F B ≫ 𝟙 _) = f ≫ hullMap P F B
  rw [Category.comp_id]

/-! ### ★保存される 11 クラスのうち、まず 3 つの成分

★`Base` / `Div` / `deg_Fr` の保存は、上の四角形と **`Remark 1.1.1`**
(合成公式)だけから出る。原文が「[cf. Remark 1.1.1]」と書くのはこれ。 -/

/-- ★**Frobenius 次数を保つ**。四角形の `deg_Fr` 成分と、
`hullMap` が pre-step(次数 1)であることから。 -/
theorem isotropification_degFr {A B : C} (f : A ⟶ B) :
    P.degFr (istrMap P F f) = P.degFr f := by
  have h := congrArg P.degFr (isotropification_square P F f)
  rw [P.degFr_comp, P.degFr_comp, (hullMap_spec P F A).2.1.1,
    (hullMap_spec P F B).2.1.1, mul_one, one_mul] at h
  exact h

/-- ★**base-isomorphism を(両向きに)保つ**。 -/
theorem isotropification_baseIso_iff {A B : C} (f : A ⟶ B) :
    IsIso (P.Base (istrMap P F f)) ↔ IsIso (P.Base f) := by
  haveI hA : IsIso (P.Base (hullMap P F A)) := (hullMap_spec P F A).2.1.2
  haveI hB : IsIso (P.Base (hullMap P F B)) := (hullMap_spec P F B).2.1.2
  have h := congrArg P.Base (isotropification_square P F f)
  rw [P.Base_comp, P.Base_comp] at h
  constructor
  · intro hi
    haveI := hi
    have hf : P.Base f = P.Base (hullMap P F A) ≫ P.Base (istrMap P F f)
        ≫ inv (P.Base (hullMap P F B)) := by
      rw [← Category.assoc, h, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
    rw [hf]
    infer_instance
  · intro hi
    haveI := hi
    have hf : P.Base (istrMap P F f) = inv (P.Base (hullMap P F A)) ≫ P.Base f
        ≫ P.Base (hullMap P F B) := by
      rw [← h, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
    rw [hf]
    infer_instance

/-- ★**isometry を(両向きに)保つ**。四角形の `Div` 成分で、
`hullMap` が isometric なので両端の寄与が消え、`Φ.map` の単射性が残る。 -/
theorem isotropification_isometric_iff {A B : C} (f : A ⟶ B) :
    IsIsometric P (istrMap P F f) ↔ IsIsometric P f := by
  have h := congrArg P.Div (isotropification_square P F f)
  rw [P.Div_comp, P.Div_comp, (hullMap_spec P F A).1, (hullMap_spec P F B).1,
    (hullMap_spec P F B).2.1.1] at h
  simp only [smul_zero, add_zero, PNat.one_coe, one_smul, map_zero, zero_add] at h
  constructor
  · intro hi
    show P.Div f = 0
    rw [← h, show P.Div (istrMap P F f) = 0 from hi, map_zero]
  · intro hi
    show P.Div (istrMap P F f) = 0
    refine Φ.map_injective (P.Base (hullMap P F A)) ?_
    rw [h, show P.Div f = 0 from hi, map_zero]

/-- ★★**isotropic な対象の isotropic hull は同型**。

`hullMap` は isometric pre-step で、`A` が isotropic ならそれは定義から同型。

★原文の「The restriction of the isotropification functor to `𝒞^istr` is
**isomorphic to the identity functor**」の核であり、**1行**である。 -/
theorem hullMap_isIso (A : C) (hA : IsIsotropic P A) : IsIso (hullMap P F A) :=
  hA _ (hullMap P F A) (hullMap_spec P F A).1 (hullMap_spec P F A).2.1

include F in
/-- ★**isotropic な対象から出る `𝒞` の射はすべて co-angular**。

`Definition 1.3, (vii), (b)` で isotropy が伝わるので、`Proposition 1.4, (i)` が使える。
★これが「`𝒞^istr` の射が `𝒞` の意味でも co-angular」を与える —— (v) の
pull-back の保存で要る（原文の「pull-back morphisms **relative to `𝒞`**」）。 -/
theorem isCoAngular_of_isotropic_dom {A B : C} (hA : IsIsotropic P A) (f : A ⟶ B) :
    IsCoAngular P f :=
  prop_1_4_i P f (fun X g => F.isotropicClosed g hA)

/-- ★**pre-step を(両向きに)保つ** —— 次数と底の同型性の合わせ技。 -/
theorem isotropification_preStep_iff {A B : C} (f : A ⟶ B) :
    IsPreStep P (istrMap P F f) ↔ IsPreStep P f := by
  constructor
  · intro h
    refine ⟨?_, (isotropification_baseIso_iff P F f).mp h.2⟩
    show P.degFr f = 1
    rw [← isotropification_degFr P F f]
    exact h.1
  · intro h
    refine ⟨?_, (isotropification_baseIso_iff P F f).mpr h.2⟩
    show P.degFr (istrMap P F f) = 1
    rw [isotropification_degFr P F f]
    exact h.1

/-- ★**pull-back を保つ**（`𝒞` の意味で）。

原文 (FrdI p.33):
> pull-back morphisms to morphisms which are pull-back morphisms relative to C,

★`Proposition 1.4, (ii)` で pull-back = co-angular ∧ isometric ∧ linear に分解し、
* co-angular は `A^istr` が isotropic だから自動（上の補題）
* isometric は保存（`isotropification_isometric_iff`）
* linear は次数の保存

の3つを合わせる。★**`Proposition 1.7` の合成補題は要らない。** -/
theorem isotropification_pullBack {A B : C} (f : A ⟶ B) (h : IsPullBack P f) :
    IsPullBack P (istrMap P F f) := by
  obtain ⟨hlb, hlin⟩ := (prop_1_4_ii P F f).mp h
  refine (prop_1_4_ii P F _).mpr ⟨⟨?_, ?_⟩, ?_⟩
  · exact isCoAngular_of_isotropic_dom P F (hullMap_spec P F A).2.2.1 _
  · exact (isotropification_isometric_iff P F f).mpr hlb.2
  · show P.degFr (istrMap P F f) = 1
    rw [isotropification_degFr P F f]
    exact hlin

/-- ★**Frobenius 型を保つ** —— co-angular（自動）＋ isometric ＋ base-isomorphism。 -/
theorem isotropification_frobType {A B : C} (f : A ⟶ B) (h : IsFrobeniusType P f) :
    IsFrobeniusType P (istrMap P F f) :=
  ⟨⟨isCoAngular_of_isotropic_dom P F (hullMap_spec P F A).2.2.1 _,
    (isotropification_isometric_iff P F f).mpr h.1.2⟩,
   (isotropification_baseIso_iff P F f).mpr h.2⟩

theorem istr_isotropic (X : Istr P) : IsIsotropic (istrPre P F) X := by
  intro Dd φ hi hs
  haveI : IsIso φ.hom := X.property Dd.obj φ.hom hi hs
  exact ⟨InducedCategory.homMk (inv φ.hom), InducedCategory.hom_ext (by simp),
    InducedCategory.hom_ext (by simp)⟩

/-- ★★**`𝒞^istr` のすべての射は co-angular**(`Proposition 1.4, (i)`)。

★原文が保存リストで「co-angular morphisms [cf. Proposition 1.4, (i)]」と
括弧書きするのはこれ —— **`𝒞^istr` では co-angular 性が自明になる**ので、
「保存する」は言うまでもない。 -/
theorem istr_coAngular {X Y : Istr P} (g : X ⟶ Y) : IsCoAngular (istrPre P F) g :=
  prop_1_4_i (istrPre P F) g (fun Z _ => istr_isotropic P F Z)

/-! ### ★★21 条の「移送」—— まず辞書と、移送しやすい条から

★`Definition 1.3` の各条を `istrPre P F` について示す。原文は
「[from the fact that `𝒞` is a Frobenioid!]」の一言だが、
**Lean では `Istr P` と `C` の間の辞書を1本ずつ引く**必要がある。
どこまでが本当に「移送」かを条ごとに測る。 -/

include F in
/-- ★辞書: `Istr P` の Frobenius 型は `C` のそれ。

★co-angular は**両側で自動**なので、実質は isometric と base-isomorphism の移送。 -/
theorem istr_frobType_iff {X Y : Istr P} (g : X ⟶ Y) :
    IsFrobeniusType (istrPre P F) g ↔ IsFrobeniusType P g.hom :=
  ⟨fun h => ⟨⟨isCoAngular_of_isotropic_dom P F X.property g.hom, h.1.2⟩, h.2⟩,
   fun h => ⟨⟨istr_coAngular P F g, h.1.2⟩, h.2⟩⟩

include F in
/-- **(iii)(a)** の移送 —— `𝒞^istr` では co-angular が自動なので**自明**。 -/
theorem istr_coAngularComp {X Y Z : Istr P} (ψ : X ⟶ Y) (φ : Y ⟶ Z) :
    IsCoAngular (istrPre P F) ψ → IsCoAngular (istrPre P F) φ →
      IsCoAngular (istrPre P F) (ψ ≫ φ) :=
  fun _ _ => istr_coAngular P F _

include F in
/-- **(iii)(b)** の移送 —— 同上。 -/
theorem istr_coAngularOfPreStep {X Y : Istr P} (α : X ⟶ Y) :
    IsCoAngular (istrPre P F) α → IsPreStep (istrPre P F) α →
      ∀ φ : X ⟶ Y, IsCoAngular (istrPre P F) φ :=
  fun _ _ φ => istr_coAngular P F φ

include F in
/-- **(vii)(b)** の移送 —— `𝒞^istr` の対象はすべて isotropic なので**自明**。 -/
theorem istr_isotropicClosed {X Y : Istr P} (_φ : X ⟶ Y) :
    IsIsotropic (istrPre P F) X → IsIsotropic (istrPre P F) Y :=
  fun _ => istr_isotropic P F Y

include F in
/-- **(vii)(a)** の移送 —— `X` 自身が isotropic なので `𝟙_X` が isotropic hull。 -/
theorem istr_isotropicHullExists (X : Istr P) :
    ∃ (Y : Istr P) (φ : X ⟶ Y), IsIsotropicHull (istrPre P F) φ :=
  ⟨X, 𝟙 X, (istrPre P F).Div_id X, isPreStep_id _ X, istr_isotropic P F X,
    fun Cc _ γ => ⟨γ, (Category.id_comp γ).symm, fun β hβ => by
      have hg : γ = β := by simpa using hβ
      exact hg.symm⟩⟩

include F in
/-- **(v)(a)** の移送 —— `C` の mono 性が充満部分圏へそのまま降りる。 -/
theorem istr_preStepMono {X Y : Istr P} (φ : X ⟶ Y) (hφ : IsPreStep (istrPre P F) φ) :
    Mono φ := by
  haveI : Mono φ.hom := F.preStepMono φ.hom hφ
  refine ⟨fun {Z} g h hgh => ?_⟩
  refine InducedCategory.hom_ext ?_
  exact (cancel_mono φ.hom).mp (congrArg InducedCategory.Hom.hom hgh)

include F in
/-- **(ii)** の本質的一意性の移送 —— `C` で得た同型を充満部分圏へ持ち上げるだけ。 -/
theorem istr_frobDegUniq (X Y Z : Istr P) (φ : X ⟶ Y) (ψ : X ⟶ Z)
    (hφ : IsFrobeniusType (istrPre P F) φ) (hψ : IsFrobeniusType (istrPre P F) ψ)
    (hd : (istrPre P F).degFr φ = (istrPre P F).degFr ψ) :
    ∃ β : Y ⟶ Z, IsIso β ∧ φ ≫ β = ψ := by
  obtain ⟨β₀, hβiso, hβ⟩ := F.frobDegUniq X.obj Y.obj Z.obj φ.hom ψ.hom
    ((istr_frobType_iff P F φ).mp hφ) ((istr_frobType_iff P F ψ).mp hψ) hd
  haveI := hβiso
  refine ⟨InducedCategory.homMk β₀, ?_, InducedCategory.hom_ext hβ⟩
  exact ⟨InducedCategory.homMk (inv β₀), InducedCategory.hom_ext (by simp),
    InducedCategory.hom_ext (by simp)⟩

/-- ★★**pull-back の移送(易しい向き)** —— `𝒞` の pull-back は `𝒞^istr` の pull-back。

原文 (FrdI p.33):
> pull-back morphisms relative to C, hence a fortiori, pull-back morphisms relative to Cistr.

★原文の「**a fortiori**」がこれ。`Definition 1.2, (ii)` の全単射は
「すべての `Z ∈ Ob(𝒞)` について」なので、**`Z` を isotropic なものに制限すれば
そのまま成り立つ**。★ただし Lean では `Istr P` と `C` の射の対応
(`InducedCategory.homEquiv`)を挟む必要がある。

★**逆向き**(「`𝒞^istr` の pull-back は `𝒞` の pull-back」)は
`Z` が isotropic でない場合を埋める必要があり、**随伴を使う**。そちらは別に扱う。 -/
theorem istr_isPullBack_of {X Y : Istr P} (g : X ⟶ Y) (h : IsPullBack P g.hom) :
    IsPullBack (istrPre P F) g := by
  intro Z
  constructor
  · intro f₁ f₂ hf
    have hp := Subtype.ext_iff.mp hf
    have h1 : (f₁ ≫ g : Z ⟶ Y) = f₂ ≫ g := congrArg Prod.fst hp
    have h2 : P.Base f₁.hom = P.Base f₂.hom := congrArg Prod.snd hp
    refine InducedCategory.hom_ext ?_
    refine (h Z.obj).1 (Subtype.ext (Prod.ext ?_ h2))
    exact congrArg InducedCategory.Hom.hom h1
  · rintro ⟨⟨a, b⟩, hab⟩
    obtain ⟨f₀, hf₀⟩ := (h Z.obj).2 ⟨(a.hom, b), hab⟩
    have hp := Subtype.ext_iff.mp hf₀
    have h1 : (f₀ ≫ g.hom : Z.obj ⟶ Y.obj) = a.hom := congrArg Prod.fst hp
    have h2 : P.Base f₀ = b := congrArg Prod.snd hp
    exact ⟨InducedCategory.homMk f₀,
      Subtype.ext (Prod.ext (InducedCategory.hom_ext h1) h2)⟩

include F in
/-- **(iv)(a)** の移送 —— `𝒞` の3分解を `𝒞^istr` へ運ぶ。

★★**中間対象が自動で isotropic になる**のが効いている ——
`Definition 1.3, (vii), (b)`(`isotropicClosed`)により、
**isotropic な対象から出る射の終域はすべて isotropic** なので、
`𝒞` の分解に現れる対象がそのまま `𝒞^istr` の対象になる。 -/
theorem istr_arbFactor {X Y : Istr P} (φ : X ⟶ Y) :
    ∃ (Z W : Istr P) (γ : X ⟶ Z) (β : Z ⟶ W) (α : W ⟶ Y),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (istrPre P F) γ ∧
        IsPreStep (istrPre P F) β ∧ IsPullBack (istrPre P F) α := by
  obtain ⟨Z₀, W₀, γ₀, β₀, α₀, heq, hγ, hβ, hα⟩ := F.arbFactor φ.hom
  have hZ : IsIsotropic P Z₀ := F.isotropicClosed γ₀ X.property
  have hW : IsIsotropic P W₀ := F.isotropicClosed β₀ hZ
  refine ⟨⟨Z₀, hZ⟩, ⟨W₀, hW⟩, InducedCategory.homMk γ₀, InducedCategory.homMk β₀,
    InducedCategory.homMk α₀, InducedCategory.hom_ext heq, ?_, hβ,
    istr_isPullBack_of P F _ hα⟩
  exact (istr_frobType_iff P F (X := X) (Y := ⟨Z₀, hZ⟩)
    (InducedCategory.homMk γ₀)).mpr hγ

include F in
/-- **(ii)** の移送 —— 各次数の Frobenius 型射。中間対象は自動で isotropic。 -/
theorem istr_frobDegSurj (X : Istr P) (n : ℕ+) :
    ∃ (Y : Istr P) (φ : X ⟶ Y), IsFrobeniusType (istrPre P F) φ ∧
      (istrPre P F).degFr φ = n := by
  obtain ⟨B₀, φ₀, hφ, hd⟩ := F.frobDegSurj X.obj n
  exact ⟨⟨B₀, F.isotropicClosed φ₀ X.property⟩, InducedCategory.homMk φ₀,
    (istr_frobType_iff P F _).mpr hφ, hd⟩

include F in
/-- **(v)(b)** の移送。 -/
theorem istr_preStepFactor {X Y : Istr P} (φ : X ⟶ Y) (hφ : IsPreStep (istrPre P F) φ) :
    ∃ (Z : Istr P) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsCoAngular (istrPre P F) β ∧ IsPreStep (istrPre P F) β ∧
        IsIsometric (istrPre P F) α ∧ IsPreStep (istrPre P F) α := by
  obtain ⟨Z₀, β₀, α₀, heq, hβc, hβs, hαi, hαs⟩ := F.preStepFactor φ.hom hφ
  refine ⟨⟨Z₀, F.isotropicClosed β₀ X.property⟩, InducedCategory.homMk β₀,
    InducedCategory.homMk α₀, InducedCategory.hom_ext heq, ?_, hβs, hαi, hαs⟩
  exact istr_coAngular P F _

include F in
/-- **(v)(c)** の移送。 -/
theorem istr_preStepFactor' {X Y : Istr P} (φ : X ⟶ Y) (hφ : IsPreStep (istrPre P F) φ) :
    ∃ (Z : Istr P) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsIsometric (istrPre P F) β ∧ IsPreStep (istrPre P F) β ∧
        IsCoAngular (istrPre P F) α ∧ IsPreStep (istrPre P F) α := by
  obtain ⟨Z₀, β₀, α₀, heq, hβi, hβs, hαc, hαs⟩ := F.preStepFactor' φ.hom hφ
  refine ⟨⟨Z₀, F.isotropicClosed β₀ X.property⟩, InducedCategory.homMk β₀,
    InducedCategory.homMk α₀, InducedCategory.hom_ext heq, hβi, hβs, ?_, hαs⟩
  exact istr_coAngular P F _

include F in
/-- ★★**isotropic な対象への pull-back の始域は isotropic**。

`Proposition 1.4, (ii)` で pull-back は **co-angular linear**、
そこに **`Proposition 1.9, (iv)`** を当てるだけ。

★★**原文が (v) の最後で「in light of Proposition 1.4, (i); assertion (iv)」と
書く理由がこれである** —— これがあるから `(𝒞^pl-bk)_{A}` の対象が
**そのまま** `((𝒞^istr)^pl-bk)_{A}` の対象になり、
`Definition 1.3, (i), (c)` の圏同値が `𝒞^istr` へ移送できる。

★**3 行**である。 -/
theorem isotropic_dom_of_pullBack {X A : C} (p : X ⟶ A) (hp : IsPullBack P p)
    (hA : IsIsotropic P A) : IsIsotropic P X := by
  obtain ⟨hlb, hlin⟩ := (prop_1_4_ii P F p).mp hp
  exact (prop_1_9_iv P F p hlb.1 hlin).mpr hA

include F in
/-- ★★**pull-back の移送(難しい向き)** —— `𝒞^istr` の pull-back は `𝒞` の pull-back。

`Definition 1.2, (ii)` の全単射は「**すべての `Z ∈ Ob(𝒞)`**」についてだが、
`istrPre` の側は「**isotropic な `Z`**」しか言わない。この差を埋めるのが

★★**isotropification が包含関手の左随伴であること**(`hullHomEquiv`)である。

`Z` の isotropic hull を取れば `Hom_𝒞(Z, X) ≅ Hom_{𝒞^istr}(Z^istr, X)` なので、
一般の `Z` についての全単射性が isotropic な `Z^istr` のそれに帰着する。
`Base` の側は `Base (hullMap Z)` が同型であることで対応する。

★★**原文はこの依存を書いていない**(「it follows immediately」で済ませている)。
★**随伴は (v) の主張の一部であるだけでなく、(v) を証明する道具でもある。** -/
theorem istr_isPullBack_to {X Y : Istr P} (g : X ⟶ Y) (h : IsPullBack (istrPre P F) g) :
    IsPullBack P g.hom := by
  intro Z
  haveI hZb : IsIso (P.Base (hullMap P F Z)) := (hullMap_spec P F Z).2.1.2
  haveI hZe : Epi (hullMap P F Z) := P.totEpiC _ _ _
  set u := hullMap P F Z with hu
  haveI hub : IsIso (P.Base u) := hZb
  constructor
  · intro f₁ f₂ hf
    have hp := Subtype.ext_iff.mp hf
    have h1 : (f₁ ≫ g.hom : Z ⟶ Y.obj) = f₂ ≫ g.hom := congrArg Prod.fst hp
    have h2 : P.Base f₁ = P.Base f₂ := congrArg Prod.snd hp
    have e₁ : u ≫ ((hullHomEquiv P F Z X).symm f₁).hom = f₁ :=
      (hullHomEquiv P F Z X).apply_symm_apply f₁
    have e₂ : u ≫ ((hullHomEquiv P F Z X).symm f₂).hom = f₂ :=
      (hullHomEquiv P F Z X).apply_symm_apply f₂
    have hgg : (hullHomEquiv P F Z X).symm f₁ = (hullHomEquiv P F Z X).symm f₂ := by
      refine (h (hullIstr P F Z)).1 (Subtype.ext (Prod.ext ?_ ?_))
      · refine InducedCategory.hom_ext ?_
        refine (cancel_epi u).mp ?_
        calc u ≫ (((hullHomEquiv P F Z X).symm f₁) ≫ g).hom
            = (u ≫ ((hullHomEquiv P F Z X).symm f₁).hom) ≫ g.hom :=
              (Category.assoc _ _ _).symm
          _ = f₁ ≫ g.hom := by rw [e₁]
          _ = f₂ ≫ g.hom := h1
          _ = (u ≫ ((hullHomEquiv P F Z X).symm f₂).hom) ≫ g.hom := by rw [e₂]
          _ = u ≫ (((hullHomEquiv P F Z X).symm f₂) ≫ g).hom := Category.assoc _ _ _
      · show P.Base ((hullHomEquiv P F Z X).symm f₁).hom
          = P.Base ((hullHomEquiv P F Z X).symm f₂).hom
        refine (cancel_epi (P.Base u)).mp ?_
        calc P.Base u ≫ P.Base ((hullHomEquiv P F Z X).symm f₁).hom
            = P.Base (u ≫ ((hullHomEquiv P F Z X).symm f₁).hom) := (P.Base_comp _ _).symm
          _ = P.Base f₁ := by rw [e₁]
          _ = P.Base f₂ := h2
          _ = P.Base (u ≫ ((hullHomEquiv P F Z X).symm f₂).hom) := by rw [e₂]
          _ = P.Base u ≫ P.Base ((hullHomEquiv P F Z X).symm f₂).hom := P.Base_comp _ _
    calc f₁ = u ≫ ((hullHomEquiv P F Z X).symm f₁).hom := e₁.symm
      _ = u ≫ ((hullHomEquiv P F Z X).symm f₂).hom := by rw [hgg]
      _ = f₂ := e₂
  · rintro ⟨⟨a, b⟩, hab⟩
    obtain ⟨w, hw1, hw2⟩ := hZb.out
    have ea : u ≫ ((hullHomEquiv P F Z Y).symm a).hom = a :=
      (hullHomEquiv P F Z Y).apply_symm_apply a
    have hcond : (istrPre P F).Base ((hullHomEquiv P F Z Y).symm a)
        = (w ≫ b) ≫ (istrPre P F).Base g := by
      show P.Base ((hullHomEquiv P F Z Y).symm a).hom
        = (w ≫ b) ≫ P.Base g.hom
      have hbase : P.Base u ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom
          = b ≫ P.Base g.hom := by
        rw [← P.Base_comp, ea, hab]
      calc P.Base ((hullHomEquiv P F Z Y).symm a).hom
          = 𝟙 _ ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom := (Category.id_comp _).symm
        _ = (w ≫ P.Base u) ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom :=
              congrArg (fun t => t ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom) hw2.symm
        _ = w ≫ (P.Base u ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom) :=
              Category.assoc _ _ _
        _ = w ≫ (b ≫ P.Base g.hom) := by rw [hbase]
        _ = (w ≫ b) ≫ P.Base g.hom := (Category.assoc _ _ _).symm
    obtain ⟨f', hf'⟩ := (h (hullIstr P F Z)).2
      ⟨((hullHomEquiv P F Z Y).symm a, w ≫ b), hcond⟩
    have hp := Subtype.ext_iff.mp hf'
    have k1 : (f' ≫ g : hullIstr P F Z ⟶ Y) = (hullHomEquiv P F Z Y).symm a :=
      congrArg Prod.fst hp
    have k2 : P.Base f'.hom = w ≫ b := congrArg Prod.snd hp
    have hk : f'.hom ≫ g.hom = ((hullHomEquiv P F Z Y).symm a).hom :=
      congrArg InducedCategory.Hom.hom k1
    refine ⟨u ≫ f'.hom, Subtype.ext (Prod.ext ?_ ?_)⟩
    · show (u ≫ f'.hom) ≫ g.hom = a
      calc (u ≫ f'.hom) ≫ g.hom = u ≫ (f'.hom ≫ g.hom) := Category.assoc _ _ _
        _ = u ≫ ((hullHomEquiv P F Z Y).symm a).hom := congrArg (fun t => u ≫ t) hk
        _ = a := ea
    · show P.Base (u ≫ f'.hom) = b
      calc P.Base (u ≫ f'.hom) = P.Base u ≫ P.Base f'.hom := P.Base_comp _ _
        _ = P.Base u ≫ (w ≫ b) := congrArg (fun t => P.Base u ≫ t) k2
        _ = (P.Base u ≫ w) ≫ b := (Category.assoc _ _ _).symm
        _ = 𝟙 _ ≫ b := by rw [hw1]
        _ = b := Category.id_comp _

/-! ### ★残りの条の移送 —— 「前向き／後ろ向き」の仕分け

★**前向き**(与えられた対象の間で射を作る・一意性を言う)は、
新しい対象を作らないので**そのまま移送できる**。
★**後ろ向き**(`𝒟` の情報から `𝒞` の対象を作る)は **isotropification が要る**。

| 条 | 向き |
|---|---|
| `pullBackLB` / `preStepFactorUniq` / `preStepFactorUniq'` / `arbFactorUniq` | **前向き** |
| `faithfulUpToUnits` / `otriFwd` / `otriBwd` / `otriBase` | **前向き** |
| `baseSurj` / `preStepSpan` | ★**後ろ向き** |
| `plBkEquiv` | 圏同値(別扱い) |
-/

include F in
/-- **(iv)(b)** の移送 —— ★**逆向きの pull-back 移送が開いた**ので通る。 -/
theorem istr_pullBackLB {X Y : Istr P} (α : X ⟶ Y) (h : IsPullBack (istrPre P F) α) :
    IsLBInvertible (istrPre P F) α ∧ IsLinear (istrPre P F) α := by
  obtain ⟨hlb, hlin⟩ := F.pullBackLB α.hom (istr_isPullBack_to P F α h)
  exact ⟨⟨istr_coAngular P F α, hlb.2⟩, hlin⟩

include F in
/-- **(v)(b)** の一意性の移送 —— `C` で得た同型を充満部分圏へ持ち上げるだけ。 -/
theorem istr_preStepFactorUniq {A B : Istr P} (X X' : Istr P)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβc : IsCoAngular (istrPre P F) β) (hβs : IsPreStep (istrPre P F) β)
    (hαi : IsIsometric (istrPre P F) α) (hαs : IsPreStep (istrPre P F) α)
    (hβc' : IsCoAngular (istrPre P F) β') (hβs' : IsPreStep (istrPre P F) β')
    (hαi' : IsIsometric (istrPre P F) α') (hαs' : IsPreStep (istrPre P F) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  obtain ⟨γ₀, h1, h2⟩ := F.preStepFactorUniq X.obj X'.obj β.hom α.hom β'.hom α'.hom
    (congrArg InducedCategory.Hom.hom heq)
    (isCoAngular_of_isotropic_dom P F A.property _) hβs hαi hαs
    (isCoAngular_of_isotropic_dom P F A.property _) hβs' hαi' hαs'
  exact ⟨InducedCategory.isoMk γ₀, InducedCategory.hom_ext h1, InducedCategory.hom_ext h2⟩

include F in
/-- **(v)(c)** の一意性の移送。 -/
theorem istr_preStepFactorUniq' {A B : Istr P} (X X' : Istr P)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβi : IsIsometric (istrPre P F) β) (hβs : IsPreStep (istrPre P F) β)
    (hαc : IsCoAngular (istrPre P F) α) (hαs : IsPreStep (istrPre P F) α)
    (hβi' : IsIsometric (istrPre P F) β') (hβs' : IsPreStep (istrPre P F) β')
    (hαc' : IsCoAngular (istrPre P F) α') (hαs' : IsPreStep (istrPre P F) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  obtain ⟨γ₀, h1, h2⟩ := F.preStepFactorUniq' X.obj X'.obj β.hom α.hom β'.hom α'.hom
    (congrArg InducedCategory.Hom.hom heq)
    hβi hβs (isCoAngular_of_isotropic_dom P F X.property _) hαs
    hβi' hβs' (isCoAngular_of_isotropic_dom P F X'.property _) hαs'
  exact ⟨InducedCategory.isoMk γ₀, InducedCategory.hom_ext h1, InducedCategory.hom_ext h2⟩

include F in
/-- **(vi)** の移送 —— `𝒪^×` の元も `C` から持ち上がる。 -/
theorem istr_faithfulUpToUnits {A B : Istr P} (φ ψ : A ⟶ B)
    (hb : BaseEquivalent (istrPre P F) φ ψ) (hm : MetricallyEquivalent (istrPre P F) φ ψ)
    (hφc : IsCoAngular (istrPre P F) φ) (hφs : IsPreStep (istrPre P F) φ)
    (hψc : IsCoAngular (istrPre P F) ψ) (hψs : IsPreStep (istrPre P F) ψ) :
    ∃ α : End B, α ∈ OTimes (istrPre P F) B ∧ φ = ψ ≫ (α : B ⟶ B) := by
  obtain ⟨α₀, hα₀, hφψ⟩ := F.faithfulUpToUnits φ.hom ψ.hom hb hm
    (isCoAngular_of_isotropic_dom P F A.property _) hφs
    (isCoAngular_of_isotropic_dom P F A.property _) hψs
  refine ⟨InducedCategory.homMk α₀, ⟨⟨hα₀.1.1, hα₀.1.2⟩, ?_⟩, InducedCategory.hom_ext hφψ⟩
  obtain ⟨v, hv⟩ := hα₀.2
  haveI : IsIso (α₀ : B.obj ⟶ B.obj) := (isUnit_iff_isIso _).mp hα₀.2
  refine (isUnit_iff_isIso _).mpr ?_
  exact ⟨InducedCategory.homMk (inv (α₀ : B.obj ⟶ B.obj)),
    InducedCategory.hom_ext (by simp), InducedCategory.hom_ext (by simp)⟩

include F in
/-- **(iv)(a)** の一意性の移送 —— 2つの同型をどちらも持ち上げる。 -/
theorem istr_arbFactorUniq {A B : Istr P} (X Y X' Y' : Istr P)
    (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B) (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B)
    (heq : (γ ≫ β ≫ α : A ⟶ B) = γ' ≫ β' ≫ α')
    (hγ : IsFrobeniusType (istrPre P F) γ) (hβ : IsPreStep (istrPre P F) β)
    (hα : IsPullBack (istrPre P F) α)
    (hγ' : IsFrobeniusType (istrPre P F) γ') (hβ' : IsPreStep (istrPre P F) β')
    (hα' : IsPullBack (istrPre P F) α') :
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  obtain ⟨δ₀, ε₀, h1, h2, h3⟩ := F.arbFactorUniq X.obj Y.obj X'.obj Y'.obj
    γ.hom β.hom α.hom γ'.hom β'.hom α'.hom
    (congrArg InducedCategory.Hom.hom heq)
    ((istr_frobType_iff P F γ).mp hγ) hβ (istr_isPullBack_to P F α hα)
    ((istr_frobType_iff P F γ').mp hγ') hβ' (istr_isPullBack_to P F α' hα')
  exact ⟨InducedCategory.isoMk δ₀, InducedCategory.isoMk ε₀,
    InducedCategory.hom_ext h1, InducedCategory.hom_ext h2, InducedCategory.hom_ext h3⟩

include F in
/-- **(iii)(c)** 順方向の移送。★`𝒪^▷` の元は `.hom` でそのまま対応する。 -/
theorem istr_otriFwd {A B : Istr P} (φ : A ⟶ B) (hst : IsPreStep (istrPre P F) φ)
    (α : End A) (hα : α ∈ OTri (istrPre P F) A) :
    ∃! β : End B, β ∈ OTri (istrPre P F) B ∧ (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  obtain ⟨β₀, ⟨hβ₀m, hβ₀e⟩, hβ₀u⟩ := F.otriFwd φ.hom
    (isCoAngular_of_isotropic_dom P F A.property _) hst α.hom hα
  refine ⟨InducedCategory.homMk β₀, ⟨hβ₀m, InducedCategory.hom_ext hβ₀e⟩, ?_⟩
  rintro β ⟨hβm, hβe⟩
  exact InducedCategory.hom_ext
    (hβ₀u β.hom ⟨hβm, congrArg InducedCategory.Hom.hom hβe⟩)

include F in
/-- **(iii)(c)** 逆方向の移送。 -/
theorem istr_otriBwd {A B : Istr P} (φ : A ⟶ B) (hst : IsPreStep (istrPre P F) φ)
    (β : End B) (hβ : β ∈ OTri (istrPre P F) B) :
    ∃! α : End A, α ∈ OTri (istrPre P F) A ∧ (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  obtain ⟨α₀, ⟨hα₀m, hα₀e⟩, hα₀u⟩ := F.otriBwd φ.hom
    (isCoAngular_of_isotropic_dom P F A.property _) hst β.hom hβ
  refine ⟨InducedCategory.homMk α₀, ⟨hα₀m, InducedCategory.hom_ext hα₀e⟩, ?_⟩
  rintro α ⟨hαm, hαe⟩
  exact InducedCategory.hom_ext
    (hα₀u α.hom ⟨hαm, congrArg InducedCategory.Hom.hom hαe⟩)

include F in
/-- **(iii)(c)** `Base` にしか依らないことの移送。 -/
theorem istr_otriBase {A B : Istr P} (φ φ' : A ⟶ B)
    (hst : IsPreStep (istrPre P F) φ) (hst' : IsPreStep (istrPre P F) φ')
    (hbase : (istrPre P F).Base φ = (istrPre P F).Base φ')
    (α : End A) (hα : α ∈ OTri (istrPre P F) A)
    (β : End B) (hβ : β ∈ OTri (istrPre P F) B)
    (h : (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ) :
    (φ' ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ' :=
  InducedCategory.hom_ext
    (F.otriBase φ.hom φ'.hom (isCoAngular_of_isotropic_dom P F A.property _) hst
      (isCoAngular_of_isotropic_dom P F A.property _) hst' hbase α.hom hα β.hom hβ
      (congrArg InducedCategory.Hom.hom h))

include F in
/-- ★**base-identity 自己射を保つ**。四角形の `Base` 成分で `Base (hullMap)` を消す。 -/
theorem isotropification_baseIdentity {A : C} (e : A ⟶ A) (h : IsBaseIdentity P e) :
    IsBaseIdentity P (istrMap P F e) := by
  haveI hb : IsIso (P.Base (hullMap P F A)) := (hullMap_spec P F A).2.1.2
  have hsq := congrArg P.Base (isotropification_square P F e)
  rw [P.Base_comp, P.Base_comp] at hsq
  have he : P.Base e = 𝟙 _ := h.trans (P.Base_id A)
  show P.Base (istrMap P F e) = P.Base (𝟙 _)
  rw [P.Base_id]
  refine (cancel_epi (P.Base (hullMap P F A))).mp ?_
  rw [hsq, he, Category.id_comp, Category.comp_id]

include F in
/-- ★**(v)** の保存リスト —— `Div-identity` 自己射。

★`Φ(Base e) = id` なので、同型で共役を取っても `id` のまま。 -/
theorem isotropification_divIdentity {A : C} (e : A ⟶ A) (h : IsDivIdentity P e) :
    IsDivIdentity P (istrMap P F e) := by
  haveI hb : IsIso (P.Base (hullMap P F A)) := (hullMap_spec P F A).2.1.2
  have hsq := congrArg P.Base (isotropification_square P F e)
  rw [P.Base_comp, P.Base_comp] at hsq
  have he : Φ.map (P.Base e) = Φ.map (𝟙 (P.toElem.obj A).base) := by
    rw [← P.Base_id A]; exact h
  have hfac : P.Base (istrMap P F e)
      = inv (P.Base (hullMap P F A)) ≫ P.Base e ≫ P.Base (hullMap P F A) := by
    rw [← hsq, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  show Φ.map (P.Base (istrMap P F e)) = Φ.map (P.Base (𝟙 _))
  rw [P.Base_id, hfac]
  ext x
  show Φ.map (inv (P.Base (hullMap P F A)) ≫ P.Base e ≫ P.Base (hullMap P F A)) x
    = Φ.map (𝟙 _) x
  have he' : ∀ y, Φ.map (P.Base e) y = y := fun y => by
    rw [DFunLike.congr_fun he y, Φ.map_id]
  rw [Φ.map_comp (P.Base e ≫ P.Base (hullMap P F A)) (inv (P.Base (hullMap P F A))),
    Φ.map_comp (P.Base (hullMap P F A)) (P.Base e), he',
    ← Φ.map_comp (P.Base (hullMap P F A)) (inv (P.Base (hullMap P F A))), IsIso.inv_hom_id]

include F in
/-- ★**(v)** の保存リスト —— `base-FSM-morphism`。

★`Base` が同型で挟まれるだけなので、`fiberwise-surjective` も `mono` も保たれる。 -/
theorem isotropification_baseFSM {A B : C} (f : A ⟶ B) (h : IsBaseFSM P f) :
    IsBaseFSM P (istrMap P F f) := by
  haveI hbA : IsIso (P.Base (hullMap P F A)) := (hullMap_spec P F A).2.1.2
  haveI hbB : IsIso (P.Base (hullMap P F B)) := (hullMap_spec P F B).2.1.2
  have hsq := congrArg P.Base (isotropification_square P F f)
  rw [P.Base_comp, P.Base_comp] at hsq
  have hfac : P.Base (istrMap P F f)
      = inv (P.Base (hullMap P F A)) ≫ P.Base f ≫ P.Base (hullMap P F B) := by
    rw [← hsq, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  obtain ⟨hfs, hmono⟩ := h
  constructor
  · intro Z γ
    obtain ⟨Dd, δB, δZ, hδ⟩ := hfs (γ ≫ inv (P.Base (hullMap P F B)))
    refine ⟨Dd, δB ≫ P.Base (hullMap P F A), δZ, ?_⟩
    rw [hfac, Category.assoc, ← Category.assoc (P.Base (hullMap P F A)), IsIso.hom_inv_id,
      Category.id_comp, ← Category.assoc, hδ, Category.assoc, Category.assoc,
      IsIso.inv_hom_id, Category.comp_id]
  · rw [hfac]
    haveI := hmono
    infer_instance

include F in
/-- ★**(v)** の保存リスト —— `LB-invertible`。

★`𝒞^istr` では co-angular が自動なので、実質は isometric の移送だけ。 -/
theorem isotropification_lbInvertible {A B : C} (f : A ⟶ B) (h : IsLBInvertible P f) :
    IsLBInvertible P (istrMap P F f) :=
  ⟨isCoAngular_of_isotropic_dom P F (hullMap_spec P F A).2.2.1 _,
   (isotropification_isometric_iff P F f).mpr h.2⟩

include F in
/-- ★**(v)** の保存リスト —— `co-angular`。★**監査で「列挙どおりの形が無い」と指摘された 1 件。**

原文 (FrdI p.32):
> morphisms, base-FSM-morphisms, base-identity en-

★★**仮定は要らない** —— `𝒞^istr` の対象は isotropic なので、
そこから出る射はすべて co-angular(`isCoAngular_of_isotropic_dom`)。
★原文の「preserves co-angular morphisms」より**強い**。 -/
theorem isotropification_coAngular {A B : C} (f : A ⟶ B) :
    IsCoAngular P (istrMap P F f) :=
  isCoAngular_of_isotropic_dom P F (hullMap_spec P F A).2.2.1 _

/-! ### ★★(v) の「through which the functor `C →FΦ` factors」

原文 (FrdI p.32):
> Bistr forms a left adjoint to the inclusion functor Cistr →C, through which

★★**原文の代名詞「which」の先行詞が曖昧である**(測定として記録する)。
直前の名詞句は「the inclusion functor `Cistr →C`」だが、
★**包含関手は `𝒞` へ**入る**ので、`𝒞 → 𝔽_Φ` がそれを経由することはできない**(向きが合わない)。
★**型が通る唯一の読みは「isotropification 関手 `𝒞 → 𝒞^istr` を経由する」**である。
★我々はこちらを採る。★**原文の曖昧さであり、我々の選択である**ことを明記しておく。

★**「経由する」の中身**: `A` とその isotropic hull `A^istr` は `𝒞` では別の対象だが、
★★**`𝔽_Φ` に落とすと同型になる。** isotropic hull は
- isometric ⟹ `Div = 0`
- pre-step ⟹ `IsLinear`(`degFr = 1`)かつ base-isomorphism
なので、`𝔽_Φ` の同型判定 `isIso_iff`(base 同型 ∧ div 可逆 ∧ deg = 1)を満たす。
-/

include F in
/-- ★★**isotropic hull は `𝔽_Φ` では同型になる**。

★これが「factors」の中身である。 -/
theorem toElem_map_hullMap_isIso (A : C) : IsIso (P.toElem.map (hullMap P F A)) := by
  obtain ⟨hisom, hstep, -, -⟩ := hullMap_spec P F A
  refine (ElemFrobCat.isIso_iff _).mpr ⟨hstep.2, ?_, hstep.1⟩
  show IsAddUnit (P.Div (hullMap P F A))
  rw [show P.Div (hullMap P F A) = 0 from hisom]
  exact isAddUnit_zero

include F in
/-- ★★**`Proposition 1.9, (v)` の「through which the functor `C →FΦ` factors」**。

原文 (FrdI p.32):
> the functor C →FΦ factors.

★**`𝒞 → 𝔽_Φ` は isotropification を経由する**(自然同型を除いて)。
★成分は `toElem_map_hullMap_isIso`、自然性は `isotropification_square` そのもの。 -/
noncomputable def isotropificationFactorIso :
    P.toElem ≅ isotropification P F ⋙ (istrPre P F).toElem :=
  NatIso.ofComponents
    (fun A => @asIso _ _ _ _ (P.toElem.map (hullMap P F A)) (toElem_map_hullMap_isIso P F A))
    (fun {A B} f => by
      show P.toElem.map f ≫ P.toElem.map (hullMap P F B)
        = P.toElem.map (hullMap P F A) ≫ P.toElem.map (istrMap P F f)
      rw [← P.toElem.map_comp, ← P.toElem.map_comp, isotropification_square P F f])

include F in
/-- ★**(v)** —— **isotropification 関手の `𝒞^istr` への制限は恒等関手と同型**。

原文 (FrdI p.32):
> morphic to the identity functor. Finally, the isotropification functor preserves

★成分は `hullMap_isIso`(isotropic な対象では hull への射が同型)、
自然性は `isotropification_square` そのもの。 -/
noncomputable def isotropificationRestrictIso :
    (isotropicProp P).ι ⋙ isotropification P F ≅ 𝟭 (Istr P) :=
  NatIso.ofComponents
    (fun X => (InducedCategory.isoMk
      (@asIso _ _ _ _ (hullMap P F X.obj) (hullMap_isIso P F X.obj X.property))).symm)
    (fun {X Y} g => by
      haveI hX : IsIso (hullMap P F X.obj) := hullMap_isIso P F X.obj X.property
      haveI hY : IsIso (hullMap P F Y.obj) := hullMap_isIso P F Y.obj Y.property
      refine InducedCategory.hom_ext ?_
      show istrMap P F g.hom ≫ inv (hullMap P F Y.obj)
        = inv (hullMap P F X.obj) ≫ g.hom
      rw [IsIso.comp_inv_eq, Category.assoc, IsIso.eq_inv_comp]
      exact isotropification_square P F g.hom)

include F in
/-- **(i)(a)** の移送(★**後ろ向き** —— isotropification が要る)。

`𝒞` で得た Frobenius-trivial な対象 `A₀` の isotropic hull を取る。
* Frobenius-trivial の `ζ` は **isotropification が関手であること**でそのまま運べる
* 底の同型は `Base (hullMap)` が同型であることで繋ぐ -/
theorem istr_baseSurj (Y : D) :
    ∃ A : Istr P, IsFrobeniusTrivial (istrPre P F) A ∧
      Nonempty (((istrPre P F).toElem.obj A).base ≅ Y) := by
  obtain ⟨A₀, ⟨ζ, hdeg, hprop⟩, ⟨e⟩⟩ := F.baseSurj Y
  haveI hb : IsIso (P.Base (hullMap P F A₀)) := (hullMap_spec P F A₀).2.1.2
  refine ⟨hullIstr P F A₀, ⟨⟨⟨fun n => (isotropification P F).map (ζ n), ?_⟩, ?_⟩, ?_, ?_⟩,
    ⟨(asIso (P.Base (hullMap P F A₀))).symm ≪≫ e⟩⟩
  · show (isotropification P F).map (ζ 1) = 𝟙 _
    rw [show ζ 1 = 𝟙 A₀ from ζ.map_one]
    exact (isotropification P F).map_id _
  · intro m n
    show (isotropification P F).map (ζ (m * n))
      = (isotropification P F).map (ζ n) ≫ (isotropification P F).map (ζ m)
    rw [ζ.map_mul]
    exact (isotropification P F).map_comp (ζ n) (ζ m)
  · intro n
    show P.degFr (istrMap P F (ζ n)) = n
    rw [isotropification_degFr]
    exact hdeg n
  · intro n
    exact ⟨isotropification_baseIdentity P F (ζ n) (hprop n).1,
      (istr_frobType_iff P F _).mpr (isotropification_frobType P F (ζ n) (hprop n).2)⟩

include F in
/-- **(i)(b)** の移送(★**後ろ向き** —— isotropification が要る)。

`𝒞` で得た span `X₀ ⟶ A.obj`, `X₀ ⟶ B.obj` を isotropification で運び、
`A`, `B` が既に isotropic なので `hullMap` の逆で戻す。

★★**手順2(`inv` を書かない)の例外条件**: **主張そのものが
`@inv _ _ _ _ (Base φ) hφ.2` を含む**ので避けられない。
`IsIso.eq_inv_comp` を**インスタンスを明示して**当てる形で扱う。
★一方、証明の中で使う逆射はすべて `IsIso.out` から明示的に取り、`inv` は書かない。
★**手順は固定するものではなく、反例が出たら条件を精密化するもの** という形の一例。 -/
theorem istr_preStepSpan (A B : Istr P)
    (α : ((istrPre P F).toElem.obj A).base ⟶ ((istrPre P F).toElem.obj B).base) (hα : IsIso α) :
    ∃ (X : Istr P) (φ : X ⟶ A) (ψ : X ⟶ B) (hφ : IsPreStep (istrPre P F) φ),
      IsPreStep (istrPre P F) ψ ∧
        α = @inv _ _ _ _ ((istrPre P F).Base φ) hφ.2 ≫ (istrPre P F).Base ψ := by
  obtain ⟨X₀, φ₀, ψ₀, hφ₀, hψ₀, heq⟩ := F.preStepSpan A.obj B.obj α hα
  haveI hmA : IsIso (hullMap P F A.obj) := hullMap_isIso P F A.obj A.property
  haveI hmB : IsIso (hullMap P F B.obj) := hullMap_isIso P F B.obj B.property
  obtain ⟨wA, hwA1, -⟩ := hmA.out
  obtain ⟨wB, hwB1, -⟩ := hmB.out
  haveI hiA : IsIso wA := ⟨hullMap P F A.obj, by
      refine (cancel_epi (hullMap P F A.obj)).mp ?_
      rw [← Category.assoc, hwA1, Category.id_comp, Category.comp_id], hwA1⟩
  haveI hiB : IsIso wB := ⟨hullMap P F B.obj, by
      refine (cancel_epi (hullMap P F B.obj)).mp ?_
      rw [← Category.assoc, hwB1, Category.id_comp, Category.comp_id], hwB1⟩
  haveI hXb : IsIso (P.Base (hullMap P F X₀)) := (hullMap_spec P F X₀).2.1.2
  obtain ⟨v, hv1, hv2⟩ := hXb.out
  -- ★`Base` の計算: hull を通しても `v` を挟むだけ
  have key : ∀ (Z : C) (f : X₀ ⟶ Z) (w : hullObj P F Z ⟶ Z),
      hullMap P F Z ≫ w = 𝟙 Z → P.Base (istrMap P F f ≫ w) = v ≫ P.Base f := by
    intro Z f w hw
    refine (cancel_epi (P.Base (hullMap P F X₀))).mp ?_
    have hsq := congrArg P.Base (isotropification_square P F f)
    rw [P.Base_comp, P.Base_comp] at hsq
    calc P.Base (hullMap P F X₀) ≫ P.Base (istrMap P F f ≫ w)
        = P.Base (hullMap P F X₀) ≫ (P.Base (istrMap P F f) ≫ P.Base w) := by
          rw [P.Base_comp]
      _ = (P.Base (hullMap P F X₀) ≫ P.Base (istrMap P F f)) ≫ P.Base w :=
          (Category.assoc _ _ _).symm
      _ = (P.Base f ≫ P.Base (hullMap P F Z)) ≫ P.Base w := by rw [hsq]
      _ = P.Base f ≫ (P.Base (hullMap P F Z) ≫ P.Base w) := Category.assoc _ _ _
      _ = P.Base f ≫ P.Base (hullMap P F Z ≫ w) := by rw [P.Base_comp]
      _ = P.Base f ≫ P.Base (𝟙 Z) := by rw [hw]
      _ = P.Base f := by rw [P.Base_id, Category.comp_id]
      _ = (P.Base (hullMap P F X₀) ≫ v) ≫ P.Base f := by rw [hv1, Category.id_comp]
      _ = P.Base (hullMap P F X₀) ≫ (v ≫ P.Base f) := Category.assoc _ _ _
  have hpsA : IsPreStep P (istrMap P F φ₀ ≫ wA) :=
    IsPreStep.comp P ((isotropification_preStep_iff P F φ₀).mpr hφ₀)
      (isPreStep_of_isIso P wA)
  have hpsB : IsPreStep P (istrMap P F ψ₀ ≫ wB) :=
    IsPreStep.comp P ((isotropification_preStep_iff P F ψ₀).mpr hψ₀)
      (isPreStep_of_isIso P wB)
  refine ⟨hullIstr P F X₀, InducedCategory.homMk (istrMap P F φ₀ ≫ wA),
    InducedCategory.homMk (istrMap P F ψ₀ ≫ wB), hpsA, hpsB, ?_⟩
  -- `heq` を `Base φ₀ ≫ α = Base ψ₀` に直す
  have heq' : P.Base φ₀ ≫ α = P.Base ψ₀ :=
    (@IsIso.eq_inv_comp _ _ _ _ _ (P.Base φ₀) hφ₀.2 _ _).mp heq
  refine (@IsIso.eq_inv_comp _ _ _ _ _ ((istrPre P F).Base
    (InducedCategory.homMk (istrMap P F φ₀ ≫ wA) : hullIstr P F X₀ ⟶ A)) hpsA.2 _ _).mpr ?_
  show P.Base (istrMap P F φ₀ ≫ wA) ≫ α = P.Base (istrMap P F ψ₀ ≫ wB)
  rw [key A.obj φ₀ wA hwA1, key B.obj ψ₀ wB hwB1, Category.assoc, heq']
  rfl

/-! ### ★**(i)(c) の圏同値の移送** —— 21 条の最後 -/

include F in
/-- 補助: `(𝒞^istr)^pl-bk` の `A` 上の対象を `𝒞^pl-bk` の `A.obj` 上へ運ぶ。

★ここで **`istr_isPullBack_to`(難しい向き)** を使う。 -/
def istrPlBkToC (A : Istr P) (Z : Over (⟨A⟩ : PlBk (istrPre P F))) :
    Over (⟨A.obj⟩ : PlBk P) :=
  Over.mk (⟨Z.hom.hom.hom, istr_isPullBack_to P F Z.hom.hom Z.hom.property⟩ :
    (⟨Z.left.obj.obj⟩ : PlBk P) ⟶ (⟨A.obj⟩ : PlBk P))

include F in
/-- 補助: 射のほうの運搬。 -/
def istrPlBkToCMap {A : Istr P} {Z W : Over (⟨A⟩ : PlBk (istrPre P F))} (h : Z ⟶ W) :
    istrPlBkToC P F A Z ⟶ istrPlBkToC P F A W :=
  Over.homMk (⟨h.left.hom.hom, istr_isPullBack_to P F h.left.hom h.left.property⟩ :
      (⟨Z.left.obj.obj⟩ : PlBk P) ⟶ (⟨W.left.obj.obj⟩ : PlBk P))
    (by
      have hw : h.left.hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w h)
      have hw2 : h.left.hom.hom ≫ W.hom.hom.hom = Z.hom.hom.hom :=
        congrArg InducedCategory.Hom.hom hw
      exact InducedWideCategory.Hom.ext hw2)

include F in
/-- **(i)(c)** の移送。★`𝒞^istr` は `𝒞` の**充満**部分圏なので、
`(𝒞^istr)^pl-bk` の `A` 上のスライスは `𝒞^pl-bk` の `A.obj` 上のスライスと
**圏として一致する** —— 対象が一致するのは
`isotropic_dom_of_pullBack`(isotropic への pull-back の始域は isotropic)による。

★充満性・忠実性・本質的全射性の**3つとも** `F.plBkEquiv A.obj` から来る。 -/
theorem istr_plBkEquiv (A : Istr P) :
    (plBkOverFunctor (istrPre P F) A).IsEquivalence := by
  haveI := F.plBkEquiv A.obj
  haveI hfaith : (plBkOverFunctor (istrPre P F) A).Faithful := by
    constructor
    intro Z W f g hfg
    have hb : (istrPre P F).Base f.left.hom = (istrPre P F).Base g.left.hom :=
      congrArg CommaMorphism.left hfg
    have hmap : (plBkOverFunctor P A.obj).map (istrPlBkToCMap P F f)
        = (plBkOverFunctor P A.obj).map (istrPlBkToCMap P F g) :=
      Over.OverMorphism.ext hb
    have h2 := (plBkOverFunctor P A.obj).map_injective hmap
    have h3 : f.left.hom.hom = g.left.hom.hom :=
      congrArg (fun x => InducedWideCategory.Hom.hom (CommaMorphism.left x)) h2
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext (InducedCategory.hom_ext h3))
  haveI hfull : (plBkOverFunctor (istrPre P F) A).Full := by
    constructor
    intro Z W h
    obtain ⟨f', hf'⟩ := (plBkOverFunctor P A.obj).map_surjective
      (show (plBkOverFunctor P A.obj).obj (istrPlBkToC P F A Z) ⟶
          (plBkOverFunctor P A.obj).obj (istrPlBkToC P F A W) from h)
    have hw : f'.left.hom ≫ W.hom.hom.hom = Z.hom.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f')
    refine ⟨Over.homMk (⟨InducedCategory.homMk f'.left.hom,
        istr_isPullBack_of P F (InducedCategory.homMk f'.left.hom) f'.left.property⟩ :
        Z.left ⟶ W.left)
      (InducedWideCategory.Hom.ext (InducedCategory.hom_ext hw)), ?_⟩
    exact Over.OverMorphism.ext (congrArg CommaMorphism.left hf')
  haveI hess : (plBkOverFunctor (istrPre P F) A).EssSurj := by
    constructor
    intro Y
    obtain ⟨Z', hZ'⟩ : ∃ Z' : Over (⟨A.obj⟩ : PlBk P),
        Z' = (plBkOverFunctor P A.obj).objPreimage Y := ⟨_, rfl⟩
    have hiZ : (plBkOverFunctor P A.obj).obj Z' ≅ Y := by
      rw [hZ']; exact (plBkOverFunctor P A.obj).objObjPreimageIso Y
    have hiso : IsIsotropic P Z'.left.obj :=
      isotropic_dom_of_pullBack P F Z'.hom.hom Z'.hom.property A.property
    refine ⟨Over.mk (⟨(InducedCategory.homMk Z'.hom.hom :
          (⟨Z'.left.obj, hiso⟩ : Istr P) ⟶ A),
        istr_isPullBack_of P F _ Z'.hom.property⟩ :
      (⟨(⟨Z'.left.obj, hiso⟩ : Istr P)⟩ : PlBk (istrPre P F)) ⟶ (⟨A⟩ : PlBk (istrPre P F))), ?_⟩
    exact ⟨hiZ⟩
  exact ⟨hfaith, hfull, hess⟩

include F in
/-- ★★★**`𝒞^istr` は `Definition 1.3` の core 21 条をすべて満たす**。

原文 (FrdI p.32):
> (v) Cistr [equipped with the restriction to C of the given functor C →FΦ] is a

★**21 条のうち 19 条は「前向き」で自動**((vii)(b) が `𝒞^istr` を閉じるから)、
**2 条だけが「後ろ向き」**で isotropification を要した(`baseSurj` / `preStepSpan`)。
`plBkEquiv` は両向きが要った(`istr_isPullBack_to` が随伴を使う)。 -/
theorem istr_frobenioidCore : FrobenioidCore (istrPre P F) where
  baseSurj := istr_baseSurj P F
  preStepSpan := istr_preStepSpan P F
  plBkEquiv := istr_plBkEquiv P F
  frobDegSurj := istr_frobDegSurj P F
  frobDegUniq := istr_frobDegUniq P F
  coAngularComp := istr_coAngularComp P F
  coAngularOfPreStep := fun α hca hst φ => istr_coAngularOfPreStep P F α hca hst φ
  otriFwd := fun φ _ hst α hα => istr_otriFwd P F φ hst α hα
  otriBwd := fun φ _ hst β hβ => istr_otriBwd P F φ hst β hβ
  otriBase := fun φ φ' _ hst _ hst' hbase α hα β hβ h =>
    istr_otriBase P F φ φ' hst hst' hbase α hα β hβ h
  arbFactor := istr_arbFactor P F
  arbFactorUniq := istr_arbFactorUniq P F
  pullBackLB := fun α h => istr_pullBackLB P F α h
  preStepMono := istr_preStepMono P F
  preStepFactor := istr_preStepFactor P F
  preStepFactorUniq := istr_preStepFactorUniq P F
  preStepFactor' := istr_preStepFactor' P F
  preStepFactorUniq' := istr_preStepFactorUniq' P F
  faithfulUpToUnits := istr_faithfulUpToUnits P F
  isotropicHullExists := istr_isotropicHullExists P F
  isotropicClosed := istr_isotropicClosed P F

/-! ### ★**(iii)(d) の圏同値2本の移送**

★`plBkEquiv` と**同じ形**である: 忠実性だけは `𝒞^istr` の中で直接出て
(全射性/単射性 —— `𝒞^istr` も totally epimorphic、pre-step も mono)、
**充満性と本質的全射性は `𝒞` 側の同値から引く**。
違うのは「対象が `𝒞^istr` に残ること」を言う補題だけ:

* コスライス側は **`isotropicClosed`**(前向き、自明)、
* スライス側は **`Proposition 1.9, (iv)`**(co-angular pre-step は co-angular linear)。
-/

section CoaPre

variable [(coaPreProp P).IsMultiplicative] [(coaPreProp (istrPre P F)).IsMultiplicative]

include F in
/-- 補助: `_A(𝒞^istr,coa-pre)` の対象を `_{A.obj}(𝒞^coa-pre)` へ。 -/
def istrCoaPreUnder (A : Istr P)
    (Z : Under (⟨A⟩ : WideSubcategory (coaPreProp (istrPre P F)))) :
    Under (⟨A.obj⟩ : WideSubcategory (coaPreProp P)) :=
  Under.mk (⟨Z.hom.hom.hom,
      isCoAngular_of_isotropic_dom P F A.property _, Z.hom.property.2⟩ :
    (⟨A.obj⟩ : WideSubcategory (coaPreProp P)) ⟶
      (⟨Z.right.obj.obj⟩ : WideSubcategory (coaPreProp P)))

include F in
/-- 補助: `(𝒞^istr,coa-pre)_A` の対象を `(𝒞^coa-pre)_{A.obj}` へ。 -/
def istrCoaPreOver (A : Istr P)
    (Z : Over (⟨A⟩ : WideSubcategory (coaPreProp (istrPre P F)))) :
    Over (⟨A.obj⟩ : WideSubcategory (coaPreProp P)) :=
  Over.mk (⟨Z.hom.hom.hom,
      isCoAngular_of_isotropic_dom P F Z.left.obj.property _, Z.hom.property.2⟩ :
    (⟨Z.left.obj.obj⟩ : WideSubcategory (coaPreProp P)) ⟶
      (⟨A.obj⟩ : WideSubcategory (coaPreProp P)))

include F in
/-- **(iii)(d)** コスライス側の移送。 -/
theorem istr_coaPreUnderEquiv
    (hC : ∀ A : C, (coaPreUnderFunctor P A).IsEquivalence) (A : Istr P) :
    (coaPreUnderFunctor (istrPre P F) A).IsEquivalence := by
  haveI := hC A.obj
  haveI hfaith : (coaPreUnderFunctor (istrPre P F) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : Z.hom.hom ≫ f.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w f)
    have h2 : Z.hom.hom ≫ g.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w g)
    haveI : Epi Z.hom.hom := (istrPre P F).totEpiC _ _ _
    exact Under.UnderMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_epi Z.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreUnderFunctor (istrPre P F) A).Full := by
    constructor
    intro Z W h
    obtain ⟨f', -⟩ := (coaPreUnderFunctor P A.obj).map_surjective
      (show (coaPreUnderFunctor P A.obj).obj (istrCoaPreUnder P F A Z) ⟶
          (coaPreUnderFunctor P A.obj).obj (istrCoaPreUnder P F A W) from h)
    have hw : Z.hom.hom.hom ≫ f'.right.hom = W.hom.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w f')
    refine ⟨Under.homMk (⟨(InducedCategory.homMk f'.right.hom :
          Z.right.obj ⟶ W.right.obj),
        ⟨istr_coAngular P F _, f'.right.property.2⟩⟩ : Z.right ⟶ W.right)
      (InducedWideCategory.Hom.ext (InducedCategory.hom_ext hw)), Subsingleton.elim _ _⟩
  haveI hess : (coaPreUnderFunctor (istrPre P F) A).EssSurj := by
    constructor
    intro Y
    obtain ⟨Z', hZ'⟩ : ∃ Z' : Under (⟨A.obj⟩ : WideSubcategory (coaPreProp P)),
        Z' = (coaPreUnderFunctor P A.obj).objPreimage Y := ⟨_, rfl⟩
    have hiZ : (coaPreUnderFunctor P A.obj).obj Z' ≅ Y := by
      rw [hZ']; exact (coaPreUnderFunctor P A.obj).objObjPreimageIso Y
    have hiso : IsIsotropic P Z'.right.obj := F.isotropicClosed Z'.hom.hom A.property
    exact ⟨Under.mk (⟨(InducedCategory.homMk Z'.hom.hom :
          A ⟶ (⟨Z'.right.obj, hiso⟩ : Istr P)),
        ⟨istr_coAngular P F _, Z'.hom.property.2⟩⟩ :
      (⟨A⟩ : WideSubcategory (coaPreProp (istrPre P F))) ⟶
        (⟨(⟨Z'.right.obj, hiso⟩ : Istr P)⟩ :
          WideSubcategory (coaPreProp (istrPre P F)))), ⟨hiZ⟩⟩
  exact ⟨hfaith, hfull, hess⟩

include F in
/-- **(iii)(d)** スライス側の移送。★対象が `𝒞^istr` に残ることに
**`Proposition 1.9, (iv)` を使う** —— co-angular pre-step は co-angular linear なので、
終域が isotropic なら始域も isotropic。 -/
theorem istr_coaPreOverEquiv
    (hC : ∀ A : C, (coaPreOverFunctor P A).IsEquivalence) (A : Istr P) :
    (coaPreOverFunctor (istrPre P F) A).IsEquivalence := by
  haveI := hC A.obj
  haveI hfaith : (coaPreOverFunctor (istrPre P F) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have h2 : g.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    haveI : Mono W.hom.hom := istr_preStepMono P F _ W.hom.property.2
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_mono W.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreOverFunctor (istrPre P F) A).Full := by
    constructor
    intro Z W h
    obtain ⟨f', -⟩ := (coaPreOverFunctor P A.obj).map_surjective
      (show (coaPreOverFunctor P A.obj).obj (istrCoaPreOver P F A Z) ⟶
          (coaPreOverFunctor P A.obj).obj (istrCoaPreOver P F A W) from h)
    have hw : f'.left.hom ≫ W.hom.hom.hom = Z.hom.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f')
    refine ⟨Over.homMk (⟨(InducedCategory.homMk f'.left.hom :
          Z.left.obj ⟶ W.left.obj),
        ⟨istr_coAngular P F _, f'.left.property.2⟩⟩ : Z.left ⟶ W.left)
      (InducedWideCategory.Hom.ext (InducedCategory.hom_ext hw)), Subsingleton.elim _ _⟩
  haveI hess : (coaPreOverFunctor (istrPre P F) A).EssSurj := by
    constructor
    intro Y
    obtain ⟨Z', hZ'⟩ : ∃ Z' : Over (⟨A.obj⟩ : WideSubcategory (coaPreProp P)),
        Z' = (coaPreOverFunctor P A.obj).objPreimage Y := ⟨_, rfl⟩
    have hiZ : (coaPreOverFunctor P A.obj).obj Z' ≅ Y := by
      rw [hZ']; exact (coaPreOverFunctor P A.obj).objObjPreimageIso Y
    have hiso : IsIsotropic P Z'.left.obj :=
      (prop_1_9_iv P F Z'.hom.hom Z'.hom.property.1 Z'.hom.property.2.1).mpr A.property
    exact ⟨Over.mk (⟨(InducedCategory.homMk Z'.hom.hom :
          (⟨Z'.left.obj, hiso⟩ : Istr P) ⟶ A),
        ⟨istr_coAngular P F _, Z'.hom.property.2⟩⟩ :
      (⟨(⟨Z'.left.obj, hiso⟩ : Istr P)⟩ : WideSubcategory (coaPreProp (istrPre P F))) ⟶
        (⟨A⟩ : WideSubcategory (coaPreProp (istrPre P F)))), ⟨hiZ⟩⟩
  exact ⟨hfaith, hfull, hess⟩

end CoaPre

/-- ★★★**`Proposition 1.9, (v)`** —— `𝒞^istr` は Frobenioid である。 -/
theorem istr_frobenioid (G : Frobenioid P) : Frobenioid (istrPre P F) := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := coaPreProp_isMultiplicative (istrPre P F)
    (istr_frobenioidCore P F).coAngularComp
  exact ⟨istr_frobenioidCore P F,
    istr_coaPreUnderEquiv P F G.coaPreUnderEquiv,
    istr_coaPreOverEquiv P F G.coaPreOverEquiv⟩

end Istr

/-! ## ★第4段 —— (ii)(iii) の梱包: `φ_*` と `φ^*` を `Functor` にする

原文 (FrdI p.31):
> (ii) Any base-isomorphism φ : A →B of C induces a functor [well-defined

原文 (FrdI p.31):
> (iii) Any pull-back morphism φ : A →B of C induces a functor [well-

★★**測定**: 関手の**対象の割り当てにだけ選択が要り**、
**射の割り当ては一意**(`imtrPre_hom_uniq`)、
**関手則(`map_id` / `map_comp`)は無料**である ——
`𝒞^imtr-pre_A` の hom 集合が**高々1元**だからである。
★原文の「[well-defined up to isomorphism]」という但し書きは
**対象の割り当てだけに掛かる**。
-/


def istr_frobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 32, item := "Proposition 1.9, (v)",
    sectionId := "frdi-prop-1-9-v" }

end ABC3.Found.FrdI
