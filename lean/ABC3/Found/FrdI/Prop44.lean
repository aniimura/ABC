import ABC3.Found.FrdI.HomColim
import ABC3.Found.FrdI.Prop111

/-!
# [FrdI] Proposition 4.4 —— `Hom^birat` と `𝒞^birat`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.82–p.83。

原文 (FrdI p.82):
> (Birationalization of a Frobenioid I) For A, B

## ★★`Hom^pf` との違い(測定、2026-08-17)

`Definition 3.1, (ii)` の `Hom^pf` と**同じ形ではない** —— **3 点**で違う:

| | `Hom^pf`(`Definition 3.1, (ii)`) | `Hom^birat`(ここ) |
|---|---|---|
| 添字 | **両側**(`A` と `B` から出る対) | **片側**(`A` へ入る) |
| 向き | **コスライス** | **スライス**(の反対圏) |
| 遷移 | `Proposition 1.10, (i)` の四角形 | **前合成** |

★共有できるのは**層 1**(`HomColim.lean`)まで —— そこは宇宙一般に書いてあるので
そのまま乗る。

## ★★★添字圏は「順序集合」である(この節の当たり)

★★**`Definition 1.3, (iii), (d)` の第 2 の圏同値**(`Frobenioid.coaPreOverEquiv`)が
```
(𝒞^{coa-pre})_A  ≃  Order(Φ(A))ᵒᵖ
```
を与える。⟹ **添字圏 `((𝒞^{coa-pre})_A)ᵒᵖ ≃ Order(Φ(A))`** は**順序圏**であり、

- **細い**(平行 2 射は高々 1 本)—— 直接示すなら
  **pre-step が mono**(`preStepMono`)で消約する。
  ★★`pf` 側が **epi**(totally epimorphic)で消約したのと**双対**である。
- **filtered** —— 上界は `Φ(A)` の**和 `a + b`** で取れる。
  ★★`pf` 側で「次数の最小公倍数」を `frobDegUniq` で作ったのに比べ、**足すだけ**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★1. `Order(M)` は有向 -/

/-- ★★**`Order(M)` は有向**(上界は**和**で取れる)。

★★これが `Hom^birat` の添字圏が filtered であることの中身である。 -/
instance orderCat_isDirected (M : Type w) [AddCommMonoid M] :
    IsDirected (OrderCat M) (· ≤ ·) := by
  refine ⟨fun a b => ⟨a + b, ?_, ?_⟩⟩
  · exact ⟨b, rfl⟩
  · exact ⟨a, add_comm b a⟩

/-! ## ★2. `𝒞^{coa-pre}` とスライス -/

/-- ★`G` を型に載せた `coaPreProp`(`instance` を回復する技、`Def31Pf.lean` と同じ)。 -/
def coaPrePropOf (_G : Frobenioid P) : MorphismProperty C := coaPreProp P

instance coaPrePropOf_isMultiplicative : (coaPrePropOf P G).IsMultiplicative :=
  coaPreProp_isMultiplicative P G.core.coAngularComp

/-- ★★**`𝒞^{coa-pre}`**(`Definition 1.3, (iii), (d)`)。 -/
abbrev CoaPre : Type u2 := WideSubcategory (coaPrePropOf P G)

/-- ★スライスの底。 -/
def coaPreObj (A : C) : CoaPre P G := ⟨A⟩

/-- ★★**スライス `(𝒞^{coa-pre})_A`** —— `A` へ入る co-angular pre-step たち。 -/
abbrev SliceA (A : C) : Type (max u2 v2) := Over (coaPreObj P G A)

/-- ★★**`Hom^birat` の添字圏** —— スライスの**反対圏**。

★★**反対圏になる理由は遷移が「前合成」だから**である ——
原文「the transition morphism induced by a co-angular pre-step `A″ → A′` … is the natural
morphism `Hom(A′,B) → Hom(A″,B)`」で、射の向きと写像の向きが**逆**になる。 -/
abbrev IdxBirat (A : C) : Type (max u2 v2) := (SliceA P G A)ᵒᵖ

/-! ## ★3. 添字圏は細い —— pre-step が mono だから -/

variable {P G} in
/-- ★★**スライスの平行 2 射は等しい** —— 構造射 `A′ → A` が **mono**(pre-step)だから。

★★`pf` 側の `idx_hom_ext` は **epi**(totally epimorphic)で消約した。**双対である。** -/
theorem slice_hom_ext {A : C} {Z W : SliceA P G A} (u v : Z ⟶ W) : u = v := by
  have hu : u.left.hom ≫ W.hom.hom = Z.hom.hom :=
    congrArg (fun t : Z.left ⟶ coaPreObj P G A => t.hom) (Over.w u)
  have hv : v.left.hom ≫ W.hom.hom = Z.hom.hom :=
    congrArg (fun t : Z.left ⟶ coaPreObj P G A => t.hom) (Over.w v)
  haveI : Mono W.hom.hom := G.core.preStepMono _ W.hom.property.2
  exact CommaMorphism.ext (WideSubcategory.hom_ext _ ((cancel_mono W.hom.hom).mp
    (hu.trans hv.symm))) (Subsingleton.elim _ _)

variable {P G} in
/-- ★添字圏(反対圏)でも同じ。 -/
theorem idxBirat_hom_ext {A : C} {Z W : IdxBirat P G A} (u v : Z ⟶ W) : u = v :=
  Quiver.Hom.unop_inj (slice_hom_ext u.unop v.unop)

/-! ## ★4. スライスは cofiltered —— 下界は `Φ(A)` の**和**で取れる -/

variable {P G} in
/-- ★スライスの対象の不変量(`Definition 1.3, (iii), (d)` の関手の値)。 -/
noncomputable def sliceInv {A : C} (Z : SliceA P G A) :
    Φ.val (P.toElem.obj A).base :=
  haveI : IsIso (P.Base Z.hom.hom) := Z.hom.property.2.2
  Φ.map (inv (P.Base Z.hom.hom)) (P.Div Z.hom.hom)

variable {P G} in
/-- ★★**スライスの共通下界** —— 不変量の**和**を実現する対象を取り、
`Definition 1.3, (iii), (d)` の**第 2 の圏同値**で両方へ落とす。

★★`pf` 側の `idx_cocone_objs` は「次数の積」を `frobDegUniq` で合わせたが、
★**こちらは単系の和 1 本**である。 -/
theorem slice_cone_objs {A : C} (Z W : SliceA P G A) :
    ∃ (V : SliceA P G A) (_ : V ⟶ Z) (_ : V ⟶ W), True := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv A
  haveI hZ : IsIso (P.Base Z.hom.hom) := Z.hom.property.2.2
  haveI hW : IsIso (P.Base W.hom.hom) := W.hom.property.2.2
  set V₀ : SliceA P G A :=
    (coaPreOverFunctor P A).objPreimage
      (Opposite.op (toOrderCat (sliceInv Z + sliceInv W))) with hV₀
  haveI hV : IsIso (P.Base V₀.hom.hom) := V₀.hom.property.2.2
  have hle0 : (toOrderCat (sliceInv Z + sliceInv W) : OrderCat (Φ.val (P.toElem.obj A).base))
      ≤ toOrderCat (sliceInv V₀) :=
    leOfHom ((coaPreOverFunctor P A).objObjPreimageIso
      (Opposite.op (toOrderCat (sliceInv Z + sliceInv W)))).hom.unop
  have hZle : MLe (sliceInv Z) (sliceInv V₀) :=
    le_trans (show (toOrderCat (sliceInv Z) : OrderCat _) ≤ _ from ⟨sliceInv W, rfl⟩) hle0
  have hWle : MLe (sliceInv W) (sliceInv V₀) :=
    le_trans (show (toOrderCat (sliceInv W) : OrderCat _) ≤ _ from
      ⟨sliceInv Z, add_comm _ _⟩) hle0
  obtain ⟨βZ, hβZc, hβZs, hβZe⟩ := coaPre_factor_of_mle P G Z.hom.hom
    Z.hom.property.1 Z.hom.property.2 V₀.hom.hom V₀.hom.property.1 V₀.hom.property.2 hZle
  obtain ⟨βW, hβWc, hβWs, hβWe⟩ := coaPre_factor_of_mle P G W.hom.hom
    W.hom.property.1 W.hom.property.2 V₀.hom.hom V₀.hom.property.1 V₀.hom.property.2 hWle
  exact ⟨V₀, Over.homMk (show V₀.left ⟶ Z.left from ⟨βZ, hβZc, hβZs⟩)
      (WideSubcategory.hom_ext _ hβZe),
    Over.homMk (show V₀.left ⟶ W.left from ⟨βW, hβWc, hβWs⟩)
      (WideSubcategory.hom_ext _ hβWe), trivial⟩

/-- ★★**スライスは cofiltered**。 -/
instance sliceA_isCofiltered {A : C} : IsCofiltered (SliceA P G A) where
  cone_objs := slice_cone_objs
  cone_maps _ _ u v := ⟨_, 𝟙 _, by rw [slice_hom_ext u v]⟩
  nonempty := ⟨Over.mk (𝟙 (coaPreObj P G A))⟩

/-- ★★★**添字圏 `((𝒞^{coa-pre})_A)ᵒᵖ` は filtered**。 -/
instance idxBirat_isFiltered {A : C} : IsFiltered (IdxBirat P G A) :=
  inferInstanceAs (IsFiltered (SliceA P G A)ᵒᵖ)

/-! ## ★5. 添字関手と `Hom^birat`

★★遷移写像は**前合成**である —— `pf` 側が `Proposition 1.10, (i)` の四角形を
解いて作ったのに比べ、**ここは合成 1 本**。 -/

/-- ★★**帰納系** —— `(A′ → A) ↦ Hom_𝒞(A′, B)`、遷移は**前合成**。

★★`ULift` を噛ませる理由は宇宙(`Def31Pf.lean` と同じ) ——
添字圏の対象は `Type (max u2 v2)` に住むので、値も同じ宇宙へ上げる。 -/
noncomputable def homFunctorBirat (A B : C) : IdxBirat P G A ⥤ Type (max v u2 v2) where
  obj Z := ULift.{max u2 v} (Z.unop.left.obj ⟶ B)
  map {Z W} u := TypeCat.ofHom fun φ : ULift.{max u2 v} (Z.unop.left.obj ⟶ B) =>
    (ULift.up (u.unop.left.hom ≫ φ.down) : ULift.{max u2 v} (W.unop.left.obj ⟶ B))
  map_id Z := by
    ext φ
    exact Category.id_comp φ.down
  map_comp u u' := by
    ext φ
    exact Category.assoc u'.unop.left.hom u.unop.left.hom φ.down

@[simp] theorem homFunctorBirat_obj {A B : C} (Z : IdxBirat P G A) :
    (homFunctorBirat P G A B).obj Z = ULift.{max u2 v} (Z.unop.left.obj ⟶ B) := rfl

variable {P G} in
@[simp] theorem homFunctorBirat_map {A B : C} {Z W : IdxBirat P G A} (u : Z ⟶ W)
    (φ : ULift.{max u2 v} (Z.unop.left.obj ⟶ B)) :
    (homFunctorBirat P G A B).map u φ = ULift.up (u.unop.left.hom ≫ φ.down) := rfl

/-- ★★★**[FrdI] Proposition 4.4** —— **`Hom^birat_𝒞(A,B)`**。 -/
noncomputable def HomBirat (A B : C) : Type (max v u2 v2) :=
  HomColim (homFunctorBirat P G A B)

variable {P G} in
/-- ★代表元 —— 「co-angular pre-step `A′ → A` と射 `A′ → B`」の組。 -/
noncomputable def HomBirat.mk {A B : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) :
    HomBirat P G A B :=
  HomColim.mk (homFunctorBirat P G A B) Z (ULift.up φ)

variable {P G} in
theorem HomBirat.exists_rep {A B : C} (z : HomBirat P G A B) :
    ∃ (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B), HomBirat.mk Z φ = z := by
  obtain ⟨Z, x, hx⟩ := HomColim.exists_rep (homFunctorBirat P G A B) z
  exact ⟨Z, x.down, hx⟩

variable {P G} in
/-- ★★**前合成で移しても同じ元**(帰納極限の要)。 -/
@[simp] theorem HomBirat.mk_map {A B : C} {Z W : IdxBirat P G A} (u : Z ⟶ W)
    (φ : Z.unop.left.obj ⟶ B) :
    HomBirat.mk W (u.unop.left.hom ≫ φ) = HomBirat.mk Z φ :=
  HomColim.mk_map (homFunctorBirat P G A B) u (ULift.up φ)

variable {P G} in
/-- ★★**共通の下界で一致すれば等しい**(filtered 性がここで効く)。 -/
theorem HomBirat.sound {A B : C} {Z W : IdxBirat P G A}
    {φ : Z.unop.left.obj ⟶ B} {ψ : W.unop.left.obj ⟶ B}
    (V : IdxBirat P G A) (u : Z ⟶ V) (v : W ⟶ V)
    (h : u.unop.left.hom ≫ φ = v.unop.left.hom ≫ ψ) :
    HomBirat.mk Z φ = HomBirat.mk W ψ :=
  HomColim.sound (homFunctorBirat P G A B) V u v (congrArg ULift.up h)

variable {P G} in
/-- ★★**逆向き**。 -/
theorem HomBirat.eq_iff {A B : C} {Z W : IdxBirat P G A}
    {φ : Z.unop.left.obj ⟶ B} {ψ : W.unop.left.obj ⟶ B} :
    HomBirat.mk Z φ = HomBirat.mk W ψ ↔
      ∃ (V : IdxBirat P G A) (u : Z ⟶ V) (v : W ⟶ V),
        u.unop.left.hom ≫ φ = v.unop.left.hom ≫ ψ := by
  refine (HomColim.eq_iff (homFunctorBirat P G A B)).trans ⟨?_, ?_⟩
  · rintro ⟨V, u, v, h⟩
    exact ⟨V, u, v, ULift.up_injective h⟩
  · rintro ⟨V, u, v, h⟩
    exact ⟨V, u, v, congrArg ULift.up h⟩

/-- ★添字圏の「頂点」`𝟙_A : A → A`。 -/
def idxBiratOne (A : C) : IdxBirat P G A :=
  Opposite.op (Over.mk (𝟙 (coaPreObj P G A)))

variable {P G} in
/-- ★★**自然な写像 `Hom_𝒞(A,B) → Hom^birat_𝒞(A,B)`**。 -/
noncomputable def toHomBirat {A B : C} (φ : A ⟶ B) : HomBirat P G A B :=
  HomBirat.mk (idxBiratOne P G A) φ

/-! ## ★6. 合成の材料

★★**合成の構図**: `f` を `(a : A′ → A, φ : A′ → B)`、`g` を `(b : B′ → B, ψ : B′ → E)` で表す。
`φ` と `b` は**同じ `B` へ**向かうが、直接は繋がらない。
★★**`b` を `φ` に沿って引き戻す**必要があり、それが
**`Proposition 1.11, (vii)`**(我々の `prop_1_11_vii` —— 任意の射について成立、取得済み)である:
```
γ ≫ φ = α ≫ b,  γ : Dd → A′ は co-angular pre-step
```
⟹ 新しい添字は `Dd → A′ → A`、射は `Dd → B′ → E`(= `α ≫ ψ`)。

★★★**well-defined である理由は `b` が mono であること**(pre-step は mono)——
下の `birat_lift_unique` がその中身である。 -/

variable {P G} in
/-- ★★★**引き戻しの一意性(共通の下界の上で)** —— **`b` が mono だから**。

★2 通りの引き戻し `(γ, α)` と `(γ′, α′)` があっても、
`Dd″` から両方へ **`A′` の上で**行けるなら、`α` たちは一致する。
★★**これが合成の well-definedness のすべて**である。 -/
theorem birat_lift_unique {A' B B' Dd Dd' Dd'' : C}
    (φ : A' ⟶ B) {b : B' ⟶ B} (hb : Mono b)
    {γ : Dd ⟶ A'} {α : Dd ⟶ B'} (hsq : γ ≫ φ = α ≫ b)
    {γ' : Dd' ⟶ A'} {α' : Dd' ⟶ B'} (hsq' : γ' ≫ φ = α' ≫ b)
    {c : Dd'' ⟶ Dd} {c' : Dd'' ⟶ Dd'} (hc : c ≫ γ = c' ≫ γ') :
    c ≫ α = c' ≫ α' := by
  haveI := hb
  refine (cancel_mono b).mp ?_
  calc (c ≫ α) ≫ b = c ≫ α ≫ b := by rw [Category.assoc]
    _ = c ≫ γ ≫ φ := by rw [hsq]
    _ = (c ≫ γ) ≫ φ := by rw [Category.assoc]
    _ = (c' ≫ γ') ≫ φ := by rw [hc]
    _ = c' ≫ γ' ≫ φ := by rw [Category.assoc]
    _ = c' ≫ α' ≫ b := by rw [hsq']
    _ = (c' ≫ α') ≫ b := by rw [Category.assoc]

/-- ★添字対象の構造射を co-angular pre-step として取り出す。 -/
theorem idxBirat_coaPre {A : C} (Z : IdxBirat P G A) :
    IsCoAngular P Z.unop.hom.hom ∧ IsPreStep P Z.unop.hom.hom :=
  ⟨Z.unop.hom.property.1, Z.unop.hom.property.2⟩

/-- ★添字対象を「co-angular pre-step」から作る。 -/
def idxBiratMk {A A' : C} (a : A' ⟶ A) (hac : IsCoAngular P a) (has : IsPreStep P a) :
    IdxBirat P G A :=
  Opposite.op (Over.mk (Y := (⟨A'⟩ : CoaPre P G))
    (show (⟨A'⟩ : CoaPre P G) ⟶ coaPreObj P G A from ⟨a, hac, has⟩))

@[simp] theorem idxBiratMk_left_obj {A A' : C} (a : A' ⟶ A) (hac : IsCoAngular P a)
    (has : IsPreStep P a) : (idxBiratMk P G a hac has).unop.left.obj = A' := rfl

@[simp] theorem idxBiratMk_hom {A A' : C} (a : A' ⟶ A) (hac : IsCoAngular P a)
    (has : IsPreStep P a) : (idxBiratMk P G a hac has).unop.hom.hom = a := rfl

variable {P G} in
/-- ★添字圏の射を「`A` の上の co-angular pre-step」から作る。 -/
def idxBiratHomMk {A : C} {Z W : IdxBirat P G A} (c : W.unop.left.obj ⟶ Z.unop.left.obj)
    (hcc : IsCoAngular P c) (hcs : IsPreStep P c)
    (hw : c ≫ Z.unop.hom.hom = W.unop.hom.hom) : Z ⟶ W :=
  Quiver.Hom.op (Over.homMk (show W.unop.left ⟶ Z.unop.left from ⟨c, hcc, hcs⟩)
    (WideSubcategory.hom_ext _ hw))

variable {P G} in
@[simp] theorem idxBiratHomMk_left {A : C} {Z W : IdxBirat P G A}
    (c : W.unop.left.obj ⟶ Z.unop.left.obj) (hcc : IsCoAngular P c) (hcs : IsPreStep P c)
    (hw : c ≫ Z.unop.hom.hom = W.unop.hom.hom) :
    (idxBiratHomMk c hcc hcs hw).unop.left.hom = c := rfl

variable {P G} in
/-- ★★**添字圏の射は「`A` 側の三角形」で決まる** —— 構造射が **mono** だから。

★★これは `slice_hom_ext`(細さ)の使いやすい言い換えであり、
合成の自然性で**繰り返し使う**。 -/
theorem idxBirat_left_ext {A : C} {Z W : IdxBirat P G A} (u v : Z ⟶ W) :
    u.unop.left.hom = v.unop.left.hom :=
  congrArg (fun t : W.unop ⟶ Z.unop => t.left.hom) (slice_hom_ext u.unop v.unop)

/-! ## ★7. 引き戻しデータ —— `Proposition 1.11, (vii)` を固定する -/

variable {P G} in
/-- ★★**引き戻しデータの存在** —— `prop_1_11_vii` の言い換え。

`Z` の代表射 `φ : A′ → B` に沿って、`B` へ入る co-angular pre-step
`b : B′ → B`(= `W` の構造射)を引き戻す。 -/
theorem birat_pull_exists (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    ∃ (Dd : C) (γ : Dd ⟶ Z.unop.left.obj) (α : Dd ⟶ W.unop.left.obj),
      IsCoAngular P γ ∧ IsPreStep P γ ∧ γ ≫ φ = α ≫ W.unop.hom.hom :=
  prop_1_11_vii P F G φ W.unop.hom.hom W.unop.hom.property.1 W.unop.hom.property.2

variable {P G} in
/-- ★引き戻しの対象。 -/
noncomputable def biratPullObj (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) : C :=
  (birat_pull_exists F Z φ W).choose

variable {P G} in
/-- ★引き戻しの co-angular pre-step `Dd → A′`。 -/
noncomputable def biratPullGamma (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    biratPullObj F Z φ W ⟶ Z.unop.left.obj :=
  (birat_pull_exists F Z φ W).choose_spec.choose

variable {P G} in
/-- ★引き戻しの射 `Dd → B′`。 -/
noncomputable def biratPullAlpha (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    biratPullObj F Z φ W ⟶ W.unop.left.obj :=
  (birat_pull_exists F Z φ W).choose_spec.choose_spec.choose

variable {P G} in
theorem biratPullGamma_coAngular (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    IsCoAngular P (biratPullGamma F Z φ W) :=
  (birat_pull_exists F Z φ W).choose_spec.choose_spec.choose_spec.1

variable {P G} in
theorem biratPullGamma_preStep (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    IsPreStep P (biratPullGamma F Z φ W) :=
  (birat_pull_exists F Z φ W).choose_spec.choose_spec.choose_spec.2.1

variable {P G} in
/-- ★★引き戻しの四角形。 -/
theorem biratPull_sq (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    biratPullGamma F Z φ W ≫ φ = biratPullAlpha F Z φ W ≫ W.unop.hom.hom :=
  (birat_pull_exists F Z φ W).choose_spec.choose_spec.choose_spec.2.2

variable {P G} in
/-- ★★引き戻しが定める `A` 側の添字対象 `Dd → A′ → A`。 -/
noncomputable def biratPullIdx (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    IdxBirat P G A :=
  idxBiratMk P G (biratPullGamma F Z φ W ≫ Z.unop.hom.hom)
    (G.core.coAngularComp _ _ (biratPullGamma_coAngular F Z φ W) Z.unop.hom.property.1)
    (IsPreStep.comp P (biratPullGamma_preStep F Z φ W) Z.unop.hom.property.2)

variable {P G} in
/-- ★添字圏の射 `Z ⟶ (引き戻した添字)`。 -/
noncomputable def biratPullHom (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    Z ⟶ biratPullIdx F Z φ W :=
  idxBiratHomMk (biratPullGamma F Z φ W) (biratPullGamma_coAngular F Z φ W)
    (biratPullGamma_preStep F Z φ W) rfl

/-! ## ★8. 第 2 変数の自然性 —— 合成の well-definedness の半分 -/

variable {P G} in
/-- ★★**第 2 変数についての自然性**。

★★**証明の骨**は 3 本:
1. `birat_lift_unique`(`W` の構造射が **mono**)
2. `idxBirat_left_ext`(添字圏が**細い**)—— 2 本の `γ` が共通の下界の上で一致する根拠
3. `HomBirat.sound`(添字圏が **filtered**)—— 共通の下界を取る -/
theorem compBirat_natural_right (F : FrobenioidCore P) {A B E : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) {W W' : IdxBirat P G B} (u : W ⟶ W')
    (ψ : W.unop.left.obj ⟶ E) :
    HomBirat.mk (biratPullIdx F Z φ W') (biratPullAlpha F Z φ W' ≫ u.unop.left.hom ≫ ψ)
      = HomBirat.mk (biratPullIdx F Z φ W) (biratPullAlpha F Z φ W ≫ ψ) := by
  haveI hb : Mono W.unop.hom.hom := G.core.preStepMono _ W.unop.hom.property.2
  set V := IsFiltered.max (biratPullIdx F Z φ W') (biratPullIdx F Z φ W) with hV
  set c := IsFiltered.leftToMax (biratPullIdx F Z φ W') (biratPullIdx F Z φ W) with hc0
  set c' := IsFiltered.rightToMax (biratPullIdx F Z φ W') (biratPullIdx F Z φ W) with hc0'
  refine HomBirat.sound V c c' ?_
  have hW' : u.unop.left.hom ≫ W.unop.hom.hom = W'.unop.hom.hom :=
    congrArg (fun t : W'.unop.left ⟶ coaPreObj P G B => t.hom) (Over.w u.unop)
  have hsq2 : biratPullGamma F Z φ W' ≫ φ
      = (biratPullAlpha F Z φ W' ≫ u.unop.left.hom) ≫ W.unop.hom.hom := by
    rw [Category.assoc, hW']
    exact biratPull_sq F Z φ W'
  have hc : c'.unop.left.hom ≫ biratPullGamma F Z φ W
      = c.unop.left.hom ≫ biratPullGamma F Z φ W' :=
    idxBirat_left_ext (biratPullHom F Z φ W ≫ c') (biratPullHom F Z φ W' ≫ c)
  have key := birat_lift_unique φ hb (biratPull_sq F Z φ W) hsq2 hc
  simp only [← Category.assoc] at key ⊢
  exact congrArg (fun t => t ≫ ψ) key.symm

/-! ## ★9. 第 2 変数についての降下 -/

variable {P G} in
/-- ★第 2 変数の余錐(`Z`, `φ` を固定)。 -/
noncomputable def compBiratCoconeRight (F : FrobenioidCore P) {A B E : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) : Cocone (homFunctorBirat P G B E) :=
  Cocone.mk (HomBirat P G A E)
    { app := fun W => TypeCat.ofHom fun ψ : ULift.{max u2 v} (W.unop.left.obj ⟶ E) =>
        HomBirat.mk (biratPullIdx F Z φ W) (biratPullAlpha F Z φ W ≫ ψ.down)
      naturality := fun W W' u => by
        ext ψ
        exact compBirat_natural_right F Z φ u ψ.down }

variable {P G} in
/-- ★★**第 2 変数についての合成**(`Z`, `φ` を固定)。 -/
noncomputable def compBiratRight (F : FrobenioidCore P) {A B E : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) :
    HomBirat P G B E → HomBirat P G A E :=
  fun g => colimit.desc (homFunctorBirat P G B E) (compBiratCoconeRight F Z φ) g

variable {P G} in
theorem compBiratRight_mk (F : FrobenioidCore P) {A B E : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B)
    (ψ : W.unop.left.obj ⟶ E) :
    compBiratRight F Z φ (HomBirat.mk W ψ)
      = HomBirat.mk (biratPullIdx F Z φ W) (biratPullAlpha F Z φ W ≫ ψ) := by
  show colimit.desc (homFunctorBirat P G B E) (compBiratCoconeRight F Z φ)
    (colimit.ι (homFunctorBirat P G B E) W (ULift.up ψ)) = _
  rw [← types_comp_apply (colimit.ι (homFunctorBirat P G B E) W)
    (colimit.desc (homFunctorBirat P G B E) (compBiratCoconeRight F Z φ)), colimit.ι_desc]
  rfl

/-! ## ★10. 第 1 変数についての降下 ⟹ **合成** -/

variable {P G} in
/-- ★★**第 1 変数についての自然性**。★骨は第 2 変数のときと同じ 3 本。 -/
theorem compBirat_natural_left (F : FrobenioidCore P) {A B E : C}
    {Z Z' : IdxBirat P G A} (u : Z ⟶ Z') (φ : Z.unop.left.obj ⟶ B) :
    compBiratRight (E := E) F Z' (u.unop.left.hom ≫ φ) = compBiratRight F Z φ := by
  funext g
  refine HomColim.induction (P := fun g =>
    compBiratRight (E := E) F Z' (u.unop.left.hom ≫ φ) g = compBiratRight F Z φ g)
    (homFunctorBirat P G B E) (fun W x => ?_) g
  show compBiratRight F Z' (u.unop.left.hom ≫ φ) (HomBirat.mk W x.down)
    = compBiratRight F Z φ (HomBirat.mk W x.down)
  rw [compBiratRight_mk, compBiratRight_mk]
  haveI hb : Mono W.unop.hom.hom := G.core.preStepMono _ W.unop.hom.property.2
  set V := IsFiltered.max (biratPullIdx F Z' (u.unop.left.hom ≫ φ) W) (biratPullIdx F Z φ W)
    with hV
  set c := IsFiltered.leftToMax (biratPullIdx F Z' (u.unop.left.hom ≫ φ) W)
    (biratPullIdx F Z φ W) with hc0
  set c' := IsFiltered.rightToMax (biratPullIdx F Z' (u.unop.left.hom ≫ φ) W)
    (biratPullIdx F Z φ W) with hc0'
  refine HomBirat.sound V c c' ?_
  have hsq1 : (biratPullGamma F Z' (u.unop.left.hom ≫ φ) W ≫ u.unop.left.hom) ≫ φ
      = biratPullAlpha F Z' (u.unop.left.hom ≫ φ) W ≫ W.unop.hom.hom := by
    rw [Category.assoc]
    exact biratPull_sq F Z' (u.unop.left.hom ≫ φ) W
  have hc : c.unop.left.hom ≫ (biratPullGamma F Z' (u.unop.left.hom ≫ φ) W ≫ u.unop.left.hom)
      = c'.unop.left.hom ≫ biratPullGamma F Z φ W := by
    have h := idxBirat_left_ext (u ≫ biratPullHom F Z' (u.unop.left.hom ≫ φ) W ≫ c)
      (biratPullHom F Z φ W ≫ c')
    exact (Category.assoc _ _ _).symm.trans h
  have key := birat_lift_unique φ hb hsq1 (biratPull_sq F Z φ W) hc
  simp only [← Category.assoc] at key ⊢
  exact congrArg (fun t => t ≫ x.down) key

/-- ★第 1 変数の余錐。 -/
noncomputable def compBiratCoconeLeft (F : FrobenioidCore P) (A B E : C) :
    Cocone (homFunctorBirat P G A B) :=
  Cocone.mk (HomBirat P G B E → HomBirat P G A E)
    { app := fun Z => TypeCat.ofHom fun φ : ULift.{max u2 v} (Z.unop.left.obj ⟶ B) =>
        compBiratRight F Z φ.down
      naturality := fun Z Z' u => by
        ext φ
        exact compBirat_natural_left F u φ.down }

/-- ★★★**[FrdI] Proposition 4.4, (i)** —— **`Hom^birat` の合成**。 -/
noncomputable def compBirat (F : FrobenioidCore P) {A B E : C}
    (f : HomBirat P G A B) (g : HomBirat P G B E) : HomBirat P G A E :=
  colimit.desc (homFunctorBirat P G A B) (compBiratCoconeLeft P G F A B E) f g

variable {P G} in
/-- ★★★**合成の計算則** —— 代表元では「引き戻して繋ぐ」。 -/
theorem compBirat_mk (F : FrobenioidCore P) {A B E : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B)
    (ψ : W.unop.left.obj ⟶ E) :
    compBirat P G F (HomBirat.mk Z φ) (HomBirat.mk W ψ)
      = HomBirat.mk (biratPullIdx F Z φ W) (biratPullAlpha F Z φ W ≫ ψ) := by
  have h : compBirat P G F (HomBirat.mk Z φ) = compBiratRight (E := E) F Z φ := by
    show colimit.desc (homFunctorBirat P G A B) (compBiratCoconeLeft P G F A B E)
      (colimit.ι (homFunctorBirat P G A B) Z (ULift.up φ)) = _
    rw [← types_comp_apply (colimit.ι (homFunctorBirat P G A B) Z)
      (colimit.desc (homFunctorBirat P G A B) (compBiratCoconeLeft P G F A B E)),
      colimit.ι_desc]
    rfl
  rw [h, compBiratRight_mk]

/-! ## ★11. 単位律 -/

variable {P G} in
/-- ★★**右単位律**。★引き戻しの四角形が `α = γ ≫ φ` を与え、
`mk_map`(前合成で不変)で `mk Z φ` に戻る。 -/
theorem compBirat_id_right (F : FrobenioidCore P) {A B : C} (f : HomBirat P G A B) :
    compBirat P G F f (toHomBirat (P := P) (G := G) (𝟙 B)) = f := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  show compBirat P G F (HomBirat.mk Z φ) (HomBirat.mk (idxBiratOne P G B) (𝟙 B)) = _
  rw [compBirat_mk]
  refine Eq.trans (congrArg (HomBirat.mk (biratPullIdx F Z φ (idxBiratOne P G B)))
    (biratPull_sq F Z φ (idxBiratOne P G B)).symm) ?_
  exact HomBirat.mk_map (biratPullHom F Z φ (idxBiratOne P G B)) φ

variable {P G} in
/-- ★★**左単位律**。★構造射が **mono** であることで、共通の下界の上で一致する。 -/
theorem compBirat_id_left (F : FrobenioidCore P) {A B : C} (f : HomBirat P G A B) :
    compBirat P G F (toHomBirat (P := P) (G := G) (𝟙 A)) f = f := by
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep f
  show compBirat P G F (HomBirat.mk (idxBiratOne P G A) (𝟙 A)) (HomBirat.mk W ψ) = _
  rw [compBirat_mk]
  haveI hb : Mono W.unop.hom.hom := G.core.preStepMono _ W.unop.hom.property.2
  set V := IsFiltered.max (biratPullIdx F (idxBiratOne P G A) (𝟙 A) W) W with hV
  set c := IsFiltered.leftToMax (biratPullIdx F (idxBiratOne P G A) (𝟙 A) W) W with hc0
  set c' := IsFiltered.rightToMax (biratPullIdx F (idxBiratOne P G A) (𝟙 A) W) W with hc0'
  refine HomBirat.sound V c c' ?_
  have hcw : c.unop.left.hom ≫ (biratPullIdx F (idxBiratOne P G A) (𝟙 A) W).unop.hom.hom
      = V.unop.hom.hom :=
    congrArg (fun t : V.unop.left ⟶ coaPreObj P G A => t.hom) (Over.w c.unop)
  have hc'w : c'.unop.left.hom ≫ W.unop.hom.hom = V.unop.hom.hom :=
    congrArg (fun t : V.unop.left ⟶ coaPreObj P G A => t.hom) (Over.w c'.unop)
  have hkey : c.unop.left.hom ≫ biratPullAlpha F (idxBiratOne P G A) (𝟙 A) W
      = c'.unop.left.hom := by
    refine (cancel_mono W.unop.hom.hom).mp ?_
    refine Eq.trans (Category.assoc _ _ _) ?_
    refine Eq.trans (congrArg (fun t => c.unop.left.hom ≫ t)
      (biratPull_sq F (idxBiratOne P G A) (𝟙 A) W).symm) ?_
    exact hcw.trans hc'w.symm
  refine Eq.trans (Category.assoc _ _ _).symm ?_
  exact congrArg (fun t => t ≫ ψ) hkey

/-! ## ★12. 結合律

★★**3 層の消約**で出る:
1. `Z` 側 —— 添字圏が**細い**(`idxBirat_left_ext`)
2. `W` 側 —— `W` の構造射が **mono**
3. `Y` 側 —— `Y` の構造射が **mono**

★★**`pf` の 4 つ組の添字圏に相当するものは、ここでは要らない。** -/

variable {P G} in
/-- ★★★**結合律**。 -/
theorem compBirat_assoc (F : FrobenioidCore P) {A B E K : C} (f : HomBirat P G A B)
    (g : HomBirat P G B E) (h : HomBirat P G E K) :
    compBirat P G F (compBirat P G F f g) h = compBirat P G F f (compBirat P G F g h) := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  obtain ⟨Y, χ, rfl⟩ := HomBirat.exists_rep h
  haveI hbW : Mono W.unop.hom.hom := G.core.preStepMono _ W.unop.hom.property.2
  haveI hbY : Mono Y.unop.hom.hom := G.core.preStepMono _ Y.unop.hom.property.2
  have hL : compBirat P G F (compBirat P G F (HomBirat.mk Z φ) (HomBirat.mk W ψ))
        (HomBirat.mk Y χ)
      = HomBirat.mk (biratPullIdx F (biratPullIdx F Z φ W) (biratPullAlpha F Z φ W ≫ ψ) Y)
        (biratPullAlpha F (biratPullIdx F Z φ W) (biratPullAlpha F Z φ W ≫ ψ) Y ≫ χ) := by
    rw [compBirat_mk, compBirat_mk]
  have hR : compBirat P G F (HomBirat.mk Z φ)
        (compBirat P G F (HomBirat.mk W ψ) (HomBirat.mk Y χ))
      = HomBirat.mk (biratPullIdx F Z φ (biratPullIdx F W ψ Y))
        (biratPullAlpha F Z φ (biratPullIdx F W ψ Y) ≫ biratPullAlpha F W ψ Y ≫ χ) := by
    rw [compBirat_mk, compBirat_mk]
  rw [hL, hR]
  set IL := biratPullIdx F (biratPullIdx F Z φ W) (biratPullAlpha F Z φ W ≫ ψ) Y with hIL
  set IR := biratPullIdx F Z φ (biratPullIdx F W ψ Y) with hIR
  set V := IsFiltered.max IL IR with hV
  set c := IsFiltered.leftToMax IL IR with hc0
  set c' := IsFiltered.rightToMax IL IR with hc0'
  refine HomBirat.sound V c c' ?_
  -- 第 1 層: `Z` 側 —— 添字圏の細さ
  have h1 : c.unop.left.hom ≫ (biratPullGamma F (biratPullIdx F Z φ W)
        (biratPullAlpha F Z φ W ≫ ψ) Y ≫ biratPullGamma F Z φ W)
      = c'.unop.left.hom ≫ biratPullGamma F Z φ (biratPullIdx F W ψ Y) := by
    have hh := idxBirat_left_ext
      (biratPullHom F Z φ W ≫ biratPullHom F (biratPullIdx F Z φ W)
        (biratPullAlpha F Z φ W ≫ ψ) Y ≫ c)
      (biratPullHom F Z φ (biratPullIdx F W ψ Y) ≫ c')
    exact (Category.assoc _ _ _).symm.trans hh
  -- 第 2 層: `W` 側 —— `W` の構造射が mono
  have h2 : c.unop.left.hom ≫ (biratPullGamma F (biratPullIdx F Z φ W)
        (biratPullAlpha F Z φ W ≫ ψ) Y ≫ biratPullAlpha F Z φ W)
      = c'.unop.left.hom ≫ (biratPullAlpha F Z φ (biratPullIdx F W ψ Y)
        ≫ biratPullGamma F W ψ Y) := by
    refine (cancel_mono W.unop.hom.hom).mp ?_
    have e1 : (c.unop.left.hom ≫ (biratPullGamma F (biratPullIdx F Z φ W)
          (biratPullAlpha F Z φ W ≫ ψ) Y ≫ biratPullAlpha F Z φ W)) ≫ W.unop.hom.hom
        = c.unop.left.hom ≫ (biratPullGamma F (biratPullIdx F Z φ W)
          (biratPullAlpha F Z φ W ≫ ψ) Y ≫ biratPullGamma F Z φ W) ≫ φ := by
      simp only [Category.assoc]
      exact congrArg (fun t => c.unop.left.hom ≫ biratPullGamma F (biratPullIdx F Z φ W)
        (biratPullAlpha F Z φ W ≫ ψ) Y ≫ t) (biratPull_sq F Z φ W).symm
    have e2 : (c'.unop.left.hom ≫ (biratPullAlpha F Z φ (biratPullIdx F W ψ Y)
          ≫ biratPullGamma F W ψ Y)) ≫ W.unop.hom.hom
        = c'.unop.left.hom ≫ biratPullGamma F Z φ (biratPullIdx F W ψ Y) ≫ φ := by
      simp only [Category.assoc]
      exact congrArg (fun t => c'.unop.left.hom ≫ t)
        (biratPull_sq F Z φ (biratPullIdx F W ψ Y)).symm
    refine e1.trans (Eq.trans ?_ e2.symm)
    simp only [← Category.assoc] at h1 ⊢
    exact congrArg (fun t => t ≫ φ) h1
  -- 第 3 層: `Y` 側 —— `Y` の構造射が mono
  have h3 : c.unop.left.hom ≫ biratPullAlpha F (biratPullIdx F Z φ W)
        (biratPullAlpha F Z φ W ≫ ψ) Y
      = c'.unop.left.hom ≫ (biratPullAlpha F Z φ (biratPullIdx F W ψ Y)
        ≫ biratPullAlpha F W ψ Y) := by
    refine (cancel_mono Y.unop.hom.hom).mp ?_
    have e1 : (c.unop.left.hom ≫ biratPullAlpha F (biratPullIdx F Z φ W)
          (biratPullAlpha F Z φ W ≫ ψ) Y) ≫ Y.unop.hom.hom
        = c.unop.left.hom ≫ (biratPullGamma F (biratPullIdx F Z φ W)
          (biratPullAlpha F Z φ W ≫ ψ) Y ≫ biratPullAlpha F Z φ W) ≫ ψ := by
      simp only [Category.assoc]
      exact congrArg (fun t => c.unop.left.hom ≫ t)
        (biratPull_sq F (biratPullIdx F Z φ W) (biratPullAlpha F Z φ W ≫ ψ) Y).symm
    have e2 : (c'.unop.left.hom ≫ (biratPullAlpha F Z φ (biratPullIdx F W ψ Y)
          ≫ biratPullAlpha F W ψ Y)) ≫ Y.unop.hom.hom
        = c'.unop.left.hom ≫ (biratPullAlpha F Z φ (biratPullIdx F W ψ Y)
          ≫ biratPullGamma F W ψ Y) ≫ ψ := by
      simp only [Category.assoc]
      exact congrArg (fun t => c'.unop.left.hom
        ≫ biratPullAlpha F Z φ (biratPullIdx F W ψ Y) ≫ t) (biratPull_sq F W ψ Y).symm
    refine e1.trans (Eq.trans ?_ e2.symm)
    simp only [← Category.assoc] at h2 ⊢
    exact congrArg (fun t => t ≫ ψ) h2
  simp only [← Category.assoc] at h3 ⊢
  exact congrArg (fun t => t ≫ χ) h3

/-! ## ★13. 圏 `𝒞^birat` と関手 `𝒞 → 𝒞^birat` -/

/-- ★★★**[FrdI] Proposition 4.4, (i)** —— **圏 `𝒞^birat`**。

原文:
> hence a category Cbirat, whose objects are the objects of

★対象は `𝒞` の対象、射は `Hom^birat`。 -/
def BiratCat (_P : PreFrobenioid C Φ) (_G : Frobenioid _P) : Type u2 := C

/-- ★`BiratCat` の対象を `𝒞` の対象として見る。 -/
def biratDown (A : BiratCat P G) : C := A

/-- ★★★**`𝒞^birat` の圏構造**。 -/
noncomputable instance biratCategory : Category.{max v u2 v2} (BiratCat P G) where
  Hom A B := HomBirat P G (biratDown P G A) (biratDown P G B)
  id A := toHomBirat (𝟙 (biratDown P G A))
  comp f g := compBirat P G G.core f g
  id_comp := compBirat_id_left G.core
  comp_id := compBirat_id_right G.core
  assoc f g h := compBirat_assoc G.core f g h

variable {P G} in
/-- ★★**`toHomBirat` は合成を保つ**。 -/
theorem toHomBirat_comp (F : FrobenioidCore P) {A B E : C} (f : A ⟶ B) (g : B ⟶ E) :
    toHomBirat (P := P) (G := G) (f ≫ g)
      = compBirat P G F (toHomBirat (P := P) (G := G) f) (toHomBirat (P := P) (G := G) g) := by
  show HomBirat.mk (idxBiratOne P G A) (f ≫ g)
    = compBirat P G F (HomBirat.mk (idxBiratOne P G A) f) (HomBirat.mk (idxBiratOne P G B) g)
  rw [compBirat_mk]
  have hα : biratPullAlpha F (idxBiratOne P G A) f (idxBiratOne P G B)
      = (biratPullHom F (idxBiratOne P G A) f (idxBiratOne P G B)).unop.left.hom ≫ f :=
    (Category.comp_id _).symm.trans
      (biratPull_sq F (idxBiratOne P G A) f (idxBiratOne P G B)).symm
  refine Eq.trans (HomBirat.mk_map
    (biratPullHom F (idxBiratOne P G A) f (idxBiratOne P G B)) (f ≫ g)).symm ?_
  refine congrArg (HomBirat.mk (biratPullIdx F (idxBiratOne P G A) f (idxBiratOne P G B))) ?_
  refine Eq.trans (Category.assoc _ _ _).symm ?_
  exact congrArg (fun t => t ≫ g) hα.symm

/-- ★★★**[FrdI] Proposition 4.4, (i)** —— **自然な関手 `𝒞 → 𝒞^birat`**。 -/
noncomputable def toBiratCat : C ⥤ BiratCat P G where
  obj A := A
  map f := toHomBirat (P := P) (G := G) f
  map_id _ := rfl
  map_comp f g := toHomBirat_comp G.core f g

/-! ## ★14. `0_𝒟` —— すべての値が 1 元単系である `𝒟` 上の単系

原文 (FrdI p.83):
> where the functors between elementary Frobenioids are those induced by the natural

★★原文は `0_𝒟` を「the monoid on `D` all of whose values on objects of `D` are equal to
the monoid with one element」と定義し、`𝔽_{0_𝒟}` を
「the product category of `D` with the one-object category determined by the monoid `N≥1`」
と説明する。★**次数だけを覚える単系**である。

★★**在庫確認(2026-08-17)**: `Φ^gp` は `ElementaryFrobenioid.lean` の `gpOn` にあるが、
**`0_𝒟` は無かった**ので、ここで置く。
★条件 (a)(b) はどちらも**台が 1 元**であることから自明
((a) の第 2 成分は `charMap_injective_of_addGroup` —— 1 元単系は群だから)。 -/

variable (D) in
/-- ★`0_𝒟` の台となる反変関手(値はすべて 1 元単系)。 -/
def trivialFunctor : Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj _ := AddCommMonCat.of (PUnit.{w + 1})
  map _ := AddCommMonCat.ofHom (0 : PUnit.{w + 1} →+ PUnit.{w + 1})
  map_id _ := AddCommMonCat.ext (fun _ => Subsingleton.elim _ _)
  map_comp _ _ := AddCommMonCat.ext (fun _ => Subsingleton.elim _ _)

variable (D) in
/-- ★★**`0_𝒟`** —— すべての値が 1 元単系である `𝒟` 上の単系。 -/
def trivialOn : MonoidOn.{v, u, w} D where
  functor := trivialFunctor D
  charInj _ :=
    ⟨fun _ _ _ => Subsingleton.elim (α := PUnit.{w + 1}) _ _,
      charMap_injective_of_addGroup (G := PUnit.{w + 1}) (H := PUnit.{w + 1}) _⟩
  fsmIso _ _ :=
    ⟨fun _ _ _ => Subsingleton.elim (α := PUnit.{w + 1}) _ _,
      fun _ => ⟨(PUnit.unit : PUnit.{w + 1}), Subsingleton.elim (α := PUnit.{w + 1}) _ _⟩⟩

/-! ## ★15. `𝒞^birat` の「次数」と「底」

★★**`𝔽_{0_𝒟}` へ落とすのに必要なのは次数と底の 2 つだけ**である
(零因子は 1 元単系なので情報を持たない)。

★★**どちらも代表元の取り替えで不変**:
- **次数** —— 遷移射は pre-step(次数 1)なので `degFr` が変わらない
- **底** —— 遷移射は base-isomorphism なので `Base(a)⁻¹ ≫ Base(φ)` が変わらない -/

/-- ★次数の余錐。 -/
noncomputable def biratDegCocone (A B : C) : Cocone (homFunctorBirat P G A B) :=
  Cocone.mk (ULift.{max v u2 v2} ℕ+)
    { app := fun Z => TypeCat.ofHom fun φ : ULift.{max u2 v} (Z.unop.left.obj ⟶ B) =>
        (ULift.up (P.degFr φ.down) : ULift.{max v u2 v2} ℕ+)
      naturality := fun Z W u => by
        ext φ
        show ULift.up (P.degFr (u.unop.left.hom ≫ φ.down)) = ULift.up (P.degFr φ.down)
        rw [P.degFr_comp,
          show P.degFr u.unop.left.hom = 1 from u.unop.left.property.2.1, mul_one] }

variable {P G} in
/-- ★★**`Hom^birat` の Frobenius 次数**。 -/
noncomputable def biratDeg {A B : C} (f : HomBirat P G A B) : ℕ+ :=
  (colimit.desc (homFunctorBirat P G A B) (biratDegCocone P G A B) f).down

variable {P G} in
@[simp] theorem biratDeg_mk {A B : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) :
    biratDeg (HomBirat.mk Z φ) = P.degFr φ := by
  show (colimit.desc (homFunctorBirat P G A B) (biratDegCocone P G A B)
    (colimit.ι (homFunctorBirat P G A B) Z (ULift.up φ))).down = _
  rw [← types_comp_apply (colimit.ι (homFunctorBirat P G A B) Z)
    (colimit.desc (homFunctorBirat P G A B) (biratDegCocone P G A B)), colimit.ι_desc]
  rfl

/-! ## ★16. `Hom^birat` の「底」

★★**測定(2026-08-17、詰まってから抜けた道)**: 底の余錐には `inv (Base a)` が要り、
`IsIso (Base a)` を入れる必要がある。★**`Z.unop.hom.hom` の形で `haveI` すると
インスタンス探索が失敗する**(スライス側の `sliceInv` では通るのに)。
★★**逃げ道は「`IsIso` を仮引数に取る形に定義を割る」** ——
`Def31Pf.lean` で `F` を型に載せて `instance` を回復したのと同じ発想である。 -/

variable {P G} in
/-- ★★**底の射**(`IsIso` を**仮引数**に取る形)。★`IsIso` は `Prop` なので値には効かない。 -/
noncomputable def sliceBaseOf {A B A' : C} (a : A' ⟶ A) (_ha : IsIso (P.Base a))
    (φ : A' ⟶ B) : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base :=
  haveI := _ha
  inv (P.Base a) ≫ P.Base φ

variable {P G} in
theorem sliceBaseOf_eq {A B A' : C} (a : A' ⟶ A) (ha : IsIso (P.Base a)) (φ : A' ⟶ B) :
    sliceBaseOf (P := P) a ha φ = haveI := ha; inv (P.Base a) ≫ P.Base φ := rfl

variable {P G} in
/-- ★`IsIso` は `Prop` なので、射が等しければ値も等しい。 -/
theorem sliceBaseOf_congr {A B A' : C} {a a' : A' ⟶ A} (h : a = a') (ha : IsIso (P.Base a))
    (ha' : IsIso (P.Base a')) (φ : A' ⟶ B) :
    sliceBaseOf (P := P) a ha φ = sliceBaseOf (P := P) a' ha' φ := by
  subst h; rfl

variable {P G} in
/-- ★★**底は前合成で不変** —— これが代表元の取り替えに対する不変性。 -/
theorem sliceBaseOf_precomp {A B A' A'' : C} (a : A' ⟶ A) (ha : IsIso (P.Base a))
    (c : A'' ⟶ A') (_hc : IsIso (P.Base c)) (hac : IsIso (P.Base (c ≫ a))) (φ : A' ⟶ B) :
    sliceBaseOf (P := P) (c ≫ a) hac (c ≫ φ) = sliceBaseOf (P := P) a ha φ := by
  haveI := ha
  haveI := hac
  show inv (P.Base (c ≫ a)) ≫ P.Base (c ≫ φ) = inv (P.Base a) ≫ P.Base φ
  rw [IsIso.inv_comp_eq, P.Base_comp, P.Base_comp, Category.assoc,
    ← Category.assoc (P.Base a), IsIso.hom_inv_id, Category.id_comp]

/-- ★底の余錐。 -/
noncomputable def biratBaseCocone (A B : C) : Cocone (homFunctorBirat P G A B) :=
  Cocone.mk (ULift.{max u2 v2} ((P.toElem.obj A).base ⟶ (P.toElem.obj B).base))
    { app := fun Z => TypeCat.ofHom fun φ : ULift.{max u2 v} (Z.unop.left.obj ⟶ B) =>
        (ULift.up (sliceBaseOf (P := P) Z.unop.hom.hom Z.unop.hom.property.2.2 φ.down) :
          ULift.{max u2 v2} ((P.toElem.obj A).base ⟶ (P.toElem.obj B).base))
      naturality := fun Z W u => by
        ext φ
        show ULift.up (sliceBaseOf (P := P) W.unop.hom.hom W.unop.hom.property.2.2
            (u.unop.left.hom ≫ φ.down))
          = ULift.up (sliceBaseOf (P := P) Z.unop.hom.hom Z.unop.hom.property.2.2 φ.down)
        refine congrArg ULift.up ?_
        have hw : u.unop.left.hom ≫ Z.unop.hom.hom = W.unop.hom.hom :=
          congrArg (fun t : W.unop.left ⟶ coaPreObj P G A => t.hom) (Over.w u.unop)
        refine Eq.trans (sliceBaseOf_congr (P := P) hw.symm W.unop.hom.property.2.2 ?_ _) ?_
        · exact hw ▸ W.unop.hom.property.2.2
        · exact sliceBaseOf_precomp (P := P) Z.unop.hom.hom Z.unop.hom.property.2.2
            u.unop.left.hom u.unop.left.property.2.2 _ φ.down }

variable {P G} in
/-- ★★**`Hom^birat` の底**。 -/
noncomputable def biratBase {A B : C} (f : HomBirat P G A B) :
    (P.toElem.obj A).base ⟶ (P.toElem.obj B).base :=
  (colimit.desc (homFunctorBirat P G A B) (biratBaseCocone P G A B) f).down

variable {P G} in
@[simp] theorem biratBase_mk {A B : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) :
    biratBase (HomBirat.mk Z φ)
      = sliceBaseOf (P := P) Z.unop.hom.hom Z.unop.hom.property.2.2 φ := by
  show (colimit.desc (homFunctorBirat P G A B) (biratBaseCocone P G A B)
    (colimit.ι (homFunctorBirat P G A B) Z (ULift.up φ))).down = _
  rw [← types_comp_apply (colimit.ι (homFunctorBirat P G A B) Z)
    (colimit.desc (homFunctorBirat P G A B) (biratBaseCocone P G A B)), colimit.ι_desc]
  rfl

/-! ## ★17. 次数と底は合成を保つ -/

variable {P G} in
/-- ★★**引き戻した射の次数は元の射の次数**。★`γ` と `b` が linear だから。 -/
theorem biratPullAlpha_degFr (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    P.degFr (biratPullAlpha F Z φ W) = P.degFr φ := by
  have hγ : P.degFr (biratPullGamma F Z φ W) = 1 := (biratPullGamma_preStep F Z φ W).1
  have hb : P.degFr W.unop.hom.hom = 1 := W.unop.hom.property.2.1
  have h := congrArg P.degFr (biratPull_sq F Z φ W)
  have hl : P.degFr (biratPullGamma F Z φ W ≫ φ) = P.degFr φ := by
    rw [P.degFr_comp, hγ, mul_one]
  have hr : P.degFr (biratPullAlpha F Z φ W ≫ W.unop.hom.hom)
      = P.degFr (biratPullAlpha F Z φ W) := by
    rw [P.degFr_comp, hb, one_mul]
  exact hr.symm.trans (h.symm.trans hl)

variable {P G} in
/-- ★★**次数は合成を保つ**(`𝔽` の合成規則 `deg (φ ≫ ψ) = deg ψ * deg φ` に一致)。 -/
theorem biratDeg_comp (F : FrobenioidCore P) {A B E : C} (f : HomBirat P G A B)
    (g : HomBirat P G B E) :
    biratDeg (compBirat P G F f g) = biratDeg g * biratDeg f := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  rw [compBirat_mk, biratDeg_mk, biratDeg_mk, biratDeg_mk]
  exact (P.degFr_comp _ _).trans
    (congrArg (fun t => P.degFr ψ * t) (biratPullAlpha_degFr F Z φ W))

variable {P G} in
/-- ★★**底の合成則**(すべて素の射で述べた形)。★引き戻しの四角形の `Base` を取って移項するだけ。 -/
theorem sliceBaseOf_comp {A B E A' B' Dd : C}
    (a : A' ⟶ A) (ha : IsIso (P.Base a)) (φ : A' ⟶ B)
    (b : B' ⟶ B) (hb : IsIso (P.Base b)) (ψ : B' ⟶ E)
    (γ : Dd ⟶ A') (_hγ : IsIso (P.Base γ)) (α : Dd ⟶ B')
    (hγa : IsIso (P.Base (γ ≫ a))) (hsq : γ ≫ φ = α ≫ b) :
    sliceBaseOf (P := P) (γ ≫ a) hγa (α ≫ ψ)
      = sliceBaseOf (P := P) a ha φ ≫ sliceBaseOf (P := P) b hb ψ := by
  haveI := ha
  haveI := hb
  haveI := hγa
  show inv (P.Base (γ ≫ a)) ≫ P.Base (α ≫ ψ)
    = (inv (P.Base a) ≫ P.Base φ) ≫ (inv (P.Base b) ≫ P.Base ψ)
  refine (IsIso.inv_comp_eq _).mpr ?_
  have hsq' := congrArg P.Base hsq
  simp only [P.Base_comp] at hsq'
  simp only [P.Base_comp, Category.assoc, IsIso.hom_inv_id_assoc]
  rw [← Category.assoc (P.Base γ), hsq', Category.assoc, IsIso.hom_inv_id_assoc]

variable {P G} in
/-- ★★**底は合成を保つ**。 -/
theorem biratBase_comp (F : FrobenioidCore P) {A B E : C} (f : HomBirat P G A B)
    (g : HomBirat P G B E) :
    biratBase (compBirat P G F f g) = biratBase f ≫ biratBase g := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  rw [compBirat_mk, biratBase_mk, biratBase_mk, biratBase_mk]
  exact sliceBaseOf_comp (P := P) Z.unop.hom.hom Z.unop.hom.property.2.2 φ
    W.unop.hom.hom W.unop.hom.property.2.2 ψ
    (biratPullGamma F Z φ W) (biratPullGamma_preStep F Z φ W).2
    (biratPullAlpha F Z φ W) _ (biratPull_sq F Z φ W)

variable {P G} in
@[simp] theorem biratDeg_toHomBirat {A B : C} (φ : A ⟶ B) :
    biratDeg (toHomBirat (P := P) (G := G) φ) = P.degFr φ :=
  biratDeg_mk (idxBiratOne P G A) φ

variable {P G} in
@[simp] theorem biratBase_toHomBirat {A B : C} (φ : A ⟶ B) :
    biratBase (toHomBirat (P := P) (G := G) φ) = P.Base φ := by
  rw [toHomBirat, biratBase_mk]
  show sliceBaseOf (P := P) (𝟙 A) _ φ = _
  haveI : IsIso (P.Base (𝟙 A)) := by rw [P.Base_id]; infer_instance
  rw [sliceBaseOf_eq]
  show inv (P.Base (𝟙 A)) ≫ P.Base φ = P.Base φ
  refine (IsIso.inv_comp_eq _).mpr ?_
  rw [P.Base_id, Category.id_comp]

/-! ## ★18. 関手 `𝒞^birat → 𝔽_{0_𝒟}`

★★原文 (i) の 1-可換図式の**右下の辺**である。
零因子は 1 元単系なので情報を持たず、**底と次数だけ**を覚える。 -/

/-- ★★★**[FrdI] Proposition 4.4, (i)** —— **関手 `𝒞^birat → 𝔽_{0_𝒟}`**。 -/
noncomputable def biratToElem : BiratCat P G ⥤ ElemFrobCat (trivialOn D) where
  obj A := ⟨(P.toElem.obj (biratDown P G A)).base⟩
  map f := ⟨biratBase f, PUnit.unit, biratDeg f⟩
  map_id A := by
    refine ElemFrobCat.Hom.ext ?_ (Subsingleton.elim (α := PUnit.{w + 1}) _ _) ?_
    · show biratBase (toHomBirat (P := P) (G := G) (𝟙 (biratDown P G A))) = 𝟙 _
      rw [biratBase_toHomBirat, P.Base_id]
    · show biratDeg (toHomBirat (P := P) (G := G) (𝟙 (biratDown P G A))) = 1
      rw [biratDeg_toHomBirat, P.degFr_id]
  map_comp f g :=
    ElemFrobCat.Hom.ext (biratBase_comp G.core f g) (Subsingleton.elim (α := PUnit.{w + 1}) _ _)
      (biratDeg_comp G.core f g)

/-! ## ★19. `𝔽_Φ → 𝔽_{Φ^gp} → 𝔽_{0_𝒟}` と 1-可換図式

原文 (FrdI p.83):
> where the functors between elementary Frobenioids are those induced by the natural

★★原文の図式は
```
𝒞        → 𝔽_Φ
↓            ↓
𝒞^birat  → 𝔽_{Φ^gp} → 𝔽_{0_𝒟}
```
であり、**単系の準同型 `Φ → Φ^gp → 0_𝒟` が誘導する関手**でつながっている。 -/

include P in
/-- ★`Φ` が divisorial なら各 `Φ(A)` は integral(`gpOn` の仮定)。 -/
theorem phi_integral (A : D) : IsIntegralMonoid (Φ.val A) := (P.divisorial A).1.1

/-- ★★**`Φ → Φ^gp`** の自然変換。 -/
noncomputable def toGpNat : Φ.functor ⟶ (Φ.gpOn (phi_integral P)).functor where
  app X := AddCommMonCat.ofHom (Algebra.GrothendieckAddGroup.of)
  naturality X Y f := by
    refine AddCommMonCat.ext (fun x => ?_)
    exact (gpMap_toGp _ (Φ.map f.unop) x).symm

/-- ★★**`Φ^gp → 0_𝒟`** の自然変換(値が 1 元なので零写像)。 -/
def gpToTrivialNat : (Φ.gpOn (phi_integral P)).functor ⟶ trivialFunctor D where
  app _ := AddCommMonCat.ofHom (0 : _ →+ PUnit.{w + 1})
  naturality _ _ _ :=
    AddCommMonCat.ext (fun _ => Subsingleton.elim (α := PUnit.{w + 1}) _ _)

/-- ★★**`𝔽_Φ → 𝔽_{Φ^gp}`**。 -/
noncomputable def elemToGp : ElemFrobCat Φ ⥤ ElemFrobCat (Φ.gpOn (phi_integral P)) :=
  ElemFrobCat.elemFrobMap (toGpNat P)

/-- ★★**`𝔽_{Φ^gp} → 𝔽_{0_𝒟}`**。 -/
noncomputable def gpToElemTrivial :
    ElemFrobCat (Φ.gpOn (phi_integral P)) ⥤ ElemFrobCat (trivialOn D) :=
  ElemFrobCat.elemFrobMap (gpToTrivialNat P)

/-- ★★★**[FrdI] Proposition 4.4, (i)** —— **1-可換図式**。

★★**「1-可換」どころか、我々の実装では図式は**恒等的に可換**である ——
`𝒞^birat → 𝔽_{0_𝒟}` は代表元の**底と次数**しか見ず、それは `𝒞 → 𝔽_Φ` の
底と次数に一致するから(`biratBase_toHomBirat` / `biratDeg_toHomBirat`)。 -/
theorem birat_square_obj (A : C) :
    (toBiratCat P G ⋙ biratToElem P G).obj A
      = (P.toElem ⋙ elemToGp P ⋙ gpToElemTrivial P).obj A := rfl

variable {P G} in
/-- ★★★**図式の可換性(射の部分)** —— 底も次数も一致する。 -/
theorem birat_square_map {A B : C} (f : A ⟶ B) :
    (toBiratCat P G ⋙ biratToElem P G).map f
      = (P.toElem ⋙ elemToGp P ⋙ gpToElemTrivial P).map f :=
  ElemFrobCat.Hom.ext (biratBase_toHomBirat f)
    (Subsingleton.elim (α := PUnit.{w + 1}) _ _) (biratDeg_toHomBirat f)

/-! ## ★20. `Proposition 4.4` の主張と実装の対応

| 原文の主張 | 実装名 |
|---|---|
| `Hom^birat_𝒞(A,B)`(帰納極限) | `HomBirat` |
| 添字圏 `𝒞^{coa-pre}_A` とその反対圏 | `SliceA` / `IdxBirat` |
| 添字圏が細い / filtered | `slice_hom_ext` / `sliceA_isCofiltered` / `idxBirat_isFiltered` |
| (i) 合成写像 | `compBirat` / `compBirat_mk` |
| (i) 圏 `𝒞^birat` | `biratCategory` |
| (i) 関手 `𝒞 → 𝒞^birat` | `toBiratCat` |
| (i) `0_𝒟` | `trivialOn` |
| (i) 関手 `𝒞^birat → 𝔽_{0_𝒟}` | `biratToElem` |
| (i) `𝔽_Φ → 𝔽_{Φ^gp} → 𝔽_{0_𝒟}` | `elemToGp` / `gpToElemTrivial` |
| (i) 図式の可換性 | `birat_square_obj` / `birat_square_map` |
| (ii) `𝒞^birat` が group-like 型の Frobenioid | ★**未** |
| (ii) `𝒞 → 𝒞^birat` が忠実 / `𝒪^▷(A)^gp ↪ 𝒪^×(A^birat)` | ★**未** |
| (iii) 部分関手 `Φ^birat ⊆ Φ^gp` | ★**未** |

★★**測定(実装上、他所でも効くもの)**

1. ★**`inv (Base a)` を扱うときは `IsIso` を仮引数に取る形に定義を割る**
   (`sliceBaseOf`)。★反対圏側(`Z.unop.hom`)では `haveI` でのインスタンス登録が
   効かないことがある。`IsIso` は `Prop` なので、仮引数にしても値は変わらない
   (`sliceBaseOf_congr`)。
2. ★**`rw` が `≫` の等式で「そこにあるのに当たらない」**ことが頻発する。
   ★`congrArg (fun t => t ≫ ψ) h` のように**項で書く**か、
   `simp only [← Category.assoc]` で正規化してから `exact` に逃がす。
3. ★**`inv` の中を `rw` すると motive が壊れる。**
   ★`IsIso.inv_comp_eq` で先に `inv` を反対側へ移してから書き換える。
4. ★**`Functor.ext` は Lean core の `Functor`(モナド)と衝突する。**
   `CategoryTheory.Functor.ext` と書く。★ただし対象が定義的に等しい場合は、
   **`obj` と `map` の成分等式で述べるほうが `eqToHom` が出ず読みやすい。**
-/

end ABC3.Found.FrdI
