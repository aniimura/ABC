import Mathlib.CategoryTheory.Filtered.Final
import ABC3.Found.FrdI.Def31
import ABC3.Found.FrdI.HomColim

/-!
# [FrdI] Definition 3.1, (ii) —— `𝒞^{Fr-tp}` / `𝒞^{bi-Fr}` と `Hom^pf`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.56。

原文 (FrdI p.56):
> (ii) Write

★★**ファイルを分けた理由**(2026-08-17 の事故から): `Definition 3.1` は
(i)(iv) と (ii)(iii) で**別のセッションが同時に書いている**。
★同じファイルを両側から書くと衝突するので、**(i)(iv) は `Def31.lean`、
(ii)(iii) はこのファイル**に置く。

## ★2 層設計(親と合意)

- **層 1**(`HomColim.lean`) —— 添字圏上の `Type` 値関手の帰納極限(汎用)。
- **層 2**(このファイル) —— `pf` 固有の添字圏と遷移写像。

★★**`Hom^pf` と `Hom^birat`(`Proposition 4.4`)は「同じ形」ではない** ——
**(1) 両側 / 片側、(2) コスライス / スライス、(3) 四角形 / 前合成**の 3 点で違う。
★共有できるのは層 1 までで、合成と圏構造は各ケース固有に書く。

## ★段取り

1. `𝒞^{Fr-tp}`(`FrTp`)—— Frobenius 型射が定める広い部分圏。
2. `𝒞^{bi-Fr}`(`BiFr`)—— **同じ Frobenius 次数**の対が定める広い部分圏。
3. コスライス `(A,B)𝒞^{bi-Fr}`(`IdxPf`)。
4. 遷移写像(`frobTransport`)—— `Proposition 1.10, (i)` の四角形。
5. 添字関手(`homFunctorPf`)。

## ★★宇宙について(測定)

★原文は添字圏を「[say, some **small skeletal** subcategory of] `(A,B)𝒞^{bi-Fr}`」
と括弧で断る。★**これは飾りではない** —— `C : Type u2`, `[Category.{v2} C]` では
コスライスの**対象は `Type (max u2 v2)`** に住み、`Hom_C(A′,B′)` は `Type v2` に住む。
★★**`Type v2` の中には、この大きさの添字圏上の余極限は一般に存在しない。**

★**我々は「小さい骨格部分圏を取る」代わりに、宇宙を上げる** ——
`ULift` で値を `Type (max u2 v2)` へ持ち上げれば `Small` から `HasColimit` が出る。
★骨格を取る選択公理的な操作が要らない分、こちらのほうが素直である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★1. `𝒞^{Fr-tp}` —— Frobenius 型射の広い部分圏 -/

include F in
/-- ★**Frobenius 型射は合成で閉じ、恒等射を含む**。

★合成は `Proposition 1.7, (i)`(`IsFrobeniusType.comp`)であり、
**`FrobenioidCore` を要求する** —— だから `instance` ではなく `theorem`。 -/
theorem frobTypeProp_isMultiplicative : (frobTypeProp P).IsMultiplicative where
  id_mem A := isFrobeniusType_of_isIso P (𝟙 A)
  comp_mem _ _ hψ hφ := IsFrobeniusType.comp P F hψ hφ

/-! ### ★★`F` を型に載せる —— `instance` を回復する技

★★**測定(2026-08-17、事故 1 件から学んだ)**: `frobTypeProp P` は `F` を含まないので
`(frobTypeProp P).IsMultiplicative` を `instance` にできない。すると
`WideSubcategory` を使うたびに `letI` が要り、★**`Under` の射を書いた瞬間に
インスタンス探索が失敗する。**

★★**逃げ道**: **`F` を引数に持つ別名**を作れば、`(frobTypePropOf P F).IsMultiplicative`
という**目標に `F` が現れる**ので `instance` にできる。
★これで `letI` が全部消え、`Under`・関手・余極限がすべて素直に通る。 -/

/-- ★`F` を型に載せた `frobTypeProp`(中身は同じ)。 -/
def frobTypePropOf (_F : FrobenioidCore P) : MorphismProperty C := frobTypeProp P

instance frobTypePropOf_isMultiplicative : (frobTypePropOf P F).IsMultiplicative :=
  frobTypeProp_isMultiplicative P F

/-- ★★**`𝒞^{Fr-tp} ⊆ 𝒞`**(`Definition 3.1, (ii)` の第 1 文)。 -/
abbrev FrTp : Type u2 := WideSubcategory (frobTypePropOf P F)

/-! ## ★2. `𝒞^{bi-Fr}` —— 次数の揃った対 -/

/-- ★★**`𝒞^{bi-Fr}` を定める `𝒞 × 𝒞` の射のクラス** ——
**両成分が Frobenius 型で、しかも Frobenius 次数が等しい**。

★★原文「pairs of morphisms of Frobenius type of the same Frobenius degree」
の直訳である。★**次数が揃っていることが `Proposition 1.10, (i)` を
呼ぶための仮定 `hd` そのもの**であり、ここが添字圏の設計の要点。 -/
def biFrProp : MorphismProperty (C × C) :=
  fun _ _ f => IsFrobeniusType P f.1 ∧ IsFrobeniusType P f.2 ∧ P.degFr f.1 = P.degFr f.2

include F in
/-- ★次数が揃う条件も合成で閉じる(`degFr` の乗法性)。 -/
theorem biFrProp_isMultiplicative : (biFrProp P).IsMultiplicative where
  id_mem X :=
    ⟨isFrobeniusType_of_isIso P (𝟙 X.1), isFrobeniusType_of_isIso P (𝟙 X.2), by
      show P.degFr (𝟙 X.1) = P.degFr (𝟙 X.2)
      rw [P.degFr_id, P.degFr_id]⟩
  comp_mem f g hf hg :=
    ⟨IsFrobeniusType.comp P F hf.1 hg.1, IsFrobeniusType.comp P F hf.2.1 hg.2.1, by
      show P.degFr (f.1 ≫ g.1) = P.degFr (f.2 ≫ g.2)
      rw [P.degFr_comp, P.degFr_comp, hf.2.2, hg.2.2]⟩

/-- ★`F` を型に載せた `biFrProp`(中身は同じ)。 -/
def biFrPropOf (_F : FrobenioidCore P) : MorphismProperty (C × C) := biFrProp P

instance biFrPropOf_isMultiplicative : (biFrPropOf P F).IsMultiplicative :=
  biFrProp_isMultiplicative P F

/-- ★★**`𝒞^{bi-Fr} ⊆ 𝒞^{Fr-tp} × 𝒞^{Fr-tp}`**(`Definition 3.1, (ii)` の第 2 文)。

★★**実装上の測定**: 原文は「`𝒞^{Fr-tp} × 𝒞^{Fr-tp}` の部分圏」と書くが、
★**`𝒞 × 𝒞` の広い部分圏として直接作るほうが等価で短い** ——
`biFrProp` が両成分の Frobenius 型を含んでいるからである。 -/
abbrev BiFr : Type u2 := WideSubcategory (biFrPropOf P F)

/-! ## ★3. コスライス `(A,B)𝒞^{bi-Fr}` —— 添字圏 -/

/-- ★添字圏の底となる対象 `(A,B) ∈ 𝒞^{bi-Fr}`。 -/
def biFrObj (A B : C) : BiFr P F := ⟨(A, B)⟩

/-- ★★**添字圏 `(A,B)𝒞^{bi-Fr}`**(`Definition 3.1, (ii)` の帰納極限の添字)。

★対象は「同次数の Frobenius 型射の対 `(A → A′, B → B′)`」、
射は「同次数の Frobenius 型射の対 `(A′ → A″, B′ → B″)` で下の三角形が可換なもの」。 -/
abbrev IdxPf (A B : C) : Type (max u2 v2) := Under (biFrObj P F A B)

/-! ## ★4. 遷移写像 —— `Proposition 1.10, (i)` の四角形 -/

variable {P F} in
/-- ★★**遷移写像**(`Definition 3.1, (ii)` の「the assignment φ → φ′」)。

`α : A′ → A″`、`β : B′ → B″` が**同次数の Frobenius 型射**のとき、
`φ : A′ → B′` に対して `φ ≫ β = α ≫ φ′` を満たす**唯一の** `φ′ : A″ → B″`。

★存在と一意性はどちらも `prop_1_10_i_exists_given`(`Proposition 1.10, (i)` の原文の形)。 -/
noncomputable def frobTransport {A' B' A'' B'' : C}
    (α : A' ⟶ A'') (hα : IsFrobeniusType P α) (β : B' ⟶ B'') (hβ : IsFrobeniusType P β)
    (hd : P.degFr α = P.degFr β) (φ : A' ⟶ B') : A'' ⟶ B'' :=
  (prop_1_10_i_exists_given P F φ α hα β hβ hd).choose

variable {P F} in
/-- ★遷移写像は四角形を可換にする。 -/
theorem frobTransport_spec {A' B' A'' B'' : C}
    (α : A' ⟶ A'') (hα : IsFrobeniusType P α) (β : B' ⟶ B'') (hβ : IsFrobeniusType P β)
    (hd : P.degFr α = P.degFr β) (φ : A' ⟶ B') :
    φ ≫ β = α ≫ frobTransport (F := F) α hα β hβ hd φ :=
  (prop_1_10_i_exists_given P F φ α hα β hβ hd).choose_spec.1

variable {P F} in
/-- ★★**一意性** —— 四角形を可換にするものは遷移写像に一致する。
★★**以下の関手性はすべてこれ 1 本から出る。** -/
theorem frobTransport_eq {A' B' A'' B'' : C}
    (α : A' ⟶ A'') (hα : IsFrobeniusType P α) (β : B' ⟶ B'') (hβ : IsFrobeniusType P β)
    (hd : P.degFr α = P.degFr β) (φ : A' ⟶ B') (ψ : A'' ⟶ B'') (h : φ ≫ β = α ≫ ψ) :
    frobTransport (F := F) α hα β hβ hd φ = ψ :=
  ((prop_1_10_i_exists_given P F φ α hα β hβ hd).choose_spec.2 ψ h).symm

variable {P F} in
/-- ★**恒等射に沿う遷移は恒等**。 -/
theorem frobTransport_id {A' B' : C} (hα : IsFrobeniusType P (𝟙 A'))
    (hβ : IsFrobeniusType P (𝟙 B')) (hd : P.degFr (𝟙 A') = P.degFr (𝟙 B')) (φ : A' ⟶ B') :
    frobTransport (F := F) (𝟙 A') hα (𝟙 B') hβ hd φ = φ :=
  frobTransport_eq _ hα _ hβ hd φ φ (by rw [Category.comp_id, Category.id_comp])

variable {P F} in
/-- ★★**合成に沿う遷移は遷移の合成**(関手性の本体)。

★証明は一意性 1 本 —— 2 つの四角形を横に並べて外側の四角形を作るだけ。 -/
theorem frobTransport_comp {A₁ B₁ A₂ B₂ A₃ B₃ : C}
    (α : A₁ ⟶ A₂) (hα : IsFrobeniusType P α) (β : B₁ ⟶ B₂) (hβ : IsFrobeniusType P β)
    (hd : P.degFr α = P.degFr β)
    (α' : A₂ ⟶ A₃) (hα' : IsFrobeniusType P α') (β' : B₂ ⟶ B₃) (hβ' : IsFrobeniusType P β')
    (hd' : P.degFr α' = P.degFr β')
    (hαc : IsFrobeniusType P (α ≫ α')) (hβc : IsFrobeniusType P (β ≫ β'))
    (hdc : P.degFr (α ≫ α') = P.degFr (β ≫ β')) (φ : A₁ ⟶ B₁) :
    frobTransport (F := F) (α ≫ α') hαc (β ≫ β') hβc hdc φ
      = frobTransport (F := F) α' hα' β' hβ' hd'
          (frobTransport (F := F) α hα β hβ hd φ) := by
  refine frobTransport_eq _ hαc _ hβc hdc φ _ ?_
  calc φ ≫ β ≫ β' = (φ ≫ β) ≫ β' := by rw [Category.assoc]
    _ = (α ≫ frobTransport (F := F) α hα β hβ hd φ) ≫ β' := by
        rw [← frobTransport_spec α hα β hβ hd φ]
    _ = α ≫ (frobTransport (F := F) α hα β hβ hd φ ≫ β') := by rw [Category.assoc]
    _ = α ≫ α' ≫ frobTransport (F := F) α' hα' β' hβ' hd'
          (frobTransport (F := F) α hα β hβ hd φ) := by
        rw [← frobTransport_spec α' hα' β' hβ' hd']
    _ = (α ≫ α') ≫ _ := by rw [Category.assoc]

/-! ## ★5. 添字関手 `(A,B)𝒞^{bi-Fr} ⥤ Type` -/

/-- ★添字圏の射が定める遷移写像(関手の `map` の中身)。 -/
noncomputable def idxTransport {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) : W.right.obj.1 ⟶ W.right.obj.2 :=
  frobTransport (F := F) u.right.hom.1 u.right.property.1 u.right.hom.2
    u.right.property.2.1 u.right.property.2.2 φ

variable {P F} in
/-- ★恒等射の遷移は恒等。 -/
theorem idxTransport_id {A B : C} (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    idxTransport P F (𝟙 Z) φ = φ :=
  frobTransport_id _ _ _ φ

variable {P F} in
/-- ★合成の遷移は遷移の合成。 -/
theorem idxTransport_comp {A B : C} {Z W V : IdxPf P F A B} (u : Z ⟶ W) (u' : W ⟶ V)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    idxTransport P F (u ≫ u') φ = idxTransport P F u' (idxTransport P F u φ) :=
  frobTransport_comp _ u.right.property.1 _ u.right.property.2.1 u.right.property.2.2
    _ u'.right.property.1 _ u'.right.property.2.1 u'.right.property.2.2 _ _ _ φ

/-- ★★**帰納系** —— `(A → A′, B → B′) ↦ Hom_𝒞(A′, B′)`、遷移写像は `frobTransport`。

★★**`ULift` を噛ませる理由**は冒頭の宇宙の測定。

★★**`TypeCat.ofHom` が要る理由(2026-08-17 の測定)**: 現行 mathlib では
**`Type u` の射は関数そのものではなく、構造体 `TypeCat.Hom` に包まれている**
(`Hom X Y` が `Fun X Y` を包む)。★したがって `X ⟶ Y` と `X → Y` は
**もはや定義的に等しくない**。★`fun φ => …` をそのまま `map` に置くと
「`→` と `⟶` の型不一致」で落ちる —— `TypeCat.ofHom`(記法 `↾`)で包む。
★逆に**関手則のほうは `ofHom` を剥がした形と定義的に等しい**ので、
`ext φ` のあと `exact` 1 本で通る。 -/
noncomputable def homFunctorPf (A B : C) : IdxPf P F A B ⥤ Type (max u2 v2) where
  obj Z := ULift.{u2} (Z.right.obj.1 ⟶ Z.right.obj.2)
  map {Z W} u := TypeCat.ofHom fun φ : ULift.{u2} (Z.right.obj.1 ⟶ Z.right.obj.2) =>
    (ULift.up (idxTransport P F u φ.down) : ULift.{u2} (W.right.obj.1 ⟶ W.right.obj.2))
  map_id Z := by
    ext φ
    exact idxTransport_id Z φ.down
  map_comp u u' := by
    ext φ
    exact idxTransport_comp u u' φ.down

@[simp] theorem homFunctorPf_obj {A B : C} (Z : IdxPf P F A B) :
    (homFunctorPf P F A B).obj Z = ULift.{u2} (Z.right.obj.1 ⟶ Z.right.obj.2) := rfl

variable {P F} in
@[simp] theorem homFunctorPf_map {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : ULift.{u2} (Z.right.obj.1 ⟶ Z.right.obj.2)) :
    (homFunctorPf P F A B).map u φ = ULift.up (idxTransport P F u φ.down) := rfl

/-! ## ★5.5. 添字圏は filtered

★★**帰納極限が「本当に帰納極限らしく振る舞う」ためには添字圏が filtered でなければ
ならない**(層 1 の `sound` / `eq_iff` / `eq_iff_same` の 3 本がそこを使う)。
★原文は「the inductive limit is parametrized by …」と述べるだけで filtered 性に
触れないが、**帰納極限と呼ぶ以上これは要る**。

★★**測定**: 2 条件のうち
- **平行 2 射の一致**は、★**`𝒞` が totally epimorphic であることから自動**
  (コスライスなので構造射が epi、それで消約できる)。**上界を取る必要すらない。**
- **共通上界**は `Definition 1.3, (ii)` の**存在と本質的一意性**の 2 つで作る。
-/

variable {P F} in
/-- ★添字圏の対象を作る補助。 -/
def idxMk {A B A' B' : C} (a : A ⟶ A') (b : B ⟶ B')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    IdxPf P F A B :=
  Under.mk (Y := (⟨(A', B')⟩ : BiFr P F))
    (show biFrObj P F A B ⟶ (⟨(A', B')⟩ : BiFr P F) from ⟨(a, b), ha, hb, hd⟩)

variable {P F} in
/-- ★★**添字圏の平行 2 射は等しい** —— `𝒞` が totally epimorphic だから。

★★**これが `cocone_maps` を自明にする。**「添字圏はほぼ順序集合」という
原文の直観(次数の可除性)の、機械的な現れである。 -/
theorem idx_hom_ext {A B : C} {Z W : IdxPf P F A B} (u v : Z ⟶ W) : u = v := by
  have h : Z.hom ≫ u.right = Z.hom ≫ v.right := by rw [Under.w u, Under.w v]
  have h1 : Z.hom.hom.1 ≫ u.right.hom.1 = Z.hom.hom.1 ≫ v.right.hom.1 :=
    congrArg (fun t : biFrObj P F A B ⟶ W.right => t.hom.1) h
  have h2 : Z.hom.hom.2 ≫ u.right.hom.2 = Z.hom.hom.2 ≫ v.right.hom.2 :=
    congrArg (fun t : biFrObj P F A B ⟶ W.right => t.hom.2) h
  haveI : Epi Z.hom.hom.1 := P.totEpiC _ _ _
  haveI : Epi Z.hom.hom.2 := P.totEpiC _ _ _
  have e1 : u.right.hom.1 = v.right.hom.1 := (cancel_epi _).mp h1
  have e2 : u.right.hom.2 = v.right.hom.2 := (cancel_epi _).mp h2
  have hr : u.right = v.right := WideSubcategory.hom_ext _ (Prod.ext e1 e2)
  exact CommaMorphism.ext (Subsingleton.elim _ _) hr

variable {P F} in
/-- ★★**共通上界** —— 2 つの添字対象は共通の上に持ち上がる。

★★**構成**: `Z` の次数を `n`、`W` の次数を `m` として、
`Z` の先から次数 `m` の Frobenius 型射 `(α, β)` を伸ばし、`V` を次数 `n·m` に取る。
★`W → V` は「`W` の先から次数 `n` を伸ばしてから、`Definition 1.3, (ii)` の
**本質的一意性**(`frobDegUniq`)の同型で `V` に合わせる」。
★★**次数の積が可換であること**が、ここでも効く(`mul_comm`)。 -/
theorem idx_cocone_objs {A B : C} (Z W : IdxPf P F A B) :
    ∃ (V : IdxPf P F A B) (_ : Z ⟶ V) (_ : W ⟶ V), True := by
  obtain ⟨ha, hb, hab⟩ := Z.hom.property
  obtain ⟨hc, hd, hcd⟩ := W.hom.property
  obtain ⟨A₃, α, hα, hαd⟩ := F.frobDegSurj Z.right.obj.1 (P.degFr W.hom.hom.1)
  obtain ⟨B₃, β, hβ, hβd⟩ := F.frobDegSurj Z.right.obj.2 (P.degFr W.hom.hom.1)
  -- `V` の構造射
  have hVd : P.degFr (Z.hom.hom.1 ≫ α) = P.degFr (Z.hom.hom.2 ≫ β) := by
    rw [P.degFr_comp, P.degFr_comp, hαd, hβd, hab]
  -- `W → V` の第 1 成分
  obtain ⟨A₄, γ₀, hγ₀, hγ₀d⟩ := F.frobDegSurj W.right.obj.1 (P.degFr Z.hom.hom.1)
  have hg1 : P.degFr (W.hom.hom.1 ≫ γ₀) = P.degFr (Z.hom.hom.1 ≫ α) := by
    rw [P.degFr_comp, P.degFr_comp, hγ₀d, hαd, mul_comm]
  obtain ⟨θ, hθ, hθe⟩ := F.frobDegUniq A A₄ A₃ (W.hom.hom.1 ≫ γ₀) (Z.hom.hom.1 ≫ α)
    (IsFrobeniusType.comp P F hc hγ₀) (IsFrobeniusType.comp P F ha hα) hg1
  haveI : IsIso θ := hθ
  -- `W → V` の第 2 成分
  obtain ⟨B₄, δ₀, hδ₀, hδ₀d⟩ := F.frobDegSurj W.right.obj.2 (P.degFr Z.hom.hom.2)
  have hg2 : P.degFr (W.hom.hom.2 ≫ δ₀) = P.degFr (Z.hom.hom.2 ≫ β) := by
    rw [P.degFr_comp, P.degFr_comp, hδ₀d, hβd, ← hab, ← hcd, mul_comm]
  obtain ⟨ζ, hζ, hζe⟩ := F.frobDegUniq B B₄ B₃ (W.hom.hom.2 ≫ δ₀) (Z.hom.hom.2 ≫ β)
    (IsFrobeniusType.comp P F hd hδ₀) (IsFrobeniusType.comp P F hb hβ) hg2
  haveI : IsIso ζ := hζ
  have hγd : P.degFr (γ₀ ≫ θ) = P.degFr Z.hom.hom.1 := by
    rw [P.degFr_comp, show P.degFr θ = 1 from isLinear_of_isIso P θ, hγ₀d, one_mul]
  have hζd : P.degFr (δ₀ ≫ ζ) = P.degFr Z.hom.hom.2 := by
    rw [P.degFr_comp, show P.degFr ζ = 1 from isLinear_of_isIso P ζ, hδ₀d, one_mul]
  refine ⟨idxMk (Z.hom.hom.1 ≫ α) (Z.hom.hom.2 ≫ β) (IsFrobeniusType.comp P F ha hα)
      (IsFrobeniusType.comp P F hb hβ) hVd,
    Under.homMk (show Z.right ⟶ _ from
        ⟨(α, β), hα, hβ, show P.degFr α = P.degFr β by rw [hαd, hβd]⟩)
      (WideSubcategory.hom_ext _ rfl),
    Under.homMk (show W.right ⟶ _ from
      ⟨(γ₀ ≫ θ, δ₀ ≫ ζ), IsFrobeniusType.comp P F hγ₀ (isFrobeniusType_of_isIso P θ),
        IsFrobeniusType.comp P F hδ₀ (isFrobeniusType_of_isIso P ζ),
        show P.degFr (γ₀ ≫ θ) = P.degFr (δ₀ ≫ ζ) by rw [hγd, hζd, hab]⟩)
      (WideSubcategory.hom_ext _ (Prod.ext ?_ ?_)), trivial⟩
  · show W.hom.hom.1 ≫ γ₀ ≫ θ = Z.hom.hom.1 ≫ α
    rw [← Category.assoc]
    exact hθe
  · show W.hom.hom.2 ≫ δ₀ ≫ ζ = Z.hom.hom.2 ≫ β
    rw [← Category.assoc]
    exact hζe

/-- ★★★**添字圏 `(A,B)𝒞^{bi-Fr}` は filtered** —— 帰納極限がまともに振る舞う根拠。 -/
instance idxPf_isFiltered {A B : C} : IsFiltered (IdxPf P F A B) where
  cocone_objs := idx_cocone_objs
  cocone_maps _ _ u v := ⟨_, 𝟙 _, by rw [idx_hom_ext u v]⟩
  nonempty := ⟨idxMk (𝟙 A) (𝟙 B) (isFrobeniusType_of_isIso P (𝟙 A))
    (isFrobeniusType_of_isIso P (𝟙 B)) (by rw [P.degFr_id, P.degFr_id])⟩

/-! ## ★6. `Hom^pf` —— 帰納極限 -/

/-- ★★★**[FrdI] Definition 3.1, (ii)** —— **`Hom^pf_𝒞(A,B)`**。

★★**余極限が存在すること**は、添字圏の対象が `Type (max u2 v2)` に住むこと
(`Small`)から自動で出る —— ★**原文が「小さい骨格部分圏を取る」で処理した点**を、
我々は宇宙を上げて処理している。 -/
noncomputable def HomPf (A B : C) : Type (max u2 v2) := HomColim (homFunctorPf P F A B)

variable {P F} in
/-- ★添字 `Z` における射 `φ : A′ ⟶ B′` が定める perfected morphism。 -/
noncomputable def HomPf.mk {A B : C} (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    HomPf P F A B :=
  HomColim.mk (homFunctorPf P F A B) Z (ULift.up φ)

variable {P F} in
/-- ★**代表元が取れる**。 -/
theorem HomPf.exists_rep {A B : C} (z : HomPf P F A B) :
    ∃ (Z : IdxPf P F A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2), HomPf.mk Z φ = z := by
  obtain ⟨Z, x, hx⟩ := HomColim.exists_rep (homFunctorPf P F A B) z
  exact ⟨Z, x.down, hx⟩

variable {P F} in
/-- ★**遷移で送っても同じ元**(帰納極限の要)。 -/
@[simp] theorem HomPf.mk_map {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    HomPf.mk W (idxTransport P F u φ) = HomPf.mk Z φ :=
  HomColim.mk_map (homFunctorPf P F A B) u (ULift.up φ)

variable {P F} in
/-- ★★**共通の上界で一致すれば等しい**(filtered 性がここで効く)。 -/
theorem HomPf.sound {A B : C} {Z W : IdxPf P F A B}
    {φ : Z.right.obj.1 ⟶ Z.right.obj.2} {ψ : W.right.obj.1 ⟶ W.right.obj.2}
    (V : IdxPf P F A B) (u : Z ⟶ V) (v : W ⟶ V)
    (h : idxTransport P F u φ = idxTransport P F v ψ) : HomPf.mk Z φ = HomPf.mk W ψ :=
  HomColim.sound (homFunctorPf P F A B) V u v (congrArg ULift.up h)

variable {P F} in
/-- ★★**逆向き** —— 等しければ共通の上界で一致する。 -/
theorem HomPf.eq_iff {A B : C} {Z W : IdxPf P F A B}
    {φ : Z.right.obj.1 ⟶ Z.right.obj.2} {ψ : W.right.obj.1 ⟶ W.right.obj.2} :
    HomPf.mk Z φ = HomPf.mk W ψ ↔
      ∃ (V : IdxPf P F A B) (u : Z ⟶ V) (v : W ⟶ V),
        idxTransport P F u φ = idxTransport P F v ψ := by
  refine (HomColim.eq_iff (homFunctorPf P F A B)).trans ⟨?_, ?_⟩
  · rintro ⟨V, u, v, h⟩
    exact ⟨V, u, v, ULift.up_injective h⟩
  · rintro ⟨V, u, v, h⟩
    exact ⟨V, u, v, congrArg ULift.up h⟩

/-! ## ★7. `𝒞 → 𝒞^pf` の出発点 —— 自然な写像 `Hom_𝒞(A,B) → Hom^pf_𝒞(A,B)` -/

/-- ★添字圏の**始対象にあたる対象** `(A →^{𝟙} A, B →^{𝟙} B)`。 -/
def idxOne (A B : C) : IdxPf P F A B :=
  idxMk (𝟙 A) (𝟙 B) (isFrobeniusType_of_isIso P (𝟙 A)) (isFrobeniusType_of_isIso P (𝟙 B))
    (by rw [P.degFr_id, P.degFr_id])

@[simp] theorem idxOne_right_obj (A B : C) : (idxOne P F A B).right.obj = (A, B) := rfl

variable {P F} in
/-- ★★**自然な写像 `Hom_𝒞(A,B) → Hom^pf_𝒞(A,B)`**(原文の関手 `𝒞 → 𝒞^pf` の射の部分)。 -/
noncomputable def toHomPf {A B : C} (φ : A ⟶ B) : HomPf P F A B :=
  HomPf.mk (idxOne P F A B) φ

/-! ## ★8. 3 つ組の添字圏 —— 合成のために「真ん中を共有」する

★★**合成の急所は「`Hom^pf(A,B)` と `Hom^pf(B,E)` の添字が別々である」こと**。
`f` の代表は `B` を `B′` に伸ばし、`g` の代表は `B` を `B″` に伸ばすので、
そのままでは繋がらない。

★★**逃げ道**: **3 脚をまとめた添字圏 `(A,B,E)𝒞^{tri-Fr}`** を作り、
そこからの 3 つの射影が**すべて cofinal(final)**であることを示す。
★すると 3 つの `Hom^pf` が**同じ添字圏の上の余極限**として書け、
合成は各段での `≫` から誘導される。
★★**射影の finality は「脚を 1 本足すだけ」で出る** —— 次数を変えずに
残りの脚を `frobDegSurj` で作ればよい。★**ここが原文の「in the evident fashion」の中身。**
-/

/-- ★**3 成分とも Frobenius 型で、次数がすべて等しい**射のクラス。 -/
def triFrProp : MorphismProperty (C × C × C) :=
  fun _ _ f => IsFrobeniusType P f.1 ∧ IsFrobeniusType P f.2.1 ∧ IsFrobeniusType P f.2.2 ∧
    P.degFr f.1 = P.degFr f.2.1 ∧ P.degFr f.2.1 = P.degFr f.2.2

include F in
theorem triFrProp_isMultiplicative : (triFrProp P).IsMultiplicative where
  id_mem X :=
    ⟨isFrobeniusType_of_isIso P (𝟙 X.1), isFrobeniusType_of_isIso P (𝟙 X.2.1),
      isFrobeniusType_of_isIso P (𝟙 X.2.2),
      show P.degFr (𝟙 X.1) = P.degFr (𝟙 X.2.1) by rw [P.degFr_id, P.degFr_id],
      show P.degFr (𝟙 X.2.1) = P.degFr (𝟙 X.2.2) by rw [P.degFr_id, P.degFr_id]⟩
  comp_mem f g hf hg :=
    ⟨IsFrobeniusType.comp P F hf.1 hg.1, IsFrobeniusType.comp P F hf.2.1 hg.2.1,
      IsFrobeniusType.comp P F hf.2.2.1 hg.2.2.1,
      show P.degFr (f.1 ≫ g.1) = P.degFr (f.2.1 ≫ g.2.1) by
        rw [P.degFr_comp, P.degFr_comp, hf.2.2.2.1, hg.2.2.2.1],
      show P.degFr (f.2.1 ≫ g.2.1) = P.degFr (f.2.2 ≫ g.2.2) by
        rw [P.degFr_comp, P.degFr_comp, hf.2.2.2.2, hg.2.2.2.2]⟩

/-- ★`F` を型に載せた `triFrProp`。 -/
def triFrPropOf (_F : FrobenioidCore P) : MorphismProperty (C × C × C) := triFrProp P

instance triFrPropOf_isMultiplicative : (triFrPropOf P F).IsMultiplicative :=
  triFrProp_isMultiplicative P F

/-- ★★**`𝒞^{tri-Fr}`** —— `𝒞^{bi-Fr}` の 3 脚版。 -/
abbrev TriFr : Type u2 := WideSubcategory (triFrPropOf P F)

/-- ★3 つ組の添字圏の底。 -/
def triFrObj (A B E : C) : TriFr P F := ⟨(A, B, E)⟩

/-- ★★**3 つ組の添字圏 `(A,B,E)𝒞^{tri-Fr}`**。 -/
abbrev IdxPf3 (A B E : C) : Type (max u2 v2) := Under (triFrObj P F A B E)

/-- ★第 1・第 2 脚への射影。 -/
def triToBi12 : TriFr P F ⥤ BiFr P F where
  obj X := ⟨(X.obj.1, X.obj.2.1)⟩
  map f := ⟨(f.hom.1, f.hom.2.1), f.property.1, f.property.2.1, f.property.2.2.2.1⟩

/-- ★第 2・第 3 脚への射影。 -/
def triToBi23 : TriFr P F ⥤ BiFr P F where
  obj X := ⟨(X.obj.2.1, X.obj.2.2)⟩
  map f := ⟨(f.hom.2.1, f.hom.2.2), f.property.2.1, f.property.2.2.1, f.property.2.2.2.2⟩

/-- ★第 1・第 3 脚への射影。★次数の一致は 2 本の等式を繋いで得る。 -/
def triToBi13 : TriFr P F ⥤ BiFr P F where
  obj X := ⟨(X.obj.1, X.obj.2.2)⟩
  map f := ⟨(f.hom.1, f.hom.2.2), f.property.1, f.property.2.2.1,
    f.property.2.2.2.1.trans f.property.2.2.2.2⟩

/-- ★★射影 `(A,B,E)𝒞^{tri-Fr} ⥤ (A,B)𝒞^{bi-Fr}`。 -/
def idx12 (A B E : C) : IdxPf3 P F A B E ⥤ IdxPf P F A B := Under.post (triToBi12 P F)

/-- ★★射影 `(A,B,E)𝒞^{tri-Fr} ⥤ (B,E)𝒞^{bi-Fr}`。 -/
def idx23 (A B E : C) : IdxPf3 P F A B E ⥤ IdxPf P F B E := Under.post (triToBi23 P F)

/-- ★★射影 `(A,B,E)𝒞^{tri-Fr} ⥤ (A,E)𝒞^{bi-Fr}`。 -/
def idx13 (A B E : C) : IdxPf3 P F A B E ⥤ IdxPf P F A E := Under.post (triToBi13 P F)

/-! ### ★3 つ組の添字圏も「細い」filtered -/

include F in
/-- ★★**2 本の Frobenius 型射に共通の上界** —— `Definition 1.3, (ii)` の
存在(`frobDegSurj`)と本質的一意性(`frobDegUniq`)から。 -/
theorem frob_common_upper {X X₁ X₂ : C} (a : X ⟶ X₁) (ha : IsFrobeniusType P a)
    (c : X ⟶ X₂) (hc : IsFrobeniusType P c) :
    ∃ (X₃ : C) (α : X₁ ⟶ X₃) (γ : X₂ ⟶ X₃), IsFrobeniusType P α ∧ IsFrobeniusType P γ ∧
      P.degFr α = P.degFr c ∧ P.degFr γ = P.degFr a ∧ a ≫ α = c ≫ γ := by
  obtain ⟨X₃, α, hα, hαd⟩ := F.frobDegSurj X₁ (P.degFr c)
  obtain ⟨X₄, γ₀, hγ₀, hγ₀d⟩ := F.frobDegSurj X₂ (P.degFr a)
  have hg : P.degFr (c ≫ γ₀) = P.degFr (a ≫ α) := by
    rw [P.degFr_comp, P.degFr_comp, hγ₀d, hαd, mul_comm]
  obtain ⟨θ, hθ, hθe⟩ := F.frobDegUniq X X₄ X₃ (c ≫ γ₀) (a ≫ α)
    (IsFrobeniusType.comp P F hc hγ₀) (IsFrobeniusType.comp P F ha hα) hg
  haveI : IsIso θ := hθ
  refine ⟨X₃, α, γ₀ ≫ θ, hα, IsFrobeniusType.comp P F hγ₀ (isFrobeniusType_of_isIso P θ), hαd,
    by rw [P.degFr_comp, show P.degFr θ = 1 from isLinear_of_isIso P θ, hγ₀d, one_mul], ?_⟩
  rw [← Category.assoc]
  exact hθe.symm

variable {P F} in
/-- ★★**3 つ組の添字圏の平行 2 射も等しい**(`idx_hom_ext` と同じ理由)。 -/
theorem idx3_hom_ext {A B E : C} {Z W : IdxPf3 P F A B E} (u v : Z ⟶ W) : u = v := by
  have h : Z.hom ≫ u.right = Z.hom ≫ v.right := by rw [Under.w u, Under.w v]
  have h1 : Z.hom.hom.1 ≫ u.right.hom.1 = Z.hom.hom.1 ≫ v.right.hom.1 :=
    congrArg (fun t : triFrObj P F A B E ⟶ W.right => t.hom.1) h
  have h2 : Z.hom.hom.2.1 ≫ u.right.hom.2.1 = Z.hom.hom.2.1 ≫ v.right.hom.2.1 :=
    congrArg (fun t : triFrObj P F A B E ⟶ W.right => t.hom.2.1) h
  have h3 : Z.hom.hom.2.2 ≫ u.right.hom.2.2 = Z.hom.hom.2.2 ≫ v.right.hom.2.2 :=
    congrArg (fun t : triFrObj P F A B E ⟶ W.right => t.hom.2.2) h
  haveI : Epi Z.hom.hom.1 := P.totEpiC _ _ _
  haveI : Epi Z.hom.hom.2.1 := P.totEpiC _ _ _
  haveI : Epi Z.hom.hom.2.2 := P.totEpiC _ _ _
  exact CommaMorphism.ext (Subsingleton.elim _ _) (WideSubcategory.hom_ext _
    (Prod.ext ((cancel_epi _).mp h1)
      (Prod.ext ((cancel_epi _).mp h2) ((cancel_epi _).mp h3))))

variable {P F} in
/-- ★★**3 つ組の添字圏の共通上界** —— 各脚に `frob_common_upper` を当てるだけ。 -/
theorem idx3_cocone_objs {A B E : C} (Z W : IdxPf3 P F A B E) :
    ∃ (V : IdxPf3 P F A B E) (_ : Z ⟶ V) (_ : W ⟶ V), True := by
  obtain ⟨ha, hb, hc, hab, hbc⟩ := Z.hom.property
  obtain ⟨ha2, hb2, hc2, hab2, hbc2⟩ := W.hom.property
  obtain ⟨A₃, α, γA, hα, hγA, hαd, hγAd, hA⟩ :=
    frob_common_upper P F Z.hom.hom.1 ha W.hom.hom.1 ha2
  obtain ⟨B₃, β, γB, hβ, hγB, hβd, hγBd, hB⟩ :=
    frob_common_upper P F Z.hom.hom.2.1 hb W.hom.hom.2.1 hb2
  obtain ⟨E₃, ε, γE, hε, hγE, hεd, hγEd, hE⟩ :=
    frob_common_upper P F Z.hom.hom.2.2 hc W.hom.hom.2.2 hc2
  refine ⟨Under.mk (Y := (⟨(A₃, B₃, E₃)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(Z.hom.hom.1 ≫ α, Z.hom.hom.2.1 ≫ β, Z.hom.hom.2.2 ≫ ε),
          IsFrobeniusType.comp P F ha hα, IsFrobeniusType.comp P F hb hβ,
          IsFrobeniusType.comp P F hc hε,
          show P.degFr (Z.hom.hom.1 ≫ α) = P.degFr (Z.hom.hom.2.1 ≫ β) by
            rw [P.degFr_comp, P.degFr_comp, hαd, hβd, hab, hab2],
          show P.degFr (Z.hom.hom.2.1 ≫ β) = P.degFr (Z.hom.hom.2.2 ≫ ε) by
            rw [P.degFr_comp, P.degFr_comp, hβd, hεd, hbc, hbc2]⟩),
    Under.homMk (show Z.right ⟶ _ from
      ⟨(α, β, ε), hα, hβ, hε, show P.degFr α = P.degFr β by rw [hαd, hβd, hab2],
        show P.degFr β = P.degFr ε by rw [hβd, hεd, hbc2]⟩)
      (WideSubcategory.hom_ext _ rfl),
    Under.homMk (show W.right ⟶ _ from
      ⟨(γA, γB, γE), hγA, hγB, hγE, show P.degFr γA = P.degFr γB by rw [hγAd, hγBd, hab],
        show P.degFr γB = P.degFr γE by rw [hγBd, hγEd, hbc]⟩)
      (WideSubcategory.hom_ext _ (Prod.ext hA.symm (Prod.ext hB.symm hE.symm))), trivial⟩

/-- ★★★**3 つ組の添字圏も filtered**。 -/
instance idxPf3_isFiltered {A B E : C} : IsFiltered (IdxPf3 P F A B E) where
  cocone_objs := idx3_cocone_objs
  cocone_maps _ _ u v := ⟨_, 𝟙 _, by rw [idx3_hom_ext u v]⟩
  nonempty := ⟨Under.mk (Y := (⟨(A, B, E)⟩ : TriFr P F))
    (show triFrObj P F A B E ⟶ _ from
      ⟨(𝟙 A, 𝟙 B, 𝟙 E), isFrobeniusType_of_isIso P (𝟙 A), isFrobeniusType_of_isIso P (𝟙 B),
        isFrobeniusType_of_isIso P (𝟙 E),
        show P.degFr (𝟙 A) = P.degFr (𝟙 B) by rw [P.degFr_id, P.degFr_id],
        show P.degFr (𝟙 B) = P.degFr (𝟙 E) by rw [P.degFr_id, P.degFr_id]⟩)⟩

/-! ### ★★★3 つの射影は cofinal

★★**これが「真ん中を共有できる」ことの正確な意味**である ——
3 つの `Hom^pf` が**同じ添字圏 `(A,B,E)𝒞^{tri-Fr}` の上の余極限**として書ける。

★★**証明は「脚を 1 本足すだけ」** —— 次数を変えずに残りの脚を `frobDegSurj` で作る。
★2 つ目の条件(平行 2 射の同一化)は**細さから自明**。 -/

/-- ★★`(A,B,E)𝒞^{tri-Fr} ⥤ (A,B)𝒞^{bi-Fr}` は cofinal。 -/
instance idx12_final (A B E : C) : (idx12 P F A B E).Final := by
  refine Functor.final_of_exists_of_isFiltered _ ?_ ?_
  · intro Z
    obtain ⟨ha, hb, hab⟩ := Z.hom.property
    obtain ⟨E', c, hc, hcd⟩ := F.frobDegSurj E (P.degFr Z.hom.hom.2)
    exact ⟨Under.mk (Y := (⟨(Z.right.obj.1, Z.right.obj.2, E')⟩ : TriFr P F))
        (show triFrObj P F A B E ⟶ _ from
          ⟨(Z.hom.hom.1, Z.hom.hom.2, c), ha, hb, hc, hab, hcd.symm⟩),
      ⟨Under.homMk (𝟙 Z.right) (WideSubcategory.hom_ext _ (Category.comp_id _))⟩⟩
  · intro _ _ s s'
    exact ⟨_, 𝟙 _, by rw [idx_hom_ext s s']⟩

/-- ★★`(A,B,E)𝒞^{tri-Fr} ⥤ (B,E)𝒞^{bi-Fr}` は cofinal。 -/
instance idx23_final (A B E : C) : (idx23 P F A B E).Final := by
  refine Functor.final_of_exists_of_isFiltered _ ?_ ?_
  · intro Y
    obtain ⟨hb, he, hbe⟩ := Y.hom.property
    obtain ⟨A', a, hafr, had⟩ := F.frobDegSurj A (P.degFr Y.hom.hom.1)
    exact ⟨Under.mk (Y := (⟨(A', Y.right.obj.1, Y.right.obj.2)⟩ : TriFr P F))
        (show triFrObj P F A B E ⟶ _ from
          ⟨(a, Y.hom.hom.1, Y.hom.hom.2), hafr, hb, he, had, hbe⟩),
      ⟨Under.homMk (𝟙 Y.right) (WideSubcategory.hom_ext _ (Category.comp_id _))⟩⟩
  · intro _ _ s s'
    exact ⟨_, 𝟙 _, by rw [idx_hom_ext s s']⟩

/-- ★★`(A,B,E)𝒞^{tri-Fr} ⥤ (A,E)𝒞^{bi-Fr}` は cofinal。 -/
instance idx13_final (A B E : C) : (idx13 P F A B E).Final := by
  refine Functor.final_of_exists_of_isFiltered _ ?_ ?_
  · intro W
    obtain ⟨ha, he, hae⟩ := W.hom.property
    obtain ⟨B', b, hbfr, hbd⟩ := F.frobDegSurj B (P.degFr W.hom.hom.1)
    exact ⟨Under.mk (Y := (⟨(W.right.obj.1, B', W.right.obj.2)⟩ : TriFr P F))
        (show triFrObj P F A B E ⟶ _ from
          ⟨(W.hom.hom.1, b, W.hom.hom.2), ha, hbfr, he, hbd.symm, hbd.trans hae⟩),
      ⟨Under.homMk (𝟙 W.right) (WideSubcategory.hom_ext _ (Category.comp_id _))⟩⟩
  · intro _ _ s s'
    exact ⟨_, 𝟙 _, by rw [idx_hom_ext s s']⟩

/-! ## ★9. `Hom^pf` の合成

★★原文は「composition of morphisms of `𝒞^pf` is defined in the evident fashion」の
一行で済ませる。★**その "evident" の中身が、上の cofinal 性と、下の
「遷移と合成の両立」の 2 つ**である。 -/

variable {P F} in
/-- ★添字圏の射に沿う遷移は四角形を可換にする(`frobTransport_spec` の言い換え)。 -/
theorem idxTransport_spec {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    φ ≫ u.right.hom.2 = u.right.hom.1 ≫ idxTransport P F u φ :=
  frobTransport_spec _ _ _ _ _ φ

variable {P F} in
/-- ★★**遷移と合成の両立** —— 3 つ組の添字圏の射に沿って遷移させると、
**合成は遷移の合成に写る**。

★証明は `frobTransport` の一意性 1 本。★2 つの四角形を横に並べて
外側の四角形を作るだけ。 -/
theorem idxTransport_comp_pair {A B E : C} {V W : IdxPf3 P F A B E} (t : V ⟶ W)
    (φ : V.right.obj.1 ⟶ V.right.obj.2.1) (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2) :
    idxTransport P F ((idx12 P F A B E).map t) φ ≫ idxTransport P F ((idx23 P F A B E).map t) ψ
      = idxTransport P F ((idx13 P F A B E).map t) (φ ≫ ψ) := by
  have h1 : φ ≫ t.right.hom.2.1
      = t.right.hom.1 ≫ idxTransport P F ((idx12 P F A B E).map t) φ :=
    idxTransport_spec ((idx12 P F A B E).map t) φ
  have h2 : ψ ≫ t.right.hom.2.2
      = t.right.hom.2.1 ≫ idxTransport P F ((idx23 P F A B E).map t) ψ :=
    idxTransport_spec ((idx23 P F A B E).map t) ψ
  refine (frobTransport_eq _ _ _ _ _ (φ ≫ ψ) _ ?_).symm
  calc (φ ≫ ψ) ≫ t.right.hom.2.2 = φ ≫ (ψ ≫ t.right.hom.2.2) := by rw [Category.assoc]
    _ = φ ≫ (t.right.hom.2.1 ≫ idxTransport P F ((idx23 P F A B E).map t) ψ) := by rw [h2]; rfl
    _ = (φ ≫ t.right.hom.2.1) ≫ idxTransport P F ((idx23 P F A B E).map t) ψ := by
        rw [Category.assoc]
    _ = (t.right.hom.1 ≫ idxTransport P F ((idx12 P F A B E).map t) φ)
          ≫ idxTransport P F ((idx23 P F A B E).map t) ψ := by rw [h1]; rfl
    _ = t.right.hom.1 ≫ (idxTransport P F ((idx12 P F A B E).map t) φ
          ≫ idxTransport P F ((idx23 P F A B E).map t) ψ) := by rw [Category.assoc]

/-- ★★**2 変数写像の材料** —— 細さ・各段の合成・自然性の 3 つ。 -/
noncomputable def compBiData (A B E : C) :
    HomColim.BiData (idx12 P F A B E ⋙ homFunctorPf P F A B)
      (idx23 P F A B E ⋙ homFunctorPf P F B E)
      (idx13 P F A B E ⋙ homFunctorPf P F A E) where
  thin := idx3_hom_ext
  m _ φ ψ := ULift.up (φ.down ≫ ψ.down)
  hm t x y := congrArg ULift.up (idxTransport_comp_pair t x.down y.down)

variable {P F} in
/-- ★`Hom^pf(A,B)` の元を、3 つ組の添字での余極限の元に読み替える。 -/
theorem colimIso12_inv_mk {A B E : C} (V : IdxPf3 P F A B E)
    (φ : V.right.obj.1 ⟶ V.right.obj.2.1) :
    (Functor.Final.colimitIso (idx12 P F A B E) (homFunctorPf P F A B)).inv
        (HomPf.mk ((idx12 P F A B E).obj V) φ)
      = HomColim.mk (idx12 P F A B E ⋙ homFunctorPf P F A B) V ⟨φ⟩ := by
  show (Functor.Final.colimitIso (idx12 P F A B E) (homFunctorPf P F A B)).inv
      (colimit.ι (homFunctorPf P F A B) ((idx12 P F A B E).obj V) (ULift.up φ)) = _
  rw [← types_comp_apply (colimit.ι (homFunctorPf P F A B) ((idx12 P F A B E).obj V))
    (Functor.Final.colimitIso (idx12 P F A B E) (homFunctorPf P F A B)).inv,
    Functor.Final.ι_colimitIso_inv]
  rfl

variable {P F} in
/-- ★同上(第 2・第 3 脚)。 -/
theorem colimIso23_inv_mk {A B E : C} (V : IdxPf3 P F A B E)
    (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2) :
    (Functor.Final.colimitIso (idx23 P F A B E) (homFunctorPf P F B E)).inv
        (HomPf.mk ((idx23 P F A B E).obj V) ψ)
      = HomColim.mk (idx23 P F A B E ⋙ homFunctorPf P F B E) V ⟨ψ⟩ := by
  show (Functor.Final.colimitIso (idx23 P F A B E) (homFunctorPf P F B E)).inv
      (colimit.ι (homFunctorPf P F B E) ((idx23 P F A B E).obj V) (ULift.up ψ)) = _
  rw [← types_comp_apply (colimit.ι (homFunctorPf P F B E) ((idx23 P F A B E).obj V))
    (Functor.Final.colimitIso (idx23 P F A B E) (homFunctorPf P F B E)).inv,
    Functor.Final.ι_colimitIso_inv]
  rfl

variable {P F} in
/-- ★同上(第 1・第 3 脚、順方向)。 -/
theorem colimIso13_hom_mk {A B E : C} (V : IdxPf3 P F A B E)
    (χ : V.right.obj.1 ⟶ V.right.obj.2.2) :
    (Functor.Final.colimitIso (idx13 P F A B E) (homFunctorPf P F A E)).hom
        (HomColim.mk (idx13 P F A B E ⋙ homFunctorPf P F A E) V ⟨χ⟩)
      = HomPf.mk ((idx13 P F A B E).obj V) χ := by
  show (Functor.Final.colimitIso (idx13 P F A B E) (homFunctorPf P F A E)).hom
      (colimit.ι (idx13 P F A B E ⋙ homFunctorPf P F A E) V (ULift.up χ)) = _
  rw [← types_comp_apply (colimit.ι (idx13 P F A B E ⋙ homFunctorPf P F A E) V)
    (Functor.Final.colimitIso (idx13 P F A B E) (homFunctorPf P F A E)).hom,
    Functor.Final.ι_colimitIso_hom]
  rfl

/-- ★★★**[FrdI] Definition 3.1, (ii)(iii)** —— **perfected morphism の合成**。

★3 つの `Hom^pf` を cofinal 性で「同じ添字圏の上の余極限」に読み替え、
各段で `≫` を取る。 -/
noncomputable def compPf {A B E : C} (f : HomPf P F A B) (g : HomPf P F B E) : HomPf P F A E :=
  (Functor.Final.colimitIso (idx13 P F A B E) (homFunctorPf P F A E)).hom
    (HomColim.bimap (compBiData P F A B E)
      ((Functor.Final.colimitIso (idx12 P F A B E) (homFunctorPf P F A B)).inv f)
      ((Functor.Final.colimitIso (idx23 P F A B E) (homFunctorPf P F B E)).inv g))

variable {P F} in
/-- ★★★**合成の計算則** —— 真ん中を共有する代表元では、合成は素直な `≫`。 -/
theorem compPf_mk {A B E : C} (V : IdxPf3 P F A B E)
    (φ : V.right.obj.1 ⟶ V.right.obj.2.1) (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2) :
    compPf P F (HomPf.mk ((idx12 P F A B E).obj V) φ) (HomPf.mk ((idx23 P F A B E).obj V) ψ)
      = HomPf.mk ((idx13 P F A B E).obj V) (φ ≫ ψ) := by
  show (Functor.Final.colimitIso (idx13 P F A B E) (homFunctorPf P F A E)).hom
      (HomColim.bimap (compBiData P F A B E) _ _) = _
  rw [colimIso12_inv_mk, colimIso23_inv_mk,
    HomColim.bimap_mk_same (compBiData P F A B E) V (ULift.up φ) (ULift.up ψ)]
  exact colimIso13_hom_mk V (φ ≫ ψ)

/-! ## ★10. 4 つ組の添字圏 —— **結合律**のために

★★**結合律には脚が 4 本要る**(`A, B, E, G`)。★3 つ組と同じ形なので、
`quadFrProp` から `IdxPf4` までを機械的にもう一段作る。

★★**要点**: 4 つ組から 3 つ組への**射影が「脚を 1 本落とすだけ」**であり、
`(idx13).obj ((q123).obj V)` と `(idx12).obj ((q134).obj V)` が
★**定義的に等しい**(どちらも脚 `(a, e)`)。★これで `compPf_mk` を 4 回当てるだけで
結合律が出る。 -/

/-- ★**4 成分とも Frobenius 型で、次数がすべて等しい**射のクラス。 -/
def quadFrProp : MorphismProperty (C × C × C × C) :=
  fun _ _ f => IsFrobeniusType P f.1 ∧ IsFrobeniusType P f.2.1 ∧ IsFrobeniusType P f.2.2.1 ∧
    IsFrobeniusType P f.2.2.2 ∧ P.degFr f.1 = P.degFr f.2.1 ∧
    P.degFr f.2.1 = P.degFr f.2.2.1 ∧ P.degFr f.2.2.1 = P.degFr f.2.2.2

include F in
theorem quadFrProp_isMultiplicative : (quadFrProp P).IsMultiplicative where
  id_mem X :=
    ⟨isFrobeniusType_of_isIso P (𝟙 X.1), isFrobeniusType_of_isIso P (𝟙 X.2.1),
      isFrobeniusType_of_isIso P (𝟙 X.2.2.1), isFrobeniusType_of_isIso P (𝟙 X.2.2.2),
      show P.degFr (𝟙 X.1) = P.degFr (𝟙 X.2.1) by rw [P.degFr_id, P.degFr_id],
      show P.degFr (𝟙 X.2.1) = P.degFr (𝟙 X.2.2.1) by rw [P.degFr_id, P.degFr_id],
      show P.degFr (𝟙 X.2.2.1) = P.degFr (𝟙 X.2.2.2) by rw [P.degFr_id, P.degFr_id]⟩
  comp_mem f g hf hg :=
    ⟨IsFrobeniusType.comp P F hf.1 hg.1, IsFrobeniusType.comp P F hf.2.1 hg.2.1,
      IsFrobeniusType.comp P F hf.2.2.1 hg.2.2.1, IsFrobeniusType.comp P F hf.2.2.2.1 hg.2.2.2.1,
      show P.degFr (f.1 ≫ g.1) = P.degFr (f.2.1 ≫ g.2.1) by
        rw [P.degFr_comp, P.degFr_comp, hf.2.2.2.2.1, hg.2.2.2.2.1],
      show P.degFr (f.2.1 ≫ g.2.1) = P.degFr (f.2.2.1 ≫ g.2.2.1) by
        rw [P.degFr_comp, P.degFr_comp, hf.2.2.2.2.2.1, hg.2.2.2.2.2.1],
      show P.degFr (f.2.2.1 ≫ g.2.2.1) = P.degFr (f.2.2.2 ≫ g.2.2.2) by
        rw [P.degFr_comp, P.degFr_comp, hf.2.2.2.2.2.2, hg.2.2.2.2.2.2]⟩

/-- ★`F` を型に載せた `quadFrProp`。 -/
def quadFrPropOf (_F : FrobenioidCore P) : MorphismProperty (C × C × C × C) := quadFrProp P

instance quadFrPropOf_isMultiplicative : (quadFrPropOf P F).IsMultiplicative :=
  quadFrProp_isMultiplicative P F

/-- ★**`𝒞^{quad-Fr}`**。 -/
abbrev QuadFr : Type u2 := WideSubcategory (quadFrPropOf P F)

/-- ★4 つ組の添字圏の底。 -/
def quadFrObj (A B E G : C) : QuadFr P F := ⟨(A, B, E, G)⟩

/-- ★★**4 つ組の添字圏**。 -/
abbrev IdxPf4 (A B E G : C) : Type (max u2 v2) := Under (quadFrObj P F A B E G)

variable {P F} in
/-- ★4 つ組の添字圏も細い。 -/
theorem idx4_hom_ext {A B E G : C} {Z W : IdxPf4 P F A B E G} (u v : Z ⟶ W) : u = v := by
  have h : Z.hom ≫ u.right = Z.hom ≫ v.right := by rw [Under.w u, Under.w v]
  have h1 : Z.hom.hom.1 ≫ u.right.hom.1 = Z.hom.hom.1 ≫ v.right.hom.1 :=
    congrArg (fun t : quadFrObj P F A B E G ⟶ W.right => t.hom.1) h
  have h2 : Z.hom.hom.2.1 ≫ u.right.hom.2.1 = Z.hom.hom.2.1 ≫ v.right.hom.2.1 :=
    congrArg (fun t : quadFrObj P F A B E G ⟶ W.right => t.hom.2.1) h
  have h3 : Z.hom.hom.2.2.1 ≫ u.right.hom.2.2.1 = Z.hom.hom.2.2.1 ≫ v.right.hom.2.2.1 :=
    congrArg (fun t : quadFrObj P F A B E G ⟶ W.right => t.hom.2.2.1) h
  have h4 : Z.hom.hom.2.2.2 ≫ u.right.hom.2.2.2 = Z.hom.hom.2.2.2 ≫ v.right.hom.2.2.2 :=
    congrArg (fun t : quadFrObj P F A B E G ⟶ W.right => t.hom.2.2.2) h
  haveI : Epi Z.hom.hom.1 := P.totEpiC _ _ _
  haveI : Epi Z.hom.hom.2.1 := P.totEpiC _ _ _
  haveI : Epi Z.hom.hom.2.2.1 := P.totEpiC _ _ _
  haveI : Epi Z.hom.hom.2.2.2 := P.totEpiC _ _ _
  exact CommaMorphism.ext (Subsingleton.elim _ _) (WideSubcategory.hom_ext _
    (Prod.ext ((cancel_epi _).mp h1) (Prod.ext ((cancel_epi _).mp h2)
      (Prod.ext ((cancel_epi _).mp h3) ((cancel_epi _).mp h4)))))

variable {P F} in
/-- ★4 つ組の添字圏の共通上界。 -/
theorem idx4_cocone_objs {A B E G : C} (Z W : IdxPf4 P F A B E G) :
    ∃ (V : IdxPf4 P F A B E G) (_ : Z ⟶ V) (_ : W ⟶ V), True := by
  obtain ⟨ha, hb, hc, hd, hab, hbc, hcd⟩ := Z.hom.property
  obtain ⟨ha2, hb2, hc2, hd2, hab2, hbc2, hcd2⟩ := W.hom.property
  obtain ⟨A₃, α, γA, hα, hγA, hαd, hγAd, hA⟩ :=
    frob_common_upper P F Z.hom.hom.1 ha W.hom.hom.1 ha2
  obtain ⟨B₃, β, γB, hβ, hγB, hβd, hγBd, hB⟩ :=
    frob_common_upper P F Z.hom.hom.2.1 hb W.hom.hom.2.1 hb2
  obtain ⟨E₃, ε, γE, hε, hγE, hεd, hγEd, hE⟩ :=
    frob_common_upper P F Z.hom.hom.2.2.1 hc W.hom.hom.2.2.1 hc2
  obtain ⟨G₃, ζ, γG, hζ, hγG, hζd, hγGd, hG⟩ :=
    frob_common_upper P F Z.hom.hom.2.2.2 hd W.hom.hom.2.2.2 hd2
  refine ⟨Under.mk (Y := (⟨(A₃, B₃, E₃, G₃)⟩ : QuadFr P F))
      (show quadFrObj P F A B E G ⟶ _ from
        ⟨(Z.hom.hom.1 ≫ α, Z.hom.hom.2.1 ≫ β, Z.hom.hom.2.2.1 ≫ ε, Z.hom.hom.2.2.2 ≫ ζ),
          IsFrobeniusType.comp P F ha hα, IsFrobeniusType.comp P F hb hβ,
          IsFrobeniusType.comp P F hc hε, IsFrobeniusType.comp P F hd hζ,
          show P.degFr (Z.hom.hom.1 ≫ α) = P.degFr (Z.hom.hom.2.1 ≫ β) by
            rw [P.degFr_comp, P.degFr_comp, hαd, hβd, hab, hab2],
          show P.degFr (Z.hom.hom.2.1 ≫ β) = P.degFr (Z.hom.hom.2.2.1 ≫ ε) by
            rw [P.degFr_comp, P.degFr_comp, hβd, hεd, hbc, hbc2],
          show P.degFr (Z.hom.hom.2.2.1 ≫ ε) = P.degFr (Z.hom.hom.2.2.2 ≫ ζ) by
            rw [P.degFr_comp, P.degFr_comp, hεd, hζd, hcd, hcd2]⟩),
    Under.homMk (show Z.right ⟶ _ from
      ⟨(α, β, ε, ζ), hα, hβ, hε, hζ,
        show P.degFr α = P.degFr β by rw [hαd, hβd, hab2],
        show P.degFr β = P.degFr ε by rw [hβd, hεd, hbc2],
        show P.degFr ε = P.degFr ζ by rw [hεd, hζd, hcd2]⟩)
      (WideSubcategory.hom_ext _ rfl),
    Under.homMk (show W.right ⟶ _ from
      ⟨(γA, γB, γE, γG), hγA, hγB, hγE, hγG,
        show P.degFr γA = P.degFr γB by rw [hγAd, hγBd, hab],
        show P.degFr γB = P.degFr γE by rw [hγBd, hγEd, hbc],
        show P.degFr γE = P.degFr γG by rw [hγEd, hγGd, hcd]⟩)
      (WideSubcategory.hom_ext _
        (Prod.ext hA.symm (Prod.ext hB.symm (Prod.ext hE.symm hG.symm)))), trivial⟩

/-- ★★4 つ組の添字圏も filtered。 -/
instance idxPf4_isFiltered {A B E G : C} : IsFiltered (IdxPf4 P F A B E G) where
  cocone_objs := idx4_cocone_objs
  cocone_maps _ _ u v := ⟨_, 𝟙 _, by rw [idx4_hom_ext u v]⟩
  nonempty := ⟨Under.mk (Y := (⟨(A, B, E, G)⟩ : QuadFr P F))
    (show quadFrObj P F A B E G ⟶ _ from
      ⟨(𝟙 A, 𝟙 B, 𝟙 E, 𝟙 G), isFrobeniusType_of_isIso P (𝟙 A),
        isFrobeniusType_of_isIso P (𝟙 B), isFrobeniusType_of_isIso P (𝟙 E),
        isFrobeniusType_of_isIso P (𝟙 G),
        show P.degFr (𝟙 A) = P.degFr (𝟙 B) by rw [P.degFr_id, P.degFr_id],
        show P.degFr (𝟙 B) = P.degFr (𝟙 E) by rw [P.degFr_id, P.degFr_id],
        show P.degFr (𝟙 E) = P.degFr (𝟙 G) by rw [P.degFr_id, P.degFr_id]⟩)⟩

/-- ★脚 1・2・3 への射影。 -/
def quadToTri123 : QuadFr P F ⥤ TriFr P F where
  obj X := ⟨(X.obj.1, X.obj.2.1, X.obj.2.2.1)⟩
  map f := ⟨(f.hom.1, f.hom.2.1, f.hom.2.2.1), f.property.1, f.property.2.1, f.property.2.2.1,
    f.property.2.2.2.2.1, f.property.2.2.2.2.2.1⟩

/-- ★脚 2・3・4 への射影。 -/
def quadToTri234 : QuadFr P F ⥤ TriFr P F where
  obj X := ⟨(X.obj.2.1, X.obj.2.2.1, X.obj.2.2.2)⟩
  map f := ⟨(f.hom.2.1, f.hom.2.2.1, f.hom.2.2.2), f.property.2.1, f.property.2.2.1,
    f.property.2.2.2.1, f.property.2.2.2.2.2.1, f.property.2.2.2.2.2.2⟩

/-- ★脚 1・2・4 への射影。 -/
def quadToTri124 : QuadFr P F ⥤ TriFr P F where
  obj X := ⟨(X.obj.1, X.obj.2.1, X.obj.2.2.2)⟩
  map f := ⟨(f.hom.1, f.hom.2.1, f.hom.2.2.2), f.property.1, f.property.2.1,
    f.property.2.2.2.1, f.property.2.2.2.2.1,
    f.property.2.2.2.2.2.1.trans f.property.2.2.2.2.2.2⟩

/-- ★脚 1・3・4 への射影。 -/
def quadToTri134 : QuadFr P F ⥤ TriFr P F where
  obj X := ⟨(X.obj.1, X.obj.2.2.1, X.obj.2.2.2)⟩
  map f := ⟨(f.hom.1, f.hom.2.2.1, f.hom.2.2.2), f.property.1, f.property.2.2.1,
    f.property.2.2.2.1, f.property.2.2.2.2.1.trans f.property.2.2.2.2.2.1,
    f.property.2.2.2.2.2.2⟩

/-- ★★`IdxPf4 ⥤ IdxPf3 A B E`。 -/
def q123 (A B E G : C) : IdxPf4 P F A B E G ⥤ IdxPf3 P F A B E :=
  Under.post (quadToTri123 P F)

/-- ★★`IdxPf4 ⥤ IdxPf3 B E G`。 -/
def q234 (A B E G : C) : IdxPf4 P F A B E G ⥤ IdxPf3 P F B E G :=
  Under.post (quadToTri234 P F)

/-- ★★`IdxPf4 ⥤ IdxPf3 A B G`。 -/
def q124 (A B E G : C) : IdxPf4 P F A B E G ⥤ IdxPf3 P F A B G :=
  Under.post (quadToTri124 P F)

/-- ★★`IdxPf4 ⥤ IdxPf3 A E G`。 -/
def q134 (A B E G : C) : IdxPf4 P F A B E G ⥤ IdxPf3 P F A E G :=
  Under.post (quadToTri134 P F)

variable {P F} in
/-- ★★**3 つの perfected morphism を 1 つの 4 つ組の添字で表す**。

★各々を 4 つ組に延ばして(残りの脚は同次数で新しく作る)、
filtered 性で共通上界を取るだけ。 -/
theorem exists_rep4 {A B E G : C} (f : HomPf P F A B) (g : HomPf P F B E) (h : HomPf P F E G) :
    ∃ (V : IdxPf4 P F A B E G)
      (φ : V.right.obj.1 ⟶ V.right.obj.2.1)
      (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2.1)
      (χ : V.right.obj.2.2.1 ⟶ V.right.obj.2.2.2),
      f = HomPf.mk ((idx12 P F A B E).obj ((q123 P F A B E G).obj V)) φ ∧
      g = HomPf.mk ((idx23 P F A B E).obj ((q123 P F A B E G).obj V)) ψ ∧
      h = HomPf.mk ((idx23 P F A E G).obj ((q134 P F A B E G).obj V)) χ := by
  obtain ⟨Zf, φ₀, rfl⟩ := HomPf.exists_rep f
  obtain ⟨Zg, ψ₀, rfl⟩ := HomPf.exists_rep g
  obtain ⟨Zh, χ₀, rfl⟩ := HomPf.exists_rep h
  obtain ⟨hfa, hfb, hfab⟩ := Zf.hom.property
  obtain ⟨hgb, hge, hgbe⟩ := Zg.hom.property
  obtain ⟨hhe, hhg, hheg⟩ := Zh.hom.property
  obtain ⟨E₁, e₁, he₁, he₁d⟩ := F.frobDegSurj E (P.degFr Zf.hom.hom.1)
  obtain ⟨G₁, g₁, hg₁, hg₁d⟩ := F.frobDegSurj G (P.degFr Zf.hom.hom.1)
  obtain ⟨A₂, a₂, ha₂, ha₂d⟩ := F.frobDegSurj A (P.degFr Zg.hom.hom.1)
  obtain ⟨G₂, g₂, hg₂, hg₂d⟩ := F.frobDegSurj G (P.degFr Zg.hom.hom.1)
  obtain ⟨A₃, a₃, ha₃, ha₃d⟩ := F.frobDegSurj A (P.degFr Zh.hom.hom.1)
  obtain ⟨B₃, b₃, hb₃, hb₃d⟩ := F.frobDegSurj B (P.degFr Zh.hom.hom.1)
  set Vf : IdxPf4 P F A B E G :=
    Under.mk (Y := (⟨(Zf.right.obj.1, Zf.right.obj.2, E₁, G₁)⟩ : QuadFr P F))
      (show quadFrObj P F A B E G ⟶ _ from
        ⟨(Zf.hom.hom.1, Zf.hom.hom.2, e₁, g₁), hfa, hfb, he₁, hg₁, hfab,
          hfab.symm.trans he₁d.symm, he₁d.trans hg₁d.symm⟩) with hVf
  set Vg : IdxPf4 P F A B E G :=
    Under.mk (Y := (⟨(A₂, Zg.right.obj.1, Zg.right.obj.2, G₂)⟩ : QuadFr P F))
      (show quadFrObj P F A B E G ⟶ _ from
        ⟨(a₂, Zg.hom.hom.1, Zg.hom.hom.2, g₂), ha₂, hgb, hge, hg₂, ha₂d, hgbe,
          hgbe.symm.trans hg₂d.symm⟩) with hVg
  set Vh : IdxPf4 P F A B E G :=
    Under.mk (Y := (⟨(A₃, B₃, Zh.right.obj.1, Zh.right.obj.2)⟩ : QuadFr P F))
      (show quadFrObj P F A B E G ⟶ _ from
        ⟨(a₃, b₃, Zh.hom.hom.1, Zh.hom.hom.2), ha₃, hb₃, hhe, hhg,
          ha₃d.trans hb₃d.symm, hb₃d, hheg⟩) with hVh
  refine ⟨IsFiltered.max (IsFiltered.max Vf Vg) Vh,
    idxTransport P F ((idx12 P F A B E).map ((q123 P F A B E G).map
      (IsFiltered.leftToMax Vf Vg ≫ IsFiltered.leftToMax (IsFiltered.max Vf Vg) Vh))) φ₀,
    idxTransport P F ((idx23 P F A B E).map ((q123 P F A B E G).map
      (IsFiltered.rightToMax Vf Vg ≫ IsFiltered.leftToMax (IsFiltered.max Vf Vg) Vh))) ψ₀,
    idxTransport P F ((idx23 P F A E G).map ((q134 P F A B E G).map
      (IsFiltered.rightToMax (IsFiltered.max Vf Vg) Vh))) χ₀,
    (HomPf.mk_map _ φ₀).symm, (HomPf.mk_map _ ψ₀).symm, (HomPf.mk_map _ χ₀).symm⟩

variable {P F} in
/-- ★★★**合成の結合律**。

★★4 つ組の添字で 3 つを表せば、`compPf_mk` を 4 回当てるだけ ——
`(idx13).obj ((q123).obj V)` と `(idx12).obj ((q134).obj V)` が
**定義的に等しい**(どちらも脚 `(a, e)`)ことが効く。 -/
theorem compPf_assoc {A B E G : C} (f : HomPf P F A B) (g : HomPf P F B E)
    (h : HomPf P F E G) :
    compPf P F (compPf P F f g) h = compPf P F f (compPf P F g h) := by
  obtain ⟨V, φ, ψ, χ, rfl, rfl, rfl⟩ := exists_rep4 f g h
  -- ★★4 本の等式。**右辺はどれも「次に `rw` したい形」で述べる**(定義的に等しい)。
  have e1 : compPf P F (HomPf.mk ((idx12 P F A B E).obj ((q123 P F A B E G).obj V)) φ)
        (HomPf.mk ((idx23 P F A B E).obj ((q123 P F A B E G).obj V)) ψ)
      = HomPf.mk ((idx12 P F A E G).obj ((q134 P F A B E G).obj V)) (φ ≫ ψ) :=
    compPf_mk ((q123 P F A B E G).obj V) φ ψ
  have e2 : compPf P F (HomPf.mk ((idx23 P F A B E).obj ((q123 P F A B E G).obj V)) ψ)
        (HomPf.mk ((idx23 P F A E G).obj ((q134 P F A B E G).obj V)) χ)
      = HomPf.mk ((idx23 P F A B G).obj ((q124 P F A B E G).obj V)) (ψ ≫ χ) :=
    compPf_mk ((q234 P F A B E G).obj V) ψ χ
  have e3 : compPf P F (HomPf.mk ((idx12 P F A B E).obj ((q123 P F A B E G).obj V)) φ)
        (HomPf.mk ((idx23 P F A B G).obj ((q124 P F A B E G).obj V)) (ψ ≫ χ))
      = HomPf.mk ((idx13 P F A B G).obj ((q124 P F A B E G).obj V)) (φ ≫ ψ ≫ χ) :=
    compPf_mk ((q124 P F A B E G).obj V) φ (ψ ≫ χ)
  have e4 : compPf P F (HomPf.mk ((idx12 P F A E G).obj ((q134 P F A B E G).obj V)) (φ ≫ ψ))
        (HomPf.mk ((idx23 P F A E G).obj ((q134 P F A B E G).obj V)) χ)
      = HomPf.mk ((idx13 P F A B G).obj ((q124 P F A B E G).obj V)) ((φ ≫ ψ) ≫ χ) :=
    compPf_mk ((q134 P F A B E G).obj V) (φ ≫ ψ) χ
  rw [e1, e2, e3, e4, Category.assoc]

/-! ## ★11. 単位律

★★**単位元は `toHomPf (𝟙)`**。★代表元の脚 `b : B → B′` に沿って `𝟙_B` を遷移させると
`𝟙_{B′}` になる(遷移の一意性)—— これが単位律の中身。 -/

variable {P F} in
/-- ★★**右単位律**。 -/
theorem compPf_id_right {A B : C} (f : HomPf P F A B) :
    compPf P F f (toHomPf (F := F) (𝟙 B)) = f := by
  obtain ⟨Z, φ, rfl⟩ := HomPf.exists_rep f
  obtain ⟨ha, hb, hab⟩ := Z.hom.property
  set V : IdxPf3 P F A B B :=
    Under.mk (Y := (⟨(Z.right.obj.1, Z.right.obj.2, Z.right.obj.2)⟩ : TriFr P F))
      (show triFrObj P F A B B ⟶ _ from
        ⟨(Z.hom.hom.1, Z.hom.hom.2, Z.hom.hom.2), ha, hb, hb, hab, rfl⟩) with hV
  have hmor : idxOne P F B B ⟶ (idx23 P F A B B).obj V :=
    Under.homMk (show (idxOne P F B B).right ⟶ _ from ⟨(Z.hom.hom.2, Z.hom.hom.2), hb, hb, rfl⟩)
      (WideSubcategory.hom_ext _ (Prod.ext (Category.id_comp _) (Category.id_comp _)))
  have hw1 : 𝟙 B ≫ hmor.right.hom.1 = Z.hom.hom.2 :=
    congrArg (fun t : biFrObj P F B B ⟶ ((idx23 P F A B B).obj V).right => t.hom.1) (Under.w hmor)
  have hw2 : 𝟙 B ≫ hmor.right.hom.2 = Z.hom.hom.2 :=
    congrArg (fun t : biFrObj P F B B ⟶ ((idx23 P F A B B).obj V).right => t.hom.2) (Under.w hmor)
  have hsq : (𝟙 B : B ⟶ B) ≫ hmor.right.hom.2
      = hmor.right.hom.1 ≫ 𝟙 Z.right.obj.2 := by
    have h1 : hmor.right.hom.1 ≫ 𝟙 Z.right.obj.2 = hmor.right.hom.1 := Category.comp_id _
    have h2 : (𝟙 B : B ⟶ B) ≫ hmor.right.hom.1 = hmor.right.hom.1 := Category.id_comp _
    rw [h1, ← h2, hw1, hw2]
  have htr : idxTransport P F hmor (𝟙 B) = 𝟙 Z.right.obj.2 :=
    frobTransport_eq _ _ _ _ _ _ _ hsq
  have hid : toHomPf (F := F) (𝟙 B) = HomPf.mk ((idx23 P F A B B).obj V) (𝟙 Z.right.obj.2) := by
    rw [← htr]
    exact (HomPf.mk_map hmor (𝟙 B)).symm
  have e : compPf P F (HomPf.mk Z φ) (HomPf.mk ((idx23 P F A B B).obj V) (𝟙 Z.right.obj.2))
      = HomPf.mk Z (φ ≫ 𝟙 Z.right.obj.2) := compPf_mk V φ (𝟙 _)
  rw [hid, e, Category.comp_id]

variable {P F} in
/-- ★★**左単位律**。 -/
theorem compPf_id_left {A B : C} (f : HomPf P F A B) :
    compPf P F (toHomPf (F := F) (𝟙 A)) f = f := by
  obtain ⟨Z, φ, rfl⟩ := HomPf.exists_rep f
  obtain ⟨ha, hb, hab⟩ := Z.hom.property
  set V : IdxPf3 P F A A B :=
    Under.mk (Y := (⟨(Z.right.obj.1, Z.right.obj.1, Z.right.obj.2)⟩ : TriFr P F))
      (show triFrObj P F A A B ⟶ _ from
        ⟨(Z.hom.hom.1, Z.hom.hom.1, Z.hom.hom.2), ha, ha, hb, rfl, hab⟩) with hV
  have hmor : idxOne P F A A ⟶ (idx12 P F A A B).obj V :=
    Under.homMk (show (idxOne P F A A).right ⟶ _ from ⟨(Z.hom.hom.1, Z.hom.hom.1), ha, ha, rfl⟩)
      (WideSubcategory.hom_ext _ (Prod.ext (Category.id_comp _) (Category.id_comp _)))
  have hw1 : 𝟙 A ≫ hmor.right.hom.1 = Z.hom.hom.1 :=
    congrArg (fun t : biFrObj P F A A ⟶ ((idx12 P F A A B).obj V).right => t.hom.1) (Under.w hmor)
  have hw2 : 𝟙 A ≫ hmor.right.hom.2 = Z.hom.hom.1 :=
    congrArg (fun t : biFrObj P F A A ⟶ ((idx12 P F A A B).obj V).right => t.hom.2) (Under.w hmor)
  have hsq : (𝟙 A : A ⟶ A) ≫ hmor.right.hom.2
      = hmor.right.hom.1 ≫ 𝟙 Z.right.obj.1 := by
    have h1 : hmor.right.hom.1 ≫ 𝟙 Z.right.obj.1 = hmor.right.hom.1 := Category.comp_id _
    have h2 : (𝟙 A : A ⟶ A) ≫ hmor.right.hom.1 = hmor.right.hom.1 := Category.id_comp _
    rw [h1, ← h2, hw1, hw2]
  have htr : idxTransport P F hmor (𝟙 A) = 𝟙 Z.right.obj.1 :=
    frobTransport_eq _ _ _ _ _ _ _ hsq
  have hid : toHomPf (F := F) (𝟙 A) = HomPf.mk ((idx12 P F A A B).obj V) (𝟙 Z.right.obj.1) := by
    rw [← htr]
    exact (HomPf.mk_map hmor (𝟙 A)).symm
  have e : compPf P F (HomPf.mk ((idx12 P F A A B).obj V) (𝟙 Z.right.obj.1)) (HomPf.mk Z φ)
      = HomPf.mk Z (𝟙 Z.right.obj.1 ≫ φ) := compPf_mk V (𝟙 _) φ
  rw [hid, e, Category.id_comp]

/-! ## ★12. `𝒞^pf`(`n = 1` の部分)と関手 `𝒞 → 𝒞^pf`

★★原文 (iii) の `𝒞^pf` の対象は**対 `(A, n)`**(`A` の `n` 乗根)である。
★ここではまず **`n = 1` の部分**、すなわち
「対象は `𝒞` の対象、射は perfected morphism」の圏を建てる。
★★**これは関手 `𝒞 → 𝒞^pf` の像**であり、(iii) の本体はこの上に `n` を足すだけになる。 -/

/-- ★★★**`𝒞^pf`(`n = 1` の部分)** —— 対象は `𝒞` の対象、射は perfected morphism。

★`P` / `F` は型に現れないが、圏構造がそれらに依存するので**引数として持たせる**。 -/
def PfCat (_P : PreFrobenioid C Φ) (_F : FrobenioidCore _P) : Type u2 := C

/-- ★`PfCat` の対象を `𝒞` の対象として見る(型は同じだが、`𝟙` の解決先を分けるために要る)。 -/
def pfDown (A : PfCat P F) : C := A

/-- ★★★**`Definition 3.1, (iii)` の圏構造** ——
結合律は `compPf_assoc`、単位律は `compPf_id_left` / `compPf_id_right`。 -/
noncomputable instance pfCatCategory : Category.{max u2 v2} (PfCat P F) where
  Hom A B := HomPf P F (pfDown P F A) (pfDown P F B)
  id A := toHomPf (F := F) (𝟙 (pfDown P F A))
  comp f g := compPf P F f g
  id_comp := compPf_id_left
  comp_id := compPf_id_right
  assoc f g h := compPf_assoc f g h

variable {P F} in
/-- ★添字 `(𝟙_A, 𝟙_B, 𝟙_E)` の 3 つ組。 -/
def idxOne3 (A B E : C) : IdxPf3 P F A B E :=
  Under.mk (Y := (⟨(A, B, E)⟩ : TriFr P F))
    (show triFrObj P F A B E ⟶ _ from
      ⟨(𝟙 A, 𝟙 B, 𝟙 E), isFrobeniusType_of_isIso P (𝟙 A), isFrobeniusType_of_isIso P (𝟙 B),
        isFrobeniusType_of_isIso P (𝟙 E),
        show P.degFr (𝟙 A) = P.degFr (𝟙 B) by rw [P.degFr_id, P.degFr_id],
        show P.degFr (𝟙 B) = P.degFr (𝟙 E) by rw [P.degFr_id, P.degFr_id]⟩)

variable {P F} in
/-- ★★**`toHomPf` は合成を保つ**。 -/
theorem toHomPf_comp {A B E : C} (f : A ⟶ B) (g : B ⟶ E) :
    toHomPf (F := F) (f ≫ g) = compPf P F (toHomPf (F := F) f) (toHomPf (F := F) g) :=
  (compPf_mk (idxOne3 A B E) f g).symm

/-- ★★★**[FrdI] Definition 3.1, (iii)** —— **自然な関手 `𝒞 → 𝒞^pf`**
(原文の「by mapping A → (A, 1)」)。 -/
noncomputable def toPfCat : C ⥤ PfCat P F where
  obj A := A
  map f := toHomPf (F := F) f
  map_id _ := rfl
  map_comp f g := toHomPf_comp f g

/-! ## ★13. **根の不変性** —— `Hom^pf(A,B) ≅ Hom^pf(A′,B′)`

★★原文 (iii) は「one verifies immediately [cf. `Definition 1.3, (ii)`] that this set of
morphisms of `𝒞^pf` from `(A, n)` to `(B, m)` is **independent** [up to uniquely determined
natural bijections] **of the choice of morphisms of Frobenius type** `A → A′`, `B → B′`」
と書く。★★**その "immediately" の中身がこれ**である。

★★**主張**: `a : A → A′`、`b : B → B′` が**同次数の Frobenius 型射**なら
**`Hom^pf(A′,B′) ≅ Hom^pf(A,B)`**。

★★**証明の形**: 添字圏の間の関手 `(A′,B′)𝒞^{bi-Fr} ⥤ (A,B)𝒞^{bi-Fr}`(構造射に `(a,b)` を
前合成する `Under.map`)が **cofinal** であることを示すだけ。
★**`Hom` 集合そのものは変わらない**(添字が動くだけ)ので、
`Under.map (a,b) ⋙ homFunctorPf(A,B) = homFunctorPf(A′,B′)` が**定義的に成り立つ**。
★★**したがって同型は `Functor.Final.colimitIso` の一言で出る。** -/

variable {P F} in
/-- ★同次数の Frobenius 型射の対が定める `𝒞^{bi-Fr}` の射。 -/
def biFrHom {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a) (b : B ⟶ B')
    (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    biFrObj P F A B ⟶ biFrObj P F A' B' := ⟨(a, b), ha, hb, hd⟩

variable {P F} in
/-- ★★**添字圏の押し出し** —— 構造射に `(a,b)` を前合成する。 -/
def pushIdx {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a) (b : B ⟶ B')
    (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    IdxPf P F A' B' ⥤ IdxPf P F A B :=
  Under.map (biFrHom a ha b hb hd)

variable {P F} in
/-- ★★★**押し出しは cofinal** —— これが根の不変性の中身。

★★**構成**: `Z ∈ (A,B)𝒞^{bi-Fr}` の次数を `t`、`(a,b)` の次数を `d` とする。
`A′`, `B′` から次数 `t` を伸ばして `V` を作ると、押し出しの次数は `t·d`。
★`Z` の先から次数 `d` を伸ばせば次数は `d·t` になり、
`Definition 1.3, (ii)` の**本質的一意性**で両者を合わせられる。
★★**ここでも効くのは `mul_comm` である。** -/
instance pushIdx_final {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a) (b : B ⟶ B')
    (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    (pushIdx (F := F) a ha b hb hd).Final := by
  refine Functor.final_of_exists_of_isFiltered _ ?_ ?_
  · intro Z
    obtain ⟨hc, he, hce⟩ := Z.hom.property
    obtain ⟨A₃, α, hα, hαd⟩ := F.frobDegSurj A' (P.degFr Z.hom.hom.1)
    obtain ⟨B₃, β, hβ, hβd⟩ := F.frobDegSurj B' (P.degFr Z.hom.hom.1)
    obtain ⟨A₄, γ₀, hγ₀, hγ₀d⟩ := F.frobDegSurj Z.right.obj.1 (P.degFr a)
    obtain ⟨B₄, δ₀, hδ₀, hδ₀d⟩ := F.frobDegSurj Z.right.obj.2 (P.degFr a)
    have h1 : P.degFr (Z.hom.hom.1 ≫ γ₀) = P.degFr (a ≫ α) := by
      rw [P.degFr_comp, P.degFr_comp, hγ₀d, hαd, mul_comm]
    obtain ⟨θ, hθ, hθe⟩ := F.frobDegUniq A A₄ A₃ (Z.hom.hom.1 ≫ γ₀) (a ≫ α)
      (IsFrobeniusType.comp P F hc hγ₀) (IsFrobeniusType.comp P F ha hα) h1
    haveI : IsIso θ := hθ
    have h2 : P.degFr (Z.hom.hom.2 ≫ δ₀) = P.degFr (b ≫ β) := by
      rw [P.degFr_comp, P.degFr_comp, hδ₀d, hβd, hd, hce, mul_comm]
    obtain ⟨ζ, hζ, hζe⟩ := F.frobDegUniq B B₄ B₃ (Z.hom.hom.2 ≫ δ₀) (b ≫ β)
      (IsFrobeniusType.comp P F he hδ₀) (IsFrobeniusType.comp P F hb hβ) h2
    haveI : IsIso ζ := hζ
    refine ⟨Under.mk (Y := (⟨(A₃, B₃)⟩ : BiFr P F))
        (show biFrObj P F A' B' ⟶ _ from
          ⟨(α, β), hα, hβ, show P.degFr α = P.degFr β by rw [hαd, hβd]⟩),
      ⟨Under.homMk (show Z.right ⟶ _ from
        ⟨(γ₀ ≫ θ, δ₀ ≫ ζ), IsFrobeniusType.comp P F hγ₀ (isFrobeniusType_of_isIso P θ),
          IsFrobeniusType.comp P F hδ₀ (isFrobeniusType_of_isIso P ζ),
          show P.degFr (γ₀ ≫ θ) = P.degFr (δ₀ ≫ ζ) by
            rw [P.degFr_comp, P.degFr_comp, hγ₀d, hδ₀d,
              show P.degFr θ = 1 from isLinear_of_isIso P θ,
              show P.degFr ζ = 1 from isLinear_of_isIso P ζ]⟩)
        (WideSubcategory.hom_ext _ (Prod.ext ?_ ?_))⟩⟩
    · show Z.hom.hom.1 ≫ γ₀ ≫ θ = a ≫ α
      rw [← Category.assoc]
      exact hθe
    · show Z.hom.hom.2 ≫ δ₀ ≫ ζ = b ≫ β
      rw [← Category.assoc]
      exact hζe
  · intro _ _ s s'
    exact ⟨_, 𝟙 _, by rw [idx_hom_ext s s']⟩

variable {P F} in
/-- ★★★**[FrdI] Definition 3.1, (iii) の「independent of the choice」** ——
**根の不変性**: 同次数の Frobenius 型射の対に沿って `Hom^pf` は変わらない。 -/
noncomputable def rootIso {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    HomPf P F A' B' ≅ HomPf P F A B :=
  Functor.Final.colimitIso (pushIdx (F := F) a ha b hb hd) (homFunctorPf P F A B)

variable {P F} in
/-- ★★**根の不変性の計算則** —— 代表元では**射そのものは変わらず、添字だけ押し出される**。 -/
theorem rootIso_hom_mk {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (V : IdxPf P F A' B') (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    (rootIso (F := F) a ha b hb hd).hom (HomPf.mk V φ)
      = HomPf.mk ((pushIdx (F := F) a ha b hb hd).obj V) φ := by
  show (Functor.Final.colimitIso (pushIdx (F := F) a ha b hb hd) (homFunctorPf P F A B)).hom
      (colimit.ι (pushIdx (F := F) a ha b hb hd ⋙ homFunctorPf P F A B) V (ULift.up φ)) = _
  rw [← types_comp_apply
    (colimit.ι (pushIdx (F := F) a ha b hb hd ⋙ homFunctorPf P F A B) V)
    (Functor.Final.colimitIso (pushIdx (F := F) a ha b hb hd) (homFunctorPf P F A B)).hom,
    Functor.Final.ι_colimitIso_hom]
  rfl

/-! ### ★3 つ組での押し出し —— 根の不変性が**合成と可換**であること -/

variable {P F} in
/-- ★同次数の Frobenius 型射の 3 つ組が定める `𝒞^{tri-Fr}` の射。 -/
def triFrHom {A B E A' B' E' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (c : E ⟶ E') (hc : IsFrobeniusType P c)
    (hab : P.degFr a = P.degFr b) (hbc : P.degFr b = P.degFr c) :
    triFrObj P F A B E ⟶ triFrObj P F A' B' E' := ⟨(a, b, c), ha, hb, hc, hab, hbc⟩

variable {P F} in
/-- ★★3 つ組の添字圏の押し出し。 -/
def pushIdx3 {A B E A' B' E' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (c : E ⟶ E') (hc : IsFrobeniusType P c)
    (hab : P.degFr a = P.degFr b) (hbc : P.degFr b = P.degFr c) :
    IdxPf3 P F A' B' E' ⥤ IdxPf3 P F A B E :=
  Under.map (triFrHom a ha b hb c hc hab hbc)

variable {P F} in
/-- ★★★**根の不変性は合成と可換**。

★★**証明は 3 行**である —— 3 つ組の押し出し `pushIdx3` の 3 つの射影が、
2 つ組の押し出し `pushIdx` と**定義的に一致する**(どちらも「脚に `(a,b,c)` を前合成」)から。 -/
theorem rootIso_comp {A B E A' B' E' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (c : E ⟶ E') (hc : IsFrobeniusType P c)
    (hab : P.degFr a = P.degFr b) (hbc : P.degFr b = P.degFr c)
    (V : IdxPf3 P F A' B' E') (φ : V.right.obj.1 ⟶ V.right.obj.2.1)
    (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2) :
    compPf P F ((rootIso (F := F) a ha b hb hab).hom
        (HomPf.mk ((idx12 P F A' B' E').obj V) φ))
      ((rootIso (F := F) b hb c hc hbc).hom (HomPf.mk ((idx23 P F A' B' E').obj V) ψ))
      = (rootIso (F := F) a ha c hc (hab.trans hbc)).hom
        (HomPf.mk ((idx13 P F A' B' E').obj V) (φ ≫ ψ)) := by
  rw [rootIso_hom_mk, rootIso_hom_mk, rootIso_hom_mk]
  exact compPf_mk ((pushIdx3 (F := F) a ha b hb c hc hab hbc).obj V) φ ψ

variable {P F} in
/-- ★★**2 つの perfected morphism を 1 つの 3 つ組の添字で表す**。 -/
theorem exists_rep3 {A B E : C} (f : HomPf P F A B) (g : HomPf P F B E) :
    ∃ (V : IdxPf3 P F A B E) (φ : V.right.obj.1 ⟶ V.right.obj.2.1)
      (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2),
      f = HomPf.mk ((idx12 P F A B E).obj V) φ ∧
      g = HomPf.mk ((idx23 P F A B E).obj V) ψ := by
  obtain ⟨Zf, φ₀, rfl⟩ := HomPf.exists_rep f
  obtain ⟨Zg, ψ₀, rfl⟩ := HomPf.exists_rep g
  obtain ⟨hfa, hfb, hfab⟩ := Zf.hom.property
  obtain ⟨hgb, hge, hgbe⟩ := Zg.hom.property
  obtain ⟨E₁, e₁, he₁, he₁d⟩ := F.frobDegSurj E (P.degFr Zf.hom.hom.1)
  obtain ⟨A₂, a₂, ha₂, ha₂d⟩ := F.frobDegSurj A (P.degFr Zg.hom.hom.1)
  set Vf : IdxPf3 P F A B E :=
    Under.mk (Y := (⟨(Zf.right.obj.1, Zf.right.obj.2, E₁)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(Zf.hom.hom.1, Zf.hom.hom.2, e₁), hfa, hfb, he₁, hfab,
          hfab.symm.trans he₁d.symm⟩) with hVf
  set Vg : IdxPf3 P F A B E :=
    Under.mk (Y := (⟨(A₂, Zg.right.obj.1, Zg.right.obj.2)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(a₂, Zg.hom.hom.1, Zg.hom.hom.2), ha₂, hgb, hge, ha₂d, hgbe⟩) with hVg
  exact ⟨IsFiltered.max Vf Vg,
    idxTransport P F ((idx12 P F A B E).map ((IsFiltered.leftToMax Vf Vg))) φ₀,
    idxTransport P F ((idx23 P F A B E).map ((IsFiltered.rightToMax Vf Vg))) ψ₀,
    (HomPf.mk_map _ φ₀).symm, (HomPf.mk_map _ ψ₀).symm⟩

variable {P F} in
/-- ★★★**根の不変性は合成と可換**(一般形)。 -/
theorem rootIso_comp' {A B E A' B' E' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (c : E ⟶ E') (hc : IsFrobeniusType P c)
    (hab : P.degFr a = P.degFr b) (hbc : P.degFr b = P.degFr c)
    (f : HomPf P F A' B') (g : HomPf P F B' E') :
    compPf P F ((rootIso (F := F) a ha b hb hab).hom f)
        ((rootIso (F := F) b hb c hc hbc).hom g)
      = (rootIso (F := F) a ha c hc (hab.trans hbc)).hom (compPf P F f g) := by
  obtain ⟨V, φ, ψ, rfl, rfl⟩ := exists_rep3 f g
  rw [compPf_mk V φ ψ]
  exact rootIso_comp a ha b hb c hc hab hbc V φ ψ

variable {P F} in
/-- ★★同上、逆向き。 -/
theorem rootIso_inv_comp {A B E A' B' E' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (c : E ⟶ E') (hc : IsFrobeniusType P c)
    (hab : P.degFr a = P.degFr b) (hbc : P.degFr b = P.degFr c)
    (f : HomPf P F A B) (g : HomPf P F B E) :
    compPf P F ((rootIso (F := F) a ha b hb hab).inv f)
        ((rootIso (F := F) b hb c hc hbc).inv g)
      = (rootIso (F := F) a ha c hc (hab.trans hbc)).inv (compPf P F f g) := by
  have h := rootIso_comp' a ha b hb c hc hab hbc
    ((rootIso (F := F) a ha b hb hab).inv f) ((rootIso (F := F) b hb c hc hbc).inv g)
  rw [Iso.inv_hom_id_apply, Iso.inv_hom_id_apply] at h
  rw [h, Iso.hom_inv_id_apply]

variable {P F} in
/-- ★★**根の不変性は恒等射を保つ**。 -/
theorem rootIso_hom_id {A A' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a) :
    (rootIso (F := F) a ha a ha rfl).hom (toHomPf (F := F) (𝟙 A'))
      = toHomPf (F := F) (𝟙 A) := by
  have hmor : idxOne P F A A ⟶ (pushIdx (F := F) a ha a ha rfl).obj (idxOne P F A' A') :=
    Under.homMk (show (idxOne P F A A).right ⟶ _ from ⟨(a, a), ha, ha, rfl⟩)
      (WideSubcategory.hom_ext _ (Prod.ext
        ((Category.id_comp a).trans (Category.comp_id a).symm)
        ((Category.id_comp a).trans (Category.comp_id a).symm)))
  have hw1 : 𝟙 A ≫ hmor.right.hom.1 = a ≫ 𝟙 A' :=
    congrArg (fun t : biFrObj P F A A ⟶
      ((pushIdx (F := F) a ha a ha rfl).obj (idxOne P F A' A')).right => t.hom.1) (Under.w hmor)
  have hw2 : 𝟙 A ≫ hmor.right.hom.2 = a ≫ 𝟙 A' :=
    congrArg (fun t : biFrObj P F A A ⟶
      ((pushIdx (F := F) a ha a ha rfl).obj (idxOne P F A' A')).right => t.hom.2) (Under.w hmor)
  have hsq : (𝟙 A : A ⟶ A) ≫ hmor.right.hom.2 = hmor.right.hom.1 ≫ 𝟙 A' := by
    have h2 : (𝟙 A : A ⟶ A) ≫ hmor.right.hom.1 = hmor.right.hom.1 := Category.id_comp _
    rw [← h2, hw1, hw2]
    exact (Category.comp_id _).symm
  have htr : idxTransport P F hmor (𝟙 A) = 𝟙 A' := frobTransport_eq _ _ _ _ _ _ _ hsq
  rw [toHomPf, toHomPf, rootIso_hom_mk, ← htr, HomPf.mk_map hmor (𝟙 A)]

variable {P F} in
/-- ★添字対象が等しければ代表元も等しい(依存型の輸送を 1 か所に閉じ込める)。 -/
theorem HomPf.mk_congr {A B : C} {Z₁ Z₂ : IdxPf P F A B} (h : Z₁ = Z₂)
    {φ : Z₁.right.obj.1 ⟶ Z₁.right.obj.2} {ψ : Z₂.right.obj.1 ⟶ Z₂.right.obj.2}
    (hφ : HEq φ ψ) : HomPf.mk Z₁ φ = HomPf.mk Z₂ ψ := by
  subst h
  rw [eq_of_heq hφ]

variable {P F} in
/-- ★★**根の不変性の 2 段合成は 1 段**。

★★添字の押し出しは「構造射に前合成する」だけなので、
2 段合成と 1 段の違いは **`Category.assoc` だけ**である。 -/
theorem rootIso_trans {A B A' B' A'' B'' : C}
    (a : A ⟶ A') (ha : IsFrobeniusType P a) (b : B ⟶ B') (hb : IsFrobeniusType P b)
    (hd : P.degFr a = P.degFr b)
    (a' : A' ⟶ A'') (ha' : IsFrobeniusType P a') (b' : B' ⟶ B'') (hb' : IsFrobeniusType P b')
    (hd' : P.degFr a' = P.degFr b')
    (hac : IsFrobeniusType P (a ≫ a')) (hbc : IsFrobeniusType P (b ≫ b'))
    (hdc : P.degFr (a ≫ a') = P.degFr (b ≫ b')) (z : HomPf P F A'' B'') :
    (rootIso (F := F) a ha b hb hd).hom ((rootIso (F := F) a' ha' b' hb' hd').hom z)
      = (rootIso (F := F) (a ≫ a') hac (b ≫ b') hbc hdc).hom z := by
  refine HomColim.induction (P := fun z =>
    (rootIso (F := F) a ha b hb hd).hom ((rootIso (F := F) a' ha' b' hb' hd').hom z)
      = (rootIso (F := F) (a ≫ a') hac (b ≫ b') hbc hdc).hom z)
    (homFunctorPf P F A'' B'') (fun V x => ?_) z
  show (rootIso (F := F) a ha b hb hd).hom
      ((rootIso (F := F) a' ha' b' hb' hd').hom (HomPf.mk V x.down))
    = (rootIso (F := F) (a ≫ a') hac (b ≫ b') hbc hdc).hom (HomPf.mk V x.down)
  rw [rootIso_hom_mk, rootIso_hom_mk, rootIso_hom_mk]
  exact HomPf.mk_congr (congrArg (Comma.mk _ _) (Category.assoc _ _ _).symm) HEq.rfl

/-! ## ★14. 選んだ根 —— `A^{(d)}` と持ち上げ

★★原文 (iii) の対象は対 `(A, n)`(「`A` の `n` 乗根」)である。
★射を書くには「次数 `d` の Frobenius 拡大」を**選ぶ**必要がある
(`Definition 1.3, (ii)` の存在から選択公理で取る)。
★★**選び方に依らないこと**は `rootIso`(根の不変性)が保証する。 -/

/-- ★★選んだ次数 `d` の Frobenius 拡大の終域 `A^{(d)}`。 -/
noncomputable def rtObj (A : C) (d : ℕ+) : C := (F.frobDegSurj A d).choose

/-- ★★選んだ次数 `d` の Frobenius 型射 `A → A^{(d)}`。 -/
noncomputable def rtExt (A : C) (d : ℕ+) : A ⟶ rtObj P F A d :=
  (F.frobDegSurj A d).choose_spec.choose

theorem rtExt_frobType (A : C) (d : ℕ+) : IsFrobeniusType P (rtExt P F A d) :=
  (F.frobDegSurj A d).choose_spec.choose_spec.1

theorem rtExt_degFr (A : C) (d : ℕ+) : P.degFr (rtExt P F A d) = d :=
  (F.frobDegSurj A d).choose_spec.choose_spec.2

/-- ★★**根の持ち上げの存在** —— `A^{(d)} → A^{(t)}`(`t = e·d`)が次数 `e` で、
しかも `A` からの構造射と両立する形で取れる。 -/
theorem rtLift_exists (A : C) {d e t : ℕ+} (h : t = e * d) :
    ∃ γ : rtObj P F A d ⟶ rtObj P F A t, IsFrobeniusType P γ ∧ P.degFr γ = e ∧
      rtExt P F A d ≫ γ = rtExt P F A t := by
  obtain ⟨X, γ₀, hγ₀, hγ₀d⟩ := F.frobDegSurj (rtObj P F A d) e
  have hdeg : P.degFr (rtExt P F A d ≫ γ₀) = P.degFr (rtExt P F A t) := by
    rw [P.degFr_comp, hγ₀d, rtExt_degFr, rtExt_degFr, h]
  obtain ⟨θ, hθ, hθe⟩ := F.frobDegUniq A X (rtObj P F A t) (rtExt P F A d ≫ γ₀) (rtExt P F A t)
    (IsFrobeniusType.comp P F (rtExt_frobType P F A d) hγ₀) (rtExt_frobType P F A t) hdeg
  haveI : IsIso θ := hθ
  refine ⟨γ₀ ≫ θ, IsFrobeniusType.comp P F hγ₀ (isFrobeniusType_of_isIso P θ),
    by rw [P.degFr_comp, show P.degFr θ = 1 from isLinear_of_isIso P θ, hγ₀d, one_mul], ?_⟩
  rw [← Category.assoc]
  exact hθe

/-- ★★**根の持ち上げ** `A^{(d)} → A^{(t)}`(`t = e·d`、次数 `e`)。 -/
noncomputable def rtLift (A : C) {d e t : ℕ+} (h : t = e * d) :
    rtObj P F A d ⟶ rtObj P F A t :=
  (rtLift_exists P F A h).choose

theorem rtLift_frobType (A : C) {d e t : ℕ+} (h : t = e * d) :
    IsFrobeniusType P (rtLift P F A h) := (rtLift_exists P F A h).choose_spec.1

theorem rtLift_degFr (A : C) {d e t : ℕ+} (h : t = e * d) :
    P.degFr (rtLift P F A h) = e := (rtLift_exists P F A h).choose_spec.2.1

theorem rtLift_ext (A : C) {d e t : ℕ+} (h : t = e * d) :
    rtExt P F A d ≫ rtLift P F A h = rtExt P F A t := (rtLift_exists P F A h).choose_spec.2.2

/-- ★★**持ち上げは一意** —— `A^{(d)}` からの射で構造射と両立するものは高々 1 本。

★理由は **`𝒞` が totally epimorphic**(構造射 `A → A^{(d)}` が epi)。 -/
theorem rtLift_uniq {A : C} {d t : ℕ+} (γ δ : rtObj P F A d ⟶ rtObj P F A t)
    (hγ : rtExt P F A d ≫ γ = rtExt P F A t) (hδ : rtExt P F A d ≫ δ = rtExt P F A t) :
    γ = δ := by
  haveI : Epi (rtExt P F A d) := P.totEpiC _ _ _
  exact (cancel_epi (rtExt P F A d)).mp (hγ.trans hδ.symm)

/-- ★★**2 段の持ち上げは 1 段の持ち上げ**(一意性から)。 -/
theorem rtLift_comp (A : C) {d e t e' s : ℕ+} (h : t = e * d) (h' : s = e' * t)
    {e'' : ℕ+} (h'' : s = e'' * d) :
    rtLift P F A h ≫ rtLift P F A h' = rtLift P F A h'' :=
  rtLift_uniq P F _ _ (by rw [← Category.assoc, rtLift_ext, rtLift_ext]) (rtLift_ext P F A h'')

/-- ★★**根の不変性(選んだ根の版)**。 -/
noncomputable def rtRootIso (A B : C) {dA dB e tA tB : ℕ+} (hA : tA = e * dA)
    (hB : tB = e * dB) :
    HomPf P F (rtObj P F A tA) (rtObj P F B tB) ≅ HomPf P F (rtObj P F A dA) (rtObj P F B dB) :=
  rootIso (F := F) (rtLift P F A hA) (rtLift_frobType P F A hA) (rtLift P F B hB)
    (rtLift_frobType P F B hB) (by rw [rtLift_degFr, rtLift_degFr])

variable {P F} in
/-- ★選んだ根の版の計算則。 -/
theorem rtRootIso_hom_mk (A B : C) {dA dB e tA tB : ℕ+} (hA : tA = e * dA) (hB : tB = e * dB)
    (V : IdxPf P F (rtObj P F A tA) (rtObj P F B tB))
    (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    (rtRootIso P F A B hA hB).hom (HomPf.mk V φ)
      = HomPf.mk ((pushIdx (F := F) (rtLift P F A hA) (rtLift_frobType P F A hA)
        (rtLift P F B hB) (rtLift_frobType P F B hB)
        (by rw [rtLift_degFr, rtLift_degFr])).obj V) φ :=
  rootIso_hom_mk _ _ _ _ _ V φ

variable {P F} in
/-- ★★**選んだ根の版でも 2 段合成は 1 段**(`rtLift_comp` から)。 -/
theorem rtRootIso_trans (A B : C) {dA dB e tA tB e' sA sB e'' : ℕ+}
    (hA : tA = e * dA) (hB : tB = e * dB) (hA' : sA = e' * tA) (hB' : sB = e' * tB)
    (hA'' : sA = e'' * dA) (hB'' : sB = e'' * dB)
    (z : HomPf P F (rtObj P F A sA) (rtObj P F B sB)) :
    (rtRootIso P F A B hA hB).hom ((rtRootIso P F A B hA' hB').hom z)
      = (rtRootIso P F A B hA'' hB'').hom z := by
  refine HomColim.induction (P := fun z =>
    (rtRootIso P F A B hA hB).hom ((rtRootIso P F A B hA' hB').hom z)
      = (rtRootIso P F A B hA'' hB'').hom z)
    (homFunctorPf P F (rtObj P F A sA) (rtObj P F B sB)) (fun V x => ?_) z
  show (rtRootIso P F A B hA hB).hom ((rtRootIso P F A B hA' hB').hom (HomPf.mk V x.down))
    = (rtRootIso P F A B hA'' hB'').hom (HomPf.mk V x.down)
  rw [rtRootIso_hom_mk, rtRootIso_hom_mk, rtRootIso_hom_mk]
  refine HomPf.mk_congr (congrArg (Comma.mk _ _) ?_) HEq.rfl
  refine Eq.trans (Category.assoc _ _ _).symm (congrArg (fun t => t ≫ V.hom) ?_)
  exact WideSubcategory.hom_ext _
    (Prod.ext (rtLift_comp P F A hA hA' hA'') (rtLift_comp P F B hB hB' hB''))

/-! ## ★15. 原文 (iii) の対象 `(A, n)`

★★原文:
「The objects of `𝒞^pf` are pairs `(A, n)`, where `A ∈ Ob(𝒞)`, `n ∈ ℕ≥1`.
The morphisms of `𝒞^pf` are given by
`Hom_{𝒞^pf}((A,n),(B,m)) := Hom^pf_𝒞(A′,B′)`
where `A → A′` is a morphism of Frobenius type of Frobenius degree `m`;
`B → B′` is a morphism of Frobenius type of Frobenius degree `n`.」

★★**次数が入れ替わる**(`(A,n)` 側に `m`、`(B,m)` 側に `n`)のは、
`(A,n)` を「`A` の `n` 乗根」と読むからである ——
どちらも `n·m` 乗根に揃えていることになる。 -/

/-- ★★★**[FrdI] Definition 3.1, (iii)** —— `𝒞^pf` の対象 `(A, n)`。 -/
structure PfRootObj (_P : PreFrobenioid C Φ) (_F : FrobenioidCore _P) where
  /-- 台となる `𝒞` の対象。 -/
  obj : C
  /-- 何乗根か。 -/
  root : ℕ+

/-- ★★★**[FrdI] Definition 3.1, (iii)** —— `𝒞^pf` の射の集合。 -/
noncomputable def HomRoot (X Y : PfRootObj P F) : Type (max u2 v2) :=
  HomPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root)

/-- ★恒等射。 -/
noncomputable def idRoot (X : PfRootObj P F) : HomRoot P F X X :=
  toHomPf (F := F) (𝟙 (rtObj P F X.obj X.root))

/-- ★★★**合成** —— 3 つを共通の根へ持ち上げてから `compPf`、最後に降ろす。

★★**根の指数の並べ方が要点**である。`X = (A,n)`, `Y = (B,m)`, `Z = (E,k)` として
- `f : Hom^pf(A^{(m)}, B^{(n)})` を次数 `k` で持ち上げ、`Hom^pf(A^{(k·m)}, B^{(n·k)})` へ
- `g : Hom^pf(B^{(k)}, E^{(m)})` を次数 `n` で持ち上げ、`Hom^pf(B^{(n·k)}, E^{(n·m)})` へ
- 合成して `Hom^pf(A^{(k·m)}, E^{(n·m)})`、これを次数 `m` の持ち上げで降ろす。

★★**真ん中が `n·k` で syntactic に一致する**ように指数を選んである
(`rtLift` が「目標の指数」を引数に取るので、`mul_comm` は**証明の側**に押し込める)。 -/
noncomputable def compRoot {X Y Z : PfRootObj P F} (f : HomRoot P F X Y) (g : HomRoot P F Y Z) :
    HomRoot P F X Z :=
  (rtRootIso P F X.obj Z.obj
      (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
      (show Y.root * X.root = Y.root * X.root from rfl)).hom
    (compPf P F
      ((rtRootIso P F X.obj Y.obj
        (show Z.root * Y.root = Z.root * Y.root from rfl)
        (show Z.root * X.root = Z.root * X.root from rfl)).inv f)
      ((rtRootIso P F Y.obj Z.obj
        (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g))

/-! ### ★単位律 -/

/-- ★★選んだ根の版でも、根の不変性は恒等射を保つ。 -/
theorem rtRootIso_hom_id (A : C) {d e t : ℕ+} (h : t = e * d) :
    (rtRootIso P F A A h h).hom (toHomPf (F := F) (𝟙 (rtObj P F A t)))
      = toHomPf (F := F) (𝟙 (rtObj P F A d)) :=
  rootIso_hom_id (rtLift P F A h) (rtLift_frobType P F A h)

theorem rtRootIso_inv_id (A : C) {d e t : ℕ+} (h : t = e * d) :
    (rtRootIso P F A A h h).inv (toHomPf (F := F) (𝟙 (rtObj P F A d)))
      = toHomPf (F := F) (𝟙 (rtObj P F A t)) := by
  rw [← rtRootIso_hom_id P F A h, Iso.hom_inv_id_apply]

variable {P F} in
/-- ★★★**右単位律**(対象 `(A,n)` の版)。 -/
theorem compRoot_id_right {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    compRoot P F f (idRoot P F Y) = f := by
  unfold compRoot idRoot
  rw [rtRootIso_inv_id, compPf_id_right, Iso.inv_hom_id_apply]

variable {P F} in
/-- ★★★**左単位律**(対象 `(A,n)` の版)。 -/
theorem compRoot_id_left {X Y : PfRootObj P F} (f : HomRoot P F X Y) :
    compRoot P F (idRoot P F X) f = f := by
  unfold compRoot idRoot
  rw [rtRootIso_inv_id, compPf_id_left, Iso.inv_hom_id_apply]

/-! ### ★★合成は「もっと大きい根」で計算してもよい

★★★**結合律の鍵**である。`compRoot` の定義は特定の指数の取り方
(`Z.root` / `X.root` / `Y.root` で持ち上げる)を使うが、
★**その共通倍数ならどの根で計算しても同じ**であることを一度示しておけば、
結合律は「両辺を `n·m·k·l` の根で計算して `compPf_assoc`」で終わる。 -/

variable {P F} in
/-- ★★★**合成は共通倍数の根で計算してよい**。

`c` は「定義の根の何倍か」、`PA` / `PB` / `PE` はその大きい根、
`ef` / `eg` / `er` は各持ち上げの次数である。★**すべて項として自由に選べる**
(等式を仮定として受け取るので、`mul_comm` / `mul_assoc` は呼び出し側に押し込める)。 -/
theorem compRoot_eq_lift {X Y Z : PfRootObj P F} (f : HomRoot P F X Y) (g : HomRoot P F Y Z)
    {c PA PB PE : ℕ+}
    (hcA : PA = c * (Z.root * Y.root)) (hcB : PB = c * (Z.root * X.root))
    (hcE : PE = c * (Y.root * X.root))
    {ef eg er : ℕ+}
    (hfA : PA = ef * Y.root) (hfB : PB = ef * X.root)
    (hgA : PB = eg * Z.root) (hgE : PE = eg * Y.root)
    (hrA : PA = er * Z.root) (hrE : PE = er * X.root) :
    compRoot P F f g
      = (rtRootIso P F X.obj Z.obj hrA hrE).hom
        (compPf P F ((rtRootIso P F X.obj Y.obj hfA hfB).inv f)
          ((rtRootIso P F Y.obj Z.obj hgA hgE).inv g)) := by
  have hf := rtRootIso_trans (F := F) X.obj Y.obj
    (show Z.root * Y.root = Z.root * Y.root from rfl)
    (show Z.root * X.root = Z.root * X.root from rfl) hcA hcB hfA hfB
    ((rtRootIso P F X.obj Y.obj hfA hfB).inv f)
  rw [Iso.inv_hom_id_apply] at hf
  have hf' : (rtRootIso P F X.obj Y.obj
        (show Z.root * Y.root = Z.root * Y.root from rfl)
        (show Z.root * X.root = Z.root * X.root from rfl)).inv f
      = (rtRootIso P F X.obj Y.obj hcA hcB).hom
        ((rtRootIso P F X.obj Y.obj hfA hfB).inv f) := by
    conv_lhs => rw [← hf]
    exact Iso.hom_inv_id_apply _ _
  have hg := rtRootIso_trans (F := F) Y.obj Z.obj
    (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
    (show Y.root * X.root = X.root * Y.root from mul_comm _ _) hcB hcE hgA hgE
    ((rtRootIso P F Y.obj Z.obj hgA hgE).inv g)
  rw [Iso.inv_hom_id_apply] at hg
  have hg' : (rtRootIso P F Y.obj Z.obj
        (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g
      = (rtRootIso P F Y.obj Z.obj hcB hcE).hom
        ((rtRootIso P F Y.obj Z.obj hgA hgE).inv g) := by
    conv_lhs => rw [← hg]
    exact Iso.hom_inv_id_apply _ _
  have hcomp : compPf P F
        ((rtRootIso P F X.obj Y.obj hcA hcB).hom
          ((rtRootIso P F X.obj Y.obj hfA hfB).inv f))
        ((rtRootIso P F Y.obj Z.obj hcB hcE).hom
          ((rtRootIso P F Y.obj Z.obj hgA hgE).inv g))
      = (rtRootIso P F X.obj Z.obj hcA hcE).hom
        (compPf P F ((rtRootIso P F X.obj Y.obj hfA hfB).inv f)
          ((rtRootIso P F Y.obj Z.obj hgA hgE).inv g)) :=
    rootIso_comp' (F := F)
      (rtLift P F X.obj hcA) (rtLift_frobType P F X.obj hcA)
      (rtLift P F Y.obj hcB) (rtLift_frobType P F Y.obj hcB)
      (rtLift P F Z.obj hcE) (rtLift_frobType P F Z.obj hcE)
      (by rw [rtLift_degFr, rtLift_degFr]) (by rw [rtLift_degFr, rtLift_degFr])
      ((rtRootIso P F X.obj Y.obj hfA hfB).inv f)
      ((rtRootIso P F Y.obj Z.obj hgA hgE).inv g)
  have hr := rtRootIso_trans (F := F) X.obj Z.obj
    (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
    (show Y.root * X.root = Y.root * X.root from rfl) hcA hcE hrA hrE
    (compPf P F ((rtRootIso P F X.obj Y.obj hfA hfB).inv f)
      ((rtRootIso P F Y.obj Z.obj hgA hgE).inv g))
  unfold compRoot
  rw [hf', hg', hcomp, hr]

variable {P F} in
set_option linter.unusedSimpArgs false in
/-- ★★★**結合律**(対象 `(A,n)` の版)。

★★**両辺を「level `n·m·k·l`」の根で計算する**と、どちらも同じ
`rtRootIso` の `.hom` を外側に持ち、中身が `compPf` の 2 通りの括弧付けになる。
★あとは `compPf_assoc` 1 本。

★★**指数の項の選び方が要点**である —— 4 回の `compRoot_eq_lift` で
「内側の結果の同型」と「外側の持ち上げの同型」が**同じ項になる**ように
`ef` / `eg` / `er` を選んである(そこが合わないと `Iso.hom_inv_id_apply` が当たらない)。 -/
theorem compRoot_assoc {X Y Z W : PfRootObj P F} (f : HomRoot P F X Y) (g : HomRoot P F Y Z)
    (h : HomRoot P F Z W) :
    compRoot P F (compRoot P F f g) h = compRoot P F f (compRoot P F g h) := by
  have hLin := compRoot_eq_lift f g (c := W.root)
    (PA := W.root * (Z.root * Y.root)) (PB := W.root * (Z.root * X.root))
    (PE := X.root * (W.root * Y.root))
    (hcA := rfl) (hcB := rfl) (hcE := by simp [mul_comm, mul_assoc, mul_left_comm])
    (ef := W.root * Z.root) (eg := W.root * X.root) (er := W.root * Y.root)
    (hfA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hfB := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hgA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hgE := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hrA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hrE := by simp [mul_comm, mul_assoc, mul_left_comm])
  have hLout := compRoot_eq_lift (compRoot P F f g) h (c := Y.root)
    (PA := W.root * (Z.root * Y.root)) (PB := X.root * (W.root * Y.root))
    (PE := X.root * (Z.root * Y.root))
    (hcA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hcB := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hcE := by simp [mul_comm, mul_assoc, mul_left_comm])
    (ef := W.root * Y.root) (eg := X.root * Y.root) (er := Z.root * Y.root)
    (hfA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hfB := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hgA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hgE := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hrA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hrE := by simp [mul_comm, mul_assoc, mul_left_comm])
  have hRin := compRoot_eq_lift g h (c := X.root)
    (PA := W.root * (Z.root * X.root)) (PB := X.root * (W.root * Y.root))
    (PE := X.root * (Z.root * Y.root))
    (hcA := by simp [mul_comm, mul_assoc, mul_left_comm]) (hcB := rfl) (hcE := rfl)
    (ef := W.root * X.root) (eg := X.root * Y.root) (er := Z.root * X.root)
    (hfA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hfB := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hgA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hgE := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hrA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hrE := by simp [mul_comm, mul_assoc, mul_left_comm])
  have hRout := compRoot_eq_lift f (compRoot P F g h) (c := Z.root)
    (PA := W.root * (Z.root * Y.root)) (PB := W.root * (Z.root * X.root))
    (PE := X.root * (Z.root * Y.root))
    (hcA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hcB := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hcE := by simp [mul_comm, mul_assoc, mul_left_comm])
    (ef := W.root * Z.root) (eg := Z.root * X.root) (er := Z.root * Y.root)
    (hfA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hfB := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hgA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hgE := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hrA := by simp [mul_comm, mul_assoc, mul_left_comm])
    (hrE := by simp [mul_comm, mul_assoc, mul_left_comm])
  rw [hLout, hRout, hLin, hRin, Iso.hom_inv_id_apply, Iso.hom_inv_id_apply, compPf_assoc]

/-- ★★★**[FrdI] Definition 3.1, (iii)** —— **`𝒞^pf`(原文どおり、対象は `(A, n)`)**。

★結合律は `compRoot_assoc`、単位律は `compRoot_id_left` / `compRoot_id_right`。 -/
noncomputable instance pfRootCategory : Category.{max u2 v2} (PfRootObj P F) where
  Hom X Y := HomRoot P F X Y
  id X := idRoot P F X
  comp f g := compRoot P F f g
  id_comp := compRoot_id_left
  comp_id := compRoot_id_right
  assoc f g h := compRoot_assoc f g h

/-- ★★**次数 1 の Frobenius 型射は同型** —— `Proposition 1.4, (iii)` から。 -/
theorem isIso_rtExt_one (A : C) : IsIso (rtExt P F A 1) :=
  prop_1_4_iii P F _ (rtExt_frobType P F A 1).1
    ⟨show P.degFr (rtExt P F A 1) = 1 from rtExt_degFr P F A 1,
      (rtExt_frobType P F A 1).2⟩

/-! ## ★16. 関手 `𝒞 → 𝒞^pf`(原文の「by mapping A → (A, 1)」) -/

/-- ★次数 1 の持ち上げは恒等射(`𝒞` の totally epimorphicity による一意性から)。 -/
theorem rtLift_eq_id (A : C) {d : ℕ+} (h : d = 1 * d) :
    rtLift P F A h = 𝟙 (rtObj P F A d) :=
  rtLift_uniq P F _ _ (rtLift_ext P F A h) (Category.comp_id _)

variable {P F} in
/-- ★★**恒等射に沿う根の不変性は恒等写像**(前合成する射が `𝟙` に等しければよい)。 -/
theorem rootIso_hom_of_eq_id {A B : C} (a : A ⟶ A) (ha : IsFrobeniusType P a) (b : B ⟶ B)
    (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (ha1 : a = 𝟙 A) (hb1 : b = 𝟙 B) (z : HomPf P F A B) :
    (rootIso (F := F) a ha b hb hd).hom z = z := by
  refine HomColim.induction (P := fun z => (rootIso (F := F) a ha b hb hd).hom z = z)
    (homFunctorPf P F A B) (fun V x => ?_) z
  show (rootIso (F := F) a ha b hb hd).hom (HomPf.mk V x.down) = HomPf.mk V x.down
  rw [rootIso_hom_mk]
  refine HomPf.mk_congr (congrArg (Comma.mk _ _) ?_) HEq.rfl
  have hid : biFrHom (F := F) a ha b hb hd = 𝟙 (biFrObj P F A B) :=
    WideSubcategory.hom_ext _ (Prod.ext ha1 hb1)
  rw [hid]
  exact Category.id_comp _

variable {P F} in
/-- ★同上、逆向き。 -/
theorem rootIso_inv_of_eq_id {A B : C} (a : A ⟶ A) (ha : IsFrobeniusType P a) (b : B ⟶ B)
    (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (ha1 : a = 𝟙 A) (hb1 : b = 𝟙 B) (z : HomPf P F A B) :
    (rootIso (F := F) a ha b hb hd).inv z = z := by
  conv_lhs => rw [← rootIso_hom_of_eq_id a ha b hb hd ha1 hb1 z]
  exact Iso.hom_inv_id_apply _ _

/-- ★選んだ根の版(次数 1 の持ち上げ)。 -/
theorem rtRootIso_hom_eq_self (A B : C) {dA dB : ℕ+} (hA : dA = 1 * dA) (hB : dB = 1 * dB)
    (z : HomPf P F (rtObj P F A dA) (rtObj P F B dB)) :
    (rtRootIso P F A B hA hB).hom z = z :=
  rootIso_hom_of_eq_id _ _ _ _ _ (rtLift_eq_id P F A hA) (rtLift_eq_id P F B hB) z

theorem rtRootIso_inv_eq_self (A B : C) {dA dB : ℕ+} (hA : dA = 1 * dA) (hB : dB = 1 * dB)
    (z : HomPf P F (rtObj P F A dA) (rtObj P F B dB)) :
    (rtRootIso P F A B hA hB).inv z = z :=
  rootIso_inv_of_eq_id _ _ _ _ _ (rtLift_eq_id P F A hA) (rtLift_eq_id P F B hB) z

variable {P F} in
/-- ★`𝒞` の射を `𝒞^pf` の射(根 1)へ移す。★`A ≅ A^{(1)}` を経由する。 -/
noncomputable def toRootHom {A B : C} (f : A ⟶ B) :
    HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩ :=
  haveI := isIso_rtExt_one P F A
  toHomPf (F := F) (inv (rtExt P F A 1) ≫ f ≫ rtExt P F B 1)

variable {P F} in
theorem toRootHom_id (A : C) : toRootHom (F := F) (𝟙 A) = idRoot P F ⟨A, 1⟩ := by
  haveI := isIso_rtExt_one P F A
  show toHomPf (F := F) (inv (rtExt P F A 1) ≫ 𝟙 A ≫ rtExt P F A 1)
    = toHomPf (F := F) (𝟙 (rtObj P F A 1))
  congr 1
  simp

variable {P F} in
theorem toRootHom_comp {A B E : C} (f : A ⟶ B) (g : B ⟶ E) :
    toRootHom (F := F) (f ≫ g)
      = compRoot P F (toRootHom (F := F) f) (toRootHom (F := F) g) := by
  have hlift := compRoot_eq_lift (toRootHom (F := F) f) (toRootHom (F := F) g)
    (c := 1) (PA := 1) (PB := 1) (PE := 1) (hcA := rfl) (hcB := rfl) (hcE := rfl)
    (ef := 1) (eg := 1) (er := 1) (hfA := rfl) (hfB := rfl) (hgA := rfl) (hgE := rfl)
    (hrA := rfl) (hrE := rfl)
  rw [hlift, rtRootIso_inv_eq_self, rtRootIso_inv_eq_self, rtRootIso_hom_eq_self]
  haveI := isIso_rtExt_one P F A
  haveI := isIso_rtExt_one P F B
  haveI := isIso_rtExt_one P F E
  show toHomPf (F := F) (inv (rtExt P F A 1) ≫ (f ≫ g) ≫ rtExt P F E 1)
    = compPf P F (toHomPf (F := F) (inv (rtExt P F A 1) ≫ f ≫ rtExt P F B 1))
      (toHomPf (F := F) (inv (rtExt P F B 1) ≫ g ≫ rtExt P F E 1))
  rw [← toHomPf_comp]
  congr 1
  simp

/-- ★★★**[FrdI] Definition 3.1, (iii)** —— **自然な関手 `𝒞 → 𝒞^pf`**
(原文の「by mapping `A → (A, 1)`」)。 -/
noncomputable def toPfRoot : C ⥤ PfRootObj P F where
  obj A := ⟨A, 1⟩
  map f := toRootHom (F := F) f
  map_id A := toRootHom_id A
  map_comp f g := toRootHom_comp f g

/-! ## ★★★`Definition 3.1, (ii)(iii)` の主張と実装の対応

★**このファイルが担当するのは (ii) と (iii)**((i)(iv) は `Def31.lean`)。

| 原文の主張 | 実装名 |
|---|---|
| (ii) `𝒞^{Fr-tp} ⊆ 𝒞` | `FrTp` |
| (ii) `𝒞^{bi-Fr} ⊆ 𝒞^{Fr-tp} × 𝒞^{Fr-tp}` | `BiFr` |
| (ii) 添字圏 `(A,B)𝒞^{bi-Fr}` | `IdxPf` |
| (ii) 遷移写像「the assignment φ → φ′」 | `frobTransport`(`Proposition 1.10, (i)`) |
| (ii) `Hom^pf_𝒞(A,B)`(帰納極限) | `HomPf` |
| (ii) perfected morphism の代表元 | `HomPf.mk` / `HomPf.exists_rep` / `HomPf.eq_iff` |
| (iii) 対象 `(A, n)` | `PfRootObj` |
| (iii) `Hom_{𝒞^pf}((A,n),(B,m)) := Hom^pf(A′,B′)` | `HomRoot` |
| (iii) 「independent of the choice of `A → A′`, `B → B′`」 | `rootIso` / `rtRootIso` |
| (iii) 「composition … in the evident fashion」 | `compPf` / `compRoot` |
| (iii) `𝒞^pf` が圏であること | `pfRootCategory`(＋ `PfCat` は `n=1` の部分) |
| (iii) 「a natural functor `𝒞 → 𝒞^pf`(`A → (A,1)`)」 | `toPfRoot`(＋ `toPfCat`) |

★★**原文が言っていないが必要だったもの(測定)**:
1. **添字圏の filtered 性**(`idxPf_isFiltered`)—— 「帰納極限」と呼ぶ以上必要。
2. **3 つ組・4 つ組の添字圏と射影の cofinal 性**(`idx12_final` ほか)——
   合成と結合律の「真ん中を共有する」を正確にするため。
3. **level(共通の根)による計算**(`compRoot_eq_lift`)—— 結合律の急所。
-/

/-! ## ★★★`Definition 3.1` 全体(条なしの `.src`)

★条なしの `.src` は「その原典項目を**完全に**実装した」という主張である。
`Definition 3.1` は **4 条・7 主張**あり、**3 ファイルに分かれて**すべて実装された。

| 条 | # | 内容 | 実装(ファイル) |
|---|---|---|---|
| (i) | 1 | `quasi-isotropic 型` | `Def31.lean` の `IsOfQuasiIsotropicType` |
| (i) | 2 | `standard 型`((a)–(e)) | `Def31.lean` の `IsOfStandardType` |
| (i) | 3 | `Frobenius-slim` な圏 | `Def31.lean` の `IsFrobeniusSlim` |
| (ii) | 4 | `𝒞^{Fr-tp}` / `𝒞^{bi-Fr}` と `Hom^pf_𝒞(A,B)` | **本ファイル** の `FrTp` / `BiFr` / `HomPf` |
| (iii) | 5 | `𝒞^pf` と関手 `𝒞 → 𝒞^pf` | **本ファイル** の `PfRootObj` / `HomRoot` / `pfRootCategory` / `toPfRoot` |
| (iv) | 6 | `unit-equivalent` と `Hom^un-tr` | `Def31.lean` の `IsUnitEquivalent` ＋ `UnTr.lean` の `HomUnTr` |
| (iv) | 7 | `𝒞^un-tr` | `UnTr.lean` の `UnTr` ＋ その `Category` インスタンス |

★★**7 主張すべてを目視で確認した(2026-08-17)**:
`Def31.lean:55/82/103`(1–3)、本ファイル(4–5)、
`Def31.lean:134` ＋ `UnTr.lean:42`(6)、`UnTr.lean:71/73`(7)。
-/

/-- ★★★**[FrdI] Definition 3.1** —— 4 条・7 主張すべて実装済み。 -/
def definition_3_1.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 56, item := "Definition 3.1",
    sectionId := "frdi-def-3-1" }

end ABC3.Found.FrdI
