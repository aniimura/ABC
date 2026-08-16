import ABC3.Found.FrdI.Def23
import Mathlib.CategoryTheory.Skeletal
import Mathlib.CategoryTheory.IsConnected

/-!
# [FrdI] Definition 2.7 —— base-section / Frobenius-section / base-Frobenius pair

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.51–p.52。

原文 (FrdI p.51):
> (i) We shall refer to as a base-section of the Frobenioid C any subcategory

原文 (FrdI p.51):
> following conditions: (a) P is a skeleton

## ★この定義の中身(測定)

**3 条、主張は 8**:

| 条 | # | 内容 |
|---|---|---|
| (i) | 1 | `base-section` = `𝒞^pl-bk` の部分圏で (a) skeleton (b) Frobenius-trivial (c) `𝒫 → 𝒟` が圏同値 |
| (i) | 2 | `P-distinguished` = `𝒫` に入る `𝒞` の射 |
| (ii) | 3 | `ϵ ∈ End(𝒫 → 𝒞)` の Frobenius 次数が**対象の取り方に依らない**(`𝒟` が連結だから) |
| (ii) | 4 | `Frobenius-section` = モノイド準同型 `F : ℕ≥1 → End(𝒫 → 𝒞)` で (a) 次数が恒等 (b) 像は base-identity な Frobenius 型 |
| (ii) | 5 | `quasi-Frobenius-section` = `ℕ≥1` の自己同型との合成を除いた `F` |
| (ii) | 6 | `F-distinguished` = `F` の像が誘導する `𝒞` の自己射 |
| (iii) | 7 | `base-Frobenius pair` / `quasi-base-Frobenius pair` |
| (iii) | 8 | `pre-model type` = base-Frobenius pair を持つこと |

## ★「部分圏」の形式化

原文は `𝒫 ⊆ 𝒞^pl-bk ⊆ 𝒞` と部分圏として述べる。
★**対象の述語 `objP` と射の述語 `homP` の組**として持ち、
そこから圏 `𝒫`(`BaseSection.Obj`)を組み立てる。
★**`𝒫 ⊆ 𝒞^pl-bk` は「`homP` の射はすべて pull-back 射」**として書く。

★**§0 の `skeleton`** は「同型な 2 対象は等しい」——
原文 (FrdI p.15):
> A category C will be called a skeleton if any two isomorphic objects of C are, in
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★★★**[FrdI] Definition 2.7, (i)** —— `𝒞` の **base-section**。

原文 (FrdI p.51):
> trivial; (c) the composite P

★`objP` / `homP` が部分圏 `𝒫 ⊆ 𝒞` を定め、
`isPullBack` が `𝒫 ⊆ 𝒞^pl-bk` を、
`skeletal` / `frobTrivial` / `essSurj`+`full`+`faithful` が (a)/(b)/(c) を与える。 -/
structure BaseSection where
  /-- `𝒫` の対象。 -/
  objP : C → Prop
  /-- `𝒫` の射(＝ **`P-distinguished`** な `𝒞` の射)。 -/
  homP : ∀ {A B : C}, (A ⟶ B) → Prop
  /-- 部分圏: 恒等射を含む。 -/
  id_mem : ∀ {A : C}, objP A → homP (𝟙 A)
  /-- 部分圏: 合成で閉じる。 -/
  comp_mem : ∀ {A B E : C} {f : A ⟶ B} {g : B ⟶ E}, homP f → homP g → homP (f ≫ g)
  /-- **`𝒫 ⊆ 𝒞^pl-bk`** —— `𝒫` の射は pull-back 射。 -/
  isPullBack : ∀ {A B : C} {f : A ⟶ B}, homP f → IsPullBack P f
  /-- **(a)** `𝒫` は §0 の意味の **skeleton** —— `𝒫` の中で同型な 2 対象は等しい。 -/
  skeletal : ∀ {A B : C} {f : A ⟶ B} {g : B ⟶ A}, homP f → homP g →
    f ≫ g = 𝟙 A → g ≫ f = 𝟙 B → A = B
  /-- **(b)** `𝒫` の対象はすべて **Frobenius-trivial**。 -/
  frobTrivial : ∀ {A : C}, objP A → IsFrobeniusTrivial P A
  /-- **(c)** `𝒫 → 𝒟` は**本質的全射**。 -/
  essSurjP : ∀ X : D, ∃ A : C, objP A ∧ Nonempty ((P.toElem.obj A).base ≅ X)
  /-- **(c)** `𝒫 → 𝒟` は**充満**。 -/
  fullP : ∀ {A B : C}, objP A → objP B →
    ∀ ψ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base, ∃ f : A ⟶ B, homP f ∧ P.Base f = ψ
  /-- **(c)** `𝒫 → 𝒟` は**忠実**。 -/
  faithfulP : ∀ {A B : C} {f g : A ⟶ B}, homP f → homP g → P.Base f = P.Base g → f = g

variable {P}

/-- ★**(i)** `P-distinguished` —— `𝒫` に入る `𝒞` の射。 -/
def IsPDistinguished (S : BaseSection P) {A B : C} (f : A ⟶ B) : Prop := S.homP f

/-! ## ★`𝒫` を圏として組み立てる -/

/-- ★`𝒫` の対象型。 -/
def BaseSection.Obj (S : BaseSection P) : Type u2 := {A : C // S.objP A}

instance (S : BaseSection P) : Category.{v2} S.Obj where
  Hom A B := {f : A.1 ⟶ B.1 // S.homP f}
  id A := ⟨𝟙 A.1, S.id_mem A.2⟩
  comp f g := ⟨f.1 ≫ g.1, S.comp_mem f.2 g.2⟩
  id_comp _ := Subtype.ext (Category.id_comp _)
  comp_id _ := Subtype.ext (Category.comp_id _)
  assoc _ _ _ := Subtype.ext (Category.assoc _ _ _)

/-- ★**(c) の関手** `𝒫 → 𝒟`。 -/
def BaseSection.toD (S : BaseSection P) : S.Obj ⥤ D where
  obj A := (P.toElem.obj A.1).base
  map f := P.Base f.1
  map_id A := P.Base_id A.1
  map_comp f g := P.Base_comp f.1 g.1

/-- ★★**(i)(c)** `𝒫 → 𝒟` は**圏同値**。 -/
theorem BaseSection.toD_isEquivalence (S : BaseSection P) : S.toD.IsEquivalence := by
  haveI : S.toD.Faithful := ⟨fun {_ _ f g} h => Subtype.ext (S.faithfulP f.2 g.2 h)⟩
  haveI : S.toD.Full := ⟨fun {A B} ψ => by
    obtain ⟨f, hf, hb⟩ := S.fullP A.2 B.2 ψ
    exact ⟨⟨f, hf⟩, hb⟩⟩
  haveI : S.toD.EssSurj := ⟨fun X => by
    obtain ⟨A, hA, ⟨e⟩⟩ := S.essSurjP X
    exact ⟨(⟨A, hA⟩ : S.Obj), ⟨e⟩⟩⟩
  exact ⟨inferInstance, inferInstance, ‹_›⟩

/-- ★**`𝒟` が連結なら `𝒫` も連結** —— 原文の
「since D, hence also P, is a connected category」。 -/
theorem BaseSection.isConnected_obj (S : BaseSection P) [IsConnected D] :
    IsConnected S.Obj := by
  haveI := S.toD_isEquivalence
  exact isConnected_of_equivalent S.toD.asEquivalence.symm

/-! ## ★(ii) —— `End(𝒫 → 𝒞)` と、その Frobenius 次数

原文 (FrdI p.51):
> sense to speak of the Frobenius degree degFr(
-/

/-- ★**`End(𝒫 → 𝒞)`** —— 包含関手の自己自然変換。 -/
@[ext]
structure SectionEnd (S : BaseSection P) where
  /-- 各対象での成分。 -/
  app : ∀ A : S.Obj, End A.1
  /-- 自然性。 -/
  naturality : ∀ {A B : S.Obj} (f : A ⟶ B), f.1 ≫ app B = app A ≫ f.1

namespace SectionEnd

variable {S : BaseSection P}

instance : Monoid (SectionEnd S) where
  one := ⟨fun _ => 𝟙 _, fun f => by simp⟩
  mul ε ε' :=
    ⟨fun A => ε.app A * ε'.app A, fun {A B} f => by
      show f.1 ≫ (ε'.app B ≫ ε.app B) = (ε'.app A ≫ ε.app A) ≫ f.1
      rw [← Category.assoc, ε'.naturality f, Category.assoc, ε.naturality f, Category.assoc]⟩
  one_mul ε := by ext A; exact Category.comp_id _
  mul_one ε := by ext A; exact Category.id_comp _
  mul_assoc ε₁ ε₂ ε₃ := by ext A; exact (Category.assoc _ _ _).symm

@[simp] theorem one_app (A : S.Obj) : (1 : SectionEnd S).app A = 𝟙 A.1 := rfl

@[simp] theorem mul_app (ε ε' : SectionEnd S) (A : S.Obj) :
    (ε * ε').app A = ε.app A * ε'.app A := rfl

/-- ★★★**(ii) の第 1 主張** —— `ϵ` の Frobenius 次数は**対象の取り方に依らない**。

★**理由は `𝒫` の連結性**(原文どおり)。自然性 `f ≫ ϵ_B = ϵ_A ≫ f` に
`degFr` を当てると `degFr ϵ_B * degFr f = degFr f * degFr ϵ_A` となり、
`ℕ≥1` は可換で簡約可能だから `degFr ϵ_A = degFr ϵ_B`。 -/
theorem degFr_const [IsConnected D] (ε : SectionEnd S) (A B : S.Obj) :
    P.degFr (ε.app A) = P.degFr (ε.app B) := by
  haveI := S.isConnected_obj
  refine constant_of_preserves_morphisms (fun X : S.Obj => P.degFr (ε.app X)) ?_ A B
  intro X Y f
  have h := congrArg P.degFr (ε.naturality f)
  rw [P.degFr_comp, P.degFr_comp, mul_comm (P.degFr f.1)] at h
  exact (mul_right_cancel h).symm

/-- ★**`ϵ` の Frobenius 次数**(対象の取り方に依らない、`degFr_const`)。 -/
noncomputable def deg [IsConnected D] (ε : SectionEnd S) : ℕ+ :=
  haveI := S.isConnected_obj
  P.degFr (ε.app (Classical.arbitrary S.Obj))

theorem deg_eq [IsConnected D] (ε : SectionEnd S) (A : S.Obj) :
    ε.deg = P.degFr (ε.app A) :=
  haveI := S.isConnected_obj
  degFr_const ε _ A

/-- ★**次数を取る写像はモノイド準同型** `End(𝒫 → 𝒞) → ℕ≥1`。 -/
noncomputable def degHom [IsConnected D] : SectionEnd S →* ℕ+ where
  toFun := deg
  map_one' :=
    haveI := S.isConnected_obj
    (deg_eq 1 (Classical.arbitrary S.Obj)).trans (P.degFr_id _)
  map_mul' ε ε' := by
    haveI := S.isConnected_obj
    let A := Classical.arbitrary S.Obj
    rw [deg_eq (ε * ε') A, deg_eq ε A, deg_eq ε' A, mul_app]
    show P.degFr (ε'.app A ≫ ε.app A) = _
    rw [P.degFr_comp]

end SectionEnd

/-! ## ★(ii) —— Frobenius-section

原文 (FrdI p.51):
> satisfying the following conditions: (a) the composite of F with the homomorphism
-/

/-- ★★★**[FrdI] Definition 2.7, (ii)** —— `𝒫`-**Frobenius-section**。

★(a) 次数を取ると `ℕ≥1` 上の恒等、★(b) 像は base-identity な Frobenius 型射。 -/
structure IsFrobeniusSection [IsConnected D] (S : BaseSection P)
    (Fs : ℕ+ →* SectionEnd S) : Prop where
  /-- **(a)** 次数との合成が恒等。 -/
  degSection : ∀ n : ℕ+, (Fs n).deg = n
  /-- **(b)** 像は base-identity。 -/
  baseIdentity : ∀ (n : ℕ+) (A : S.Obj), IsBaseIdentity P ((Fs n).app A)
  /-- **(b)** 像は Frobenius 型。 -/
  frobType : ∀ (n : ℕ+) (A : S.Obj), IsFrobeniusType P ((Fs n).app A)

/-- ★**(ii)** `F-distinguished` —— `F` の像が誘導する `𝒞` の自己射。 -/
def IsFDistinguished [IsConnected D] {S : BaseSection P} (Fs : ℕ+ →* SectionEnd S)
    {A : C} (hA : S.objP A) (f : End A) : Prop :=
  ∃ n : ℕ+, (Fs n).app ⟨A, hA⟩ = f

/-- ★**(ii)** `quasi-Frobenius-section` —— `ℕ≥1` の自己同型との合成を除いて等しいこと。

★**同値関係**であり(`quasiFrobRel_equivalence`)、その商が
原文の「regarded as being known only up to composition with automorphisms」である。 -/
def QuasiFrobRel [IsConnected D] {S : BaseSection P} (F₁ F₂ : ℕ+ →* SectionEnd S) : Prop :=
  ∃ σ : ℕ+ ≃* ℕ+, F₂ = F₁.comp (σ : ℕ+ →* ℕ+)

theorem quasiFrobRel_equivalence [IsConnected D] (S : BaseSection P) :
    Equivalence (QuasiFrobRel (S := S)) where
  refl F := ⟨MulEquiv.refl _, by ext n; rfl⟩
  symm := fun {F₁ F₂} ⟨σ, h⟩ => ⟨σ.symm, by
    subst h; ext n; simp⟩
  trans := fun {F₁ F₂ F₃} ⟨σ, h⟩ ⟨σ', h'⟩ => ⟨σ'.trans σ, by
    subst h; subst h'; ext n; rfl⟩

/-! ## ★(iii) —— base-Frobenius pair と pre-model type

原文 (FrdI p.52):
> (iii) We shall refer to a pair (P, F), where P is a base-section of C, and F is
-/

variable (P) in
/-- ★★★**[FrdI] Definition 2.7, (iii)** —— **base-Frobenius pair**。 -/
structure BaseFrobeniusPair [IsConnected D] where
  /-- base-section `𝒫`。 -/
  sec : BaseSection P
  /-- `𝒫`-Frobenius-section `F`。 -/
  frob : ℕ+ →* SectionEnd sec
  /-- `F` が Frobenius-section であること。 -/
  isFrobSection : IsFrobeniusSection sec frob

variable (P) in
/-- ★★**[FrdI] Definition 2.7, (iii)** —— `𝒞` が **pre-model type** である。

★原文が「or, equivalently, a quasi-base-Frobenius pair」と言うとおり、
quasi 版の存在と同値である(`isPreModelType_iff_quasi`)。 -/
def IsPreModelType [IsConnected D] : Prop := Nonempty (BaseFrobeniusPair P)

/-- ★**(iii)** `quasi-base-Frobenius pair` —— `F` を `ℕ≥1` の自己同型との合成を
除いてしか知らないものと見た組。 -/
def QuasiBaseFrobRel [IsConnected D] (p q : BaseFrobeniusPair P) : Prop :=
  ∃ h : p.sec = q.sec, QuasiFrobRel (h ▸ p.frob) q.frob

/-- ★★**「base-Frobenius pair を持つ ⟺ quasi-base-Frobenius pair を持つ」** ——
原文の「or, equivalently」。

★quasi 版は base-Frobenius pair の**同値類**なので、存在は同値である。 -/
theorem isPreModelType_iff_quasi [IsConnected D] :
    IsPreModelType P ↔ Nonempty (Quot (QuasiBaseFrobRel (P := P))) :=
  ⟨fun ⟨p⟩ => ⟨Quot.mk _ p⟩, fun ⟨c⟩ => Quot.recOnSubsingleton c fun p => ⟨p⟩⟩

/-! ## ★Remark 2.7.2 の第 1 主張

原文 (FrdI p.52):
> pair of C. Then the only arrows of C which are both F- and P-distinguished [hence

★**F- かつ P-distinguished な射は恒等射に限る。**
-/

/-- ★★★**[FrdI] Remark 2.7.2 の第 1 主張** ——
`F-distinguished` かつ `P-distinguished` な射は**恒等射**である。

★**理由は次数の消去**: `P-distinguished` なら pull-back 射ゆえ linear(次数 1)、
`F-distinguished` なら次数は `n`。よって `n = 1` で、`F` はモノイド準同型だから
`F(1) = 1`、すなわち恒等射。

★**base-trivial 型も `𝒞` が skeleton であることも要らない**
(原文はそれらを「Suppose further」として第 2 主張のためにだけ課している)。 -/
theorem fdist_and_pdist_eq_id [IsConnected D] (Fc : FrobenioidCore P) {S : BaseSection P}
    {Fs : ℕ+ →* SectionEnd S} (hFs : IsFrobeniusSection S Fs)
    {A : C} (hA : S.objP A) (f : End A)
    (hF : IsFDistinguished Fs hA f) (hP : S.homP f) : f = 𝟙 A := by
  obtain ⟨n, hn⟩ := hF
  have hlin : P.degFr f = 1 := (Fc.pullBackLB f (S.isPullBack hP)).2
  have hdeg : P.degFr f = n := by
    rw [← hn, ← SectionEnd.deg_eq (Fs n) ⟨A, hA⟩, hFs.degSection n]
  have hn1 : n = 1 := by rw [← hdeg, hlin]
  rw [← hn, hn1, map_one]
  rfl

/-! ## ★Remark 2.7.2 の第 2・第 3 主張 —— 3 分解とその一意性

原文 (FrdI p.52):
> is a base-identity pre-step endomorphism;

★**`𝒞` が base-trivial 型かつ skeleton** のとき、どの射も
`φ = γ ≫ β ≫ α`(原文の `α ∘ β ∘ γ`)と**一意に**分解する。
-/

section Rem272

variable [IsConnected D] {S : BaseSection P}

/-- ★**base-trivial 型かつ `𝒞` が skeleton なら、`𝒞` のどの対象も `𝒫` の対象**。

★`𝒫 → 𝒟` の本質的全射性で底の同型な `𝒫`-対象を取り、
base-trivial 型でそれが `𝒞` で同型になり、skeleton で**等しく**なる。 -/
theorem objP_of_baseTrivial_of_skeletal (S : BaseSection P)
    (hbt : ∀ A : C, IsBaseTrivial P A) (hskel : Skeletal C) (A : C) : S.objP A := by
  obtain ⟨A', hA', ⟨e⟩⟩ := S.essSurjP ((P.toElem.obj A).base)
  obtain ⟨g⟩ := hbt A A' ⟨e.symm⟩
  exact (hskel ⟨g⟩ : A' = A) ▸ hA'

/-- ★★★**[FrdI] Remark 2.7.2 の第 2 主張** —— 3 分解の**存在**。

★**手順**:
1. `𝒫 → 𝒟` の充満性で `Base α = Base φ` なる `P-distinguished` な `α` を取る
2. `α` は pull-back 射なので、`f ≫ α = φ` かつ `Base f = 𝟙` なる `f` が**一意に**取れる
3. `Definition 1.3, (iv), (a)` で `f = γ₁ ≫ β₁ ≫ α₁` と分け、
   base-trivial ＋ skeleton で中間対象が `A` に**潰れる**
4. `Base f = 𝟙` から `Base α₁` が同型、よって `α₁` 自身が同型
   (`Remark 1.2.1`、`isIso_of_isPullBack_of_isBaseIso`)
5. `Definition 1.3, (ii)` の一意性(`frobDegUniq`)で `γ₁` を
   `F` の像 `γ = F(n)_A` に取り替える —— ずれた同型は `β` に吸収する -/
theorem rem_2_7_2_factor (Fc : FrobenioidCore P) {Fs : ℕ+ →* SectionEnd S}
    (hFs : IsFrobeniusSection S Fs)
    (hbt : ∀ A : C, IsBaseTrivial P A) (hskel : Skeletal C)
    {A B : C} (hA : S.objP A) (hB : S.objP B) (φ : A ⟶ B) :
    ∃ (n : ℕ+) (β : A ⟶ A) (α : A ⟶ B),
      φ = ((Fs n).app ⟨A, hA⟩ : A ⟶ A) ≫ β ≫ α ∧
      IsBaseIdentity P β ∧ IsPreStep P β ∧ S.homP α := by
  -- 段 1
  obtain ⟨α, hαP, hαb⟩ := S.fullP hA hB (P.Base φ)
  have hαpb : IsPullBack P α := S.isPullBack hαP
  -- 段 2
  obtain ⟨-, hsurj⟩ := hαpb A
  obtain ⟨f, hf⟩ := hsurj ⟨(φ, 𝟙 _), by rw [hαb, Category.id_comp]⟩
  have hf' := congrArg Subtype.val hf
  have hfα : f ≫ α = φ := congrArg Prod.fst hf'
  have hfb : P.Base f = 𝟙 _ := congrArg Prod.snd hf'
  -- 段 3
  obtain ⟨X, Y, γ₁, β₁, α₁, hfac, hγ₁, hβ₁, hα₁⟩ := Fc.arbFactor f
  haveI : IsIso (P.Base γ₁) := hγ₁.2
  obtain ⟨gX⟩ := hbt A X ⟨asIso (P.Base γ₁)⟩
  have hXA : A = X := (hskel ⟨gX⟩).symm
  subst hXA
  haveI : IsIso (P.Base β₁) := hβ₁.2
  obtain ⟨gY⟩ := hbt A Y ⟨asIso (P.Base β₁)⟩
  have hYA : A = Y := (hskel ⟨gY⟩).symm
  subst hYA
  -- 段 4
  have hbα₁ : IsIso (P.Base α₁) := by
    have h : P.Base γ₁ ≫ P.Base β₁ ≫ P.Base α₁ = 𝟙 _ := by
      rw [← P.Base_comp, ← P.Base_comp, ← hfac, hfb]
    haveI : IsIso (P.Base γ₁ ≫ P.Base β₁ ≫ P.Base α₁) := by rw [h]; infer_instance
    haveI := IsIso.of_isIso_comp_left (P.Base γ₁) (P.Base β₁ ≫ P.Base α₁)
    exact IsIso.of_isIso_comp_left (P.Base β₁) (P.Base α₁)
  haveI : IsIso α₁ := isIso_of_isPullBack_of_isBaseIso P Fc α₁ hα₁ hbα₁
  have hβ₂ : IsPreStep P (β₁ ≫ α₁) := IsPreStep.comp P hβ₁ (isPreStep_of_isIso P α₁)
  -- 段 5
  set n : ℕ+ := P.degFr γ₁ with hn
  have hγdeg : P.degFr ((Fs n).app ⟨A, hA⟩) = n :=
    (SectionEnd.deg_eq (Fs n) ⟨A, hA⟩).symm.trans (hFs.degSection n)
  obtain ⟨e, hei, hee⟩ := Fc.frobDegUniq A A A γ₁ ((Fs n).app ⟨A, hA⟩) hγ₁
    (hFs.frobType n ⟨A, hA⟩) (by rw [hγdeg])
  haveI := hei
  have hγ₁e : γ₁ = ((Fs n).app ⟨A, hA⟩ : A ⟶ A) ≫ inv e := by rw [← hee]; simp
  refine ⟨n, inv e ≫ β₁ ≫ α₁, α, ?_, ?_, ?_, hαP⟩
  · rw [← hfα, hfac, hγ₁e]
    simp
  · -- base-identity
    show P.Base (inv e ≫ β₁ ≫ α₁) = P.Base (𝟙 A)
    have hγb : P.Base ((Fs n).app ⟨A, hA⟩) = P.Base (𝟙 A) := hFs.baseIdentity n ⟨A, hA⟩
    have h : P.Base ((Fs n).app ⟨A, hA⟩) ≫ P.Base (inv e ≫ β₁ ≫ α₁) = P.Base (𝟙 A) := by
      rw [← P.Base_comp, ← Category.assoc, ← hγ₁e, ← hfac, hfb, P.Base_id]
    rw [hγb, P.Base_id, Category.id_comp] at h
    rw [P.Base_id]
    exact h
  · exact IsPreStep.comp P (isPreStep_of_isIso P (inv e)) hβ₂

/-- ★★★**[FrdI] Remark 2.7.2 の第 3 主張** —— 3 分解の**一意性**(strict)。

★**手順**: `γ`, `β` が base-identity なので `Base α = Base φ = Base α'`、
`𝒫 → 𝒟` の忠実性で `α = α'`。次に `α` が pull-back 射だから
`γ ≫ β = γ' ≫ β'`。次数を取ると `β` が linear なので `n = m`。
最後は `𝒞` の totally epimorphicity(`γ` が epi)で `β = β'`。 -/
theorem rem_2_7_2_uniq {Fs : ℕ+ →* SectionEnd S} (hFs : IsFrobeniusSection S Fs)
    {A B : C} (hA : S.objP A) {n m : ℕ+} {β β' : A ⟶ A} {α α' : A ⟶ B}
    (hβ : IsBaseIdentity P β) (hβs : IsPreStep P β) (hαP : S.homP α)
    (hβ' : IsBaseIdentity P β') (hβ's : IsPreStep P β') (hα'P : S.homP α')
    (heq : ((Fs n).app ⟨A, hA⟩ : A ⟶ A) ≫ β ≫ α
        = ((Fs m).app ⟨A, hA⟩ : A ⟶ A) ≫ β' ≫ α') :
    n = m ∧ β = β' ∧ α = α' := by
  have hγb : ∀ k : ℕ+, P.Base ((Fs k).app ⟨A, hA⟩) = 𝟙 _ := by
    intro k
    have h : P.Base ((Fs k).app ⟨A, hA⟩) = P.Base (𝟙 A) := hFs.baseIdentity k ⟨A, hA⟩
    rwa [P.Base_id] at h
  have hβb : P.Base β = 𝟙 _ := by have h : P.Base β = P.Base (𝟙 A) := hβ; rwa [P.Base_id] at h
  have hβ'b : P.Base β' = 𝟙 _ := by have h : P.Base β' = P.Base (𝟙 A) := hβ'; rwa [P.Base_id] at h
  -- `α = α'`
  have hbase : P.Base α = P.Base α' := by
    have h := congrArg P.Base heq
    rw [P.Base_comp, P.Base_comp, P.Base_comp, P.Base_comp, hγb n, hγb m, hβb, hβ'b] at h
    simpa using h
  have hαα : α = α' := S.faithfulP hαP hα'P hbase
  subst hαα
  -- `γ ≫ β = γ' ≫ β'`
  have hcomp : ((Fs n).app ⟨A, hA⟩ : A ⟶ A) ≫ β = ((Fs m).app ⟨A, hA⟩ : A ⟶ A) ≫ β' := by
    obtain ⟨hinj, -⟩ := S.isPullBack hαP A
    refine hinj (Subtype.ext (Prod.ext ?_ ?_))
    · show (((Fs n).app ⟨A, hA⟩ : A ⟶ A) ≫ β) ≫ α = (((Fs m).app ⟨A, hA⟩ : A ⟶ A) ≫ β') ≫ α
      rw [Category.assoc, Category.assoc]; exact heq
    · show P.Base (((Fs n).app ⟨A, hA⟩ : A ⟶ A) ≫ β)
        = P.Base (((Fs m).app ⟨A, hA⟩ : A ⟶ A) ≫ β')
      rw [P.Base_comp, P.Base_comp, hγb n, hγb m, hβb, hβ'b]
  -- `n = m`
  have hnm : n = m := by
    have h := congrArg P.degFr hcomp
    rw [P.degFr_comp, P.degFr_comp,
      show P.degFr β = 1 from hβs.1, show P.degFr β' = 1 from hβ's.1,
      ← SectionEnd.deg_eq (Fs n) ⟨A, hA⟩, ← SectionEnd.deg_eq (Fs m) ⟨A, hA⟩,
      hFs.degSection n, hFs.degSection m] at h
    simpa using h
  subst hnm
  -- `β = β'`
  haveI : Epi ((Fs n).app ⟨A, hA⟩ : A ⟶ A) := P.totEpiC _ _ _
  exact ⟨rfl, (cancel_epi ((Fs n).app ⟨A, hA⟩ : A ⟶ A)).mp hcomp, rfl⟩

end Rem272

/-! ## ★★★出典の紐付け(`.src`) -/

variable (P) in
/-- ★★★**[FrdI] Definition 2.7** —— 3 条・8 主張すべてが実装された。

| # | 主張 | 実装 |
|---|---|---|
| 1 | (i) base-section | `BaseSection` |
| 2 | (i) `P-distinguished` | `IsPDistinguished` |
| 3 | (ii) 次数が対象に依らない | `SectionEnd.degFr_const` / `SectionEnd.deg` |
| 4 | (ii) Frobenius-section | `IsFrobeniusSection` |
| 5 | (ii) quasi-Frobenius-section | `QuasiFrobRel` / `quasiFrobRel_equivalence` |
| 6 | (ii) `F-distinguished` | `IsFDistinguished` |
| 7 | (iii) base-Frobenius pair / quasi 版 | `BaseFrobeniusPair` / `QuasiBaseFrobRel` |
| 8 | (iii) pre-model type(＋ quasi との同値) | `IsPreModelType` / `isPreModelType_iff_quasi` |

★**補助として `𝒫 → 𝒟` の圏同値**(`BaseSection.toD_isEquivalence`)と
★**`𝒫` の連結性**(`BaseSection.isConnected_obj`)を出した——
後者が (ii) の「次数が対象に依らない」の根拠である。 -/
def BaseSection.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 51, item := "Definition 2.7",
    sectionId := "frdi-def-2-7" }

/-- ★★★**[FrdI] Remark 2.7.2** —— 3 主張すべてが実装された。

| # | 主張 | 実装 |
|---|---|---|
| 1 | F- かつ P-distinguished な射は恒等射 | `fdist_and_pdist_eq_id` |
| 2 | base-trivial ＋ skeleton なら 3 分解が**存在** | `rem_2_7_2_factor` |
| 3 | その分解は(strict に)**一意** | `rem_2_7_2_uniq` |

★補助として「どの対象も `𝒫` の対象」(`objP_of_baseTrivial_of_skeletal`)を出した。 -/
def rem_2_7_2_factor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 52, item := "Remark 2.7.2",
    sectionId := "frdi-remark-2-7-2" }

end ABC3.Found.FrdI
