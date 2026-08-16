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

end ABC3.Found.FrdI
