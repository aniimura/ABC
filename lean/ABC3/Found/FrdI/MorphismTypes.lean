import ABC3.Found.FrdI.ElementaryFrobenioid
import Mathlib.CategoryTheory.Endomorphism
import Mathlib.Data.Nat.Prime.Basic
import ABC3.Meta.Claim

/-!
# [FrdI] Definition 1.2 —— pre-Frobenioid の射と対象の類型

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.21–p.23(**400 dpi 目視確認 2026-08-15**)。

原文 (FrdI p.21):
> Definition 1.2. Let Φ be a divisorial monoid on a connected, totally epimorphic

原文 (FrdI p.21):
> category D; C →FΦ a pre-Frobenioid; φ ∈Arr(C). Write φ : A →B [where

## ★この定義の規模と、ここで写した範囲

`Definition 1.2` は (i)–(v) からなり、**40 語を超える**類型を一度に定める。
このファイルが写したのは:

* **(i) 全部** —— `linear` / `isometric` / `metrically equivalent`
* **(ii) 全部** —— `base-isomorphism` / `base-FSM-morphism` / `base-isomorphic` /
  `pull-back morphism` / `base-equivalent` / `Div-equivalent` /
  `base-identity` / `Div-identity` / `𝒪^×(A)` / `𝒪^▷(A)`
* **(iii) 全部** —— `pre-step` / `step` / `primary pre-step` / `co-angular` /
  `LB-invertible` / `morphism of Frobenius type` / `prime-Frobenius morphism`
* **(iv) のうち `Definition 1.3` が引くもの** —— `Frobenius-trivial` / `isotropic` /
  `isotropic hull` / `unit-trivial` / `group-like` / `base-trivial` /
  `metrically trivial` / `Div-Frobenius-trivial` /
  `universally Div-Frobenius-trivial` / `Frobenius-ample` / `quasi-Frobenius-trivial` /
  `sub-quasi-Frobenius-trivial`

★**写していないもの(明示する)**: (iv) の
`Aut-ample` / `Aut^sub-ample` / `End-ample` / `perfect` /
`Frobenius-compact` / `Frobenius-normalized` / `Frobenius-isotropic`、
および (v)(pre-Frobenioid に「of X type」を付ける言い換え)。
`Frobenius-compact` は `𝒪^×(C)^pf`(**perfection**)を要求し、
`perfect object` は「base-isomorphic な対象すべて」に渡る条件を含む。
どちらも `Definition 1.3` は引かない。

## ★mathlib の実測(2026-08-15)

| 必要なもの | mathlib | 判定 |
|---|---|---|
| `End A` のモノイド構造 | `CategoryTheory.End`(`x * y = y ≫ x`) | ★**使う** |
| 部分モノイド | `Submonoid` | ★**使う** |
| `𝔓rimes` | `Nat.Prime` | ★**使う** |
| `Aut A` | `CategoryTheory.Aut` | ★**使う** |
| `isotropic hull` に相当する普遍性 | 無い | 自前 |
| `co-angular` | 無い | 自前 |

★`co-angular` は「任意の3分解について β が同型」という**全称条件**であり、
mathlib の直交因子系(`MorphismProperty.factorization`)とは形が違う
(あちらは分解の**存在と一意性**、こちらは分解が存在したときの**帰結**)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w

variable {D : Type u} [Category.{v} D] {C : Type u} [Category.{v} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ### (i) —— `linear` / `isometric` / `metrically equivalent`

原文 (FrdI p.21):
> (i) We shall say that φ is linear if degFr(φ) = 1. We shall say that φ is isometric,

原文 (FrdI p.21):
> or, alternatively, an isometry, if Div(φ) = 0 [cf. Definition 1.1, (iii)]. If ψ ∈Arr(C)

原文 (FrdI p.21):
> is co-objective with φ [cf. §0], then we shall say that φ, ψ are metrically equivalent

原文 (FrdI p.21):
> if Div(φ) = Div(ψ).
-/

/-- **(i)** `linear` —— Frobenius 次数が 1。 -/
def IsLinear {A B : C} (φ : A ⟶ B) : Prop := P.degFr φ = 1

/-- **(i)** `isometric`(= `isometry`)—— 零因子が 0。 -/
def IsIsometric {A B : C} (φ : A ⟶ B) : Prop := P.Div φ = 0

/-- **(i)** `metrically equivalent`。

★原文は `co-objective`(始域と終域が一致)を仮定するが、Lean では
`φ ψ : A ⟶ B` と型で書いた時点でそれが満たされている。 -/
def MetricallyEquivalent {A B : C} (φ ψ : A ⟶ B) : Prop := P.Div φ = P.Div ψ

/-! ### (ii) —— base 系の類型

原文 (FrdI p.21):
> (ii) We shall refer to φ as a base-isomorphism (respectively, base-FSM-morphism)

原文 (FrdI p.21):
> if Base(φ) is an isomorphism (respectively, FSM-morphism [cf. §0]) in D. We shall

原文 (FrdI p.21):
> refer to two objects of C that map to isomorphic objects of D as base-isomorphic.

原文 (FrdI p.21):
> We shall refer to φ as a pull-back morphism if the natural transformation of con-

原文 (FrdI p.21):
> travariant functors on C

原文 (FrdI p.21):
> projection functor C →D] induced by φ is an isomorphism. If ψ ∈Arr(C) is co-

原文 (FrdI p.21):
> objective with φ [cf. §0], then we shall say that φ, ψ are base-equivalent (respectively,

原文 (FrdI p.21):
> Div-equivalent) if Base(φ) = Base(ψ) (respectively, Φ(φ) = Φ(ψ)). If A = B [i.e.,

原文 (FrdI p.22):
> Div-identity) endomorphism if it is base-equivalent (respectively, Div-equivalent)

原文 (FrdI p.22):
> to the identity endomorphism of A. Write

原文 (FrdI p.22):
> for the submonoids of base-identity linear endomorphisms.
-/

/-- **(ii)** `base-isomorphism`。 -/
def IsBaseIsomorphism {A B : C} (φ : A ⟶ B) : Prop := IsIso (P.Base φ)

/-- **(ii)** `base-FSM-morphism`。 -/
def IsBaseFSM {A B : C} (φ : A ⟶ B) : Prop := IsFSMMorphism (P.Base φ)

/-- **(ii)** `base-isomorphic`(対象について)。 -/
def BaseIsomorphic (A B : C) : Prop :=
  Nonempty ((P.toElem.obj A).base ≅ (P.toElem.obj B).base)

/-- **(ii)** `base-equivalent`。 -/
def BaseEquivalent {A B : C} (φ ψ : A ⟶ B) : Prop := P.Base φ = P.Base ψ

/-- **(ii)** `Div-equivalent` —— 原文は `Φ(φ) = Φ(ψ)`、すなわち
**誘導される単項式射 `Φ(B_𝒟) → Φ(A_𝒟)` が一致すること**である。

★`metrically equivalent`(`Div(φ) = Div(ψ)`)とは**別の条件**であることに注意。
原文は (i) と (ii) で別々に定めている。 -/
def DivEquivalent {A B : C} (φ ψ : A ⟶ B) : Prop :=
  Φ.map (P.Base φ) = Φ.map (P.Base ψ)

/-- **(ii)** `base-identity` 自己射。 -/
def IsBaseIdentity {A : C} (φ : A ⟶ A) : Prop := BaseEquivalent P φ (𝟙 A)

/-- **(ii)** `Div-identity` 自己射。 -/
def IsDivIdentity {A : C} (φ : A ⟶ A) : Prop := DivEquivalent P φ (𝟙 A)

/-- **(ii)** `𝒪^▷(A)` —— base-identity かつ linear な自己射のなす部分モノイド。

★原文は「submonoids」と言うだけで、部分モノイドであることを示していない。
**それを示すのがこの定義の中身**である(単位元と合成で閉じること)。 -/
def OTri (A : C) : Submonoid (End A) where
  carrier := {φ : End A | IsBaseIdentity P φ ∧ IsLinear P φ}
  one_mem' := by
    constructor
    · show P.Base (𝟙 A) = P.Base (𝟙 A)
      rfl
    · show P.degFr (𝟙 A) = 1
      simp
  mul_mem' := by
    rintro x y ⟨hx1, hx2⟩ ⟨hy1, hy2⟩
    constructor
    · show P.Base (y ≫ x) = P.Base (𝟙 A)
      rw [P.Base_comp]
      show P.Base y ≫ P.Base x = _
      rw [show P.Base y = P.Base (𝟙 A) from hy1, show P.Base x = P.Base (𝟙 A) from hx1]
      simp
    · show P.degFr (y ≫ x) = 1
      rw [P.degFr_comp]
      rw [show P.degFr x = 1 from hx2, show P.degFr y = 1 from hy2, mul_one]

/-- **(ii)** `𝒪^×(A)` —— `𝒪^▷(A)` のうち可逆なもの(原文は `Aut_𝒞(A)` の部分モノイド)。 -/
def OTimes (A : C) : Submonoid (End A) where
  carrier := {φ : End A | (IsBaseIdentity P φ ∧ IsLinear P φ) ∧ IsUnit φ}
  one_mem' := ⟨(OTri P A).one_mem, isUnit_one⟩
  mul_mem' hx hy := ⟨(OTri P A).mul_mem hx.1 hy.1, hx.2.mul hy.2⟩

theorem OTimes_le_OTri (A : C) : OTimes P A ≤ OTri P A := fun _ hx => hx.1

/-- **(ii)** `pull-back morphism`。

原文 (FrdI p.21):
> Hom C(−, A) →Hom C(−, B) ×Hom D(−,B D)|C (Hom D(−, A D)|C)

★原文のファイバー積を展開すると、各 `X ∈ Ob(𝒞)` について

`Hom_𝒞(X, A) → { (g, h) | g : X ⟶ B, h : X_𝒟 ⟶ A_𝒟, Base(g) = h ≫ Base(φ) }`、
`f ↦ (f ≫ φ, Base f)`

が**全単射**である、という条件になる。★`Base(f ≫ φ) = Base f ≫ Base φ`
(Remark 1.1.1)なので、この写像は確かにファイバー積へ落ちる。 -/
def IsPullBack {A B : C} (φ : A ⟶ B) : Prop :=
  ∀ X : C, Function.Bijective
    (fun f : X ⟶ A =>
      (⟨(f ≫ φ, P.Base f), by rw [P.Base_comp]⟩ :
        {p : (X ⟶ B) × ((P.toElem.obj X).base ⟶ (P.toElem.obj A).base) //
          P.Base p.1 = p.2 ≫ P.Base φ}))

/-! ### (iii) —— step / co-angular / Frobenius type

原文 (FrdI p.22):
> (iii) We shall say that φ is a pre-step [a term motivated by the point of view

原文 (FrdI p.22):
> constituted by a non-zero zero divisor] if it is a linear base-isomorphism. If φ is

原文 (FrdI p.22):
> a pre-step, then we shall say that it is a step (respectively, a primary pre-step) if

原文 (FrdI p.22):
> φ is not an isomorphism (respectively, if the zero divisor Div(φ) ∈Φ(A) of φ is a

原文 (FrdI p.22):
> primary [cf. §0] element of the monoid Φ(A)). We shall say that φ is co-angular [a

原文 (FrdI p.22):
> [Mzk15], Definition 3.1, (iii)] if, for any factorization φ = α ◦β ◦γ in C, where α is

原文 (FrdI p.22):
> linear, β is an isometric pre-step, and either α or γ is a base-isomorphism, it follows

原文 (FrdI p.22):
> that β is an isomorphism. We shall say that φ is LB-invertible [i.e., “line bundle-

原文 (FrdI p.22):
> if it is co-angular and isometric. We shall say that φ is a morphism of Frobenius

原文 (FrdI p.22):
> tensor power” for some n ∈N≥1] if φ is an LB-invertible base-isomorphism. We shall

原文 (FrdI p.22):
> say that φ is a prime-Frobenius morphism, or, alternatively, a degFr(φ)-Frobenius
-/

/-- **(iii)** `pre-step` —— linear な base-isomorphism。 -/
def IsPreStep {A B : C} (φ : A ⟶ B) : Prop := IsLinear P φ ∧ IsBaseIsomorphism P φ

/-- **(iii)** `step` —— 同型でない pre-step。 -/
def IsStep {A B : C} (φ : A ⟶ B) : Prop := IsPreStep P φ ∧ ¬ IsIso φ

/-- **(iii)** `primary pre-step` —— `Div(φ)` が `Φ(A)` の primary 元である pre-step。 -/
def IsPrimaryPreStep {A B : C} (φ : A ⟶ B) : Prop :=
  IsPreStep P φ ∧ IsPrimaryElt (P.Div φ)

/-- **(iii)** `co-angular`。

原文の `φ = α ◦ β ◦ γ` は「先に `γ`」なので、Lean では `φ = γ ≫ β ≫ α`。 -/
def IsCoAngular {A B : C} (φ : A ⟶ B) : Prop :=
  ∀ (X Y : C) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
    φ = γ ≫ β ≫ α → IsLinear P α → IsIsometric P β → IsPreStep P β →
    (IsBaseIsomorphism P α ∨ IsBaseIsomorphism P γ) → IsIso β

/-- **(iii)** `LB-invertible` —— co-angular かつ isometric。 -/
def IsLBInvertible {A B : C} (φ : A ⟶ B) : Prop := IsCoAngular P φ ∧ IsIsometric P φ

/-- **(iii)** `morphism of Frobenius type` —— LB-invertible な base-isomorphism。 -/
def IsFrobeniusType {A B : C} (φ : A ⟶ B) : Prop :=
  IsLBInvertible P φ ∧ IsBaseIsomorphism P φ

/-- **(iii)** `prime-Frobenius morphism` —— Frobenius 次数が素数である Frobenius 型の射。 -/
def IsPrimeFrobenius {A B : C} (φ : A ⟶ B) : Prop :=
  IsFrobeniusType P φ ∧ Nat.Prime ((P.degFr φ : ℕ+) : ℕ)

/-! ### (iv) —— 対象の類型(`Definition 1.3` が引くもの)

原文 (FrdI p.22):
> (iv) A Frobenius-ample object of C is defined to be an object C such that for any

原文 (FrdI p.22):
> object of C is defined to be an object C such that there exists a homomorphism

原文 (FrdI p.22):
> of monoids ζ : N≥1 →EndC(C) which satisfies the following properties: (a) the

原文 (FrdI p.22):
> composite of ζ with the map to N≥1 given by the Frobenius degree is the identity

原文 (FrdI p.22):
> of Frobenius type. A Div-Frobenius-trivial object of C is defined to be an object

原文 (FrdI p.22):
> the image of ζ are Div-identity endomorphisms of Frobenius type. A universally

原文 (FrdI p.22):
> Div-Frobenius-trivial object of C is defined to be an object C such that for every

原文 (FrdI p.22):
> for any n ∈N≥1, C admits a base-identity endomorphism [which is not necessarily

原文 (FrdI p.22):
> is defined to be an object C such that there exists a co-angular pre-step D →C in

原文 (FrdI p.23):
> such that any object D ∈Ob(C) such that Base(C) ∼= Base(D) [in D] is, in fact,

原文 (FrdI p.23):
> an object C such that O×(C) = {1}. An isotropic object [a term motivated by the

原文 (FrdI p.23):
> C such that any isometric pre-step C →D in C is, in fact, an isomorphism. We

原文 (FrdI p.23):
> [of A] if φ is an isometric pre-step, B is isotropic, and for every morphism γ : A →C,

原文 (FrdI p.23):
> where C is isotropic, there exists a unique morphism β : B →C such that γ = β ◦φ.
-/

/-- **(iv)** `Frobenius-ample object` —— 任意の次数の自己射を持つ。 -/
def IsFrobeniusAmple (A : C) : Prop := ∀ n : ℕ+, ∃ φ : A ⟶ A, P.degFr φ = n

/-- **(iv)** `Frobenius-trivial object` —— 条件 (a)(b) を満たすモノイド準同型
`ζ : ℕ≥1 → End_𝒞(C)` が存在する。 -/
def IsFrobeniusTrivial (A : C) : Prop :=
  ∃ ζ : ℕ+ →* End A,
    (∀ n : ℕ+, P.degFr (ζ n) = n) ∧
    (∀ n : ℕ+, IsBaseIdentity P (ζ n) ∧ IsFrobeniusType P (ζ n))

/-- **(iv)** `Div-Frobenius-trivial object` —— 上の (b) を `Div-identity` に替えたもの。 -/
def IsDivFrobeniusTrivial (A : C) : Prop :=
  ∃ ζ : ℕ+ →* End A,
    (∀ n : ℕ+, P.degFr (ζ n) = n) ∧
    (∀ n : ℕ+, IsDivIdentity P (ζ n) ∧ IsFrobeniusType P (ζ n))

/-- **(iv)** `universally Div-Frobenius-trivial object`。 -/
def IsUnivDivFrobeniusTrivial (A : C) : Prop :=
  ∀ (A' : C) (φ : A' ⟶ A), IsPullBack P φ → IsDivFrobeniusTrivial P A'

/-- **(iv)** `quasi-Frobenius-trivial object` —— 任意の次数の
**base-identity** 自己射を持つ(Frobenius 型でなくてよい)。 -/
def IsQuasiFrobeniusTrivial (A : C) : Prop :=
  ∀ n : ℕ+, ∃ φ : A ⟶ A, IsBaseIdentity P φ ∧ P.degFr φ = n

/-- **(iv)** `sub-quasi-Frobenius-trivial object`。 -/
def IsSubQuasiFrobeniusTrivial (A : C) : Prop :=
  ∃ (Dd : C) (α : Dd ⟶ A), IsCoAngular P α ∧ IsPreStep P α ∧ IsQuasiFrobeniusTrivial P Dd

/-- **(iv)** `metrically trivial object`。 -/
def IsMetricallyTrivial (A : C) : Prop :=
  ∀ (Dd : C) (φ : A ⟶ Dd), IsCoAngular P φ → IsPreStep P φ → Nonempty (Dd ≅ A)

/-- **(iv)** `base-trivial object`。 -/
def IsBaseTrivial (A : C) : Prop := ∀ Dd : C, BaseIsomorphic P A Dd → Nonempty (Dd ≅ A)

/-- **(iv)** `group-like object` —— `Φ(C) = 0`。

★`Definition 1.1, (i)` の `IsGroupLike`(`M^char` が零)を、
`isGroupLike_iff` で「全元が可逆」に翻訳した形で使う。 -/
def IsGroupLikeObj (A : C) : Prop := IsGroupLike (Φ.val (P.toElem.obj A).base)

/-- **(iv)** `unit-trivial object` —— `𝒪^×(C) = {1}`。 -/
def IsUnitTrivial (A : C) : Prop := OTimes P A = ⊥

/-- **(iv)** `isotropic object` —— 任意の isometric pre-step `C → D` が同型。 -/
def IsIsotropic (A : C) : Prop :=
  ∀ (Dd : C) (φ : A ⟶ Dd), IsIsometric P φ → IsPreStep P φ → IsIso φ

/-- **(iv)** `isotropic hull` —— isometric pre-step であって、終域が isotropic で、
isotropic な対象への射を**一意に**通す(普遍性)。 -/
def IsIsotropicHull {A B : C} (φ : A ⟶ B) : Prop :=
  IsIsometric P φ ∧ IsPreStep P φ ∧ IsIsotropic P B ∧
    ∀ (Cc : C), IsIsotropic P Cc → ∀ γ : A ⟶ Cc, ∃! β : B ⟶ Cc, γ = φ ≫ β

/-! ### ★Remark 1.2.1 —— 原文が「定義から形式的に従う」と述べる含意

原文 (FrdI p.24):
> Remark 1.2.1. The following implications follow formally from the definitions:

原文 (FrdI p.24):
> pull-back morphism which is a base-isomorphism ⇐⇒ isomorphism

原文 (FrdI p.24):
> base-trivial =⇒metrically trivial

原文 (FrdI p.24):
> base-identity =⇒ Div-identity

原文 (FrdI p.24):
> universally Div-Frobenius-trivial =⇒ Div-Frobenius-trivial

★原文は「follow formally」と書くだけで証明を置かない。**それを実際に示す。**

★**含意は4本である**(2026-08-15 訂正)。一度この docstring は3本しか引用しておらず、
2本目(`base-trivial ⟹ metrically trivial`)を写し落としていた。
下の `isMetricallyTrivial_of_isBaseTrivial` がその4本目である。
-/

/-- ★**`base-identity ⟹ Div-identity`**(Remark 1.2.1 の第3行)。

`Base(φ) = Base(𝟙)` なら、`Φ` を当てた `Φ(φ) = Φ(𝟙)` も従う。 -/
theorem isDivIdentity_of_isBaseIdentity {A : C} {φ : A ⟶ A}
    (h : IsBaseIdentity P φ) : IsDivIdentity P φ := by
  show Φ.map (P.Base φ) = Φ.map (P.Base (𝟙 A))
  rw [show P.Base φ = P.Base (𝟙 A) from h]

/-- ★**`base-trivial ⟹ metrically trivial`**(Remark 1.2.1 の第2行)。

co-angular pre-step は特に pre-step なので base-isomorphism、よって始域と終域は
base-isomorphic。`base-trivial` はそれを同型に格上げする。 -/
theorem isMetricallyTrivial_of_isBaseTrivial {A : C}
    (h : IsBaseTrivial P A) : IsMetricallyTrivial P A := by
  intro Dd φ _ hstep
  exact h Dd ⟨@asIso _ _ _ _ (P.Base φ) hstep.2⟩

/-- ★**`base-equivalent ⟹ Div-equivalent`**(上の一般形)。 -/
theorem divEquivalent_of_baseEquivalent {A B : C} {φ ψ : A ⟶ B}
    (h : BaseEquivalent P φ ψ) : DivEquivalent P φ ψ := by
  show Φ.map (P.Base φ) = Φ.map (P.Base ψ)
  rw [show P.Base φ = P.Base ψ from h]

/-- ★**`universally Div-Frobenius-trivial ⟹ Div-Frobenius-trivial`**
(Remark 1.2.1 の第4行)。

恒等射が pull-back morphism であることを使う。 -/
theorem isDivFrobeniusTrivial_of_univ {A : C}
    (h : IsUnivDivFrobeniusTrivial P A) : IsDivFrobeniusTrivial P A := by
  refine h A (𝟙 A) ?_
  intro X
  constructor
  · intro f g hfg
    have := congrArg (fun p => (p : (X ⟶ A) × _).1) (congrArg Subtype.val hfg)
    simpa using this
  · rintro ⟨⟨g, h'⟩, hgh⟩
    have hb : P.Base g = h' := by simpa using hgh
    exact ⟨g, Subtype.ext (by simp [hb])⟩

/-- ★**`pull-back morphism` かつ `base-isomorphism` ⟸ `isomorphism`**
(Remark 1.2.1 の第1行のうち、示せる向き)。

★もう一方の向き(⟹)は `Definition 1.3`(Frobenioid の条件)を要る——
`Proposition 1.4, (ii)` がそれを扱う。ここでは**片方だけ**示す。 -/
theorem isPullBack_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] : IsPullBack P φ := by
  intro X
  constructor
  · intro f g hfg
    have h1 : f ≫ φ = g ≫ φ :=
      congrArg (fun p => (p : (X ⟶ B) × _).1) (congrArg Subtype.val hfg)
    exact (cancel_mono φ).mp h1
  · rintro ⟨⟨g, h'⟩, hgh⟩
    refine ⟨g ≫ inv φ, Subtype.ext ?_⟩
    have h2 : P.Base φ ≫ P.Base (inv φ) = 𝟙 _ := by
      rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id]
    have hb : P.Base (g ≫ inv φ) = h' := by
      rw [P.Base_comp, hgh, Category.assoc, h2, Category.comp_id]
    simp [hb]

theorem isBaseIsomorphism_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] :
    IsBaseIsomorphism P φ := by
  show IsIso (P.Base φ)
  refine ⟨P.Base (inv φ), ?_, ?_⟩
  · rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id]
  · rw [← P.Base_comp, IsIso.inv_hom_id, P.Base_id]

/-! ### ★形式的に従う小さな含意(定義が空でないことの確認を兼ねる) -/

/-- 恒等射は linear。 -/
theorem isLinear_id (A : C) : IsLinear P (𝟙 A) := by simp [IsLinear]

/-- 恒等射は isometric。 -/
theorem isIsometric_id (A : C) : IsIsometric P (𝟙 A) := by simp [IsIsometric]

/-- 恒等射は pre-step。 -/
theorem isPreStep_id (A : C) : IsPreStep P (𝟙 A) :=
  ⟨isLinear_id P A, isBaseIsomorphism_of_isIso P (𝟙 A)⟩

/-- ★恒等射は **step ではない**(step は「同型でない」pre-step だから)。
`pre-step` と `step` が別物であることの確認。 -/
theorem not_isStep_id (A : C) : ¬ IsStep P (𝟙 A) := fun h => h.2 inferInstance

/-- `linear` は合成で閉じる。 -/
theorem IsLinear.comp {A B E : C} {ψ : A ⟶ B} {φ : B ⟶ E}
    (hψ : IsLinear P ψ) (hφ : IsLinear P φ) : IsLinear P (ψ ≫ φ) := by
  show P.degFr (ψ ≫ φ) = 1
  rw [P.degFr_comp, show P.degFr φ = 1 from hφ, show P.degFr ψ = 1 from hψ, mul_one]

/-- `isometric` は合成で閉じる。

`Div(ψ ≫ φ) = Φ(Base ψ)(Div φ) + deg_Fr(φ) · Div ψ` なので、
両方が isometric なら `Div = 0 + deg · 0 = 0`。★`Φ` が加法を保つこと
(`map_zero`)が効いており、**`linear` は要らない**。 -/
theorem IsIsometric.comp {A B E : C} {ψ : A ⟶ B} {φ : B ⟶ E}
    (hψ : IsIsometric P ψ) (hφ : IsIsometric P φ) : IsIsometric P (ψ ≫ φ) := by
  show P.Div (ψ ≫ φ) = 0
  rw [P.Div_comp, show P.Div φ = 0 from hφ, show P.Div ψ = 0 from hψ, map_zero, smul_zero,
    add_zero]

/-! ### ★出典の紐付け(`.src`) -/

def isDivIdentity_of_isBaseIdentity.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Remark 1.2.1",
    sectionId := "frdi-remark-1-2-1" }

def IsLinear.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (i) — linear",
    sectionId := "frdi-def-1-2-i" }

def IsIsometric.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (i) — isometric",
    sectionId := "frdi-def-1-2-i" }

def MetricallyEquivalent.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (i) — metrically equivalent",
    sectionId := "frdi-def-1-2-i" }

def IsBaseIsomorphism.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (ii) — base-isomorphism",
    sectionId := "frdi-def-1-2-ii" }

def IsBaseFSM.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (ii) — base-FSM-morphism",
    sectionId := "frdi-def-1-2-ii" }

def IsPullBack.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (ii) — pull-back morphism",
    sectionId := "frdi-def-1-2-ii" }

def DivEquivalent.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (ii) — Div-equivalent",
    sectionId := "frdi-def-1-2-ii" }

def OTri.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (ii) — 𝒪^▷(A)",
    sectionId := "frdi-def-1-2-ii" }

def OTimes.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (ii) — 𝒪^×(A)",
    sectionId := "frdi-def-1-2-ii" }

def IsPreStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iii) — pre-step",
    sectionId := "frdi-def-1-2-iii" }

def IsStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iii) — step",
    sectionId := "frdi-def-1-2-iii" }

def IsPrimaryPreStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iii) — primary pre-step",
    sectionId := "frdi-def-1-2-iii" }

def IsCoAngular.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iii) — co-angular",
    sectionId := "frdi-def-1-2-iii" }

def IsLBInvertible.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iii) — LB-invertible",
    sectionId := "frdi-def-1-2-iii" }

def IsFrobeniusType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iii) — morphism of Frobenius type",
    sectionId := "frdi-def-1-2-iii" }

def IsPrimeFrobenius.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iii) — prime-Frobenius morphism",
    sectionId := "frdi-def-1-2-iii" }

def IsFrobeniusAmple.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — Frobenius-ample object",
    sectionId := "frdi-def-1-2-iv" }

def IsFrobeniusTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — Frobenius-trivial object",
    sectionId := "frdi-def-1-2-iv" }

def IsDivFrobeniusTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — Div-Frobenius-trivial object",
    sectionId := "frdi-def-1-2-iv" }

def IsUnivDivFrobeniusTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22,
    item := "Definition 1.2, (iv) — universally Div-Frobenius-trivial object",
    sectionId := "frdi-def-1-2-iv" }

def IsQuasiFrobeniusTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — quasi-Frobenius-trivial object",
    sectionId := "frdi-def-1-2-iv" }

def IsSubQuasiFrobeniusTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22,
    item := "Definition 1.2, (iv) — sub-quasi-Frobenius-trivial object",
    sectionId := "frdi-def-1-2-iv" }

def IsMetricallyTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — metrically trivial object",
    sectionId := "frdi-def-1-2-iv" }

def IsBaseTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — base-trivial object",
    sectionId := "frdi-def-1-2-iv" }

def IsGroupLikeObj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — group-like object",
    sectionId := "frdi-def-1-2-iv" }

def IsUnitTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — unit-trivial object",
    sectionId := "frdi-def-1-2-iv" }

def IsIsotropic.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — isotropic object",
    sectionId := "frdi-def-1-2-iv" }

def IsIsotropicHull.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22, item := "Definition 1.2, (iv) — isotropic hull",
    sectionId := "frdi-def-1-2-iv" }

end ABC3.Found.FrdI
