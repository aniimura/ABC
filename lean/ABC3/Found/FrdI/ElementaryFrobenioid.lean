import ABC3.Found.FrdI.MonoidVocabulary
import ABC3.Found.FrdI.CategoryVocabulary
import Mathlib.CategoryTheory.IsConnected
import Mathlib.Data.PNat.Basic
import Mathlib.Algebra.Category.MonCat.Basic

/-!
# [FrdI] Definition 1.1, (ii)(iii)(iv) —— elementary Frobenioid と pre-Frobenioid

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.19–p.20(**400 dpi 目視確認 2026-08-15**)。

## ★何が既にあって、何が無かったか

前回 `MonoidVocabulary.lean` で `Definition 1.1, (i)` を閉じた。
それは **モノイドの側**(`pre-divisorial` / `divisorial` / `group-like`)である。
`Definition 1.2` が依拠するのは **(ii)(iii)(iv)、すなわち圏の側**であり、
そこはまだ何も無かった。このファイルがそこを作る。

## この節の構成

* `(iii)` **`𝔽_M`** = モノイド `M` に付随する elementary Frobenioid(モノイド。台は `M × ℕ≥1`)
* `(iii)` **`𝔽 = 𝔽_{ℤ≥0}`** = standard Frobenioid
* `(ii)` **`monoid on 𝒟`** = 反変関手 `Φ : 𝒟 → 𝔐𝔬𝔫` で2条件を満たすもの
* `(iii)` **`𝔽_Φ`** = `Φ` に付随する elementary Frobenioid(圏)
* `(iv)` **pre-Frobenioid structure** `𝒞 → 𝔽_Φ`

## ★mathlib の実測(2026-08-15)

`𝔽_M` の積 `(a₁,n₁)·(a₂,n₂) = (a₁ + n₁·a₂, n₁·n₂)` は**半直積の形**をしているが、
mathlib の `SemidirectProduct (φ : G →* MulAut N)` は

* `G`、`N` が**群**であることを要求し
* 作用が**自己同型** (`MulAut N`) であることを要求する

のに対し、ここでの `ℕ≥1` の `M` への作用 `a ↦ n · a` は
**自己同型ではない**(`M = ℕ` で `2 · (−)` は全射でない)。★したがって
mathlib の `SemidirectProduct` は使えない。自前で書く理由はこれである。

使うもの:

* `ℕ≥1` = mathlib の **`ℕ+`**(`PNat`、`CancelCommMonoid ℕ+`)
* `𝔐𝔬𝔫`(可換モノイドの圏) = mathlib の **`AddCommMonCat`**
  (原文 p.11: 「the category of commutative monoids」、演算は加法で書く)
* 反変関手 = `Dᵒᵖ ⥤ AddCommMonCat`
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

/-! ### ★`𝔽_M` —— モノイドに付随する elementary Frobenioid

原文 (FrdI p.20):
> set of FM is the product

原文 (FrdI p.20):
> equipped with the monoid structure is given as follows: if a1, a2 ∈M, n1, n2 ∈N≥1,

原文 (FrdI p.20):
> then (a1, n1) · (a2, n2) = (a1 + n1 · a2, n1 · n2). Also, we shall write F def = FZ≥0 and
-/

/-- **[FrdI] Definition 1.1, (iii)** —— `𝔽_M`。台は `M × ℕ≥1`。

`div` は零因子 `Div(−)`、`deg` は Frobenius 次数 `deg_Fr(−)` に対応する。 -/
@[ext]
structure ElemFrob (M : Type u) [AddCommMonoid M] where
  /-- 零因子(原文の `Div(−)`) -/
  div : M
  /-- Frobenius 次数(原文の `deg_Fr(−)`) -/
  deg : ℕ+

namespace ElemFrob

variable {M : Type u} [AddCommMonoid M]

instance : Mul (ElemFrob M) where
  mul x y := ⟨x.div + (x.deg : ℕ) • y.div, x.deg * y.deg⟩

instance : One (ElemFrob M) where
  one := ⟨0, 1⟩

@[simp] theorem mul_div (x y : ElemFrob M) :
    (x * y).div = x.div + (x.deg : ℕ) • y.div := rfl
@[simp] theorem mul_deg (x y : ElemFrob M) : (x * y).deg = x.deg * y.deg := rfl
@[simp] theorem one_div : (1 : ElemFrob M).div = 0 := rfl
@[simp] theorem one_deg : (1 : ElemFrob M).deg = 1 := rfl

instance : Monoid (ElemFrob M) where
  mul_assoc x y z := by
    ext
    · simp [smul_add, mul_smul, add_assoc]
    · simp [mul_assoc]
  one_mul x := by ext <;> simp
  mul_one x := by ext <;> simp

/-- **Frobenius 次数は乗法的** —— 原文 Remark 1.1.1 の
`deg_Fr(φ ∘ ψ) = deg_Fr(φ) · deg_Fr(ψ)`。

原文 (FrdI p.21):
> Indeed, this follows immediately from the definition of an elementary Frobenioid in
-/
def degHom : ElemFrob M →* ℕ+ where
  toFun := deg
  map_one' := rfl
  map_mul' _ _ := rfl

/-! ### 非退化 -/

/-- ★`𝔽_M` は**可換でない** —— `M = ℕ` での反例。

`(0,2)·(1,1) = (2,2)` だが `(1,1)·(0,2) = (1,2)`。

★これが無いと「`M × ℕ≥1` に何か積を入れた」だけで済んでしまう。
原文の積は**半直積であって直積ではない**。 -/
theorem not_commutative_elemFrob :
    ¬ ∀ x y : ElemFrob ℕ, x * y = y * x := by
  intro h
  have := congrArg div (h ⟨0, 2⟩ ⟨1, 1⟩)
  simp at this

/-- 可逆元の Frobenius 次数は `1` —— 単位群は「次数 1 の層」に閉じ込められる。

★「`𝔽_M` の単位群が全体ではない」ことの具体形。 -/
theorem deg_eq_one_of_isUnit {x : ElemFrob M} (hx : IsUnit x) : x.deg = 1 := by
  obtain ⟨u, rfl⟩ := hx
  have h1 : ((u : ElemFrob M) * (↑u⁻¹ : ElemFrob M)).deg = 1 := by
    rw [u.mul_inv]; rfl
  rw [mul_deg] at h1
  have h2 : ((u : ElemFrob M).deg : ℕ) * ((↑u⁻¹ : ElemFrob M).deg : ℕ) = 1 := by
    exact_mod_cast congrArg (fun n : ℕ+ => (n : ℕ)) h1
  exact PNat.coe_eq_one_iff.mp (Nat.eq_one_of_mul_eq_one_right h2)

/-- ★`𝔽_ℕ` には**可逆でない元がある**(例: `(0, 2)`)。 -/
theorem not_isUnit_deg_two : ¬ IsUnit (⟨0, 2⟩ : ElemFrob ℕ) := by
  intro h
  have := deg_eq_one_of_isUnit h
  exact absurd this (by decide)

/-- **[FrdI] Definition 1.1, (iii)** —— **standard Frobenioid** `𝔽 = 𝔽_{ℤ≥0}`。

原文 (FrdI p.20):
> then (a1, n1) · (a2, n2) = (a1 + n1 · a2, n1 · n2). Also, we shall write F def = FZ≥0 and

原文 (FrdI p.20):
> refer to F as the standard Frobenioid.

★`ℤ_{≥0}` は `ℕ` である。 -/
abbrev Standard : Type := ElemFrob ℕ

/-- standard Frobenioid も非可換(上の反例そのもの)。 -/
theorem standard_not_commutative : ¬ ∀ x y : Standard, x * y = y * x :=
  not_commutative_elemFrob

end ElemFrob

/-! ### ★`monoid on 𝒟` —— [FrdI] Definition 1.1, (ii)

原文 (FrdI p.19):
> as a monoid on D if the following conditions are satisfied: (a) every morphism of

原文 (FrdI p.19):
> monoids α∗ : Φ(A) →Φ(B) induced by a morphism α : B →A of D is char-

原文 (FrdI p.19):
> acteristically injective [cf. §0]; (b) if α is an FSM-morphism [cf. §0] of D, then

原文 (FrdI p.19):
> α∗ : Φ(A) →Φ(B) is an isomorphism of monoids.
-/

variable (D : Type u) [Category.{v} D]

/-- **[FrdI] Definition 1.1, (ii)** —— `monoid on 𝒟`。

`𝒟` 上の**反変**関手 `Φ : 𝒟 → 𝔐𝔬𝔫`(= `Dᵒᵖ ⥤ AddCommMonCat`)であって、
原文の (a)(b) を満たすもの。

★(a) は `MonoidVocabulary.lean` の `IsCharacteristicallyInjective`、
(b) は `CategoryVocabulary.lean` の `IsFSMMorphism` を使う——
**§0 の2語がここで初めて合流する。** -/
structure MonoidOn where
  /-- 反変関手 `Φ : 𝒟 → 𝔐𝔬𝔫` -/
  functor : Dᵒᵖ ⥤ AddCommMonCat.{w}
  /-- **(a)** `α*` は characteristically injective -/
  charInj : ∀ {A B : D} (α : B ⟶ A),
    IsCharacteristicallyInjective (functor.map α.op).hom
  /-- **(b)** `α` が FSM-morphism なら `α*` は同型 -/
  fsmIso : ∀ {A B : D} (α : B ⟶ A), IsFSMMorphism α →
    Function.Bijective (functor.map α.op).hom

variable {D}

namespace MonoidOn

variable (Φ : MonoidOn.{v, u, w} D)

/-- `Φ(A)` の台。 -/
abbrev val (A : D) : Type w := (Φ.functor.obj (op A) : AddCommMonCat.{w})

/-- **`α*`** —— 原文の誘導射 `Φ(A) → Φ(B)`(`α : B ⟶ A` に対し**逆向き**)。 -/
def map {A B : D} (α : B ⟶ A) : Φ.val A →+ Φ.val B := (Φ.functor.map α.op).hom

@[simp] theorem map_id (A : D) (x : Φ.val A) : Φ.map (𝟙 A) x = x := by
  simp [map]

@[simp] theorem map_comp {A B E : D} (α : B ⟶ A) (β : E ⟶ B) (x : Φ.val A) :
    Φ.map (β ≫ α) x = Φ.map β (Φ.map α x) := by
  simp [map, op_comp, Φ.functor.map_comp]

/-- **(a)** は特に `α*` が単射であることを含む。 -/
theorem map_injective {A B : D} (α : B ⟶ A) : Function.Injective (Φ.map α) :=
  (Φ.charInj α).1

/-- `Φ` が **divisorial** であること。

原文 (FrdI p.19):
> Φ(A) [as A ranges over the objects of D] (respectively, some monoid Φ(A) [where
-/
def IsDivisorialOn : Prop := ∀ A : D, IsDivisorial (Φ.val A)

end MonoidOn

/-! ### ★`Φ^char` —— [FrdI] Proposition 1.5 が使うモノイド

原文 (FrdI p.27):
> (i) FΦ, equipped with the natural functor FΦ →FΦchar, is a Frobenioid of Aut-

★`Proposition 1.5` は `Φ` が **pre-divisorial**(sharp とは限らない)のとき、
pre-Frobenioid 構造を **`Φ^char` 経由**で入れる。`PreFrobenioid` は `Φ` が
**divisorial**(= pre-divisorial + sharp)を要求するので、
`Φ^char` が無条件に sharp であること(`isSharp_mChar`)が鍵である。
-/

namespace MonoidOn

/-- `Φ^char` の台となる反変関手。 -/
def charFunctor (Φ : MonoidOn.{v, u, w} D) : Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj X := AddCommMonCat.of (MChar (Φ.val X.unop))
  map {X Y} f := AddCommMonCat.ofHom (charMap (Φ.map f.unop))
  map_id X := by
    refine AddCommMonCat.ext (fun x => ?_)
    obtain ⟨a, rfl⟩ := toChar_surjective (Φ.val X.unop) x
    show charMap (Φ.map (𝟙 X.unop)) (toChar a) = _
    rw [charMap_toChar, Φ.map_id]
    rfl
  map_comp {X Y Z} f g := by
    refine AddCommMonCat.ext (fun x => ?_)
    obtain ⟨a, rfl⟩ := toChar_surjective (Φ.val X.unop) x
    show charMap (Φ.map (g.unop ≫ f.unop)) (toChar a) = _
    rw [charMap_toChar, Φ.map_comp]
    show _ = charMap (Φ.map g.unop) (charMap (Φ.map f.unop) (toChar a))
    rw [charMap_toChar, charMap_toChar]

/-- ★**`Φ^char`** —— `Φ` の characteristic を取った `𝒟` 上のモノイド。

★(a) は `Φ.charInj` の**第2成分そのもの**(もう一段は `isSharp_mChar` で出る)、
(b) は `Φ.fsmIso` から `charMap` の全単射性で出る。 -/
def charOn (Φ : MonoidOn.{v, u, w} D) : MonoidOn.{v, u, w} D where
  functor := Φ.charFunctor
  charInj α := by
    show IsCharacteristicallyInjective (charMap (Φ.map α))
    exact isCharacteristicallyInjective_of_injective_mChar (Φ.charInj α).2
  fsmIso α hα := by
    show Function.Bijective (charMap (Φ.map α))
    exact ⟨(Φ.charInj α).2, charMap_surjective (Φ.fsmIso α hα).2⟩

/-! ### ★★`Φ^gp`（2026-08-16 追加）

原文 (FrdI p.19):
> by assigning A →Φ(A)char, A →Φ(A)gp, A →Φ(A)pf], which we shall refer to,

★★**条件 (a) の第2成分は無料である** —— 値が群なので
`M^char` が自明になる（`charMap_injective_of_addGroup`）。
★第1成分は `gpMap_injective` で、**各 `Φ(A)` の簡約性**を要求する。
★`Φ` が pre-divisorial なら integral であり、`isCancelAdd_of_isIntegralMonoid` が与える。 -/

/-- `Φ^gp` の台となる反変関手。 -/
noncomputable def gpFunctor (Φ : MonoidOn.{v, u, w} D) : Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj X := AddCommMonCat.of (Gp (Φ.val X.unop))
  map {X Y} f := AddCommMonCat.ofHom (gpMap _ (Φ.map f.unop))
  map_id X := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (𝟙 X.unop) = AddMonoidHom.id _ := by
      ext a; exact Φ.map_id _ a
    show gpMap _ (Φ.map (𝟙 X.unop)) x = _
    rw [h, gpMap_id]
    rfl
  map_comp {X Y Z} f g := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (g.unop ≫ f.unop) = (Φ.map g.unop).comp (Φ.map f.unop) := by
      ext a; exact Φ.map_comp _ _ a
    show gpMap _ (Φ.map (g.unop ≫ f.unop)) x = _
    rw [h, gpMap_comp]
    rfl

/-! ★★**原文の一文についての測定（2026-08-16）**

原文 (FrdI p.19) は仮定なしに「Φ determines monoids Φchar, Φgp, Φpf on D」と
述べるが、★★**`Φ^gp` についてはそのままでは成り立たないらしい**。

★**検証役が示した障害**: `ℕ ↪ ℕ∪{∞}` は characteristically injective だが、
`Gp(ℕ) = ℤ → Gp(ℕ∪{∞}) = 0` は単射でない。
つまり**整域性を落とすと条件 (a) が壊れる**。

★**ただしこれはモノイドの層での障害であって、
そういう`Φ` と `D` を実際に構成したわけではない** ——
★★**「整域性なしでは危ない」までは言えるが、「原文が誤っている」とは断定しない。**

★**下流では失われない**: `Definition 1.1, (iv)` 以降 `Φ` は常に divisorial
（integral ∧ sharp を含む）なので、`gpOn` の `hint` も `pfOn` の `hsh` も自動で満たされる。

★**`pfOn` の `hsh` が本当に必要かは UNVERIFIED** である（落とせる可能性がある）。 -/

/-- ★★**`Φ^gp`** —— `Φ` の groupification。

★`hint` は各 `Φ(A)` が integral であること（`Φ` が pre-divisorial なら従う）。 -/
noncomputable def gpOn (Φ : MonoidOn.{v, u, w} D)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) : MonoidOn.{v, u, w} D where
  functor := gpFunctor Φ
  charInj {A B} α := by
    letI := isCancelAdd_of_isIntegralMonoid (Φ.val A) (hint A)
    letI := isCancelAdd_of_isIntegralMonoid (Φ.val B) (hint B)
    show IsCharacteristicallyInjective (gpMap _ (Φ.map α))
    exact ⟨gpMap_injective _ (Φ.map_injective α), charMap_injective_of_addGroup _⟩
  fsmIso {A B} α hα := by
    letI := isCancelAdd_of_isIntegralMonoid (Φ.val A) (hint A)
    letI := isCancelAdd_of_isIntegralMonoid (Φ.val B) (hint B)
    show Function.Bijective (gpMap _ (Φ.map α))
    exact ⟨gpMap_injective _ (Φ.fsmIso α hα).1, gpMap_surjective _ (Φ.fsmIso α hα).2⟩

/-! ### ★★`Φ^pf`（2026-08-16 追加）

★★**`Φ^gp` と違って条件 (a) の第2成分が無料にならない**が、
★**`Φ` が divisorial なら各 `Φ(A)` は sharp で、
`isSharp_pf` で `Pf (Φ(A))` も sharp になる** ——
すると `toChar` が単射になり、`charMap` の単射性が
`Pf.map` の単射性に帰着する。 -/

/-- `Φ^pf` の台となる反変関手。 -/
noncomputable def pfFunctor (Φ : MonoidOn.{v, u, w} D) : Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj X := AddCommMonCat.of (Pf (Φ.val X.unop))
  map {X Y} f := AddCommMonCat.ofHom (Pf.map (Φ.map f.unop))
  map_id X := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (𝟙 X.unop) = AddMonoidHom.id _ := by
      ext a; exact Φ.map_id _ a
    show Pf.map (Φ.map (𝟙 X.unop)) x = _
    rw [h, Pf.map_id]
    rfl
  map_comp {X Y Z} f g := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (g.unop ≫ f.unop) = (Φ.map g.unop).comp (Φ.map f.unop) := by
      ext a; exact Φ.map_comp _ _ a
    show Pf.map (Φ.map (g.unop ≫ f.unop)) x = _
    rw [h, Pf.map_comp]
    rfl

/-- ★★**`Φ^pf`** —— `Φ` の perfection。

★`hsh` は各 `Φ(A)` が sharp であること（`Φ` が divisorial なら従う）。 -/
noncomputable def pfOn (Φ : MonoidOn.{v, u, w} D)
    (hsh : ∀ A : D, IsSharp (Φ.val A)) : MonoidOn.{v, u, w} D where
  functor := pfFunctor Φ
  charInj {A B} α := by
    show IsCharacteristicallyInjective (Pf.map (Φ.map α))
    exact ⟨Pf.map_injective (Φ.map_injective α),
      charMap_pfMap_injective (hsh B) (Φ.map_injective α)⟩
  fsmIso {A B} α hα := by
    show Function.Bijective (Pf.map (Φ.map α))
    exact ⟨Pf.map_injective (Φ.fsmIso α hα).1, Pf.map_surjective (Φ.fsmIso α hα).2⟩

/-! ### ★★`Φ` が `non-dilating`（2026-08-16 追加）

原文 (FrdI p.19):
> divisorial, then we shall say that Φ is non-dilating if the endomorphisms of Φ(A),

★**対象ごとの `IsNonDilating`（`MonoidVocabulary`）を、
`End_𝒟(A)` が誘導する自己準同型全体に課したもの**である。 -/

/-- **[FrdI] Definition 1.1, (ii)** `Φ` が `non-dilating`。 -/
def IsNonDilatingOn (Φ : MonoidOn.{v, u, w} D) : Prop :=
  ∀ (A : D) (e : A ⟶ A), IsNonDilating (Φ.map e)


@[simp] theorem charOn_val (Φ : MonoidOn.{v, u, w} D) (A : D) :
    (Φ.charOn).val A = MChar (Φ.val A) := rfl

/-- ★**`Φ` が pre-divisorial なら `Φ^char` は divisorial** —— `isDivisorial_mChar` の系。

★これが `PreFrobenioid` の `divisorial` フィールドを埋める。 -/
theorem charOn_isDivisorialOn (Φ : MonoidOn.{v, u, w} D)
    (h : ∀ A : D, IsPreDivisorial (Φ.val A)) : (Φ.charOn).IsDivisorialOn :=
  fun A => isDivisorial_mChar (Φ.val A) (h A)

end MonoidOn

/-! ### ★`𝔽_Φ` —— [FrdI] Definition 1.1, (iii) の一般形

原文 (FrdI p.20):
> defined as follows: The objects of FΦ are the objects of D. If A, B ∈Ob(FΦ), whose

原文 (FrdI p.20):
> of FΦ is defined to be a collection of data

原文 (FrdI p.20):
> where φD : AD →BD is a morphism of D; Zφ ∈Φ(AD); nφ ∈N≥1. Here,

原文 (FrdI p.20):
> is given as follows:
-/

/-- `𝔽_Φ` の対象 = `𝒟` の対象。★型としては `𝒟` そのものだが、
圏の構造が違うので新しい型として包む。 -/
structure ElemFrobCat (Φ : MonoidOn.{v, u, w} D) where
  /-- `𝒟` の対象 -/
  base : D

namespace ElemFrobCat

variable {Φ : MonoidOn.{v, u, w} D}

/-- `𝔽_Φ` の射 `A → B` = 3つ組 `(φ_𝒟, Z_φ, n_φ)`。

原文 (FrdI p.20):
> where φD : AD →BD is a morphism of D; Zφ ∈Φ(AD); nφ ∈N≥1. Here,
-/
@[ext]
structure Hom (A B : ElemFrobCat Φ) where
  /-- `Base(φ)` -/
  base : A.base ⟶ B.base
  /-- `Div(φ) ∈ Φ(A_𝒟)` -/
  div : Φ.val A.base
  /-- `deg_Fr(φ) ∈ ℕ≥1` -/
  deg : ℕ+

/-- 合成。

原文 (FrdI p.20):
> ψ ◦φ = (ψD ◦φD, φ∗ D(Zψ) + nψ · Zφ, nψ · nφ) : A →C
-/
def Hom.comp {A B E : ElemFrobCat Φ} (φ : Hom A B) (ψ : Hom B E) : Hom A E where
  base := φ.base ≫ ψ.base
  div := Φ.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div
  deg := ψ.deg * φ.deg

instance : Category (ElemFrobCat Φ) where
  Hom := Hom
  id A := ⟨𝟙 A.base, 0, 1⟩
  comp := Hom.comp
  id_comp φ := by
    refine Hom.ext ?_ ?_ ?_
    · exact Category.id_comp _
    · show Φ.map (𝟙 _) φ.div + (φ.deg : ℕ) • (0 : Φ.val _) = φ.div
      simp
    · show φ.deg * 1 = φ.deg
      simp
  comp_id φ := by
    refine Hom.ext ?_ ?_ ?_
    · exact Category.comp_id _
    · show Φ.map φ.base (0 : Φ.val _) + ((1 : ℕ+) : ℕ) • φ.div = φ.div
      simp
    · show (1 : ℕ+) * φ.deg = φ.deg
      simp
  assoc φ ψ χ := by
    refine Hom.ext ?_ ?_ ?_
    · exact Category.assoc _ _ _
    · show Φ.map (φ.base ≫ ψ.base) χ.div
          + (χ.deg : ℕ) • (Φ.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div)
        = Φ.map φ.base (Φ.map ψ.base χ.div + (χ.deg : ℕ) • ψ.div)
          + ((χ.deg * ψ.deg : ℕ+) : ℕ) • φ.div
      rw [Φ.map_comp, map_add, map_nsmul, smul_add, ← add_assoc, ← mul_smul]
      norm_cast
    · simp [Hom.comp, mul_assoc]

@[simp] theorem comp_base {A B E : ElemFrobCat Φ} (φ : A ⟶ B) (ψ : B ⟶ E) :
    (φ ≫ ψ).base = φ.base ≫ ψ.base := rfl
@[simp] theorem comp_deg {A B E : ElemFrobCat Φ} (φ : A ⟶ B) (ψ : B ⟶ E) :
    (φ ≫ ψ).deg = ψ.deg * φ.deg := rfl
/-- ★★**合成の `div` 成分**(2026-08-15 追加)。

★**これが無かった。** `comp_base` / `comp_deg` / `id_base` / `id_deg` / `id_div` は
`@[simp]` で揃っていたのに、**`comp_div` だけが欠けていた。**

★子は「構造体の射影・合成が、定義的には計算できるのに `simp`/`rw` が触れない」と
報告し、それを **Lean 技術的な未完 4 件に共通する壁**だと見立てた。見立ては正しく、
★**壁の正体は「`simp` 集合に穴が 1 つ開いていた」ことだった。**
`.div` を含む目標だけが閉じなかったのは、そこだけ補題が無かったからである。

★**測定**: 症状は 4 箇所に別々に現れたが、原因は 1 つで、しかも
**「規則が難しい」のではなく「規則が 1 本足りない」**だった。
★我々は「症状を数えるのではなく原因を特定する」を繰り返してきたが、
**今回は原因が「欠落」であって「機構」ではなかった。** 欠落は、
**揃っているはずのものを並べて見ないと見えない。** -/
@[simp] theorem comp_div {A B E : ElemFrobCat Φ} (φ : A ⟶ B) (ψ : B ⟶ E) :
    (φ ≫ ψ).div = Φ.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div := rfl
@[simp] theorem id_base (A : ElemFrobCat Φ) : Hom.base (𝟙 A) = 𝟙 A.base := rfl
@[simp] theorem id_deg (A : ElemFrobCat Φ) : Hom.deg (𝟙 A) = 1 := rfl
@[simp] theorem id_div (A : ElemFrobCat Φ) : Hom.div (𝟙 A) = 0 := rfl

/-- **自然な射影関手 `𝔽_Φ → 𝒟`**。

原文 (FrdI p.20):
> functor:
-/
def proj : ElemFrobCat Φ ⥤ D where
  obj A := A.base
  map φ := φ.base

/-! ### ★★`Φ ↦ 𝔽_Φ` の関手性（2026-08-16 追加）

原文 (FrdI p.20):
> Observe that the assignment Φ →FΦ is functorial with respect to homomorphisms

★★**監査で「指摘したのはこちらだ」と区別された項目**である ——
`PreFrobenioid.divisorFunctor`（原文 p.21）は `Φ` を**固定して** `𝒞` 上へ持ち上げるもので、
★こちらは **`Φ` を動かす**。★**量化している対象が違う。**

★原文の「homomorphisms of functors [on D] valued in monoids」は
★**単なる自然変換** `Φ.functor ⟶ Φ'.functor` である（(a)(b) は各々の `Φ` の側の条件）。

★★**関手性の非自明な部分は `Div` 成分だけで、
そこで使うのが `η` の自然性である** —— `η (Φ.map f x) = Φ'.map f (η x)`。 -/

/-- ★★**`Φ ↦ 𝔽_Φ` の関手性**。自然変換 `η : Φ ⟶ Φ'` が
関手 `𝔽_Φ ⥤ 𝔽_{Φ'}` を誘導する。 -/
def elemFrobMap {Φ Ψ : MonoidOn.{v, u, w} D} (η : Φ.functor ⟶ Ψ.functor) :
    ElemFrobCat Φ ⥤ ElemFrobCat Ψ where
  obj A := ⟨A.base⟩
  map {A B} φ := ⟨φ.base, (η.app (op A.base)).hom φ.div, φ.deg⟩
  map_id A := by
    refine ElemFrobCat.Hom.ext rfl ?_ rfl
    show (η.app (op A.base)).hom 0 = 0
    exact map_zero _
  map_comp {A B E} φ ψ := by
    refine ElemFrobCat.Hom.ext rfl ?_ rfl
    show (η.app (op A.base)).hom (Φ.map φ.base ψ.div + ((ψ.deg : ℕ+) : ℕ) • φ.div)
      = Ψ.map φ.base ((η.app (op B.base)).hom ψ.div)
        + ((ψ.deg : ℕ+) : ℕ) • (η.app (op A.base)).hom φ.div
    rw [map_add, map_nsmul]
    congr 1
    exact congrArg (fun t => (AddCommMonCat.Hom.hom t) ψ.div) (η.naturality φ.base.op)

/-- ★★**割り当てが関手的であること（恒等）**。

★★**監査で「各 `η` が関手を誘導する」までは示せているが、
「割り当てが**関手的である**」は示せていないと指摘された**。
★`pushFunctor` に `map_id`/`map_comp` がなければ関手と呼べないのと同じ形である。 -/
@[simp] theorem elemFrobMap_id (Φ : MonoidOn.{v, u, w} D) :
    elemFrobMap (𝟙 Φ.functor) = 𝟭 (ElemFrobCat Φ) := rfl

/-- ★★**割り当てが関手的であること（合成）**。 -/
theorem elemFrobMap_comp {Φ Ψ Θ : MonoidOn.{v, u, w} D}
    (η : Φ.functor ⟶ Ψ.functor) (θ : Ψ.functor ⟶ Θ.functor) :
    elemFrobMap (η ≫ θ) = elemFrobMap η ⋙ elemFrobMap θ := rfl

/-! ★**`elemFrobToChar`（`Prop15`）はこの特殊例である** —— 監査の指摘。
★実際の写像は下流にあるので、ここでは形の一致だけを記録する。 -/

/-! ### ★Remark 1.1.1 の3式

原文 (FrdI p.21):
> Base(φ ◦ψ) = Base(φ) ◦Base(ψ)

原文 (FrdI p.21):
> Div(φ ◦ψ) = (Base(ψ))∗(Div(φ)) + degFr(φ) · Div(ψ)

原文 (FrdI p.21):
> degFr(φ ◦ψ) = degFr(φ) · degFr(ψ)

原文 (FrdI p.21):
> Indeed, this follows immediately from the definition of an elementary Frobenioid in

★原文は `φ ◦ ψ`(先に `ψ`)で書く。Lean の `ψ ≫ φ` が同じものである。
★「follows immediately」と書かれた3式を**実際に示す**。
-/

/-- `Base(φ ∘ ψ) = Base(φ) ∘ Base(ψ)`。 -/
theorem base_comp {A B E : ElemFrobCat Φ} (ψ : A ⟶ B) (φ : B ⟶ E) :
    (ψ ≫ φ).base = ψ.base ≫ φ.base := rfl

/-- `Div(φ ∘ ψ) = (Base(ψ))*(Div(φ)) + deg_Fr(φ) · Div(ψ)`。 -/
theorem div_comp {A B E : ElemFrobCat Φ} (ψ : A ⟶ B) (φ : B ⟶ E) :
    (ψ ≫ φ).div = Φ.map ψ.base φ.div + (φ.deg : ℕ) • ψ.div := rfl

def div_comp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Remark 1.1.1",
    sectionId := "frdi-remark-1-1-1" }

/-- `deg_Fr(φ ∘ ψ) = deg_Fr(φ) · deg_Fr(ψ)`。 -/
theorem degFr_comp {A B E : ElemFrobCat Φ} (ψ : A ⟶ B) (φ : B ⟶ E) :
    (ψ ≫ φ).deg = φ.deg * ψ.deg := rfl

/-- ★**divisor monoid が sharp なら、`𝔽_Φ` の同型は零因子 0 を持つ**。

`φ ≫ φ⁻¹ = 𝟙` の零因子成分は
`Φ(Base φ)(Div φ⁻¹) + deg_Fr(φ⁻¹) · Div φ = 0`。
sharp なモノイドでは `x + y = 0` から `y` が可逆、よって `y = 0`。
さらに `deg_Fr(φ⁻¹) ≥ 1` なので、`MonoidVocabulary.lean` の
`isTorsionFreeNaive_of_isSharp` により `Div φ = 0`。

★これは `Definition 1.2` の `isometric` と `isomorphism` を結ぶ含意であり、
`Definition 1.2` の非退化 witness(`Witness.lean`)で使う。 -/
theorem div_eq_zero_of_isIso (hsharp : ∀ A : D, IsSharp (Φ.val A))
    {A B : ElemFrobCat Φ} (φ : A ⟶ B) [IsIso φ] : Hom.div φ = 0 := by
  have hd := congrArg Hom.div (IsIso.hom_inv_id φ)
  rw [div_comp, id_div] at hd
  have hunit : IsAddUnit (((Hom.deg (inv φ) : ℕ+) : ℕ) • Hom.div φ) :=
    ⟨⟨_, Φ.map (Hom.base φ) (Hom.div (inv φ)), by rw [add_comm]; exact hd, hd⟩, rfl⟩
  have hz := hsharp A.base _ hunit
  exact isTorsionFreeNaive_of_isSharp (hsharp A.base) _ _ (Hom.deg (inv φ)).pos hz

/-- ★★**`𝔽_Φ` の同型の完全な判定** ——
`φ` が同型 ⟺ 「底が同型」「零因子が**可逆**」「Frobenius 次数が 1」。

★`div_eq_zero_of_isIso` は `Φ(A)` が **sharp** な場合の系である
(sharp なら可逆元は `0` だけ)。**一般には「可逆」であって「零」ではない。**
`Proposition 1.5` が `𝔽_Φ → 𝔽_{Φ^char}` を経由する理由がここにある。 -/
theorem isIso_iff {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    IsIso φ ↔ IsIso (Hom.base φ) ∧ IsAddUnit (Hom.div φ) ∧ Hom.deg φ = 1 := by
  constructor
  · intro h
    have hd1 := congrArg Hom.deg (IsIso.hom_inv_id φ)
    rw [comp_deg, id_deg] at hd1
    have hd1' : ((Hom.deg (inv φ) : ℕ+) : ℕ) * ((Hom.deg φ : ℕ+) : ℕ) = 1 := by
      exact_mod_cast congrArg (fun n : ℕ+ => (n : ℕ)) hd1
    have hdegi : Hom.deg (inv φ) = 1 :=
      PNat.coe_eq_one_iff.mp (Nat.dvd_one.mp ⟨(Hom.deg φ : ℕ+), hd1'.symm⟩)
    have hdeg : Hom.deg φ = 1 :=
      PNat.coe_eq_one_iff.mp
        (Nat.dvd_one.mp ⟨(Hom.deg (inv φ) : ℕ+), by rw [mul_comm]; exact hd1'.symm⟩)
    have hdv := congrArg Hom.div (IsIso.hom_inv_id φ)
    rw [div_comp, id_div, hdegi] at hdv
    simp only [PNat.one_coe, one_smul] at hdv
    refine ⟨⟨Hom.base (inv φ), ?_, ?_⟩, ⟨⟨Hom.div φ, Φ.map (Hom.base φ) (Hom.div (inv φ)),
      by rw [add_comm]; exact hdv, hdv⟩, rfl⟩, hdeg⟩
    · rw [← comp_base, IsIso.hom_inv_id, id_base]
    · rw [← comp_base, IsIso.inv_hom_id, id_base]
  · rintro ⟨hbase, ⟨U, hU⟩, hdeg⟩
    refine ⟨⟨Hom.mk (inv (Hom.base φ)) (Φ.map (inv (Hom.base φ)) U.neg) 1, ?_, ?_⟩⟩
    · refine Hom.ext (IsIso.hom_inv_id _) ?_ ?_
      · show Φ.map (Hom.base φ) (Φ.map (inv (Hom.base φ)) U.neg) + ((1 : ℕ+) : ℕ) • Hom.div φ
          = (0 : Φ.val A.base)
        rw [← Φ.map_comp, IsIso.hom_inv_id, Φ.map_id]
        simp only [PNat.one_coe, one_smul, ← hU]
        rw [add_comm]
        exact U.val_neg
      · show (1 : ℕ+) * Hom.deg φ = 1
        simp [hdeg]
    · refine Hom.ext (IsIso.inv_hom_id _) ?_ ?_
      · show Φ.map (inv (Hom.base φ)) (Hom.div φ)
            + ((Hom.deg φ : ℕ+) : ℕ) • Φ.map (inv (Hom.base φ)) U.neg
          = (0 : Φ.val B.base)
        rw [hdeg]
        simp only [PNat.one_coe, one_smul, ← map_add, ← hU]
        rw [U.val_neg, map_zero]
      · show Hom.deg φ * 1 = 1
        simp [hdeg]

/-- ★**負の対照(`sub-automorphism`)** —— Frobenius 次数が `1` でない自己射は
**sub-automorphism ではない**。

`β ≫ φ = φ ≫ α` の次数成分は `deg φ · deg β = deg α · deg φ` であり、
`β` が同型なら `deg β = 1` なので `ℕ+` の簡約で `deg α = 1` が出る。

★`CategoryVocabulary.lean` の `isSubAutomorphism_of_isInitial` が示すとおり、
始対象を持つ圏ではこの反例は作れない。**`𝔽_Φ` を要した**のはそのためである。 -/
theorem not_isSubAutomorphism_of_deg_ne_one {A : ElemFrobCat Φ} (α : End A)
    (h : Hom.deg α ≠ 1) : ¬ IsSubAutomorphism α := by
  rintro ⟨B, φ, β, hβ, hcomm⟩
  have hβdeg : Hom.deg β = 1 := ((isIso_iff β).mp hβ).2.2
  have hh := congrArg Hom.deg hcomm
  rw [comp_deg, comp_deg, hβdeg, mul_one] at hh
  have hc : (1 : ℕ+) * Hom.deg φ = Hom.deg α * Hom.deg φ := by rw [one_mul]; exact hh
  exact h (mul_right_cancel hc).symm

end ElemFrobCat

/-- ★**`𝔽_Φ` は totally epimorphic である**。

原文 (FrdI p.27):
> Proof. Since D is a connected, totally epimorphic category, the fact that FΦ is

原文 (FrdI p.27):
> tion 1.1, (iii); the fact that a pre-divisorial monoid is integral [cf. Definition 1.1,

★これは `Proposition 1.5` の証明の**第1段**である。原文が挙げる入力は3つ:

1. **`𝒟` が totally epimorphic**(`hD`)
2. **`Φ(A)` が integral**(`hint`)—— pre-divisorial なら integral
3. **`Definition 1.1, (ii), (a)` の単射性** —— `MonoidOn.charInj` として構造に入っている

★**3つがそのまま3成分に対応する**: `base` は 1 で、`div` は 2 と 3 で、`deg` は `ℕ+` の簡約で消える。

★**2026-08-15 訂正の記録**: 一度この補題を
「`𝒟` の hom が subsingleton」+「`Φ(A)` が `IsCancelAdd`」という**原文より強い仮定**で書いていた。
模型(`Vee` は前順序)では通るが**一般には使えない**形で、`Proposition 1.5` を阻む。
★「模型で通る」と「一般に使える」は違う、という実例。 -/
theorem isTotallyEpimorphic_elemFrobCat {Φ : MonoidOn.{v, u, w} D}
    (hD : IsTotallyEpimorphic D) (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) :
    IsTotallyEpimorphic (ElemFrobCat Φ) := by
  intro X Y f
  refine ⟨fun {Z} g h hgh => ?_⟩
  -- `deg`: `ℕ+` は簡約的
  have hdeg : ElemFrobCat.Hom.deg g = ElemFrobCat.Hom.deg h := by
    have hh := congrArg ElemFrobCat.Hom.deg hgh
    simp only [ElemFrobCat.comp_deg] at hh
    exact mul_right_cancel hh
  -- `base`: ★入力1 —— `𝒟` が totally epimorphic なので `Base f` は epi
  have hbase : ElemFrobCat.Hom.base g = ElemFrobCat.Hom.base h := by
    haveI : Epi (ElemFrobCat.Hom.base f) := hD _ _ _
    have hh := congrArg ElemFrobCat.Hom.base hgh
    simp only [ElemFrobCat.comp_base] at hh
    exact (cancel_epi (ElemFrobCat.Hom.base f)).mp hh
  -- `div`: ★入力2 と 3
  have hdiv : ElemFrobCat.Hom.div g = ElemFrobCat.Hom.div h := by
    have hh := congrArg ElemFrobCat.Hom.div hgh
    simp only [ElemFrobCat.div_comp, hdeg] at hh
    letI := isCancelAdd_of_isIntegralMonoid _ (hint X.base)
    exact Φ.map_injective (ElemFrobCat.Hom.base f) (add_right_cancel hh)
  exact ElemFrobCat.Hom.ext hbase hdiv hdeg

/-! ### ★`pre-Frobenioid structure` —— [FrdI] Definition 1.1, (iv)

原文 (FrdI p.20):
> (iv) Let D, Φ, FΦ be as in (iii); C a category. Assume further that Φ is divisorial,

原文 (FrdI p.20):
> and that C, D are connected, totally epimorphic categories [cf. §0]. Then we shall

原文 (FrdI p.20):
> as a pre-Frobenioid structure on C. The natural projection functor FΦ →D thus
-/

/-- ★**`𝒟` が connected なら `𝔽_Φ` も connected**(2026-08-16 に追加)。

★`𝔽_Φ` の対象は `𝒟` の対象と 1 対 1 で、`𝒟` の射 `f` は `⟨f, 0, 1⟩` として持ち上がる。
★したがって `𝒟` の連結成分がそのまま `𝔽_Φ` の連結成分になる。 -/
theorem isConnected_elemFrobCat (Φ : MonoidOn.{v, u, w} D) [IsConnected D] :
    IsConnected (ElemFrobCat Φ) := by
  obtain ⟨d₀⟩ := (inferInstance : Nonempty D)
  refine IsConnected.of_induct (j₀ := (⟨d₀⟩ : ElemFrobCat Φ)) ?_
  intro p hp0 hstep A
  have key : ∀ d : D, (⟨d⟩ : ElemFrobCat Φ) ∈ p :=
    induct_on_objects (J := D) {d | (⟨d⟩ : ElemFrobCat Φ) ∈ p} hp0
      (fun {d₁ d₂} f => hstep (⟨f, 0, 1⟩ : (⟨d₁⟩ : ElemFrobCat Φ) ⟶ (⟨d₂⟩ : ElemFrobCat Φ)))
  exact key A.base

/-- **[FrdI] Definition 1.1, (iv)** —— `𝒞` 上の pre-Frobenioid structure。

★★**この docstring は誤っていた。2026-08-16 に訂正した。**

以前ここには「`connected`(**圏**について)は §0 が**対象**についてしか定めていないため、
ここには入れていない。原文の曖昧さであり…」と書かれていた。
★★**これは事実に反する。** 原文 §0(物理 p.16)は圏について定義している:

原文 (FrdI p.16):
> as a connected component of C. In particular, we shall say that C is connected if

★しかも**同じリポジトリの `CategoryVocabulary.lean` が同じ行を引用して
「§0 は圏についての `connected` を定義している」と正しく書いていた** ——
★★**2 つのファイルが矛盾しており、誤った側の理由で条件が落ちていた。**
文脈を持たない検証役の監査で判明した。

★**現状(測定)**: 原文 `Definition 1.1, (iv)` は

原文 (FrdI p.20):
> and that C, D are connected, totally epimorphic categories [cf. §0]. Then we shall

と述べるが、下のフィールドには `connected` が **2 つとも無い**。
したがって我々の `PreFrobenioid` は**原文より広い**(定理は強くなるが、
★**定義そのものは原文と一致していない**)。

★**足すときの手**(測定済み): mathlib の `CategoryTheory.IsConnected` を使う。
構成サイトは 6 つ —— `wP` / `nP` / `cP` / `elemPreFrobenioid` /
`cfpPreFrobenioid` / `istrPre`。★`Istr P` の連結性は `isotropification` が
`C` の zigzag を送ることから出る(isotropic な `A` では `hullMap` が同型)。

★現在フィールドとして持つのは:
`Φ` が divisorial、`𝒞` と `𝒟` が totally epimorphic、そして関手 `𝒞 → 𝔽_Φ`。 -/
structure PreFrobenioid (C : Type u2) [Category.{v2} C] (Φ : MonoidOn.{v, u, w} D) where
  /-- 原文の `𝒞 → 𝔽_Φ` -/
  toElem : C ⥤ ElemFrobCat Φ
  /-- `Φ` は divisorial -/
  divisorial : Φ.IsDivisorialOn
  /-- `𝒞` は totally epimorphic -/
  totEpiC : IsTotallyEpimorphic C
  /-- `𝒟` は totally epimorphic -/
  totEpiD : IsTotallyEpimorphic D
  /-- ★`𝒞` は connected(2026-08-16 に追加。原文 p.20 の仮定) -/
  connectedC : IsConnected C
  /-- ★`𝒟` は connected(同上) -/
  connectedD : IsConnected D

namespace PreFrobenioid

variable {C : Type u2} [Category.{v2} C] {Φ : MonoidOn.{v, u, w} D}

/-- **自然な射影関手 `𝒞 → 𝒟`**。

原文 (FrdI p.20):
> restricts to a natural projection functor
-/
def proj (P : PreFrobenioid C Φ) : C ⥤ D := P.toElem ⋙ ElemFrobCat.proj

/-- ★★**「`Φ` を `𝒞` 上の関手とみなす」**（2026-08-16 追加）。

原文 (FrdI p.21):
> often regard Φ as a functor on C [i.e., by composing the original functor Φ with

★**原文はこれを abuse of notation として導入している**が、
★★**Lean では同じ名前を使い回せないので、別の定義として置く必要がある**。
★原文が「記法の乱用」で済ませているものは、形式化すると
**別の射になる** —— この項目での 1 例目。 -/
def divisorFunctor (P : PreFrobenioid C Φ) : Cᵒᵖ ⥤ AddCommMonCat.{w} :=
  P.proj.op ⋙ Φ.functor

@[simp] theorem divisorFunctor_obj (P : PreFrobenioid C Φ) (A : C) :
    ((divisorFunctor P).obj (op A) : AddCommMonCat.{w}) = Φ.functor.obj (op (P.proj.obj A)) := rfl

/-- **`Base(−)`** —— `𝒞` の射に制限した操作。 -/
abbrev Base (P : PreFrobenioid C Φ) {A B : C} (φ : A ⟶ B) :
    (P.toElem.obj A).base ⟶ (P.toElem.obj B).base := (P.toElem.map φ).base

/-- **`Div(−)`** —— `𝒞` の射に制限した操作。 -/
abbrev Div (P : PreFrobenioid C Φ) {A B : C} (φ : A ⟶ B) : Φ.val (P.toElem.obj A).base :=
  (P.toElem.map φ).div

/-- **`deg_Fr(−)`** —— `𝒞` の射に制限した操作。 -/
abbrev degFr (P : PreFrobenioid C Φ) {A B : C} (φ : A ⟶ B) : ℕ+ := (P.toElem.map φ).deg

@[simp] theorem Base_id (P : PreFrobenioid C Φ) (A : C) : P.Base (𝟙 A) = 𝟙 _ := by
  simp [PreFrobenioid.Base]

@[simp] theorem degFr_id (P : PreFrobenioid C Φ) (A : C) : P.degFr (𝟙 A) = 1 := by
  simp [PreFrobenioid.degFr]

@[simp] theorem Div_id (P : PreFrobenioid C Φ) (A : C) : P.Div (𝟙 A) = 0 := by
  simp [PreFrobenioid.Div]

/-- ★**Remark 1.1.1 の3式は `𝒞` の側でも成り立つ**(関手性から)。 -/
theorem degFr_comp (P : PreFrobenioid C Φ) {A B E : C} (ψ : A ⟶ B) (φ : B ⟶ E) :
    P.degFr (ψ ≫ φ) = P.degFr φ * P.degFr ψ := by
  simp [degFr, P.toElem.map_comp]

/-- `Base(−)` の合成則(`𝒞` の側)。 -/
theorem Base_comp (P : PreFrobenioid C Φ) {A B E : C} (ψ : A ⟶ B) (φ : B ⟶ E) :
    P.Base (ψ ≫ φ) = P.Base ψ ≫ P.Base φ := by
  simp [PreFrobenioid.Base, P.toElem.map_comp]

/-- `Div(−)` の合成則(`𝒞` の側)。 -/
theorem Div_comp (P : PreFrobenioid C Φ) {A B E : C} (ψ : A ⟶ B) (φ : B ⟶ E) :
    P.Div (ψ ≫ φ) = Φ.map (P.Base ψ) (P.Div φ) + (P.degFr φ : ℕ) • P.Div ψ := by
  simp [PreFrobenioid.Div, PreFrobenioid.Base, PreFrobenioid.degFr, P.toElem.map_comp,
    ElemFrobCat.div_comp]

end PreFrobenioid

/-! ### ★出典の紐付け(`.src`) -/

def ElemFrob.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 20, item := "Definition 1.1, (iii) — 𝔽_M = M × ℕ≥1",
    sectionId := "frdi-def-1-1-iii" }

def MonoidOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 19, item := "Definition 1.1, (ii) — monoid on 𝒟",
    sectionId := "frdi-def-1-1-ii" }

def ElemFrobCat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 20, item := "Definition 1.1, (iii) — 𝔽_Φ",
    sectionId := "frdi-def-1-1-iii" }

def PreFrobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 20, item := "Definition 1.1, (iv) — pre-Frobenioid structure",
    sectionId := "frdi-def-1-1-iv" }

/-! ### ★★★`Definition 1.1` 全体(条なしの `.src`)

★条なしの `.src` は「その原典項目を**完全に**実装した」という主張である。
`Definition 1.1` は (i)〜(iv) の 18 項目がすべて実装されたので、ここで付ける。

* **(i)** 単元・可換性の語彙: `IsPreDivisorial` / `IsDivisorial` / `IsGroupLike` /
  `isDivisorial_mChar` / `IsNonDilating`
* **(ii)** `𝒟` 上の単系: `MonoidOn` / `charOn` / `gpOn` / `pfOn` / `IsNonDilatingOn`
* **(iii)** `𝔽_Φ` と `𝔽_M`: `ElemFrobCat` / `proj` / `elemFrobMap` ＋ 関手則 2 本
  (`elemFrobMap_id` / `elemFrobMap_comp`) / `ElemFrob` / `Standard` /
  `constPhi` / `constObj` / `constEndEquiv`
* **(iv)** `PreFrobenioid` / `PreFrobenioid.divisorFunctor`

★★**最後に埋まったのは (iii) の `Φ ↦ 𝔽_Φ` の関手性**である。
`elemFrobMap_id` / `elemFrobMap_comp` は結局 `rfl` で通った ——
恒等側は射の構造体の eta、合成側は関手圏の合成が成分ごとであることによる。
★「`rfl` で通るから中身が無い」ではない: `ElemFrobCat Φ` は一般に
subsingleton ではなく、文そのものは原文が要求する関手則である
(検証役の監査で確認した)。 -/
def definition_1_1.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 19, item := "Definition 1.1",
    sectionId := "frdi-def-1-1-i" }

end ABC3.Found.FrdI
