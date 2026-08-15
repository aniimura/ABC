import ABC3.Found.FrdI.Prop111

/-!
# [FrdI] Proposition 1.13 —— Rigidity and Slimness

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.39–40。§0 の語彙は物理 p.14(**600 dpi 拡大確認 2026-08-16**)。

## ★この命題の規模(着手前の測定)

**3 条 (i)–(iii)、主張は 3**:

| 条 | 主張 | 内容 |
|---|---|---|
| (i) | 1 | `𝒞_A → 𝒟` が **rigid**(`𝒟` が slim のとき) |
| (ii) | 1 | `𝒞_A → 𝔽_Φ` が **rigid** |
| (iii) | 1 | 条件 (a) または (b) の下で `𝒞` が **slim** |

★★**残る項目の中で最小である**(`1.10` は 21、`1.11` は 15)。

## ★不足していた語彙

`slim` と `rigid`(§0、物理 p.14)が未実装だった。ここで足す。

原文 (FrdI p.14):
> that φ is rigid if Aut(φ) is trivial. A category C will be called slim if the natural

★**`Aut(φ)` は関手の自己同型**(自然同型 `φ ≅ φ`)であって、対象の自己同型ではない。
★**`𝒞_A` はスライス圏**で、「natural functor `𝒞_A → 𝒞`」は忘却関手である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

/-! ## ★§0 の語彙 —— `rigid` と `slim`

原文 (FrdI p.14):
> — the monoid of natural transformations from the functor φ to itself. We shall say

原文 (FrdI p.14):
> functor CA →C is rigid, for every A ∈Ob(C). We recall that if Π is a profinite
-/

/-- **[FrdI] §0** 関手が `rigid` —— `Aut(φ)` が自明。

★**`Aut(φ)` は関手の自己同型**(自然同型 `φ ≅ φ`)である。 -/
def IsRigidFunctor {C₁ : Type u} [Category.{v} C₁] {C₂ : Type u2} [Category.{v2} C₂]
    (F : C₁ ⥤ C₂) : Prop := ∀ η : F ≅ F, η = Iso.refl F

/-- **[FrdI] §0** 圏が `slim` —— すべての `A` についてスライス `𝒞_A → 𝒞` が rigid。

★**`𝒞_A` はスライス圏**(`A` へ入る射のなす圏)で、
「natural functor」は忘却関手 `Over.forget A` である。 -/
def IsSlimCat (C₁ : Type u) [Category.{v} C₁] : Prop :=
  ∀ A : C₁, IsRigidFunctor (Over.forget A)

def IsRigidFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — rigid",
    sectionId := "frdi-s0-rigid-slim" }

def IsSlimCat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — slim",
    sectionId := "frdi-s0-rigid-slim" }

/-! ### ★非退化を両側から

★**定義を置いたら、空でも全体でもないことを確かめる**(我々の作法)。
-/

/-- ★**非退化(下)**: 終域が薄い圏(射が高々1本)なら、どんな関手も rigid。

★`PUnit` を終域に取れば、自然変換の成分がすべて `𝟙` しかないので `Aut` は自明。 -/
theorem isRigidFunctor_to_punit {C₁ : Type u} [Category.{v} C₁]
    (F : C₁ ⥤ Discrete PUnit) : IsRigidFunctor F := by
  intro η
  apply Iso.ext
  apply NatTrans.ext
  funext X
  apply Subsingleton.elim

/-- ★**非退化(上)の形** —— 非自明な自己同型が 1 つでもあれば rigid でない。

★**定義の言い換えだが、「何を見つければ rigid でないと言えるか」を明示する。** -/
theorem not_isRigidFunctor_of_ne {C₁ : Type u} [Category.{v} C₁]
    {C₂ : Type u2} [Category.{v2} C₂] {F : C₁ ⥤ C₂}
    (η : F ≅ F) (h : η ≠ Iso.refl F) : ¬ IsRigidFunctor F := fun hr => h (hr η)

/-! ### ★具体的な非退化(上)は切った —— 記録

★**`Bool` への定数関手 `Discrete PUnit ⥤ Type` に `not` による自己同型を入れる**、
という具体例を試したが、Lean で止まった:
```
type mismatch: not has type Bool → Bool
but is expected to have type constBool.obj x ⟶ constBool.obj x
```
★**`Type` の hom が `ConcreteCategory.hom` を経由して簡約しない**(mathlib の
この版では `Type` の圏構造がそう組まれている)。

★**これは分類表 #1(2 つの綴り)の新しい現れ**である ——
`Bool → Bool` と `Bool ⟶ Bool` が「同じ型の 2 つの綴り」で、
★**関手の成分として書くと後者が要求され、`not` は前者として推論される。**

★**`sorry` は置かない**(`Found/` の規律)。
★**次に当たるときの手**: `show (Bool ⟶ Bool) from Bool.not` のように
**綴りを先に決めてから書く**(表 #3)。★**今回は追わない** ——
`not_isRigidFunctor_of_ne` で「何を見つければよいか」は明示されており、
★**具体例は `Proposition 1.13` の主張には要らない。**
-/

/-! ## ★(ii) の核心 —— `𝔽_Φ` では base-identity な自己同型は恒等射

★**原文**(p.40、目視):
> By assertion (i), it follows that the automorphisms of objects of FΦ … induced by α
> are **base-identity automorphisms**. Since Φ is divisorial, hence, in particular,
> **sharp** [cf. Definition 1.1, (i), (ii)], it thus follows that all of these automorphisms
> are **trivial**, hence that α is trivial.

★★**「sharp だから自明」の中身**:
`𝔽_Φ` の射は `⟨base, div, deg⟩` の三つ組で、恒等射は `⟨𝟙, 0, 1⟩`。
base-identity な自己同型 `f` について
- `f ≫ inv f = 𝟙` の `deg` 成分から **`deg f = 1`**(`ℕ+` の消去)
- 同じ式の `div` 成分から **`(inv f).div + f.div = 0`**、すなわち `f.div` は可逆
- ★**sharp から `f.div = 0`**

したがって `f = ⟨𝟙, 0, 1⟩ = 𝟙`。

★**`divisorial` の `sharp` の部分が、ここでも単独で効いている**
(`Proposition 1.10` で 4 箇所、`1.11` で 0 箇所、ここで 1 箇所)。
-/

/-- ★★**`𝔽_Φ` では base-identity な自己同型は恒等射**。

★**`Φ` が sharp であることだけが効く。** -/
theorem elemFrob_baseIdentity_aut_eq_id {D : Type u} [Category.{v} D]
    {Φ : MonoidOn.{v, u, w} D} (hsh : ∀ A : D, IsSharp (Φ.val A))
    {A : ElemFrobCat Φ} (f : A ⟶ A) [hiso : IsIso f]
    (hb : ElemFrobCat.Hom.base f = 𝟙 A.base) : f = 𝟙 A := by
  -- `deg f = 1`
  have hcomp : f ≫ inv f = 𝟙 A := IsIso.hom_inv_id f
  have hdegs : ElemFrobCat.Hom.deg f = 1 ∧ ElemFrobCat.Hom.deg (inv f) = 1 := by
    have h := congrArg ElemFrobCat.Hom.deg hcomp
    rw [ElemFrobCat.comp_deg] at h
    have h1 : ElemFrobCat.Hom.deg (inv f) * ElemFrobCat.Hom.deg f = 1 := by simpa using h
    exact ⟨pnat_right_eq_one h1, pnat_left_eq_one h1⟩
  have hdeg := hdegs.1
  -- `div f` は可逆 ⟹ sharp から 0
  have hdiv : ElemFrobCat.Hom.div f = 0 := by
    have h := congrArg ElemFrobCat.Hom.div hcomp
    rw [ElemFrobCat.comp_div, hb, Φ.map_id, hdegs.2] at h
    refine hsh A.base _ ?_
    refine ⟨⟨ElemFrobCat.Hom.div f, ElemFrobCat.Hom.div (inv f), ?_, ?_⟩, rfl⟩
    · simpa [add_comm] using h
    · simpa using h
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · simpa using hb
  · simpa using hdiv
  · simpa using hdeg

/-! ## ★(i) —— 原文が最後の一歩を書いていない

★**原文の証明**(p.40、目視):
> Any automorphism α of the functor CA →D determines an automorphism of the
> composite functor C^pl-bk_A →CA →D. On the other hand, this composite functor
> factors as a composite C^pl-bk_A →D_{A_D} →D, where the first functor is
> [by Definition 1.3, (i), (c)] an equivalence of categories. Thus, we conclude that α
> determines an automorphism of the natural functor D_{A_D} →D, which is
> **necessarily trivial, since D is slim**.

★★**ここで原文は終わっている。しかし示したのは
「`α` が `𝒞^pl-bk_A` の上で自明」だけである。**
★**`𝒞_A` 全体で自明であることは、まだ言えていない。**

★**残りの一歩**(原文が書いていない):
`𝒞_A` の任意の対象 `(X, f : X ⟶ A)` は、`Definition 1.3, (iv), (a)` の分解で
`f = γ ≫ β ≫ α₀`(`α₀` は pull-back)と書ける。
`(X, f) ⟶ (Y, α₀)` は `𝒞_A` の射なので、★**`α` の自然性**が
```
α_{(X,f)} ≫ Base(γ ≫ β) = Base(γ ≫ β) ≫ α_{(Y,α₀)} = Base(γ ≫ β)
```
を与える。★**`Base(γ ≫ β)` は同型**(`γ` は Frobenius 型、`β` は pre-step、
どちらも base-isomorphism)なので **`α_{(X,f)} = 𝟙`**。

★★**「pull-back の部分で自明 ⟹ 全体で自明」は、自然性と分解の合わせ技である。**
★**原文はこの一歩を書いていない** —— 我々が数えてきた
「原文が書かない一歩」の、この項目での 1 例目。
-/

/-- ★**(i) の残りの一歩** —— 自然変換の成分が、同型に沿って自明性を伝える。

★`h : a ≫ u = u` で `u` が同型なら `a = 𝟙`。
★**自然性が与える式がこの形になる**(上の docstring を見よ)。 -/
theorem eq_id_of_comp_eq_of_isIso {E : Type u} [Category.{v} E] {X Y : E}
    (a : X ⟶ X) (u : X ⟶ Y) [IsIso u]
    (h : a ≫ u = u) : a = 𝟙 X := by
  have : a ≫ u = 𝟙 X ≫ u := by rw [Category.id_comp]; exact h
  exact (cancel_mono u).mp this

end ABC3.Found.FrdI
