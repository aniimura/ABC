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

universe v u w v2 u2 v3 u3

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

/-- ★★**「Since A is arbitrary」の中身** —— すべてのスライスで rigid なら全体で rigid。

★原文 (FrdI p.40) は (i)(ii) の各条を
> we thus conclude that both CA →D and C →D are rigid. This completes the proof

で締める。★**この一歩は `Over.mk (𝟙 X)` を取れば閉じる**(`X` 自身をスライスの底に使う)。
★**原文はこの選び方を書いていない**(本項目で 3 例目)。 -/
theorem isRigidFunctor_of_forall_over {C₁ : Type u} [Category.{v} C₁]
    {C₂ : Type u2} [Category.{v2} C₂] (G : C₁ ⥤ C₂)
    (h : ∀ A : C₁, IsRigidFunctor (Over.forget A ⋙ G)) : IsRigidFunctor G := by
  intro η
  apply Iso.ext
  apply NatTrans.ext
  funext X
  exact congrArg (fun t => t.hom.app (Over.mk (𝟙 X)))
    (h X (Functor.isoWhiskerLeft (Over.forget X) η))

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
    (a : X ⟶ X) (u : X ⟶ Y) (hu : IsIso u)
    (h : a ≫ u = u) : a = 𝟙 X := by
  haveI := hu
  have h2 : a ≫ u = 𝟙 X ≫ u := by rw [Category.id_comp]; exact h
  exact (cancel_mono u).mp h2

/-- ★`eq_id_of_comp_eq_of_isIso` の**鏡像** —— `u ≫ a = u` で `u` が同型なら `a = 𝟙`。

★**分類表 #5 を避けるため `hu` は明示引数**(インスタンスにしない)。 -/
theorem eq_id_of_eq_comp_of_isIso {E : Type u} [Category.{v} E] {X Y : E}
    (a : Y ⟶ Y) (u : X ⟶ Y) (hu : IsIso u)
    (h : u ≫ a = u) : a = 𝟙 Y := by
  haveI := hu
  have h2 : u ≫ a = u ≫ 𝟙 Y := by rw [Category.comp_id]; exact h
  exact (cancel_epi u).mp h2

/-! ### ★★原文の「the first functor is an equivalence of categories」の一歩

★原文 (FrdI p.40) は
> factors as a composite C^pl-bk_A →D_{A_D} →D, where the first functor is
> [by Definition 1.3, (i), (c)] an equivalence of categories. Thus, we conclude that α
> determines an automorphism of the natural functor D_{A_D} →D

と書く。★★**「圏同値だから自己同型が移る」を、実際に移すのがこの補題である。**

★**向きに注意**: 欲しいのは「`G` が rigid ⟹ `F ⋙ G` が rigid」(`F` は圏同値)。
原文の言い回しは逆向き(「`α` が `D_{A_D}→D` の自己同型を**定める**」)だが、
★**使うのは対偶側の含意**である —— slim から来るのは `G` の rigid 性で、
結論したいのは合成の rigid 性。

★**証明の骨**: `η : F ⋙ G ≅ F ⋙ G` を `e := F.asEquivalence` で
`G ≅ G` に共役し、`hG` で潰す。戻すとき **`e.inverse.obj b` の上でしか成分が取れない**ので、
最後に **`e.unit` に沿った自然性**で任意の対象へ運ぶ。
★**この最後の運搬も、原文には書かれていない**(本項目で 2 例目)。 -/
theorem isRigidFunctor_comp_of_equivalence {A₁ : Type u} [Category.{v} A₁]
    {B₁ : Type u2} [Category.{v2} B₁] {E₁ : Type u3} [Category.{v3} E₁]
    (e : A₁ ≌ B₁) {G : B₁ ⥤ E₁} (hG : IsRigidFunctor G) :
    IsRigidFunctor (e.functor ⋙ G) := by
  intro η
  set X : e.inverse ⋙ e.functor ⋙ G ≅ G := e.invFunIdAssoc G with hX
  -- 共役して `G ≅ G` にし、`hG` で潰す
  have hθ : (X.symm ≪≫ Functor.isoWhiskerLeft e.inverse η ≪≫ X) = Iso.refl G :=
    hG (X.symm ≪≫ Functor.isoWhiskerLeft e.inverse η ≪≫ X)
  -- `e.inverse` の像の上で `η` は自明
  have key : ∀ b : B₁,
      η.hom.app (e.inverse.obj b) = 𝟙 ((e.functor ⋙ G).obj (e.inverse.obj b)) := by
    intro b
    have h : X.inv.app b ≫ η.hom.app (e.inverse.obj b) ≫ X.hom.app b = 𝟙 (G.obj b) :=
      congrArg (fun t : G ≅ G => t.hom.app b) hθ
    -- ★**左から `X.hom.app b` を掛け戻す。`𝟙` を作らない `_assoc` 版が要点** ——
    -- ★`𝟙` を経由すると「対象の 2 つの綴り」(表 #1)で `Category.comp_id` が発火しない。
    have h2 : η.hom.app (e.inverse.obj b) ≫ X.hom.app b = X.hom.app b := by
      have h3 := congrArg (fun t => X.hom.app b ≫ t) h
      simp only [Iso.hom_inv_id_app_assoc, Category.comp_id] at h3
      -- ★`simpa` ではなく `exact` —— **既定透明度でないと 2 つの綴りが同一視されない**(表 #1)
      exact h3
    -- ★**可逆性はインスタンス探索に任せず逆射を手で書く**(表 #5)。
    have hXi : IsIso (X.hom.app b) :=
      ⟨⟨X.inv.app b, X.hom_inv_id_app b, X.inv_hom_id_app b⟩⟩
    exact eq_id_of_comp_eq_of_isIso _ (X.hom.app b) hXi h2
  -- ★`e.unit` に沿って任意の対象へ運ぶ(原文が書いていない段)
  apply Iso.ext
  apply NatTrans.ext
  funext a
  show η.hom.app a = 𝟙 ((e.functor ⋙ G).obj a)
  set u : (e.functor ⋙ G).obj (e.inverse.obj (e.functor.obj a)) ⟶ (e.functor ⋙ G).obj a :=
    (e.functor ⋙ G).map (e.unitInv.app a) with hu
  set w : (e.functor ⋙ G).obj a ⟶ (e.functor ⋙ G).obj (e.inverse.obj (e.functor.obj a)) :=
    (e.functor ⋙ G).map (e.unit.app a) with hw
  have hwu : w ≫ u = 𝟙 ((e.functor ⋙ G).obj a) := by
    rw [hw, hu, ← Functor.map_comp]
    -- ★最後に残るのは `𝟙 (G.obj (e.functor.obj a)) = 𝟙 ((e.functor ⋙ G).obj a)`、
    -- ★すなわち**同じ対象の 2 つの綴り**。`simp` は届かないので `rfl` で閉じる。
    simp
    rfl
  have hnat : u ≫ η.hom.app a = η.hom.app (e.inverse.obj (e.functor.obj a)) ≫ u :=
    η.hom.naturality (e.unitInv.app a)
  rw [key (e.functor.obj a), Category.id_comp] at hnat
  calc η.hom.app a = (w ≫ u) ≫ η.hom.app a := by rw [hwu, Category.id_comp]
    _ = w ≫ u ≫ η.hom.app a := by rw [Category.assoc]
    _ = w ≫ u := by rw [hnat]
    _ = 𝟙 _ := hwu

/-- ★`isRigidFunctor_comp_of_equivalence` の `[F.IsEquivalence]` 版。

★`F.asEquivalence.functor` は `F` に**定義的に等しい**ので、そのまま渡せる。 -/
theorem isRigidFunctor_of_isEquivalence_comp {A₁ : Type u} [Category.{v} A₁]
    {B₁ : Type u2} [Category.{v2} B₁] {E₁ : Type u3} [Category.{v3} E₁]
    (F : A₁ ⥤ B₁) [F.IsEquivalence] {G : B₁ ⥤ E₁} (hG : IsRigidFunctor G) :
    IsRigidFunctor (F ⋙ G) :=
  isRigidFunctor_comp_of_equivalence F.asEquivalence hG

section PropI

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

include P in
/-- ★★**`Proposition 1.13, (i)` の本体** ——
「pull-back の対象で自明」から「`𝒞_A` 全体で自明」へ。

★★**原文はこの段を書いていない**(上の docstring を見よ)。
仮定 `htriv` が原文の示した部分(圏同値 ＋ `𝒟` の slim 性)であり、
★**ここで実装するのはそこから先である。**

★**構成**: `𝒞_A` の対象 `X` の構造射を `F.arbFactor` で `γ ≫ β ≫ α₀` に分解し、
`X ⟶ Over.mk α₀` を作って**自然性**を当てる。
`Base (γ ≫ β)` が同型なので消せる。 -/
theorem prop_1_13_i_from_pullBack (F : FrobenioidCore P) (A : C)
    (η : (Over.forget A ⋙ P.proj) ≅ (Over.forget A ⋙ P.proj))
    (htriv : ∀ Y : Over A, IsPullBack P Y.hom →
      η.hom.app Y = 𝟙 ((Over.forget A ⋙ P.proj).obj Y)) :
    η = Iso.refl _ := by
  apply Iso.ext
  apply NatTrans.ext
  funext X
  obtain ⟨Y₀, Z₀, γ, β, α₀, hfac, hγ, hβ, hα₀⟩ := F.arbFactor X.hom
  -- `X ⟶ Over.mk α₀`
  set Y : Over A := Over.mk α₀ with hY
  have hw : (γ ≫ β) ≫ α₀ = X.hom := by rw [Category.assoc, ← hfac]
  set f : X ⟶ Y := Over.homMk (γ ≫ β) hw with hf
  -- 自然性
  have hnat := η.hom.naturality f
  rw [htriv Y hα₀, Category.comp_id] at hnat
  -- `Base (γ ≫ β)` は同型
  have hbi : IsBaseIsomorphism P (γ ≫ β) := isBaseIsomorphism_comp P hγ.2 hβ.2
  have hmap : (Over.forget A ⋙ P.proj).map f = P.Base (γ ≫ β) := rfl
  rw [hmap] at hnat
  show η.hom.app X = 𝟙 _
  exact eq_id_of_comp_eq_of_isIso (η.hom.app X) (P.Base (γ ≫ β)) hbi hnat.symm

/-! ### ★原文の合成 `C^pl-bk_A → C_A → D` を Lean に置く -/

/-- ★**`(𝒞^pl-bk)_A → 𝒞_A` の包含関手** —— 原文の合成の最初の矢。 -/
def plBkToOver (A : C) : Over (⟨A⟩ : PlBk P) ⥤ Over A :=
  Over.post (X := (⟨A⟩ : PlBk P)) (wideSubcategoryInclusion (pullBackProp P))

/-- ★★**原文の「factors as a composite」は Lean では `rfl` である。**

原文 (FrdI p.40):
> the other hand, this composite functor factors as a composite Cpl-bk

★**2 通りの合成が定義的に等しい** —— `P.proj = P.toElem ⋙ ElemFrobCat.proj` と
`P.Base φ = (P.toElem.map φ).base` が同じものを指しているため。 -/
theorem plBkToOver_comp (A : C) :
    plBkToOver P A ⋙ Over.forget A ⋙ P.proj
      = plBkOverFunctor P A ⋙ Over.forget ((P.toElem.obj A).base) := rfl

include P in
/-- ★★**原文が示した部分** —— `𝒟` が slim なら `α` は pull-back 対象の上で自明。

★**圏同値**(`Definition 1.3, (i), (c)`)と **slim 性**の合わせ技。
`isRigidFunctor_comp_of_equivalence` がその「合わせ」を担う。 -/
theorem prop_1_13_i_pullBack_trivial (F : FrobenioidCore P) (hslim : IsSlimCat D) (A : C)
    (η : (Over.forget A ⋙ P.proj) ≅ (Over.forget A ⋙ P.proj))
    (Y : Over A) (hpb : IsPullBack P Y.hom) :
    η.hom.app Y = 𝟙 ((Over.forget A ⋙ P.proj).obj Y) := by
  haveI := F.plBkEquiv A
  -- `α` を `(𝒞^pl-bk)_A` に制限する
  have hrig : IsRigidFunctor (plBkOverFunctor P A ⋙ Over.forget ((P.toElem.obj A).base)) :=
    isRigidFunctor_of_isEquivalence_comp (plBkOverFunctor P A) (hslim _)
  have hη : (Functor.isoWhiskerLeft (plBkToOver P A) η :
      (plBkOverFunctor P A ⋙ Over.forget ((P.toElem.obj A).base)) ≅ _) = Iso.refl _ :=
    hrig _
  -- `Y` を `(𝒞^pl-bk)_A` の対象として持ち上げる(`Y.hom` が pull-back だから可能)
  exact congrArg (fun t => t.hom.app (Over.mk (⟨Y.hom, hpb⟩ : (⟨Y.left⟩ : PlBk P) ⟶ ⟨A⟩))) hη

include P in
/-- ★★**`Proposition 1.13, (i)`** —— `𝒟` が slim なら `𝒞_A → 𝒟` は rigid。

原文 (FrdI p.40):
> (i) The composite CA →D of the natural functor CA →C with the natural
-/
theorem prop_1_13_i (F : FrobenioidCore P) (hslim : IsSlimCat D) (A : C) :
    IsRigidFunctor (Over.forget A ⋙ P.proj) := fun η =>
  prop_1_13_i_from_pullBack P F A η
    (fun Y hpb => prop_1_13_i_pullBack_trivial P F hslim A η Y hpb)

include P in
/-- ★**`Proposition 1.13, (i)` の「In particular」**。

原文 (FrdI p.40):
> projection functor C →D is rigid [cf. §0]. In particular, the functor C →D is
-/
theorem prop_1_13_i_global (F : FrobenioidCore P) (hslim : IsSlimCat D) :
    IsRigidFunctor P.proj :=
  isRigidFunctor_of_forall_over P.proj fun A => prop_1_13_i P F hslim A

/-! ## ★(ii) —— `𝒞_A → 𝔽_Φ`

★**組み立ては 2 段**:
1. `ElemFrobCat.proj` と合成すれば (i) が使えて、**成分は base-identity**。
2. `elemFrob_baseIdentity_aut_eq_id`(sharp から)で**恒等射**。

★**原文の順序どおりだが、1 段目の「合成すれば (i)」は原文が
「By assertion (i), it follows that …」と一言で済ませている部分**である。 -/

include P in
/-- ★★**`Proposition 1.13, (ii)`** —— `𝒞_A → 𝔽_Φ` は rigid。

原文 (FrdI p.40):
> (ii) The composite CA →FΦ of the natural functor CA →C with the functor
-/
theorem prop_1_13_ii (F : FrobenioidCore P) (hslim : IsSlimCat D) (A : C) :
    IsRigidFunctor (Over.forget A ⋙ P.toElem) := by
  intro η
  apply Iso.ext
  apply NatTrans.ext
  funext Y
  -- 1 段目: `ElemFrobCat.proj` と合成して (i) を当てる ⟹ base-identity
  have hbase : (Functor.isoWhiskerRight η ElemFrobCat.proj :
      (Over.forget A ⋙ P.proj) ≅ (Over.forget A ⋙ P.proj)) = Iso.refl _ :=
    prop_1_13_i P F hslim A _
  have hb : ElemFrobCat.Hom.base (η.hom.app Y) = 𝟙 (P.toElem.obj Y.left).base :=
    congrArg (fun t => t.hom.app Y) hbase
  -- 2 段目: sharp から恒等射
  haveI : IsIso (η.hom.app Y) := ⟨⟨η.inv.app Y, η.hom_inv_id_app Y, η.inv_hom_id_app Y⟩⟩
  show η.hom.app Y = 𝟙 _
  exact elemFrob_baseIdentity_aut_eq_id (fun X => (P.divisorial X).2) _ hb

include P in
/-- ★**`Proposition 1.13, (ii)` の「In particular」**。

原文 (FrdI p.40):
> C →FΦ is rigid. In particular, the functor C →FΦ is rigid.
-/
theorem prop_1_13_ii_global (F : FrobenioidCore P) (hslim : IsSlimCat D) :
    IsRigidFunctor P.toElem :=
  isRigidFunctor_of_forall_over P.toElem fun A => prop_1_13_ii P F hslim A

/-! ## ★(iii) —— `𝒞` が slim

原文 (FrdI p.40):
> (iii) Suppose, moreover, that every object A ∈Ob(C) satisfies [at least] one of

★**原文の証明は 3 段**:
1. `α` の成分は base-identity 自己同型、すなわち `𝒪^×(−)` に属する。
2. **isometric pre-step についての自然性**から `𝒪^×(−)^{imtr-pre}` に属する。
3. 条件 (a) または (b) から自明。

★★**1 段目は原文より短く済む** —— 原文は (i) を使って base-identity を出すが、
★**(ii) を使えば `𝔽_Φ` の像が `𝟙` になり、base-identity・`Div = 0`・次数 1 が
一度に出る**。`𝒪^×` の membership に必要なものがそのまま揃う。 -/

include P in
/-- ★★**(iii) の第1段** —— `α` の成分は `𝒪^×` に属する。

★**(ii) から直ちに出る**(原文は (i) を経由するが、(ii) のほうが強い)。 -/
theorem prop_1_13_iii_mem_otimes (F : FrobenioidCore P) (hslim : IsSlimCat D) (A : C)
    (α : (Over.forget A : Over A ⥤ C) ≅ Over.forget A) (X : Over A) :
    (α.hom.app X : End X.left) ∈ OTimes P X.left := by
  have h := prop_1_13_ii P F hslim A (Functor.isoWhiskerRight α P.toElem)
  have hc : P.toElem.map (α.hom.app X) = 𝟙 (P.toElem.obj X.left) :=
    congrArg (fun t => t.hom.app X) h
  refine ⟨⟨?_, ?_⟩, ?_⟩
  · show P.Base (α.hom.app X) = P.Base (𝟙 X.left)
    show (P.toElem.map (α.hom.app X)).base = _
    rw [hc]; simp
  · show P.degFr (α.hom.app X) = 1
    show (P.toElem.map (α.hom.app X)).deg = 1
    rw [hc]; rfl
  · exact (CategoryTheory.isUnit_iff_isIso (α.hom.app X : End X.left)).mpr inferInstance

include P in
/-- ★★**(iii) の第2段** —— `α` の成分は `𝒪^×(−)^{imtr-pre}` に属する。

★★**自然性がそのまま「同型を除いて可換」を与える。**
`γ : Cc ⟶ X.left` が isometric pre-step なら、`Cc` は `γ ≫ X.hom` によって
`𝒞_A` の対象になり、`Over.homMk γ` が射になる。`α` の自然性は
```
γ ≫ α_X = α_{(Cc, γ ≫ X.hom)} ≫ γ
```
を与え、右の `α_{...}` は同型である。★**これが
`pushFunctorIsoOfCommutes` の仮定そのもの。** -/
theorem prop_1_13_iii_mem_imtrPre (F : FrobenioidCore P) (hslim : IsSlimCat D) (A : C)
    (α : (Over.forget A : Over A ⥤ C) ≅ Over.forget A) (X : Over A) :
    (α.hom.app X : End X.left) ∈ OTimesImtrPre P F X.left := by
  refine ⟨prop_1_13_iii_mem_otimes P F hslim A α X,
    isBaseIsomorphism_of_isIso P (α.hom.app X), ⟨?_⟩⟩
  refine pushFunctorIsoOfCommutes P F (α.hom.app X)
    (isBaseIsomorphism_of_isIso P (α.hom.app X)) (fun Z => ?_)
  -- `Z.hom.hom : Z.left.obj ⟶ X.left` は isometric pre-step
  refine ⟨α.hom.app (Over.mk (Z.hom.hom ≫ X.hom)),
    ⟨⟨⟨α.inv.app (Over.mk (Z.hom.hom ≫ X.hom)),
      α.hom_inv_id_app (Over.mk (Z.hom.hom ≫ X.hom)),
      α.inv_hom_id_app (Over.mk (Z.hom.hom ≫ X.hom))⟩⟩⟩, ?_⟩
  exact α.hom.naturality (Over.homMk Z.hom.hom rfl : Over.mk (Z.hom.hom ≫ X.hom) ⟶ X)

include P in
/-- ★★**`Proposition 1.13, (iii)` の条件 (a) の場合** —— `𝒞` は slim。

原文 (FrdI p.40):
> the following two conditions: (a) O×(A)imtr-pre = {1} [cf. Proposition 1.9, (ii)];

★第 1 段と第 2 段を繋ぐだけで閉じる。 -/
theorem prop_1_13_iii_a (F : FrobenioidCore P) (hslim : IsSlimCat D)
    (ha : ∀ A : C, OTimesImtrPre P F A = {1}) : IsSlimCat C := by
  intro A α
  apply Iso.ext
  apply NatTrans.ext
  funext X
  have h := prop_1_13_iii_mem_imtrPre P F hslim A α X
  rw [ha X.left] at h
  exact h

include P in
/-- ★★**(iii)(b) の核心** —— Frobenius-normalized な対象では、
`α` の成分は**任意の次数 `n` について `n` 乗**である。

原文 (FrdI p.40):
> functoriality of the automorphisms induced by α with respect to base-identity en-

★★**自然性と Frobenius-normalized の式を合わせるだけで出る**:
- 自然性: `ζ ≫ α_X = α_Z ≫ ζ`（`Z` は `ζ ≫ X.hom` を構造射に持つ対象）
- normalized: `ζ ≫ (α_Z)^n = α_Z ≫ ζ`
- ★**`𝒞` は totally epimorphic なので `ζ` は epi**、左から消せる。 -/
theorem prop_1_13_iii_pow (F : FrobenioidCore P) (hslim : IsSlimCat D) (A : C)
    (α : (Over.forget A : Over A ⥤ C) ≅ Over.forget A) (X : Over A)
    (hfn : IsFrobeniusNormalized P X.left)
    (ζ : End X.left) (hζb : IsBaseIdentity P (ζ : X.left ⟶ X.left)) (n : ℕ+)
    (hζd : P.degFr (ζ : X.left ⟶ X.left) = n) :
    ∃ β : End X.left, β ∈ OTimes P X.left ∧
      ((α.hom.app X : End X.left) : X.left ⟶ X.left)
        = ((β ^ (n : ℕ) : End X.left) : X.left ⟶ X.left) := by
  set Z : Over A := Over.mk ((ζ : X.left ⟶ X.left) ≫ X.hom) with hZ
  have hnat := α.hom.naturality
    (Over.homMk (ζ : X.left ⟶ X.left) rfl : Z ⟶ X)
  -- 自然性: `ζ ≫ α_X = α_Z ≫ ζ`
  have hn : (ζ : X.left ⟶ X.left) ≫ α.hom.app X
      = α.hom.app Z ≫ (ζ : X.left ⟶ X.left) := hnat
  have hmemZ : (α.hom.app Z : End X.left) ∈ OTimes P X.left :=
    prop_1_13_iii_mem_otimes P F hslim A α Z
  -- normalized: `ζ ≫ (α_Z)^{degFr ζ} = α_Z ≫ ζ`
  have hfn' := hfn ζ hζb (α.hom.app Z : End X.left) (OTimes_le_OTri P X.left hmemZ)
  rw [hζd] at hfn'
  refine ⟨(α.hom.app Z : End X.left), hmemZ, ?_⟩
  haveI : Epi (ζ : X.left ⟶ X.left) := P.totEpiC _ _ _
  refine (cancel_epi (ζ : X.left ⟶ X.left)).mp ?_
  rw [hn]
  exact hfn'.symm

/-- ★`⋂_{n∈ℕ≥1} {𝒪^×(A)}^n` —— 原文 (b) の第1条件が `{1}` だと言っている集合。 -/
def OTimesPowInter (A : C) : Set (End A) :=
  ⋂ n : ℕ+, (fun x : End A => x ^ (n : ℕ)) '' (OTimes P A : Set (End A))

include P in
/-- ★★**quasi-Frobenius-trivial かつ Frobenius-normalized な対象では、
`α` の成分は `⋂_n {𝒪^×}^n` に属する**。

★`IsQuasiFrobeniusTrivial` は「各次数の base-identity 自己射がある」であり、
★**`prop_1_13_iii_pow` が要求するものそのもの**である。 -/
theorem prop_1_13_iii_mem_powInter (F : FrobenioidCore P) (hslim : IsSlimCat D) (A : C)
    (α : (Over.forget A : Over A ⥤ C) ≅ Over.forget A) (X : Over A)
    (hq : IsQuasiFrobeniusTrivial P X.left) (hfn : IsFrobeniusNormalized P X.left) :
    (α.hom.app X : End X.left) ∈ OTimesPowInter P X.left := by
  refine Set.mem_iInter.mpr (fun n => ?_)
  obtain ⟨φ, hφb, hφd⟩ := hq n
  obtain ⟨β, hβ, hEq⟩ := prop_1_13_iii_pow P F hslim A α X hfn φ hφb n hφd
  exact ⟨β, hβ, hEq.symm⟩

include P in
/-- ★★**`Proposition 1.13, (iii)` の条件 (b) の場合**。

原文 (FrdI p.40):
> such that B is quasi-Frobenius-trivial and Frobenius-normalized. Then the

★★**原文より強い仮定であることを明記する**。原文は
「co-angular pre-step `B → A` で `B` が quasi-Frobenius-trivial かつ
Frobenius-normalized なものがある」という形で、`Definition 1.3, (iii), (c)` の
全単射 `𝒪^×(B) ≅ 𝒪^×(A)` で `A` へ送る。
★**ここでは `A` 自身がその性質を持つ場合を実装している**。 -/
theorem prop_1_13_iii_b (F : FrobenioidCore P) (hslim : IsSlimCat D)
    (hq : ∀ A : C, IsQuasiFrobeniusTrivial P A)
    (hfn : ∀ A : C, IsFrobeniusNormalized P A)
    (hb : ∀ A : C, OTimesPowInter P A = {1}) : IsSlimCat C := by
  intro A α
  apply Iso.ext
  apply NatTrans.ext
  funext X
  have h := prop_1_13_iii_mem_powInter P F hslim A α X (hq X.left) (hfn X.left)
  rw [hb X.left] at h
  exact h

/-- ★**共役はべきに伸びる** —— `ψ ≫ γ = β ≫ ψ` なら `ψ ≫ γ^n = β^n ≫ ψ`。

★**`End` の積は `x * y = y ≫ x`** なので `γ^(n+1) = γ^n * γ = γ ≫ γ^n` である。 -/
theorem comp_pow_of_comp {A B : C} (ψ : B ⟶ A) {β : End B} {γ : End A}
    (h : ψ ≫ (γ : A ⟶ A) = (β : B ⟶ B) ≫ ψ) (n : ℕ) :
    ψ ≫ ((γ ^ n : End A) : A ⟶ A) = ((β ^ n : End B) : B ⟶ B) ≫ ψ := by
  induction n with
  | zero => simp
  | succ n ih =>
    have hg : ((γ ^ (n + 1) : End A) : A ⟶ A)
        = (γ : A ⟶ A) ≫ ((γ ^ n : End A) : A ⟶ A) := by rw [pow_succ]; rfl
    have hb : ((β ^ (n + 1) : End B) : B ⟶ B)
        = (β : B ⟶ B) ≫ ((β ^ n : End B) : B ⟶ B) := by rw [pow_succ]; rfl
    rw [hg, hb, ← Category.assoc, h, Category.assoc, ih, ← Category.assoc]

include P in
/-- ★**単元の逆元も `𝒪^▷` に入る** —— `Base` と `degFr` の関手性から。 -/
theorem otri_inv_mem {B : C} (β β' : End B) (hβ : β ∈ OTimes P B)
    (hcomp : (β : B ⟶ B) ≫ (β' : B ⟶ B) = 𝟙 B) : β' ∈ OTri P B := by
  constructor
  · show P.Base (β' : B ⟶ B) = P.Base (𝟙 B)
    have hb := congrArg P.Base hcomp
    rw [P.Base_comp, show P.Base (β : B ⟶ B) = P.Base (𝟙 B) from hβ.1.1,
      P.Base_id, Category.id_comp] at hb
    simpa using hb
  · show P.degFr (β' : B ⟶ B) = 1
    have hd := congrArg P.degFr hcomp
    rw [P.degFr_comp, show P.degFr (β : B ⟶ B) = 1 from hβ.1.2, P.degFr_id, mul_one] at hd
    exact hd

include P in
/-- ★★**(iii)(b) の移送** —— co-angular pre-step `ψ : B ⟶ X.left` に沿って
`B` 側の情報を `X.left` へ送る。

原文 (FrdI p.40):
> Frobenius-normalized objects — hence also [cf. Definition 1.3, (iii), (c)] objects as in (b)

★★**`Definition 1.3, (iii), (c)` の全単射がここで使われる。**

★**Lean の手**(2026-08-16 の記録どおり): `Over.mk (ψ ≫ X.hom)` を `set` で置くと
`End (Over.mk …).left` と `End B` が**同じ型の 2 つの綴り**になり、
インスタンス探索が通らない(表 #1)。
★**`have` の型注釈で一度だけ定義的等しさを使い、以後は `B` の綴りに固定する。** -/
theorem prop_1_13_iii_transfer (F : FrobenioidCore P) (hslim : IsSlimCat D)
    (A₀ : C) (α : (Over.forget A₀ : Over A₀ ⥤ C) ≅ Over.forget A₀) (X : Over A₀)
    {B : C} (ψ : B ⟶ X.left) (hψc : IsCoAngular P ψ) (hψs : IsPreStep P ψ)
    (hq : IsQuasiFrobeniusTrivial P B) (hfn : IsFrobeniusNormalized P B) :
    (α.hom.app X : End X.left) ∈ OTimesPowInter P X.left := by
  refine Set.mem_iInter.mpr (fun n => ?_)
  haveI : Epi ψ := P.totEpiC _ _ _
  obtain ⟨φ, hφb, hφd⟩ := hq n
  -- ★型注釈で `B` の綴りに固定する
  have hnat : ψ ≫ (α.hom.app X : X.left ⟶ X.left)
      = ((α.hom.app (Over.mk (ψ ≫ X.hom)) : End B) : B ⟶ B) ≫ ψ :=
    α.hom.naturality (Over.homMk ψ rfl : Over.mk (ψ ≫ X.hom) ⟶ X)
  have key : ∃ β : End B, β ∈ OTimes P B ∧
      ((α.hom.app (Over.mk (ψ ≫ X.hom)) : End B) : B ⟶ B)
        = ((β ^ (n : ℕ) : End B) : B ⟶ B) :=
    prop_1_13_iii_pow P F hslim A₀ α (Over.mk (ψ ≫ X.hom)) hfn φ hφb n hφd
  obtain ⟨β₀, hβ₀, hEq⟩ := key
  obtain ⟨u, hu⟩ := hβ₀.2
  have hiu : (β₀ : B ⟶ B) ≫ ((↑u⁻¹ : End B) : B ⟶ B) = 𝟙 B := by
    rw [← hu]; exact u.inv_val
  have hui : ((↑u⁻¹ : End B) : B ⟶ B) ≫ (β₀ : B ⟶ B) = 𝟙 B := by
    rw [← hu]; exact u.val_inv
  obtain ⟨γ, ⟨hγo, hγ⟩, -⟩ := F.otriFwd ψ hψc hψs β₀ (OTimes_le_OTri P B hβ₀)
  obtain ⟨γ', ⟨hγo', hγ'⟩, -⟩ := F.otriFwd ψ hψc hψs (↑u⁻¹ : End B)
    (otri_inv_mem P β₀ (↑u⁻¹ : End B) hβ₀ hiu)
  have hgg' : (γ : X.left ⟶ X.left) ≫ (γ' : X.left ⟶ X.left) = 𝟙 X.left := by
    refine (cancel_epi ψ).mp ?_
    rw [← Category.assoc, hγ, Category.assoc, hγ', ← Category.assoc, hiu,
      Category.id_comp, Category.comp_id]
  have hg'g : (γ' : X.left ⟶ X.left) ≫ (γ : X.left ⟶ X.left) = 𝟙 X.left := by
    refine (cancel_epi ψ).mp ?_
    rw [← Category.assoc, hγ', Category.assoc, hγ, ← Category.assoc, hui,
      Category.id_comp, Category.comp_id]
  refine ⟨γ, ⟨hγo, (CategoryTheory.isUnit_iff_isIso γ).mpr
    ⟨⟨(γ' : X.left ⟶ X.left), hgg', hg'g⟩⟩⟩, ?_⟩
  refine (cancel_epi ψ).mp ?_
  have e1 : ψ ≫ ((γ ^ (n : ℕ) : End X.left) : X.left ⟶ X.left)
      = ((β₀ ^ (n : ℕ) : End B) : B ⟶ B) ≫ ψ := comp_pow_of_comp ψ hγ (n : ℕ)
  rw [e1, ← hEq]
  exact hnat.symm

include P in
/-- ★★★**`Proposition 1.13, (iii)`** —— 各対象が (a) または (b) を満たせば `𝒞` は slim。

原文 (FrdI p.40):
> (iii) Suppose, moreover, that every object A ∈Ob(C) satisfies [at least] one of

★**原文の「at least one of」を `∨` でそのまま書く。** -/
theorem prop_1_13_iii (F : FrobenioidCore P) (hslim : IsSlimCat D)
    (h : ∀ A : C,
      OTimesImtrPre P F A = {1} ∨
      (OTimesPowInter P A = {1} ∧
        ∃ (B : C) (ψ : B ⟶ A), IsCoAngular P ψ ∧ IsPreStep P ψ ∧
          IsQuasiFrobeniusTrivial P B ∧ IsFrobeniusNormalized P B)) :
    IsSlimCat C := by
  intro A α
  apply Iso.ext
  apply NatTrans.ext
  funext X
  rcases h X.left with ha | ⟨hb, B, ψ, hψc, hψs, hq, hfn⟩
  · have hm := prop_1_13_iii_mem_imtrPre P F hslim A α X
    rw [ha] at hm
    exact hm
  · have hm := prop_1_13_iii_transfer P F hslim A α X ψ hψc hψs hq hfn
    rw [hb] at hm
    exact hm

/-! ## ★★命題全体の `.src`（2026-08-16）

★**原文の主張は 3 ではなく 5 であった** —— (i)(ii) にそれぞれ
「In particular」の第2主張がある。ファイル冒頭の見積もり表は誤っていた。

| 条 | 主張 | 宣言 |
|---|---|---|
| (i) | `𝒞_A → 𝒟` が rigid | `prop_1_13_i` |
| (i) | **In particular**: `𝒞 → 𝒟` が rigid | `prop_1_13_i_global` |
| (ii) | `𝒞_A → 𝔽_Φ` が rigid | `prop_1_13_ii` |
| (ii) | **In particular**: `𝒞 → 𝔽_Φ` が rigid | `prop_1_13_ii_global` |
| (iii) | (a) または (b) の下で `𝒞` が slim | `prop_1_13_iii` |

★**原文が一言で済ませていた段をすべて実体にした**:
- 「the first functor is an equivalence of categories」→ `isRigidFunctor_comp_of_equivalence`
- 「Since A is arbitrary」→ `isRigidFunctor_of_forall_over`
- 「By assertion (i), it follows that …」→ (ii) の 2 段
- 「the functoriality … with respect to isometric pre-steps」→ `prop_1_13_iii_mem_imtrPre`
- 「the functoriality … with respect to base-identity endomorphisms」→ `prop_1_13_iii_pow`
- 「hence also [cf. Definition 1.3, (iii), (c)] objects as in (b)」→ `prop_1_13_iii_transfer`

★★**この項目で繰り返し現れた骨格**: 原文が一言で済ませる段は、
たいてい**自然性が式を与え、`epi` / `mono` / `hom の subsingleton` のいずれかが
消去を与える**という形であった。 -/

def prop_1_13.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 39, item := "Proposition 1.13",
    sectionId := "frdi-prop-1-13-i" }

end PropI

end ABC3.Found.FrdI
