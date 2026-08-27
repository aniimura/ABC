/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AHeightWitness
import ABC3.Found.Arakelov.ADegBase

/-!
# [GenEll] Definition 1.2, (i) —— **高さ関数 `ht_M̄ : X(ℚ̄) → ℝ`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★何が足りていなかったか

`Found/Arakelov/AHeightWitness.lean` は `htOf X L x = deg_F(x_F^* L̄)` を作っていたが、
そこの `x : NFPoint X` は**数体 `F` と射 `x_F` の対**であって、`X(ℚ̄)` の点ではない。
★同じ幾何的な点が別の `(F, x_F)` で 2 通りに表せるので、
**「`X(ℚ̄)` 上の関数」と言うには代表元の取り替えで値が変わらないことが要る**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★「any morphism」——**これが well-defined 性の要求そのもの**である。

## ★★★埋めた段

| 段 | 宣言 |
|---|---|
| 引き戻しの関手性（アルキメデス側） | ★`arcCMPullback_comp`（本ファイル） |
| 引き戻しの関手性（`APic` 全体） | ★`aPicPullback_comp`（本ファイル） |
| 高さの底変換不変性 | ★★`htOf_baseChange`（本ファイル） |
| 商としての `X(ℚ̄)` | ★★`NFAlgPoint`（本ファイル） |
| ★**高さ関数そのもの** | ★★★`ht`（本ファイル） |

★次数側の `degFOf_baseChange`（`ADegBase.lean`、`[K:ℚ] = [F:ℚ]·[K:F]` が
正規化の分母をちょうど約分する）が既にあったので、残りは関手性だけだった。

## ★★★★★ついでに `Proposition 1.4, (i)` の但し書きが外れた

`Found/GenEll/HeightAdditive.lean` の `htArith_tensor_unconditional` は
**因子表示**なので「`x` が `D` を通らない」を仮定していた
（原文は可逆層で `X(ℚ̄)` 全体）。
★★**線束表示では `ht` が `X(ℚ̄)` 全体で定義されている**ので、
`ht_mul` は原文どおり**全体で**成り立つ。

## ★★★★★★G3 負の対照 —— 何が降りないか

`finrank_not_descends` —— **体の次数 `[F:ℚ]` は点の関数にならない**。
★底変換で関係づけられた 2 つの代表元は次数が違う。

★★★**これが `Definition 1.5` が「最小定義体」を要る理由と同じ現象**である
（`LogDiffPoint.lean` の `logDiffAt_le_baseChange`：`log-diff` も底変換で増える。
★`ht` だけが不変で、`log-diff` は不変でない——だから最小定義体で測る）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField ABC3.Found.Arakelov ABC3.Interface.Arakelov

/-! ## ★引き戻しの関手性 -/

/-- ★**アルキメデス側の引き戻しは関手的**。 -/
theorem arcCMPullback_comp {X Y Z : Scheme.{0}} (g : X ⟶ Y) (f : Y ⟶ Z) (c : arcCM Z) :
    arcCMPullback (g ≫ f) c = arcCMPullback g (arcCMPullback f c) :=
  @ContinuousMap.ext _ ℝ (arcTopology X) _ _ _ (fun p => by
    simp only [arcCMPullback, arcCM.mk, arcCM.fn, ContinuousMap.coe_mk, Category.assoc])

/-- ★★**算術直線束の引き戻しは関手的**——可逆層側は `picardDataWitness.pullback_comp`、
アルキメデス側は `arcCMPullback_comp`。 -/
theorem aPicPullback_comp {X Y Z : Scheme.{0}} (g : X ⟶ Y) (f : Y ⟶ Z)
    (L : picardDataWitness.Pic Z × Multiplicative (arcCM Z)) :
    aPicDataWitness.pullback (g ≫ f) L
      = aPicDataWitness.pullback g (aPicDataWitness.pullback f L) :=
  Prod.ext (picardDataWitness.pullback_comp g f L.1) (arcCMPullback_comp g f L.2)

/-! ## ★★★★★★高さは代表元の取り替えで変わらない -/

/-- ★★★★★★**高さの底変換不変性**——`Definition 1.2, (i)` の実質。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★★★「any morphism」を成り立たせているのは
`degFOf_baseChange`（`[K:ℚ] = [F:ℚ]·[K:F]` が正規化の分母を約分する）と
`aPicPullback_comp`（引き戻しの関手性）の 2 本である。 -/
theorem htOf_baseChange {X : Scheme.{0}} (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] [IsScalarTower (𝓞 F) (𝓞 K) K]
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X)
    (L : picardDataWitness.Pic X × Multiplicative (arcCM X)) :
    htOf X L ⟨K, Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF⟩
      = htOf X L ⟨F, xF⟩ := by
  show degFOf K (aPicDataWitness.pullback
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) L)
    = degFOf F (aPicDataWitness.pullback xF L)
  rw [aPicPullback_comp]
  exact degFOf_baseChange F K (aPicDataWitness.pullback xF L)

/-! ## ★★★★★★★`X(ℚ̄)` と高さ関数 -/

/-- ★**底変換の関係**——`(K, x_F の底変換)` は `(F, x_F)` と同じ `X(ℚ̄)` の点である。 -/
inductive NFBaseChangeRel (X : Scheme.{0}) : NFPoint X → NFPoint X → Prop
  | up (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
      [Algebra F K] [IsScalarTower (𝓞 F) (𝓞 K) K]
      (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
      NFBaseChangeRel X ⟨F, xF⟩
        ⟨K, Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF⟩

/-- ★★★★**`X(ℚ̄)`** —— 数体つきの射を底変換で同一視したもの。

★★宇宙は `Type 1`（`NFPoint` が数体 `F : Type` にわたるため）。 -/
def NFAlgPoint (X : Scheme.{0}) : Type 1 := Quot (NFBaseChangeRel X)

/-- ★代表元から `X(ℚ̄)` の点を作る。 -/
def NFAlgPoint.mk {X : Scheme.{0}} (x : NFPoint X) : NFAlgPoint X := Quot.mk _ x

/-- ★★★★★★★**[GenEll] Definition 1.2, (i)** —— 高さ関数 `ht_M̄ : X(ℚ̄) → ℝ`。

原文 (GenEll p.5):
> (i) We shall refer to the function

★★★**`htOf_baseChange` があるから商へ降りる**——それが本定義の内容である。 -/
noncomputable def ht (X : Scheme.{0})
    (L : picardDataWitness.Pic X × Multiplicative (arcCM X)) : NFAlgPoint X → ℝ :=
  Quot.lift (htOf X L) (fun a b hab => by
    cases hab with
    | up F K xF => exact (htOf_baseChange F K xF L).symm)

@[simp] theorem ht_mk {X : Scheme.{0}}
    (L : picardDataWitness.Pic X × Multiplicative (arcCM X)) (x : NFPoint X) :
    ht X L (NFAlgPoint.mk x) = htOf X L x := rfl

/-- ★★**G2 非空虚**——`X` に ℚ-点があれば `X(ℚ̄)` は空でない。 -/
theorem nfAlgPoint_nonempty {X : Scheme.{0}}
    (xQ : Spec (CommRingCat.of (𝓞 ℚ)) ⟶ X) : Nonempty (NFAlgPoint X) :=
  ⟨NFAlgPoint.mk ⟨ℚ, xQ⟩⟩

/-! ## ★★★★★★★`Proposition 1.4, (i)` —— 但し書きなし -/

attribute [local instance] aPicGroup

/-- ★★★**代表元の水準での加法性**。 -/
theorem htOf_mul (X : Scheme.{0})
    (L M : picardDataWitness.Pic X × Multiplicative (arcCM X)) (x : NFPoint X) :
    htOf X (L * M) x = htOf X L x + htOf X M x :=
  (congrArg (degFOf x.F) (aPicDataWitness.pullback_mul x.xF L M)).trans
    (degFOf_mul x.F _ _)

/-- ★★★★★★★**[GenEll] Proposition 1.4, (i)** —— `X(ℚ̄)` **全体**での加法性。

原文 (GenEll p.6):
> htL⊗M(x) = htL(x) + htM(x)

★★★`Found/GenEll/HeightAdditive.lean` の `htArith_tensor_unconditional` は
**因子表示**なので「`x` が `D` を通らない」を仮定していた。
★★**線束表示では `ht` が `X(ℚ̄)` 全体で定義されている**ので但し書きが要らない
——原文の形そのものである。 -/
theorem ht_mul (X : Scheme.{0})
    (L M : picardDataWitness.Pic X × Multiplicative (arcCM X)) (x : NFAlgPoint X) :
    ht X (L * M) x = ht X L x + ht X M x := by
  induction x using Quot.inductionOn with
  | _ y => exact htOf_mul X L M y

/-! ## ★★★★★★G3 負の対照 -/

/-- ★★★★**負の対照** —— **体の次数 `[F:ℚ]` は点の関数にならない**。

底変換で関係づけられた 2 つの代表元は、同じ `X(ℚ̄)` の点を表しながら次数が違う。
★したがって `Quot.lift` に `Module.finrank ℚ (·).F` を流すことは**できない**
——`ht` が流せるのは `htOf_baseChange` があるからであって、自明ではない。

★★★**これが `Definition 1.5` が「最小定義体」を要る理由と同じ現象**である。
`LogDiffPoint.lean` の `logDiffAt_le_baseChange` が示すとおり `log-diff` も
底変換で増える。★`ht` だけが不変で、`log-diff` は不変でない。 -/
theorem finrank_not_descends {X : Scheme.{0}} (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] [IsScalarTower (𝓞 F) (𝓞 K) K]
    (hne : Module.finrank F K ≠ 1)
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    NFBaseChangeRel X ⟨F, xF⟩
        ⟨K, Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF⟩
      ∧ Module.finrank ℚ F ≠ Module.finrank ℚ K := by
  refine ⟨NFBaseChangeRel.up F K xF, ?_⟩
  haveI : FiniteDimensional ℚ F := inferInstance
  haveI : FiniteDimensional F K := Module.Finite.of_restrictScalars_finite ℚ F K
  have htower : Module.finrank ℚ F * Module.finrank F K = Module.finrank ℚ K :=
    Module.finrank_mul_finrank ℚ F K
  intro h
  rw [← htower] at h
  have hF : 0 < Module.finrank ℚ F := Module.finrank_pos
  have : Module.finrank F K = 1 := by
    nlinarith [Module.finrank_pos (R := F) (M := K)]
  exact hne this

/-! ### ★★★★★★★★項目全体の `.src`

★`.src` は「その原典項目を**完全に**実装した」という主張である
（`tools/genell-progress.mjs` の規則）。下の 1 つは、
`Definition 1.2` の (i)(ii) がともに `Found/` に揃ったので置く。 -/

/-- ★★★★★★★★**[GenEll] Definition 1.2** —— (i)(ii) がともに実装された。

原文 (GenEll p.5):
> Definition 1.2.

## ★主張

| 原文 | 宣言 |
|---|---|
| (i) 高さ関数 `ht_M̄ : X(ℚ̄) → ℝ` | ★`ht`（本ファイル） |
| (i) 代表元の取り替えで不変（"any morphism"） | ★★`htOf_baseChange`（本ファイル） |
| (ii) `≾ / ≿ / ≈` の定義 | `BDle` / `BDge` / `BDeq`（`BDClass.lean`） |
| (ii) `≈` は同値関係 | `bdeq_equivalence`（`BDClass.lean`） |
| (ii) BD-類 | `BDClass` / `bdSetoid`（`BDClass.lean`） |
| (ii) `≾ ∧ ≿ ↔ ≈` | `BDClass.eq_of_le_of_ge`（`BDClass.lean`） |
| (ii) BD-類への `≾` の持ち上げ | `BDClass.le` / `BDClass.ge`（`BDClass.lean`） |

## ★★G2 非空虚 / G3 負の対照

* G2: `nfAlgPoint_nonempty`（ℚ-点があれば `X(ℚ̄)` は空でない）
* G3: ★`finrank_not_descends`（**体の次数は点の関数にならない**）
* G3: `bdle_ne_bdge`（`≾` と `≿` は別物）

## ★★★★逸脱の記録（CLAUDE.md の「逸脱」）

### 1. `X(ℚ̄)` を**底変換の商**として作った

原文は `X(ℚ̄)` を代数閉包 `ℚ̄` 上の点の集合として扱い、
`x_F : Spec 𝓞_F → X` を「`x` を生じさせる任意の射」と言う。
★本実装は逆向きに、**対 `(F, x_F)` を底変換で同一視した商**として `X(ℚ̄)` を作る
（`NFAlgPoint`）。★★これは同じものの別の作り方であって、
「any morphism がすべて同じ値を与える」という原文の要求は
`htOf_baseChange` が**定理として**満たしている。

★★★`ℚ̄` 上の点との同一視（`X(ℚ̄) ≃ NFAlgPoint X`）は本項目には要らない
——原文が使うのは高さの well-defined 性だけである。

### 2. 算術直線束は `picardDataWitness.Pic X × Multiplicative (arcCM X)` である

`Definition 1.1` の `APic(X)` の実装（`Found/Arakelov/APicWitness.lean`、
可逆層の類と `X^arc` 上の連続関数の対）を使っている。
★`Definition 1.1` そのものの項目全体はまだ閉じていない
（層の側の欄が残っている）が、本項目が使うのは
`pullback` / `pullback_mul` / `pullback_comp` と `deg_F` だけであり、
それらはすべて `Found/` にある。 -/
def definition_1_2.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5, item := "Definition 1.2",
    sectionId := "genell-def-1-2" }

def definition_1_2.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ht(高さ関数 ht_M̄ : X(ℚ̄) → ℝ ——(i))"
      (.inProject "ABC3" "ABC3.Found.GenEll.ht") 5,
    .citation "[ABC3]" "htOf_baseChange(代表元の取り替えで不変——原文の「any morphism」)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htOf_baseChange") 4,
    .citation "[ABC3]" "degFOf_baseChange([K:ℚ] = [F:ℚ]·[K:F] が正規化の分母を約分する)"
      (.inProject "ABC3" "ABC3.Interface.Arakelov.degFOf_baseChange") 4,
    .citation "[ABC3]" "bdeq_equivalence(≈ は同値関係——(ii))"
      (.inProject "ABC3" "ABC3.Found.GenEll.bdeq_equivalence") 5,
    .citation "[ABC3]" "BDClass(BD-類——(ii))"
      (.inProject "ABC3" "ABC3.Found.GenEll.BDClass") 5,
    .implicitStep
      ("★逸脱 1: X(ℚ̄) を底変換の商 NFAlgPoint として作った。原文は ℚ̄ 上の点の集合として" ++
       "扱い x_F を「x を生じさせる任意の射」と言うが、その要求(any morphism が" ++
       "すべて同じ値を与える)は htOf_baseChange が定理として満たしている") 4,
    .implicitStep
      ("★逸脱 2: 算術直線束は APicWitness.lean の実装(可逆層の類と X^arc 上の連続関数の対)" ++
       "を使う。Definition 1.1 の項目全体はまだ閉じていないが、本項目が使うのは" ++
       "pullback / pullback_mul / pullback_comp と deg_F だけで、すべて Found/ にある。" ++
       "★★Definition 1.1 が閉じない理由は 2 点で測れた(2026-08-27): " ++
       "(a) APicData に ι_X(複素共役)との両立を課す欄が無い —— arcCM X = C(X^arc, ℝ) は" ++
       "無条件なので、実装の APic は原文の APic より大きい。因子表示側には" ++
       "IsConjInvariant があり ArchADivBase.lean で実際に使っている; " ++
       "(b) APicData に pullback_comp の欄が無い —— 本ファイルの aPicPullback_comp は" ++
       "witness についての定理であって Interface の保証ではない。" ++
       "★どちらも本項目の証明には要らない(deg_F は全埋め込みにわたる和なので" ++
       "共役対称性が自動で出る)") 3,
    .implicitStep
      ("★★G3 負の対照: finrank_not_descends —— 体の次数 [F:ℚ] は点の関数にならない。" ++
       "したがって ht が商へ降りるのは自明ではなく htOf_baseChange のおかげである。" ++
       "★これは Definition 1.5 が最小定義体を要る理由と同じ現象") 5 ]

/-! ### ★出典の紐付け(`.src`) -/

def arcCMPullback_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.2, (i)(アルキメデス側の引き戻しの関手性)",
    sectionId := "genell-def-1-2-i" }

def aPicPullback_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.2, (i)(算術直線束の引き戻しの関手性)",
    sectionId := "genell-def-1-2-i" }

def htOf_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.2, (i)(高さの底変換不変性——原文の「any morphism」)",
    sectionId := "genell-def-1-2-i" }

def NFAlgPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(X(ℚ̄)——底変換の商として)",
    sectionId := "genell-def-1-2-i" }

def ht.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(高さ関数 ht_M̄ : X(ℚ̄) → ℝ)",
    sectionId := "genell-def-1-2-i" }

def ht_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(X(ℚ̄) 全体——線束表示なので但し書きなし)",
    sectionId := "genell-prop-1-4" }

def ht_mul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ht(高さが X(ℚ̄) 全体で定義されていること)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ht") 5,
    .implicitStep
      ("★★HeightAdditive.lean の htArith_tensor_unconditional は因子表示なので" ++
       "「x が D を通らない」を仮定していた(原文は可逆層で X(ℚ̄) 全体)。" ++
       "線束表示では ht が X(ℚ̄) 全体で定義されているので但し書きが要らない") 6 ]

def finrank_not_descends.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(★G3 負の対照——体の次数は点の関数にならない)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
