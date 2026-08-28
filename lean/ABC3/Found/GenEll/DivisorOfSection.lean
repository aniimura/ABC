/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.VeryAmpleDef
import ABC3.Meta.Claim

/-!
# ★★★★★★★**切断の零点因子** `div(s)` —— 加群層と因子を繋ぐ橋（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★これは台帳の段 E0 である

`ResearchPaper/mathlib-gap.json` の `ample-and-projective-embedding` は、
Serre の定理（段 E）を 4 つに割った（2026-08-28）:

| 段 | 内容 |
|---|---|
| **E0** | **加群層 ⟷ 因子の橋**（切断の零点因子 `div(s)` の構成） ← ★本ファイル |
| E1 | 大域切断が定める `ℙᴺ` への射（`Proj` の普遍性） |
| E2 | ample ⟹ 有限個の切断で覆える＋次数を揃える |
| E3 | その射が immersion であること |

★`IsAmple`（`AmpleDef.lean`）は**加群層**の述語、
`IsVeryAmpleDiv`（`VeryAmpleDef.lean`）は**因子**の述語なので、
段 E を述べるにはこの 2 つを繋ぐ必要がある。

## ★★★★★機構 —— 自明化ごとの `span` を `⨅` で束ねる

自明化 `e : M|_V ≅ 𝟙_` のもとで `s` は関数 `trivValue M V e s ∈ Γ(X,V)` になる。
★2 つの自明化は**単元倍**で違う（`trivValue_congr`）ので、
生成するイデアル `span {trivValue}` は**自明化に依らない**（`span_trivValue_congr`）。

★★そこで

    `divIdeal M s U ≝ ⨅_{e} span {trivValue M U e s}`

と置く。★★★`⨅` は**空なら `⊤`**なので、自明化が無いアフィン開では `⊤` になる
——これが素直な定義である。

## ★★★★★★`IdealSheafData` へ載せる

mathlib の `Scheme.IdealSheafData` は「**すべての**アフィン開に対するイデアル」と
`map_ideal_basicOpen`（`I(D(f)) = I(U)_f`）を要求する。

★**自明化がある開の上ではその両立が成り立つ**（`divIdeal_map_basicOpen`）
——機構は `trivValue_restrict`（段 D2）だけである。

★★しかし**自明化が無い開では成り立たない**（`⊤` の像は `⊤` だが、
基本開集合の側には自明化があって `⊤` でないことがある）。
★★★そこで `IdealSheafData.ofIdeals`（族に含まれる**最大の**イデアル層）で受ける。

## ★残っている段（明示）

★**`(divisorOfSection M s).ideal U = divIdeal M s U`（`U` が自明化する場合）は
まだ示していない**。`ofIdeals` は `≤` しか与えない（`divisorOfSection_ideal_le`）。
★★等号には「自明化する開が基底をなす」ことを使う議論が要る。
★★★それが済むまで `div(s)` は**上からの近似**である——本ファイルはそう記録する。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★★(1) 生成するイデアルは自明化に依らない -/

/-- ★★★★★**切断が生成するイデアルは自明化の取り方に依らない**。

★機構は `trivValue_congr`（2 つの自明化は単元倍で違う）だけである。 -/
theorem span_trivValue_congr (M : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) :
    Ideal.span {trivValue M V e' s} = Ideal.span {trivValue M V e s} := by
  obtain ⟨u, hu, huv⟩ := trivValue_congr M V e e' s
  rw [huv]
  exact Ideal.span_singleton_mul_right_unit hu _

/-! ## ★★★★★★(2) アフィン開ごとのイデアル -/

open scoped Classical in
/-- ★★★★★★**切断 `s` の零点因子のイデアル**（アフィン開 `U` ごと）。

★自明化が無い `U` では `⊤`（`⨅` が空だから）。 -/
noncomputable def divIdeal (M : X.PresheafOfModules) (s : M.obj (op ⊤))
    (U : X.affineOpens) : Ideal Γ(X, U.1) :=
  ⨅ e : ((restrictPresheafFunctor X U.1).obj M ≅ 𝟙_ (PresheafModulesOn X U.1)),
    Ideal.span {trivValue M U.1 e s}

/-- ★★**自明化があれば、どれで測ってもよい**。 -/
theorem divIdeal_eq (M : X.PresheafOfModules) (s : M.obj (op ⊤)) (U : X.affineOpens)
    (e : (restrictPresheafFunctor X U.1).obj M ≅ 𝟙_ (PresheafModulesOn X U.1)) :
    divIdeal M s U = Ideal.span {trivValue M U.1 e s} := by
  refine le_antisymm (iInf_le _ e) (le_iInf fun e' => ?_)
  rw [span_trivValue_congr M U.1 e e' s]

/-- ★自明化が無ければ `⊤` である。 -/
theorem divIdeal_eq_top (M : X.PresheafOfModules) (s : M.obj (op ⊤)) (U : X.affineOpens)
    (h : IsEmpty ((restrictPresheafFunctor X U.1).obj M ≅ 𝟙_ (PresheafModulesOn X U.1))) :
    divIdeal M s U = ⊤ :=
  iInf_of_empty _

/-! ## ★★★★★★★(3) 基本開集合との両立（自明化がある側） -/

/-- ★★★★★★★**自明化がある開の上では `map_ideal_basicOpen` が成り立つ**。

    `(divIdeal U).map (制限) = divIdeal (D(f))`

★機構は `trivValue_restrict`（段 D2）だけである
——自明化は `trivialOfLe` で `D(f)` へ降り、`trivValue` は制限と可換だからである。 -/
theorem divIdeal_map_basicOpen (M : X.PresheafOfModules) (s : M.obj (op ⊤))
    (U : X.affineOpens) (f : Γ(X, U.1))
    (e : (restrictPresheafFunctor X U.1).obj M ≅ 𝟙_ (PresheafModulesOn X U.1)) :
    (divIdeal M s U).map (X.presheaf.map (homOfLE (X.basicOpen_le f)).op).hom
      = divIdeal M s (X.affineBasicOpen f) := by
  have hle : (X.affineBasicOpen f).1 ≤ U.1 := X.basicOpen_le f
  rw [divIdeal_eq M s U e, Ideal.map_span,
    divIdeal_eq M s (X.affineBasicOpen f) (trivialOfLe M hle e),
    trivValue_restrict M hle e s, Set.image_singleton]
  rfl

/-! ## ★★★★★★(4) イデアル層として受ける -/

/-- ★★★★★★**切断の零点因子** `div(s)`。

★★`IdealSheafData.ofIdeals`（族に含まれる**最大の**イデアル層）で受ける
——自明化が無いアフィン開があるので、族そのものはイデアル層とは限らない。 -/
noncomputable def divisorOfSection (M : X.PresheafOfModules) (s : M.obj (op ⊤)) :
    X.IdealSheafData :=
  Scheme.IdealSheafData.ofIdeals (divIdeal M s)

/-- ★**`ofIdeals` が与えるのは上からの近似である**（等号はまだ示していない）。 -/
theorem divisorOfSection_ideal_le (M : X.PresheafOfModules) (s : M.obj (op ⊤)) :
    (divisorOfSection M s).ideal ≤ divIdeal M s :=
  Scheme.IdealSheafData.ideal_ofIdeals_le _

/-- ★★★**健全性（空虚封じ）**——構造層 `𝒪_X` では `div(s)` の欄は
`s` が生成するイデアルそのものである。

★これで `divIdeal` が「知っている `span {s}`」の一般化だと型で言える。 -/
theorem divIdeal_unit (s : ((𝟙_ X.PresheafOfModules).obj (op (⊤ : X.Opens)) : Type))
    (U : X.affineOpens) :
    divIdeal (𝟙_ X.PresheafOfModules) s U
      = Ideal.span {(X.presheaf.map (homOfLE (le_top : U.1 ≤ ⊤)).op).hom s} := by
  rw [divIdeal_eq _ s U (trivialOfLe _ (le_top : U.1 ≤ ⊤) restrictPresheafUnit.symm),
    trivValue_restrict _ (le_top : U.1 ≤ ⊤) restrictPresheafUnit.symm s,
    trivValue_unit_top]

/-! ## ★出典の紐付け(`.src`) -/

def span_trivValue_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(切断が生成するイデアルは自明化に依らない)",
    sectionId := "genell-prop-1-4" }

def divIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——切断の零点因子のイデアル、アフィン開ごと)",
    sectionId := "genell-prop-1-4" }

def divIdeal_map_basicOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(自明化がある開の上では基本開集合と両立する)",
    sectionId := "genell-prop-1-4" }

def divisorOfSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(切断の零点因子——ofIdeals で受けた上からの近似)",
    sectionId := "genell-prop-1-4" }

def divIdeal_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(健全性——構造層なら div(s) の欄は span {s})",
    sectionId := "genell-prop-1-4" }

def divisorOfSection.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "trivValue_congr(2 つの自明化は単元倍で違う)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_congr") 6,
    .citation "[ABC3]" "trivValue_restrict(自明化の制限と trivValue は可換)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_restrict") 6,
    .citation "[mathlib]" "Scheme.IdealSheafData.ofIdeals / ideal_ofIdeals_le"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.ofIdeals") 6,
    .implicitStep
      ("★**`(divisorOfSection M s).ideal U = divIdeal M s U`（`U` が自明化する場合）は" ++
       "まだ示していない**——`ofIdeals` は `≤` しか与えない。" ++
       "★★等号には「自明化する開が基底をなす」ことを使う議論が要る。" ++
       "★★★それが済むまで div(s) は**上からの近似**である") 6,
    .implicitStep
      ("★これは mathlib-gap.json の ample-and-projective-embedding の**段 E0** である。" ++
       "★★段 E1(大域切断が定める ℙᴺ への射)・E2・E3(immersion 性)が残る") 6 ]

end ABC3.Found.GenEll
