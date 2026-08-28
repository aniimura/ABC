/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.ArcConjInvol
import ABC3.Found.Arakelov.ArchDegSmul
import ABC3.Found.Arakelov.AMetricNorm
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★算術直線束の計量が**複素共役 `ι_X` と両立する**こと（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX.

## ★★★★★★★★★★これは何か——Definition 1.1 (i) の**最後の条項**

原文の算術直線束は `(L, |−|_L)` であって、その計量に

> compatible [in the evident sense] with the complex conjugation automorphism `ι_X`

が課されている。★★いまの設計（`AMetric` ＝ 前層＋局所自明性＋`LocalMetric`）には
**この条項が無かった**（2026-08-28 の実測）——`IsConjCompatible` は
`GreenMetric`（`ArcGreenConj.lean`）と `TorsorMetric`（`ArcTorsorMetric.lean`）にしか無い。

★★★本ファイルはそれを `LocalMetric` / `AMetric` に入れる。

## ★★★★★「evident sense」を読む

`ι_X` は複素点を `conjPoint p = conjSpec ≫ p` に送る（`ArchConj.lean`）。
★アフィンでは **`ι_X` は関数の値を複素共役にする**（`evalAffine_conjPoint`、`ArcConjInvol.lean`）。

    `a(ι_X p) = conj (a(p))`

★★したがって計量の側で `ι_X` と両立するとは

    `h_{V,e}(ι_X p) = h_{V,e}(p)`   （基準ノルムが共役で不変）

ということである。★★★このとき**切断のノルムも共役で不変**になる:

    `|s|(ι_X p) = |s|(p)`

——`|conj z| = |z|` だからである。それが `normOf_conjPoint` / `norm_conjPoint` である。

## ★★★★★★機構（3 段）

| 段 | 内容 |
|---|---|
| `conjPoint_preimage_eq_top` | ★`ι_X p` は `p` と同じ開集合に落ちる |
| `evalOn_conjPoint` | ★★★**開集合上でも `ι_X` は値を複素共役にする**（`evalOn_specMap` ＋ `appLE_comp_appLE`） |
| `normOf_conjPoint` | ★★★★★切断のノルムが共役で不変になる |

★★二段目が本体である。在庫の `evalOn_specMap`（`§9-775`）が
`conjSpec = Spec.map (starRingEnd ℂ)` にそのまま効き、
`Scheme.Hom.appLE_comp_appLE` で合成を割る——**新しい数学は要らない**。

## ★測定の記録

`conjSpec.appLE ⊤ ⊤ _` は `(Spec.map (ofHom (starRingEnd ℂ))).appLE ⊤ ⊤ _` と
**定義的に等しい**ので、`evalOn_specMap` を `rw [evalOn_appLE]` で剥がすだけで渡れる
（2026-08-28 実測）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

/-! ## ★★★(1) `ι_X` は開集合を跨がない -/

/-- ★**`ι_X p` は `p` と同じ開集合に落ちる**。

★`conjPoint p = conjSpec ≫ p` であり、`conjSpec ⁻¹ᵁ ⊤ = ⊤` だからである。 -/
theorem conjPoint_preimage_eq_top {X : Scheme.{0}} {V : X.Opens}
    {p : Spec (CommRingCat.of ℂ) ⟶ X} (h : p ⁻¹ᵁ V = ⊤) : conjPoint p ⁻¹ᵁ V = ⊤ := by
  show (conjSpec ≫ p) ⁻¹ᵁ V = ⊤
  rw [Scheme.Hom.comp_preimage, h]
  simp

/-! ## ★★★★★(2) 開集合上でも `ι_X` は値を複素共役にする -/

/-- ★★★**`conjSpec` の大域切断への作用は複素共役そのもの**。

★在庫の `evalOn_specMap`（`§9-775`）を `conjSpec = Spec.map (starRingEnd ℂ)` に当てる。 -/
theorem evalOn_conjSpec (y : (Γ(Spec (CommRingCat.of ℂ),
    (⊤ : (Spec (CommRingCat.of ℂ)).Opens)) : Type)) (hc) :
    (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((conjSpec.appLE ⊤ ⊤ hc).hom y)
      = starRingEnd ℂ ((Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom y) := by
  have h2 := evalOn_specMap (CommRingCat.of ℂ) (CommRingCat.ofHom (starRingEnd ℂ))
    (by simp) y
  rw [evalOn_appLE] at h2
  exact h2

/-- ★★★★★**開集合上でも `ι_X` は関数の値を複素共役にする**。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX.

    `c(ι_X p) = conj (c(p))`

★`ArcConjInvol.lean` の `evalAffine_conjPoint` を任意の開集合 `V` に広げたものである。 -/
theorem evalOn_conjPoint {X : Scheme.{0}} (V : X.Opens) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (h : p ⁻¹ᵁ V = ⊤) (c : ((X.presheaf.obj (op V)) : Type)) :
    evalOn (conjPoint p) V (conjPoint_preimage_eq_top h) c
      = starRingEnd ℂ (evalOn p V h c) := by
  have hcomp : (p.appLE V ⊤ (le_of_eq h.symm)) ≫ (conjSpec.appLE ⊤ ⊤ (by simp))
      = (conjSpec ≫ p).appLE V ⊤ (le_of_eq (conjPoint_preimage_eq_top h).symm) :=
    Scheme.Hom.appLE_comp_appLE conjSpec p V ⊤ ⊤ _ _
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      (((conjSpec ≫ p).appLE V ⊤ (le_of_eq (conjPoint_preimage_eq_top h).symm)).hom c) = _
  rw [← hcomp]
  exact evalOn_conjSpec ((p.appLE V ⊤ (le_of_eq h.symm)).hom c) (by simp)

/-! ## ★★★★★★★(3) 計量が `ι_X` と両立するということ -/

variable {X : Scheme.{0}} {F : X.PresheafOfModules}

/-- ★★★★★★★★★★**[GenEll] Definition 1.1, (i)** ——
計量が**複素共役 `ι_X` と両立する**こと。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX.

★「evident sense」の中身は「基準ノルムが `ι_X` で不変」である
——切断の側の値は共役になるが、絶対値は変わらないからである（`normOf_conjPoint`）。

★★`p ⁻¹ᵁ V = ⊤`（`p` が `V` に落ちる）に限って課す。★★★そこを外した `h V e p` は
`dite` の `else` 枝が返す**意味の無い値**であり、等長同型（`IsIsometry`）もそこでは
何も言わないからである——この制限があって初めて条件が等長で移る。 -/
def LocalMetric.IsConjCompatible (m : LocalMetric X F) : Prop :=
  ∀ (V : X.Opens) (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X), p ⁻¹ᵁ V = ⊤ → m.h V e (conjPoint p) = m.h V e p

/-- ★★★★★★★★**`ι_X` と両立する計量では切断のノルムが共役で不変**。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX.

    `|s|(ι_X p) = |s|(p)`

★機構は `evalOn_conjPoint`（値が共役になる）＋ `Complex.abs_conj`（絶対値は変わらない）。 -/
theorem LocalMetric.normOf_conjPoint {m : LocalMetric X F} (hm : m.IsConjCompatible)
    (V : X.Opens) (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : F.obj (op ⊤)) :
    m.normOf V e (conjPoint p) (conjPoint_preimage_eq_top hp) s = m.normOf V e p hp s := by
  show trivSecNorm F V e (conjPoint p) (conjPoint_preimage_eq_top hp) s
      * m.h V e (conjPoint p) = trivSecNorm F V e p hp s * m.h V e p
  rw [hm V e p hp]
  congr 1
  show ‖evalOn (conjPoint p) V (conjPoint_preimage_eq_top hp) (trivValue F V e s)‖
    = ‖evalOn p V hp (trivValue F V e s)‖
  rw [evalOn_conjPoint V p hp (trivValue F V e s)]
  exact RCLike.norm_conj _

/-! ## ★★★★★★★★★(4) 算術直線束の側 -/

/-- ★★★★★★★★★★**[GenEll] Definition 1.1, (i)** ——
算術直線束 `L̄ = (L, |−|_L)` の計量が `ι_X` と両立すること。 -/
def AMetric.IsConjCompatible (L : AMetric X) : Prop :=
  L.metric.IsConjCompatible

/-- ★★★★★★★★★★**大域ノルムは `ι_X` で不変である**。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX.

    `|s|_L̄(ι_X p) = |s|_L̄(p)`

★★これが原文の「compatible with `ι_X`」が**切断の水準で言っていること**である。 -/
theorem AMetric.norm_conjPoint {L : AMetric X} (hL : L.IsConjCompatible)
    (s : L.sheaf.obj (op ⊤)) (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    L.norm s (conjPoint p) = L.norm s p := by
  obtain ⟨c⟩ := nonempty_normChart L.triv p
  rw [AMetric.norm_eq L s (conjPoint p) c.V c.e (conjPoint_preimage_eq_top c.hp),
    AMetric.norm_eq L s p c.V c.e c.hp]
  exact LocalMetric.normOf_conjPoint hL c.V c.e p c.hp s

/-! ## ★出典の紐付け(`.src`) -/

def evalOn_conjPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(開集合上でも ι_X は関数の値を複素共役にする)",
    sectionId := "genell-def-1-1-i" }

def LocalMetric.IsConjCompatible.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(計量が複素共役 ι_X と両立すること——「evident sense」の中身)",
    sectionId := "genell-def-1-1-i" }

def AMetric.norm_conjPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(ι_X と両立する計量では |s|(ι_X p) = |s|(p))",
    sectionId := "genell-def-1-1-i" }

def AMetric.norm_conjPoint.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "evalOn_specMap(Spec.map f が定める点での値は f そのもの、§9-775)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.evalOn_specMap") 4,
    .citation "[mathlib]" "Scheme.Hom.appLE_comp_appLE"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.appLE_comp_appLE") 3,
    .implicitStep
      ("★★残るのは閉性である: 1（構造層の標準計量）・⊗・引き戻し・等長で " ++
       "IsConjCompatible が保たれることを示すと、APicM の共役両立な部分群の上で " ++
       "deg_F / ht がそのまま使える") 3 ]

end ABC3.Found.Arakelov
