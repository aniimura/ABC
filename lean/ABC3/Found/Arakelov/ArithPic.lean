/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.APicWitness
import ABC3.Found.Arakelov.ArcConjInvol
import ABC3.Found.Arakelov.ArithSections

/-!
# Arakelov —— **`APic(X)`（`ι_X` 両立を型に入れた版）と引き戻し**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★`Definition 1.1` の組み上げ（あと 1 段）

`Found/Arakelov/APicWitness.lean` の `aPicDataWitness` は
`APic(X) = Pic(X) × Multiplicative C(X^arc, ℝ)` を作っているが、
★**原文が要求する「計量は `ι_X` と両立する」を型に入れていない**
（`Interface` の `APicData` にその欄が無いため）。

★★本ファイルはそれを型に入れた版 `ArithPic X` を作る:

    ArithPic X ≝ Pic(X) × Multiplicative (共役不変な連続関数のなす部分群)

| 原文の語 | 本ファイル |
|---|---|
| `L̄ = (L, |−|_L)`、`|−|_L` は `ι_X` 両立 | `arithPicOfMetric L m hm` |
| テンソル積 | ★**群の乗法がそのままテンソル積**（`arithPicOfMetric_mul` は `rfl`） |
| `APic(X)` が群であること | `arithPicGroup` |
| (ii) 引き戻し `φ^*L̄` | `arithPicPullback` / `arithPicPullback_mul` |
| 射・`Γ(L̄)` | `Found/Arakelov/ArithSections.lean` |

## ★★★★★捻れ集合の表示が「余計な元を持たない」こと

`ArithPic` は Green 関数の側だけを持つので、
「本当に計量から来ているのか」を別に言う必要がある。それが

* `arithPicOfMetric_surjective` —— すべての類は `ι_X` 両立な `TorsorMetric` から来る
* `arithPicForget_ofMetric` —— 忘れると元の直線束に戻る

の 2 本である。★`TorsorMetric` は `(green, green_cont, triv)` の 3 つ組なので
どちらも `rfl` で出る。

## ★★★★★★★引き戻しが共役不変性を保つ理由

`conjPoint (p ≫ f) = conjPoint p ≫ f`（`ArcConjInvol.lean` の `conjPoint_comp`）
——★`ι_X` は射に沿って自然だからである。

## ★残っている 1 段（明示）

★★**`Γ(L̄) = Hom(Ō_X, L̄)` の同一視**が残っている。
原文は `Γ(L̄)` を「`Ō_X → L̄` の射の集合」と**定義**するが、
本プロジェクトは `Γ(L̄)` を「`|s|_L ≤ 1` なる大域切断の集合」として定義した
（`ArithSections.lean`）。★★★両者が一致するには
`|φ(s)|_L = |s|·|φ(1)|_L`（計量の `𝒪_X`-乗法性）が要る。
★これが `Definition 1.1` を**条なしで**閉じるための最後の 1 本である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite ABC3.Interface.Arakelov ABC3.Found.GenEll

/-! ## ★★共役不変な Green 関数のなす部分群 -/

/-- ★★**共役不変な Green 関数のなす部分群**——原文『compatible with `ι_X`』。 -/
noncomputable def conjArcCM (X : Scheme.{0}) : AddSubgroup (arcCM X) where
  carrier := {c | IsConjInvariant c.fn}
  add_mem' := by
    intro a b ha hb p
    show (a.fn (conjPoint p) + b.fn (conjPoint p)) = a.fn p + b.fn p
    rw [ha p, hb p]
  zero_mem' := by
    intro p; rfl
  neg_mem' := by
    intro a ha p
    show -(a.fn (conjPoint p)) = -(a.fn p)
    rw [ha p]

theorem mem_conjArcCM_iff {X : Scheme.{0}} (c : arcCM X) :
    c ∈ conjArcCM X ↔ ∀ p : complexPoints X, c.fn (conjPoint p) = c.fn p := Iff.rfl

/-- ★★**引き戻しは共役不変性を保つ**——`conjPoint_comp`（`ι_X` の自然性）による。 -/
theorem arcCMPullback_mem_conjArcCM {X Y : Scheme.{0}} (f : X ⟶ Y) {c : arcCM Y}
    (hc : c ∈ conjArcCM Y) : arcCMPullback f c ∈ conjArcCM X := by
  intro p
  show c.fn (conjPoint p ≫ f) = c.fn (p ≫ f)
  rw [← conjPoint_comp f p]
  exact hc (p ≫ f)

/-! ## ★★★★★★★★`APic(X)` -/

/-- ★★★★★★★★**算術 Picard 群 `APic(X)`**（`ι_X` 両立を型に入れた版）。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★捻れ集合の表示——直線束と、`ι_X` 両立な Green 関数の対。 -/
@[reducible] noncomputable def ArithPic (X : Scheme.{0}) : Type 1 :=
  picardDataWitness.Pic X × Multiplicative ↥(conjArcCM X)

/-- ★★**群構造**——`Pic` の積と Green 関数の和。 -/
@[reducible] noncomputable def arithPicGroup (X : Scheme.{0}) : CommGroup (ArithPic X) :=
  letI := picardDataWitness.group X
  inferInstance

/-- ★★**[GenEll] Definition 1.1, (ii)** —— 引き戻し `φ^*L̄`。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z -proper, Z -flat schemes. Then -/
noncomputable def arithPicPullback {X Y : Scheme.{0}} (f : X ⟶ Y) (L : ArithPic Y) :
    ArithPic X :=
  (picardDataWitness.pullback f L.1,
    Multiplicative.ofAdd (⟨arcCMPullback f ((Multiplicative.toAdd L.2 : ↥(conjArcCM Y)) : arcCM Y),
      arcCMPullback_mem_conjArcCM f (Multiplicative.toAdd L.2).2⟩ : ↥(conjArcCM X)))

/-- ★引き戻しは群準同型である。 -/
theorem arithPicPullback_mul {X Y : Scheme.{0}} (f : X ⟶ Y) (L M : ArithPic Y) :
    arithPicPullback f
        (@HMul.hMul _ _ _
          (@instHMul _ (arithPicGroup Y).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
      = @HMul.hMul _ _ _
          (@instHMul _ (arithPicGroup X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
          (arithPicPullback f L) (arithPicPullback f M) := by
  refine Prod.ext (picardDataWitness.pullback_mul f L.1 M.1) ?_
  rfl

/-! ## ★★★★計量との対応 -/

/-- ★計量を忘れる。 -/
noncomputable def arithPicForget {X : Scheme.{0}} (L : ArithPic X) :
    picardDataWitness.Pic X := L.1

/-- ★対から作る（Green 関数の側を直に与える形）。 -/
noncomputable def arithPicMk {X : Scheme.{0}} (L : picardDataWitness.Pic X)
    (c : arcCM X) (hc : c ∈ conjArcCM X) : ArithPic X :=
  (L, Multiplicative.ofAdd (⟨c, hc⟩ : ↥(conjArcCM X)))

@[simp] theorem arithPicForget_mk {X : Scheme.{0}} (L : picardDataWitness.Pic X)
    (c : arcCM X) (hc : c ∈ conjArcCM X) :
    arithPicForget (arithPicMk L c hc) = L := rfl

theorem arithPicMk_surjective {X : Scheme.{0}} (x : ArithPic X) :
    ∃ (L : picardDataWitness.Pic X) (c : arcCM X) (hc : c ∈ conjArcCM X),
      arithPicMk L c hc = x :=
  ⟨x.1, ((Multiplicative.toAdd x.2 : ↥(conjArcCM X)) : arcCM X),
    (Multiplicative.toAdd x.2).2, rfl⟩

/-- ★★★★**`ι_X` 両立な計量から算術直線束の類を作る**——原文の `L̄ = (L, |−|_L)` そのもの。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
noncomputable def arithPicOfMetric {X : Scheme.{0}} (L : picardDataWitness.Pic X)
    (m : TorsorMetric X (picardDataWitness.sheafOf X L)) (hm : m.IsConjCompatible) :
    ArithPic X :=
  arithPicMk L (arcCM.mk X m.green m.green_cont) hm

@[simp] theorem arithPicForget_ofMetric {X : Scheme.{0}} (L : picardDataWitness.Pic X)
    (m : TorsorMetric X (picardDataWitness.sheafOf X L)) (hm : m.IsConjCompatible) :
    arithPicForget (arithPicOfMetric L m hm) = L := rfl

/-- ★★★★★**すべての類は `ι_X` 両立な計量から来る**——`ArithPic` に余計な元は無い。 -/
theorem arithPicOfMetric_surjective {X : Scheme.{0}} (x : ArithPic X) :
    ∃ (L : picardDataWitness.Pic X) (m : TorsorMetric X (picardDataWitness.sheafOf X L))
      (hm : m.IsConjCompatible), arithPicOfMetric L m hm = x :=
  ⟨x.1, ⟨((Multiplicative.toAdd x.2 : ↥(conjArcCM X)) : arcCM X).fn,
    arcCM.fn_cont _, picSheaf_locallyTrivial X x.1⟩, (Multiplicative.toAdd x.2).2, rfl⟩

/-- ★★★★★★**群の乗法がそのままテンソル積である**。

原文 (GenEll p.3):
> of tensor product, thus determine a group APic(X).

★**`rfl` で出る**——`TorsorMetric.tensor` の Green 関数が和だからである。 -/
theorem arithPicOfMetric_mul {X : Scheme.{0}}
    (L L' : picardDataWitness.Pic X)
    (m : TorsorMetric X (picardDataWitness.sheafOf X L)) (hm : m.IsConjCompatible)
    (m' : TorsorMetric X (picardDataWitness.sheafOf X L')) (hm' : m'.IsConjCompatible) :
    @HMul.hMul _ _ _
        (@instHMul _ (arithPicGroup X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
        (arithPicOfMetric L m hm) (arithPicOfMetric L' m' hm')
      = arithPicOfMetric
          (@HMul.hMul _ _ _
            (@instHMul _ (picardDataWitness.group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
            L L')
          (TorsorMetric.tensor m m' (picSheaf_locallyTrivial X _))
          (IsConjInvariant.add hm hm') := rfl

/-! ## ★★★`Γ(L̄)`（`TorsorMetric` 版） -/

variable {X : Scheme.{0}} {F : X.Modules}

/-- ★★★**`Γ(L̄)`** ——`|s|_L ≤ 1` なる大域切断の集合（`TorsorMetric` 版）。 -/
def TorsorMetric.arithSections (m : TorsorMetric X F) : Set (F.val.obj (op ⊤) : Type) :=
  {s | ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, m.norm s p ≤ 1}

theorem TorsorMetric.zero_mem_arithSections (m : TorsorMetric X F) :
    (0 : (F.val.obj (op ⊤) : Type)) ∈ m.arithSections := by
  intro p
  have h : m.norm (0 : (F.val.obj (op ⊤) : Type)) p = 0 := by
    show Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F 0) = 0
    rw [arcEval_zero, (m.base.eq_zero_iff p 0).2 rfl, mul_zero]
  rw [h]; norm_num

/-! ### ★出典の紐付け(`.src`)

★★★**項目全体の `.src` はまだ置かない。** 残っているのは
`Γ(L̄) = Hom(Ō_X, L̄)` の同一視 1 本である(上の docstring)。 -/

def ArithPic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(APic(X)——ι_X 両立を型に入れた版)",
    sectionId := "genell-def-1-1-i" }

def arithPicOfMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(ι_X 両立な計量から算術直線束の類を作る)",
    sectionId := "genell-def-1-1-i" }

def arithPicOfMetric_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(群の乗法がテンソル積であること)",
    sectionId := "genell-def-1-1-i" }

def arithPicPullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻し φ^*L̄——ι_X 両立を保つこと込み)",
    sectionId := "genell-def-1-1-i" }

def ArithPic.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "picardDataWitness(スキーム上の直線束の群 Pic(X))"
      (.inProject "ABC3" "ABC3.Interface.Arakelov.picardDataWitness") 3,
    .citation "[ABC3]" "TorsorMetric(計量 = 基準計量 × exp(-Green))"
      (.inProject "ABC3" "ABC3.Found.Arakelov.TorsorMetric") 3,
    .citation "[ABC3]" "conjPoint_comp(ι_X は射に沿って自然)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.conjPoint_comp") 3,
    .implicitStep
      ("★aPicDataWitness は APic(X) = Pic(X) × Multiplicative C(X^arc, ℝ) を作るが、" ++
       "Interface の APicData に「計量は ι_X と両立する」欄が無いので" ++
       "その条件を型に持っていない。本ファイルはそれを型に入れた版である") 3,
    .implicitStep
      ("★★残っている 1 段: Γ(L̄) = Hom(Ō_X, L̄) の同一視。" ++
       "原文は Γ(L̄) を「Ō_X → L̄ の射の集合」と定義するが、" ++
       "本プロジェクトは「|s|_L ≤ 1 なる大域切断の集合」として定義した。" ++
       "★★★一致には計量の 𝒪_X-乗法性 |φ(s)|_L = |s|·|φ(1)|_L が要る") 3 ]

end ABC3.Found.Arakelov
