/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.APicWitness
import ABC3.Found.Arakelov.ArcConjInvol
import ABC3.Found.Arakelov.ArithSections
import ABC3.Found.Arakelov.ArcUnitEval
import ABC3.Found.Arakelov.ArcEvalGlobal

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

## ★★★★★★計量の `𝒪_X`-乗法性は本ファイルに入っている

    `|c · s|_L (p) = ‖c(p)‖ · |s|_L (p)`（`TorsorMetric.norm_smul`）

★これが `Γ(L̄) = Hom(Ō_X, L̄)` の同一視の**数学的な核**である
——`φ(s) = s · φ(1)` なので `|φ(s)| = ‖s(p)‖ · |φ(1)|` となり、
`|φ(s)| ≤ |s|_{Ō}` は `|φ(1)| ≤ 1` と同値になる。

## ★残っている 1 段（明示）

★★**正規化した自明な算術直線束 `Ō_X`（`|1| = 1`）**だけが残っている。
内訳は 3 つで、**2 つはもうある**:

| 段 | 状態 |
|---|---|
| 正性 `u(p) ≠ 0` | ★**閉じた**——`arcEval_unit_one_ne_zero`（2026-08-28） |
| 連続性 | ★**ある**——`continuous_evalGlobal`（`ArcEvalGlobal.lean`） |
| ファイバーの**正規な ℂ-ノルム** | ★★残っている |

★★★三つ目の中身は、`Scheme.ΓSpecIso` と
`moduleSpecΓFunctor` が `Γ(Spec ℂ, ⊤)` に入れる ℂ-加群構造の**両立**:

    `(ΓSpecIso ℂ).hom (c • y) = c · (ΓSpecIso ℂ).hom y`

★mathlib にこの補題は無い（2026-08-28 実測、`exact?` ・ `simp` とも失敗）。
★★★これが `Definition 1.1` を**条なしで**閉じるための最後の 1 本である。
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

/-! ## ★★★★★★ノルムは `Γ(X, ⊤)` の作用について乗法的 -/

/-- ★前層加群の恒等射での移送は元を動かさない。

★★mathlib の `PresheafOfModules.map_id` は `restrictScalarsId'` を経由するので、
**元の水準の形**を別に用意する。 -/
theorem PresheafOfModules.map_id_apply' {C : Type} [Category C] {R : Cᵒᵖ ⥤ RingCat}
    (M : PresheafOfModules R) (Y : Cᵒᵖ) (x : (M.obj Y : Type)) :
    (M.map (𝟙 Y)).hom x = x := by
  rw [PresheafOfModules.map_id]; rfl

variable (F) (p : Spec (CommRingCat.of ℂ) ⟶ X)

theorem preimage_top_eq_top : p ⁻¹ᵁ (⊤ : X.Opens) = ⊤ := rfl

/-- ★`⊤` 上では `arcEvalOnTop` は `arcEval` そのものである。 -/
theorem arcEvalOnTop_top (s : (F.val.obj (op (⊤ : X.Opens)) : Type)) :
    arcEvalOnTop F p ⊤ (preimage_top_eq_top p) s = arcEval p F s := by
  show ((((Scheme.Modules.pullback p).obj F).val.map
    (homOfLE (le_of_eq (preimage_top_eq_top p).symm)).op)).hom (arcEvalOn F p ⊤ s) = _
  have h1 : ((homOfLE (le_of_eq (preimage_top_eq_top p).symm)
      : (⊤ : (Spec (CommRingCat.of ℂ)).Opens) ⟶ p ⁻¹ᵁ (⊤ : X.Opens))).op
      = 𝟙 (op (p ⁻¹ᵁ (⊤ : X.Opens))) := rfl
  rw [h1]
  exact PresheafOfModules.map_id_apply' _ _ _

/-- ★★**評価は `Γ(X, ⊤)` の作用について半線形**——第 257 の `arcEvalOn_smul` を `⊤` で読む。 -/
theorem arcEval_smul (c : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type))
    (s : (F.val.obj (op (⊤ : X.Opens)) : Type)) :
    arcEval p F (c • s) = (evalOn p ⊤ (preimage_top_eq_top p) c) • arcEval p F s := by
  have h := arcEvalOnTop_smul F p ⊤ (preimage_top_eq_top p) c s
  rw [arcEvalOnTop_top, arcEvalOnTop_top] at h
  exact h

/-- ★★★★★★**ノルムは `Γ(X, ⊤)` の作用について乗法的**:

    `|c · s|_L (p) = ‖c(p)‖ · |s|_L (p)`

★★★これが `Γ(L̄) = Hom(Ō_X, L̄)` の同一視の**数学的な核**である
——`φ(s) = s · φ(1)` なので `|φ(s)| = ‖s(p)‖ · |φ(1)|` となり、
`|φ(s)| ≤ |s|_{Ō}` は `|φ(1)| ≤ 1`（`Ō_X` を `|1| = 1` で正規化したとき）と同値になる。 -/
theorem TorsorMetric.norm_smul (m : TorsorMetric X F)
    (c : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type))
    (s : (F.val.obj (op (⊤ : X.Opens)) : Type)) :
    m.norm (c • s) p = ‖evalOn p ⊤ (preimage_top_eq_top p) c‖ * m.norm s p := by
  show Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F (c • s)) = _
  rw [arcEval_smul, m.base.smul]
  show _ = ‖evalOn p ⊤ (preimage_top_eq_top p) c‖
    * (Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F s))
  ring

/-! ## ★★★★構造層のファイバーで `1` は消えない -/

/-- ★★★★**構造層のファイバーで `1` は消えない**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これは**正規化した自明な算術直線束 `Ō_X`（`|1| = 1`）を作るための正性**である
——`|1| = exp(-green)·u` を `1` にするには `u ≝ base.nrm p (arcEval p 1) > 0` が要り、
`ArcMetric.eq_zero_iff` によりそれは `arcEval p 1 ≠ 0` に等しい。

★機構は `evalUnit_eq`（`ArcUnitEval.lean`）＋「環準同型は `1` を `1` へ送る」＋
`Scheme.ΓSpecIso` で `Γ(Spec ℂ, ⊤) ≅ ℂ` に落として `0 ≠ 1` を使う、の 3 段。 -/
theorem arcEval_unit_one_ne_zero (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    arcEval p (unitModules X) (1 : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)) ≠ 0 := by
  intro h
  have h1 := evalUnit_eq p (1 : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type))
  rw [h, map_zero] at h1
  have h2 : ((Scheme.Hom.toRingCatSheafHom p).hom.app (op ⊤)) 1 = 1 := map_one _
  have h3 : (0 : ((Spec (CommRingCat.of ℂ)).presheaf.obj
        (op (⊤ : (Spec (CommRingCat.of ℂ)).Opens)) : Type)) = 1 := h1.trans h2
  have h4 := congrArg (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom h3
  rw [map_zero, map_one] at h4
  exact zero_ne_one h4

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

def TorsorMetric.norm_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(ノルムは Γ(X, ⊤) の作用について乗法的)",
    sectionId := "genell-def-1-1-i" }

def arcEval_unit_one_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(構造層のファイバーで 1 は消えない)",
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
       "★★★計量の 𝒪_X-乗法性は 2026-08-27 に閉じた" ++
       "——TorsorMetric.norm_smul。" ++
       "★★★★残るのは **正規化した自明な算術直線束 Ō_X（|1| = 1）**だけであり、" ++
       "それには arcEval p 1 ≠ 0 と HasContMetrics X が要る") 3 ]

end ABC3.Found.Arakelov
