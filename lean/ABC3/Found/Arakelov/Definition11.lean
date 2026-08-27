/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.ArithPic
import ABC3.Found.Arakelov.APicQuot
import ABC3.Meta.Claim

/-!
# [GenEll] Definition 1.1 —— **算術直線束・`APic(X)`・引き戻し**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★★到達点 —— 項目 `Definition 1.1` の組み上げ

本ファイルは残っていた最後の 1 本、**正規化した自明な算術直線束 `Ō_X`（`|1| = 1`）**を
作り、それによって原文の

> `Γ(L̄)` は `Ō_X → L̄` の射の集合

を**そのままの形で**定義する（`arithGamma`）。

### ★★原文の語と、それを作った場所

| 原文の語 | どこ |
|---|---|
| `X^arc` はコンパクト正規複素解析空間、`ι_X` つき | `ArcTopology.lean` / `ArcConjInvol.lean` / `UltraCompact.lean` |
| 直線束 `L` on `X` | `PicWitness.lean` の `picardDataWitness.Pic` |
| `L^arc` 上の hermitian 計量、`ι_X` 両立 | `ArcTorsorMetric.lean` の `TorsorMetric` ＋ `IsConjCompatible` |
| `L̄ = (L, |−|_L)` | `ArithPic.lean` の `arithPicOfMetric` |
| **射**（`|−|_L ≤ 1` の切断が `|−|_M ≤ 1` へ） | `ArithSections.lean` の `ArithHom` |
| **`Γ(L̄)` は `Ō_X → L̄` の射の集合** | ★**本ファイル**の `arithGamma` |
| テンソル積で群 `APic(X)` | `ArithPic.lean` の `arithPicGroup` ＋ `arithPicOfMetric_mul`（`rfl`） |
| (ii) 引き戻し `φ^*L̄` | `ArithPic.lean` の `arithPicPullback` ＋ `arithPicPullback_mul` |

## ★★★★★★`Ō_X` の作り方 —— `TorsorMetric` では作れない

`TorsorMetric` の `base` は `Classical.choice` なので `|1|` を直接は選べない。
★そこで `GreenMetric`（`base` を**フィールドに持つ**）を使い、
基準計量そのものを**正規なもの**として作る:

    `|v|_p ≝ ‖(ΓSpecIso ℂ)(unitFiberIso p (v))‖`

★★`unitFiberIso`（`ArcUnitEval.lean`）が構造層のファイバーを
`Γ(Spec ℂ, ⊤)` へ移し、`Scheme.ΓSpecIso` がそれを `ℂ` へ移す。

### ★★★★★3 つの公理はいずれも既存の在庫から出た

| 公理 | 根拠 |
|---|---|
| `eq_zero_iff` | 同型は単射（`hom_inv_id`） |
| `smul` | ★`Scheme.Modules.smul_Spec_def` が `rfl`——`ℂ`-作用は `ΓSpecIso` 経由 |
| `cont`（連続性） | ★★`continuous_evalGlobal`（`ArcEvalGlobal.lean`） |

★★★`|1| = 1` は `evalGlobal p 1 = 1`（環準同型）から出る。

## ★逸脱（明示）

★**射の条件を強い側で取った**。原文は「**局所的に** `|−|_L ≤ 1` の切断が
`|−|_M ≤ 1` へ移る」だが、`ArithHom` は `|φ(s)|_M ≤ |s|_L`（作用素ノルム `≤ 1`）を要求する。
★★直線束の上では両者は同値（局所的に `φ` は切断 `u` 倍で、
「`≤ 1` を保つ」⟺ `|u| ≤ 1` ⟺ `|φ(s)| ≤ |s|`）だが、**その同値は証明していない**。
★★★向きは**強い側**なので下流は弱くならない（`ArithHom.mapsTo_arithSections`）。

## ★★★★★逸脱 3 は消えた（2026-08-28）——`APic(X)` は**同型類**の群である

原文は「テンソル積で**同型類**が群 `APic(X)` をなす」と書く。
`ArithPic X = Pic(X) × Multiplicative (共役不変な連続関数)` は**対の群**であり、
差は `Γ(X, 𝒪_X^×)` の作用（単元 `u` は計量を `|u|` 倍する）である。

★**`Found/Arakelov/APicQuot.lean` がその部分群で商を取る**（`APicOf`）。
★★鍵は `evalGlobal` の共役両立
（`evalGlobal (ι_X p) g = conj (evalGlobal p g)`）で、
これは `Scheme.ΓSpecIso_naturality` を `starRingEnd ℂ` に当てるだけで出た。
★★★引き戻しも降りる（`APicOfPullback`）ので (ii) も同型類の側で取れている。

## ★★★★★★★★★訂正（2026-08-28）——項目全体の `.src` を下げた

★**当初ここに `item := "Definition 1.1"`（条なし）を置いていた。過剰な主張であった。**

原文は「**テンソル積で**同型類が群 `APic(X)` をなす」と書く。
ところが `TorsorMetric` の設計では

* 計量は `base_F · exp(-green)` で表され、`base_F` は **`Classical.choice`**
* `TorsorMetric.tensor` は **`green` を足すだけ**

なので、群法則が表す計量は `base_{G⊗H} · e^{-(g+h)}` であり、
真のテンソル積 `(base_G ⊗ base_H) · e^{-(g+h)}` と**一致しない**。

★★差は 2-コサイクル `c(L,M) = log(base_{[LM]} / (base_{[L]} ⊗ base_{[M]}))` である。
★★★したがって `ArithPic` / `APicOf` は **`APic(X)` と集合としては対応するが、
群としては基準計量を整合的に選ばない限り一致しない**。

★★★★**失われていないもの**: 対の群・同型類への商・引き戻し・
射と `Γ(L̄)`・`Ō_X`（`|1| = 1`）はすべて**真のまま**である。
欠けているのは「群法則がテンソル積である」の 1 本だけである。

★★★★★これは `Interface/Arakelov/APic.lean` の `APicData` の穴と同種であり
（`logMetric_tensor` は Green の和しか要求しない）、
`Def12Height.lean` が 2026-08-27 に `Definition 1.2` を下げたのと**同じ根**である。

★`X` に正規性・`ℤ`-固有性・`ℤ`-平坦性を**課していない**。それらは §1 の地の文の
標準仮定であって、定義そのものには要らない——要る所（`X^arc` のコンパクト性）では
`compactSpace_arc` が固有性から出す。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite ABC3.Interface.Arakelov ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-! ## ★★★★★★構造層のファイバーの正規な `ℂ`-ノルム -/

/-- ★**構造層のファイバーの正規な `ℂ`-ノルム**。 -/
noncomputable def unitArcNrm (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (v : (arcFiber p (unitModules X) : Type)) : ℝ :=
  ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((unitFiberIso p).hom.hom v)‖

/-- ★★`0` になるのは `0` のときだけ——同型は単射だからである。 -/
theorem unitArcNrm_eq_zero_iff (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (v : (arcFiber p (unitModules X) : Type)) :
    unitArcNrm p v = 0 ↔ v = 0 := by
  constructor
  · intro h
    have h1 : (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((unitFiberIso p).hom.hom v) = 0 :=
      norm_eq_zero.1 h
    have h2 : (unitFiberIso p).hom.hom v = 0 := by
      have h3 := congrArg (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom h1
      rw [map_zero] at h3
      rwa [← CommRingCat.comp_apply, (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom_inv_id] at h3
    have h4 := congrArg (unitFiberIso p).inv.hom h2
    rw [map_zero] at h4
    rwa [← ModuleCat.comp_apply, (unitFiberIso p).hom_inv_id] at h4
  · intro h
    subst h
    have e1 : (unitFiberIso p).hom.hom (0 : (arcFiber p (unitModules X) : Type)) = 0 := map_zero _
    have e2 : (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
        ((unitFiberIso p).hom.hom (0 : (arcFiber p (unitModules X) : Type))) = 0 := by
      rw [e1]; exact map_zero _
    show ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((unitFiberIso p).hom.hom (0 : (arcFiber p (unitModules X) : Type)))‖ = 0
    rw [e2, norm_zero]

/-- ★★★★★**`ΓSpecIso` は `ℂ`-作用と両立する**。

★機構は mathlib の `Scheme.Modules.smul_Spec_def`（`rfl`）——
`Γ(Spec ℂ, ⊤)` の `ℂ`-加群構造は `(ΓSpecIso ℂ).inv` 倍で入っている。

★★`rw [Scheme.Modules.smul_Spec_def]` は具体化すると instance が出ないので通らない
（`tools\lean-idioms.md`）。**前層の綴りで補題を立てて `exact` で当てる**のが直し方である。 -/
theorem gammaSpecIso_smul (c : ℂ) (y : ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).obj
        (unitModules (Spec (CommRingCat.of ℂ))) : Type)) :
    (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (c • y)
      = c * (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom y := by
  have key : ∀ (z : ((Spec (CommRingCat.of ℂ)).presheaf.obj
        (op (⊤ : (Spec (CommRingCat.of ℂ)).Opens)) : Type)),
      (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
          (((Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom c) * z)
        = c * (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom z := by
    intro z
    rw [map_mul]
    congr 1
    exact congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) c)
      (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv_hom_id
  exact key y

theorem unitArcNrm_smul (p : Spec (CommRingCat.of ℂ) ⟶ X) (c : ℂ)
    (v : (arcFiber p (unitModules X) : Type)) :
    unitArcNrm p (c • v) = ‖c‖ * unitArcNrm p v := by
  show ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((unitFiberIso p).hom.hom (c • v))‖ = _
  rw [map_smul, gammaSpecIso_smul, norm_mul]
  rfl

/-- ★★★**大域切断のノルムはその値の絶対値である**——`evalUnit_eq` による。 -/
theorem unitArcNrm_arcEval (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (s : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)) :
    unitArcNrm p (arcEval p (unitModules X) s) = ‖evalGlobal p s‖ := by
  show ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
    ((unitFiberIso p).hom.hom (arcEval p (unitModules X) s))‖ = _
  congr 1
  rw [evalUnit_eq]
  rfl

/-- ★★★★★★**構造層の正規な計量**。 -/
noncomputable def unitArcMetric (X : Scheme.{0}) : ArcMetric X (unitModules X) where
  nrm := unitArcNrm
  nonneg _ _ := norm_nonneg _
  eq_zero_iff p v := unitArcNrm_eq_zero_iff p v
  smul p c v := unitArcNrm_smul p c v

/-- ★★★連続性は `continuous_evalGlobal`（大域の正則関数は `X^arc` 上連続）から出る。 -/
noncomputable def unitContArcMetric (X : Scheme.{0}) : ContArcMetric X (unitModules X) where
  toArcMetric := unitArcMetric X
  cont s := by
    letI := arcTopology X
    have h : (fun p : Spec (CommRingCat.of ℂ) ⟶ X =>
        (unitArcMetric X).nrm p (arcEval p (unitModules X) s))
        = fun p => ‖evalGlobal p s‖ := funext fun p => unitArcNrm_arcEval p s
    rw [h]
    exact (continuous_evalGlobal s).norm

/-! ## ★★★★★★★★自明な算術直線束 `Ō_X` -/

/-- ★構造層の単位切断 `1 ∈ Γ(X, 𝒪_X)`。 -/
noncomputable def structOne (X : Scheme.{0}) :
    ((unitModules X).val.obj (op (⊤ : X.Opens)) : Type) :=
  (1 : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type))

/-- ★★★★★★★★**自明な算術直線束 `Ō_X`**——構造層と正規な計量（Green 関数は `0`）。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
noncomputable def trivialGreenMetric (X : Scheme.{0}) : GreenMetric X (unitModules X) where
  base := unitContArcMetric X
  green := fun _ => 0
  green_cont := by letI := arcTopology X; exact continuous_const

/-- ★`Ō_X` の計量は `ι_X` と両立する（Green 関数が定数 `0` だから）。 -/
theorem trivialGreenMetric_isConjCompatible :
    (trivialGreenMetric X).IsConjCompatible := fun _ => rfl

/-- ★★★★★★**`|1|_{Ō} = 1`**——正規化の中身。 -/
theorem trivialGreenMetric_norm_one (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    (trivialGreenMetric X).norm (structOne X) p = 1 := by
  show Real.exp (-(0:ℝ)) * unitArcNrm p
    (arcEval p (unitModules X) (1 : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type))) = 1
  rw [unitArcNrm_arcEval, evalGlobal_one, neg_zero, Real.exp_zero, norm_one, one_mul]

/-! ## ★★★★★★★★`Γ(L̄) ≝ Hom(Ō_X, L̄)` -/

variable {F : X.Modules}

/-- ★★★★★★★★**`Γ(L̄) ≝ Hom(Ō_X, L̄)`** —— 原文の定義そのもの。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
def arithGamma (m : GreenMetric X F) : Type :=
  ArithHom (trivialGreenMetric X) m

/-- ★★`Γ(L̄)` の元を切断と見る（`1` の像）。 -/
noncomputable def arithGammaToSection (m : GreenMetric X F) (φ : arithGamma m) :
    (F.val.obj (op ⊤) : Type) :=
  (φ.hom.val.app (op ⊤)).hom (structOne X)

/-- ★★★★★★**`Γ(L̄)` の元は単位球に入る**——`|1|_{Ō} = 1` がそのまま効く。

★★これが「原文の `Γ(L̄)`」と「`|s|_L ≤ 1` なる大域切断」を繋ぐ橋である。 -/
theorem arithGamma_apply_one_mem (m : GreenMetric X F) (φ : arithGamma m) :
    arithGammaToSection m φ ∈ m.arithSections := by
  intro p
  have h := φ.norm_le (structOne X) p
  rw [trivialGreenMetric_norm_one] at h
  exact h

/-! ### ★★★★★★★★★★★出典の紐付け(`.src`) —— **項目全体**

★★これで `Definition 1.1` の (i)(ii) が揃った。上の docstring の表がその内訳で、
逸脱は 2 つ（射の条件を強い側で取ったこと、`X` に §1 の標準仮定を課していないこと）。 -/

def definition_1_1.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(語彙は揃った——ただし群法則が計量のテンソル積を表すには基準計量の整合性が要る)",
    sectionId := "genell-def-1-1-i" }

def trivialGreenMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(自明な算術直線束 Ō_X——|1| = 1)",
    sectionId := "genell-def-1-1-i" }

def arithGamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(Γ(L̄) ≝ Hom(Ō_X, L̄))",
    sectionId := "genell-def-1-1-i" }

def gammaSpecIso_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(ΓSpecIso は ℂ-作用と両立する)",
    sectionId := "genell-def-1-1-i" }

def definition_1_1.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "picardDataWitness(スキーム上の直線束の群 Pic(X))"
      (.inProject "ABC3" "ABC3.Interface.Arakelov.picardDataWitness") 3,
    .citation "[ABC3]" "TorsorMetric / IsConjCompatible(ι_X 両立な hermitian 計量)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.TorsorMetric") 3,
    .citation "[ABC3]" "ArithPic / arithPicOfMetric / arithPicGroup(算術直線束と APic(X))"
      (.inProject "ABC3" "ABC3.Found.Arakelov.ArithPic") 3,
    .citation "[ABC3]" "arithPicPullback((ii) 引き戻し φ^*L̄)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.arithPicPullback") 3,
    .citation "[ABC3]" "ArithHom(算術直線束の射)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.ArithHom") 3,
    .citation "[ABC3]" "compactSpace_arc(固有性から X^arc がコンパクト)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.compactSpace_arc") 3,
    .citation "[mathlib]" "Scheme.Modules.smul_Spec_def(Γ(Spec R, U) の R-作用は ΓSpecIso 経由)"
      (.inMathlib "AlgebraicGeometry.Scheme.Modules.smul_Spec_def") 3,
    .implicitStep
      ("★逸脱 1: 射の条件を強い側で取った。原文は「**局所的に** |−|_L ≤ 1 の切断が " ++
       "|−|_M ≤ 1 へ移る」だが、ArithHom は |φ(s)|_M ≤ |s|_L(作用素ノルム ≤ 1)を要求する。" ++
       "★直線束の上では両者は同値だが、その同値は証明していない。" ++
       "★★向きは強い側なので下流は弱くならない(ArithHom.mapsTo_arithSections)") 3,
    .implicitStep
      ("★★逸脱 2: X に正規性・ℤ-固有性・ℤ-平坦性を課していない。" ++
       "それらは §1 の地の文の標準仮定であって定義そのものには要らず、" ++
       "要る所(X^arc のコンパクト性)では compactSpace_arc が固有性から出す") 3,
    .citation "[ABC3]" "APicOf(APic(X)は同型類の群——逸脱 3 は 2026-08-28 に消えた)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.APicOf") 3,
    .citation "[ABC3]" "APicOfPullback((ii) は同型類の側でも取れている)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.APicOfPullback") 3,
    .implicitStep
      ("★★★★計量は「基準計量 × exp(-Green)」の捻れ集合表示で持つ。" ++
       "これは Arakelov 理論の標準的な表示であり、" ++
       "arithPicOfMetric_surjective が「余計な元が無い」ことを保証する") 3 ]

end ABC3.Found.Arakelov
