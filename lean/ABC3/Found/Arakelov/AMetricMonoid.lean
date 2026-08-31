/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.TensorMetric
import ABC3.Meta.Claim

/-!
# **算術直線束**とそのテンソル積（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★原文の `L = (L, |−|_L)` をそのまま型にする

    `AMetric X ≝ (可逆層 L) × (局所自明性) × (計量 |−|_L)`

★テンソル積は `LocalMetric.tensor`（`TensorMetric.lean`）で入る。
★★単位元は `Ō_X`（`structLocalMetric`）である。

## ★★★★★★これで **`|s ⊗ t| = |s| · |t|` が語彙の側で言える**

    `(L * M).normOf V (tensorTriv eL eM) p hp (s ⊗ₜ t)`
      `= L.normOf V eL p hp s * M.normOf V eM p hp t`

★2026-08-28 に `Definition 1.1` の項目全体の `.src` を下げた理由
（`TorsorMetric.base` が `Classical.choice` で群法則が計量のテンソル積を表さない）は
**本ファイルの水準では解消している**——`mul` は**本物のテンソル積**である。

## ★残っている段（明示）

★★**同型類による商**——原文の `APic(X)` は**同型類**の群である。
そこには (1) 等長同型の同値関係、(2) 積が同値関係で降りること、
(3) 結合律・単位律、(4) **逆元（双対計量）** が要る。
★★★本ファイルはその手前——**モノイダルな演算そのもの**——までである。

★★★★もう 1 つ: `Definition 1.1` の項目全体には (ii) の `deg_F`
（台帳 `arakelov-degF-finite-places`）も要る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-- ★★★★★★★★**算術直線束** `L̄ = (L, |−|_L)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
structure AMetric (X : Scheme.{0}) where
  /-- 可逆層（前層の水準）。 -/
  sheaf : X.PresheafOfModules
  /-- ★局所自明性——これが「直線束」の中身である。 -/
  triv : IsLocallyTrivial X sheaf
  /-- ★★計量 `|−|_L`。 -/
  metric : LocalMetric X sheaf

namespace AMetric

/-- ★★**単位元 `Ō_X`**——構造層に標準計量を載せたもの。 -/
noncomputable def one (X : Scheme.{0}) : AMetric X :=
  ⟨𝟙_ X.PresheafOfModules, isLocallyTrivial_unit, structLocalMetric X⟩

/-- ★★★★★★**テンソル積** `L̄ ⊗ M̄`。 -/
noncomputable def mul (L M : AMetric X) : AMetric X :=
  ⟨L.sheaf ⊗ M.sheaf, L.triv.tensor M.triv,
    LocalMetric.tensor L.triv M.triv L.metric M.metric⟩

noncomputable instance : Mul (AMetric X) := ⟨mul⟩
noncomputable instance : One (AMetric X) := ⟨one X⟩

@[simp] theorem mul_sheaf (L M : AMetric X) : (L * M).sheaf = L.sheaf ⊗ M.sheaf := rfl

@[simp] theorem one_sheaf (X : Scheme.{0}) :
    (1 : AMetric X).sheaf = 𝟙_ X.PresheafOfModules := rfl

/-- ★★★**切断のノルム**。 -/
noncomputable def normOf (L : AMetric X) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : L.sheaf.obj (op ⊤)) : ℝ :=
  L.metric.normOf V e p hp s

theorem normOf_nonneg (L : AMetric X) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : L.sheaf.obj (op ⊤)) :
    0 ≤ L.normOf V e p hp s :=
  L.metric.normOf_nonneg V e p hp s

/-- ★★★★**ノルムはチャートに依らない**。 -/
theorem normOf_chart_indep (L : AMetric X) {V V' : X.Opens}
    (e : (restrictPresheafFunctor X V).obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (e' : (restrictPresheafFunctor X V').obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V'))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (hp' : p ⁻¹ᵁ V' = ⊤)
    (s : L.sheaf.obj (op ⊤)) :
    L.normOf V e p hp s = L.normOf V' e' p hp' s :=
  L.metric.normOf_chart_indep e e' p hp hp' s

/-- ★★★★★★★★★**`|s ⊗ t|_{L̄⊗M̄} = |s|_L̄ · |t|_M̄`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★これが原文の「テンソル積」が**本物のテンソル積である**ことの中身である。 -/
theorem normOf_mul (L M : AMetric X) (V : X.Opens)
    (eL : (restrictPresheafFunctor X V).obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (eM : (restrictPresheafFunctor X V).obj M.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : L.sheaf.obj (op ⊤)) (t : M.sheaf.obj (op ⊤)) :
    (L * M).normOf V (tensorTriv eL eM) p hp
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
      = L.normOf V eL p hp s * M.normOf V eM p hp t :=
  normOf_tensor_metric L.triv M.triv L.metric M.metric V eL eM p hp s t

/-- ★★★★★★★**テンソル自明化を選ばない形**。 -/
theorem normOf_mul' (L M : AMetric X) (V : X.Opens)
    (eL : (restrictPresheafFunctor X V).obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (eM : (restrictPresheafFunctor X V).obj M.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (f : (restrictPresheafFunctor X V).obj (L * M).sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : L.sheaf.obj (op ⊤)) (t : M.sheaf.obj (op ⊤)) :
    (L * M).normOf V f p hp (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
      = L.normOf V eL p hp s * M.normOf V eM p hp t :=
  normOf_tensor_metric' L.triv M.triv L.metric M.metric V eL eM f p hp s t

/-- ★★**対数は足し算になる**——Green 関数の側の形。 -/
theorem log_normOf_mul (L M : AMetric X) (V : X.Opens)
    (eL : (restrictPresheafFunctor X V).obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (eM : (restrictPresheafFunctor X V).obj M.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : L.sheaf.obj (op ⊤)) (t : M.sheaf.obj (op ⊤))
    (hs : L.normOf V eL p hp s ≠ 0) (ht : M.normOf V eM p hp t ≠ 0) :
    Real.log ((L * M).normOf V (tensorTriv eL eM) p hp
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t))
      = Real.log (L.normOf V eL p hp s) + Real.log (M.normOf V eM p hp t) := by
  rw [normOf_mul, Real.log_mul hs ht]

/-- ★★**単位元では `|1| = 1`**——基準の自明化で見たとき。 -/
theorem one_normOf_baseTriv (X : Scheme.{0}) (V : X.Opens)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : (1 : AMetric X).sheaf.obj (op ⊤)) :
    (1 : AMetric X).normOf V (baseTriv X V) p hp s
      = ‖evalOn p V hp (trivValue (𝟙_ X.PresheafOfModules) V (baseTriv X V) s)‖ :=
  structLocalMetric_normOf_baseTriv X V p hp s

end AMetric

/-! ### ★出典の紐付け(`.src`) -/

def AMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(算術直線束 L̄ = (L, |−|_L) の型)",
    sectionId := "genell-def-1-1-i" }

def AMetric.mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(算術直線束のテンソル積——本物のテンソル積)",
    sectionId := "genell-def-1-1-i" }

def AMetric.normOf_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(|s ⊗ t| = |s| · |t|)",
    sectionId := "genell-def-1-1-i" }

def AMetric.mul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "LocalMetric.tensor(テンソル積の計量の構成)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.LocalMetric.tensor") 3,
    .citation "[ABC3]" "IsLocallyTrivial.tensor(局所自明性はテンソル積で閉じる)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.IsLocallyTrivial.tensor") 3,
    .citation "[ABC3]" "structLocalMetric(Ō_X の標準計量——単位元)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.structLocalMetric") 3,
    .implicitStep
      ("★★★★★2026-08-28 に Definition 1.1 の項目全体の .src を下げた理由" ++
       "(TorsorMetric.base が Classical.choice で群法則が計量のテンソル積を表さない)は" ++
       "**本ファイルの水準では解消している**——mul は本物のテンソル積である") 3,
    .implicitStep
      ("★★残っているのは**同型類による商**である。原文の APic(X) は同型類の群であり、" ++
       "(1) 等長同型の同値関係、(2) 積が同値関係で降りること、" ++
       "(3) 結合律・単位律、(4) 逆元(双対計量)が要る。" ++
       "★本ファイルはその手前——モノイダルな演算そのもの——までである") 3,
    .implicitStep
      ("★★★もう 1 つ: Definition 1.1 の項目全体には (ii) の deg_F" ++
       "(台帳 arakelov-degF-finite-places)も要る") 3 ]

end ABC3.Found.Arakelov
