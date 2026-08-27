/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AMetricPic
import ABC3.Meta.Claim

/-!
# **`|s|_L̄` —— 切断の大域ノルム**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★原文の `| − |_L` そのもの

`LocalMetric.normOf` はチャート `(V, e)` を引数に取るが、
`normOf_chart_indep` によって**値はチャートに依らない**。
★したがって局所自明性を使ってチャートを選べば、

    `|s|_L̄ : X^arc → ℝ`

という**引数がチャートを含まない**関数が定まる。★★これが原文の `| − |_L` である。

## ★★★★★★★そして掛け算になる

    `|s ⊗ t|_{L̄⊗M̄}(p) = |s|_L̄(p) · |t|_M̄(p)`   （`norm_mul`）

★選択（`Classical.choice`）は使うが、**値は選択に依らない**
（`norm_eq` が任意のチャートでの値と一致することを言う）。

## ★残っている段（明示）

★★これで台帳 `arakelov-degF-finite-places` の**段 C（アルキメデス部分）の材料**が揃う
——複素点は埋め込み `σ : F ↪ ℂ` に対応するので、`Σ_σ log |s|(σ)` が書ける。
★★★残るのは段 A（アフィンの橋）・段 B（有限部分）・段 D（積公式）・段 E（加法性）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-- ★`p` の周りで自明になるチャート。 -/
structure NormChart (X : Scheme.{0}) (F : X.PresheafOfModules)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) where
  /-- 開集合。 -/
  V : X.Opens
  /-- `p` を含む。 -/
  hp : p ⁻¹ᵁ V = ⊤
  /-- 自明化。 -/
  e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)

/-- ★★**局所自明ならどの複素点にもチャートが取れる**。 -/
theorem nonempty_normChart {F : X.PresheafOfModules} (h : IsLocallyTrivial X F)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : Nonempty (NormChart X F p) := by
  obtain ⟨S, hS, htriv⟩ := h ⊤
  obtain ⟨W, g, hg, hW⟩ := hS (p.base default) trivial
  have hgt : S.arrows (homOfLE (le_top : W ≤ ⊤)) :=
    (Subsingleton.elim g (homOfLE (le_top : W ≤ ⊤))) ▸ hg
  have htop : p ⁻¹ᵁ W = ⊤ :=
    preimage_eq_top_of_mem p W (fun z => by rw [Subsingleton.elim z default]; exact hW)
  exact ⟨⟨W, htop, (htriv (homOfLE (le_top : W ≤ ⊤)) hgt).some⟩⟩

namespace AMetric

open scoped Classical

/-- ★★★★★★★★★**切断の大域ノルム** `|s|_L̄ : X^arc → ℝ`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★チャートは選ぶが、`norm_eq` が「値はどのチャートでも同じ」を言う。 -/
noncomputable def norm (L : AMetric X) (s : L.sheaf.obj (op ⊤))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : ℝ :=
  if hc : Nonempty (NormChart X L.sheaf p) then
    L.metric.normOf hc.some.V hc.some.e p hc.some.hp s
  else 0

/-- ★★★★★★**値はどのチャートで測っても同じ**——選択に依らない。 -/
theorem norm_eq (L : AMetric X) (s : L.sheaf.obj (op ⊤))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (hp : p ⁻¹ᵁ V = ⊤) :
    L.norm s p = L.metric.normOf V e p hp s := by
  have hc : Nonempty (NormChart X L.sheaf p) := ⟨⟨V, hp, e⟩⟩
  show (if hc : Nonempty (NormChart X L.sheaf p) then
      L.metric.normOf hc.some.V hc.some.e p hc.some.hp s else 0) = _
  rw [dif_pos hc]
  exact L.metric.normOf_chart_indep hc.some.e e p hc.some.hp hp s

theorem norm_nonneg (L : AMetric X) (s : L.sheaf.obj (op ⊤))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : 0 ≤ L.norm s p := by
  obtain ⟨c⟩ := nonempty_normChart L.triv p
  rw [norm_eq L s p c.V c.e c.hp]
  exact L.metric.normOf_nonneg c.V c.e p c.hp s

/-- ★★★★★★★★★**`|s ⊗ t| = |s| · |t|`**——原文の「テンソル積」の中身。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
theorem norm_mul (L M : AMetric X) (s : L.sheaf.obj (op ⊤)) (t : M.sheaf.obj (op ⊤))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    (L * M).norm (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t) p
      = L.norm s p * M.norm t p := by
  obtain ⟨c⟩ := nonempty_normChart L.triv p
  obtain ⟨c'⟩ := nonempty_normChart M.triv p
  have hpW : p ⁻¹ᵁ (c.V ⊓ c'.V) = ⊤ := by
    show p ⁻¹ᵁ c.V ⊓ p ⁻¹ᵁ c'.V = ⊤
    rw [c.hp, c'.hp, inf_idem]
  set eL := trivialOfLe L.sheaf (inf_le_left : c.V ⊓ c'.V ≤ c.V) c.e with heL
  set eM := trivialOfLe M.sheaf (inf_le_right : c.V ⊓ c'.V ≤ c'.V) c'.e with heM
  rw [norm_eq (L * M) _ p (c.V ⊓ c'.V) (tensorTriv eL eM) hpW,
    norm_eq L s p (c.V ⊓ c'.V) eL hpW, norm_eq M t p (c.V ⊓ c'.V) eM hpW]
  exact normOf_tensor_metric L.triv M.triv L.metric M.metric (c.V ⊓ c'.V) eL eM p hpW s t

/-- ★★**対数は足し算になる**——Green 関数の側の形。 -/
theorem log_norm_mul (L M : AMetric X) (s : L.sheaf.obj (op ⊤)) (t : M.sheaf.obj (op ⊤))
    (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (hs : L.norm s p ≠ 0) (ht : M.norm t p ≠ 0) :
    Real.log ((L * M).norm (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t) p)
      = Real.log (L.norm s p) + Real.log (M.norm t p) := by
  rw [norm_mul, Real.log_mul hs ht]

end AMetric

/-! ### ★出典の紐付け(`.src`) -/

def AMetric.norm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(切断の大域ノルム |s|_L̄ : X^arc → ℝ)",
    sectionId := "genell-def-1-1-i" }

def AMetric.norm_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(大域ノルムの値はチャートに依らないこと)",
    sectionId := "genell-def-1-1-i" }

def AMetric.norm_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(|s ⊗ t| = |s| · |t|——大域ノルムの形)",
    sectionId := "genell-def-1-1-i" }

def AMetric.norm.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "LocalMetric.normOf_chart_indep(ノルムはチャートに依らない)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.LocalMetric.normOf_chart_indep") 3,
    .citation "[ABC3]" "normOf_tensor_metric(|s ⊗ t| = |s| · |t|)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.normOf_tensor_metric") 3,
    .implicitStep
      ("★選択(Classical.choice)は使うが、**値は選択に依らない**" ++
       "——norm_eq が任意のチャートでの値と一致することを言う") 3,
    .implicitStep
      ("★★これで台帳 arakelov-degF-finite-places の**段 C(アルキメデス部分)の材料**が揃う" ++
       "——複素点は埋め込み σ : F ↪ ℂ に対応するので Σ_σ log |s|(σ) が書ける。" ++
       "★★★残るのは段 A(アフィンの橋)・段 B(有限部分)・段 D(積公式)・段 E(加法性)") 3 ]

end ABC3.Found.Arakelov
