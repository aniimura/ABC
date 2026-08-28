/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DivisorOfSection
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★段 E0 —— `div(s)` は近似でなく等号である（大域自明の場合）（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★これは何か —— 段 E0 の残りを 2 つに割る

`§9-809`（`DivisorOfSection.lean`）は `div(s) ≝ ofIdeals (divIdeal)` と置き、
`ofIdeals` が `≤` しか与えないので

    `(divisorOfSection M s).ideal U = divIdeal M s U`

を**未証明のまま残していた**（`div(s)` は上からの近似）。

★本ファイルはそれを 2 段に割る:

| 段 | 内容 | 状態 |
|---|---|---|
| (a) | **証人があれば等号**（`divisorOfSection_ideal_eq_of_witness`） | ★★済 |
| (b) | 証人の構成（一般の局所自明な `M`） | ★開 |

★★そして **`M` が大域的に自明な場合**は (b) が直ちに済む——
証人は `ofIdealTop (span {trivValue ⊤ e s})` である。

## ★★★測定の記録 —— 何が本当の障害だったか

★★★★障害は「`ofIdeals` が `≤` しか与えない」ことではなく、
**`divIdeal` が自明化の無いアフィン開で `⊤` になること**である。
★大域的に自明なら**すべての**アフィン開に自明化が降りるので `⊤` の場合が消え、
`divIdeal` 自身が `IdealSheafData` の族になる（2026-08-28 実測）。

★★したがって一般の場合に要るのは「自明化の無い開でも正しい値を与える定義」であり、
それは `⨅_{V ⊆ U 自明化する} (制限の comap)` の形になる
——`map_ideal_basicOpen` の `⊆` 向きに**分母を払う議論**（`§9-826`・`§9-831`）が要る。

## ★これで足りる場面（明示）

★★★★★消費側（段 E3 の `X_{s_i}` 上の議論）では **`M` は `X_{s_i}` の上で自明**である
——`s_i` が消えない開だから。したがって本ファイルの大域自明版でちょうど足りる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★(a) 証人があれば等号 -/

/-- ★★★★★★★★**証人があれば `div(s)` は等号である**。

★`J.ideal ≤ divIdeal` なら `J ≤ ofIdeals divIdeal`（`le_ofIdeals_iff`）なので、
`divIdeal U ≤ J.ideal U ≤ (ofIdeals divIdeal).ideal U ≤ divIdeal U` で挟める。 -/
theorem divisorOfSection_ideal_eq_of_witness (M : X.PresheafOfModules)
    (s : M.obj (op ⊤)) (U : X.affineOpens)
    (J : X.IdealSheafData) (hJle : J.ideal ≤ divIdeal M s)
    (hJU : divIdeal M s U ≤ J.ideal U) :
    (divisorOfSection M s).ideal U = divIdeal M s U := by
  refine le_antisymm (divisorOfSection_ideal_le M s U) (le_trans hJU ?_)
  exact (Scheme.IdealSheafData.le_ofIdeals_iff.2 hJle) U

/-! ## ★★★★★★★★★★(b) 大域的に自明なら証人がある -/

/-- ★★★★★★★★★★**`M` が大域的に自明なら `div(s)` は等号である**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★証人は `ofIdealTop (span {trivValue ⊤ e s})` である。
★★大域的に自明なら**すべての**アフィン開に自明化が降りるので
`divIdeal` の `⊤` の場合が消え、`divIdeal` 自身が `IdealSheafData` の族になる。
★★★消費側（段 E3 の `X_{s_i}` 上の議論）ではこれでちょうど足りる
——`s_i` が消えない開だから `M` はそこで自明である。 -/
theorem divisorOfSection_ideal_eq_of_globalTriv (M : X.PresheafOfModules)
    (s : M.obj (op ⊤))
    (e : (restrictPresheafFunctor X ⊤).obj M ≅ 𝟙_ (PresheafModulesOn X ⊤))
    (U : X.affineOpens) :
    (divisorOfSection M s).ideal U = divIdeal M s U := by
  have key : ∀ V : X.affineOpens,
      (Scheme.IdealSheafData.ofIdealTop (Ideal.span {trivValue M ⊤ e s})).ideal V
        = divIdeal M s V := by
    intro V
    rw [Scheme.IdealSheafData.ofIdealTop_ideal,
      divIdeal_eq M s V (trivialOfLe M (le_top : V.1 ≤ ⊤) e), Ideal.map_span,
      trivValue_restrict M (le_top : V.1 ≤ ⊤) e s, Set.image_singleton]
  exact divisorOfSection_ideal_eq_of_witness M s U
    (Scheme.IdealSheafData.ofIdealTop (Ideal.span {trivValue M ⊤ e s}))
    (fun V => le_of_eq (key V)) (le_of_eq (key U).symm)

/-- ★★**大域自明なら `div(s)` は `ofIdealTop` そのものである**。 -/
theorem divisorOfSection_eq_ofIdealTop (M : X.PresheafOfModules) (s : M.obj (op ⊤))
    (e : (restrictPresheafFunctor X ⊤).obj M ≅ 𝟙_ (PresheafModulesOn X ⊤)) :
    divisorOfSection M s
      = Scheme.IdealSheafData.ofIdealTop (Ideal.span {trivValue M ⊤ e s}) := by
  refine Scheme.IdealSheafData.ext (funext fun V => ?_)
  rw [divisorOfSection_ideal_eq_of_globalTriv M s e V,
    Scheme.IdealSheafData.ofIdealTop_ideal,
    divIdeal_eq M s V (trivialOfLe M (le_top : V.1 ≤ ⊤) e), Ideal.map_span,
    trivValue_restrict M (le_top : V.1 ≤ ⊤) e s, Set.image_singleton]

/-! ## ★出典の紐付け(`.src`) -/

def divisorOfSection_ideal_eq_of_witness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2(証人があれば div(s) は等号である)",
    sectionId := "genell-def-1-2" }

def divisorOfSection_ideal_eq_of_globalTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2(M が大域的に自明なら div(s) は等号である)",
    sectionId := "genell-def-1-2" }

def divisorOfSection_eq_ofIdealTop.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2(大域自明なら div(s) は ofIdealTop そのもの)",
    sectionId := "genell-def-1-2" }

def divisorOfSection_ideal_eq_of_globalTriv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "divisorOfSection / divIdeal(段 E0 の入口、§9-809)"
      (.inProject "ABC3" "ABC3.Found.GenEll.divisorOfSection") 2,
    .citation "[mathlib]" "IdealSheafData.le_ofIdeals_iff / ofIdealTop_ideal"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.le_ofIdeals_iff") 2,
    .implicitStep
      ("★★★★測定: 障害は「ofIdeals が ≤ しか与えない」ことではなく、" ++
       "**divIdeal が自明化の無いアフィン開で ⊤ になること**である。" ++
       "★大域的に自明なら**すべての**アフィン開に自明化が降りるので ⊤ の場合が消え、" ++
       "divIdeal 自身が IdealSheafData の族になる(2026-08-28 実測)") 3,
    .implicitStep
      ("★★一般の場合に要るのは「自明化の無い開でも正しい値を与える定義」であり、" ++
       "それは ⨅_{V ⊆ U 自明化する} (制限の comap) の形になる" ++
       "——map_ideal_basicOpen の ⊆ 向きに**分母を払う議論**(§9-826・§9-831)が要る") 5,
    .implicitStep
      ("★★★消費側(段 E3 の X_{s_i} 上の議論)では M は X_{s_i} の上で**自明**である" ++
       "(s_i が消えない開だから)。したがって本ファイルの大域自明版でちょうど足りる") 3 ]

end ABC3.Found.GenEll
