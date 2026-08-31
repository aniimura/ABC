/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.ArcGreenConj
import ABC3.Meta.Claim

/-!
# Arakelov —— **`Γ(L̄)` と算術直線束の射**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★`Definition 1.1, (i)` の**残っていた 2 語**

`Definition 1.1, (i)` は算術直線束 `L̄ = (L, |−|_L)` を導入したうえで、
続けて**射**と `Γ(L̄)` を導入する:

> 射は「局所的に `|−|_L ≤ 1` の切断が `|−|_M ≤ 1` へ移る」もの。
> `Γ(L̄)` は `Ō_X → L̄` の射の集合。

★2026-08-27 の在庫確認で、`Found/Arakelov`（317 ファイル）には
**この 2 語だけが無かった**。本ファイルがそれを埋める。

| 語 | 本ファイル |
|---|---|
| `Γ(L̄)` | `GreenMetric.arithSections` ——`|s| ≤ 1` なる大域切断の集合 |
| 射 | `ArithHom` ——`|φ(s)|_M ≤ |s|_L` |

## ★★★★射は「`≤ 1` を保つ」より**強い形**で入れた

原文は「`≤ 1` の切断が `≤ 1` へ移る」と書くが、本ファイルは
**`|φ(s)|_M ≤ |s|_L`（作用素ノルム `≤ 1`）**を要求する。

★これは**逸脱である**（記録）。理由は 2 つ:
* 前者は「局所的に」が要るが、`Γ` を大域切断で書いている本ファイルでは
  局所化の段が別に要る
* 後者は前者を**含意する**（`mapsTo_arithSections`）ので、**下流は弱くならない**

★★★逆向き（`≤ 1` 保存 ⟹ 作用素ノルム `≤ 1`）は局所自明化と
スカラー倍（`ArcMetric.smul`）から出るはずだが、本ファイルには入れていない。

## ★★★★★★`X^arc` がコンパクトなら、どんな切断も計量をずらせば `Γ` に入る

原文 (GenEll p.7) の `Proposition 1.4, (ii)` の証明は

> to a section of L ⊗M over X such that |t|L⊗M ≤1 on Xarc [where we recall that

と書く。★その「そういう切断が取れる」の根拠が
`exists_scale_mem_arithSections` である——`X^arc` がコンパクトなので
`|s|` は最大値を持ち、`c ≝ log(C+1)` だけ計量をずらせばよい。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

variable {X : Scheme.{0}} {F G H : X.Modules}

/-! ## ★零切断のノルム -/

/-- ★評価は `0` を `0` へ送る。 -/
theorem arcEval_zero (p : Spec (CommRingCat.of ℂ) ⟶ X) (L : X.Modules) :
    arcEval p L 0 = 0 := map_zero _

/-- ★★零切断のノルムは `0`。 -/
theorem GreenMetric.norm_zero (m : GreenMetric X F) (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    m.norm (0 : (F.val.obj (op ⊤) : Type)) p = 0 := by
  show Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F 0) = 0
  rw [arcEval_zero, (m.base.eq_zero_iff p 0).2 rfl, mul_zero]

/-! ## ★★★★★★★★`Γ(L̄)` -/

/-- ★★★★★★★★**`Γ(L̄)`** —— 原文『`Γ(L̄)` は `Ō_X → L̄` の射の集合』の中身。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★`Ō_X → L̄` の射は「`1` の像」で決まり、その像が `|−|_L ≤ 1` を満たすことが
射であることの条件である。★★したがって集合としては
**`|s|_L ≤ 1` なる大域切断の全体**に等しい。本ファイルはそちらを定義とする。 -/
def GreenMetric.arithSections (m : GreenMetric X F) : Set (F.val.obj (op ⊤) : Type) :=
  {s | ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, m.norm s p ≤ 1}

theorem GreenMetric.mem_arithSections_iff (m : GreenMetric X F)
    (s : (F.val.obj (op ⊤) : Type)) :
    s ∈ m.arithSections ↔ ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, m.norm s p ≤ 1 := Iff.rfl

/-- ★**`Γ(L̄)` は空でない**——零切断が入る。 -/
theorem GreenMetric.zero_mem_arithSections (m : GreenMetric X F) :
    (0 : (F.val.obj (op ⊤) : Type)) ∈ m.arithSections := by
  intro p
  rw [m.norm_zero p]
  norm_num

/-! ## ★★★★★★★算術直線束の射 -/

/-- ★★★★★★★**算術直線束の射** `L̄ → M̄`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★層の射 `φ : L → M` であって、**ノルムを増やさない**もの。
★★原文の「`|−|_L ≤ 1` の切断が `|−|_M ≤ 1` へ移る」は
これから出る（`mapsTo_arithSections`）——**逸脱の向きは強い側**である。 -/
structure ArithHom (mF : GreenMetric X F) (mG : GreenMetric X G) where
  /-- 下にある層の射。 -/
  hom : F ⟶ G
  /-- ★ノルムを増やさない。 -/
  norm_le : ∀ (s : (F.val.obj (op ⊤) : Type)) (p : Spec (CommRingCat.of ℂ) ⟶ X),
    mG.norm ((hom.val.app (op ⊤)).hom s) p ≤ mF.norm s p

namespace ArithHom

/-- ★恒等射。 -/
def id (mF : GreenMetric X F) : ArithHom mF mF where
  hom := 𝟙 F
  norm_le _ _ := le_of_eq (by rfl)

/-- ★★合成。 -/
def comp {mF : GreenMetric X F} {mG : GreenMetric X G} {mH : GreenMetric X H}
    (φ : ArithHom mF mG) (ψ : ArithHom mG mH) : ArithHom mF mH where
  hom := φ.hom ≫ ψ.hom
  norm_le s p := by
    have h1 := ψ.norm_le ((φ.hom.val.app (op ⊤)).hom s) p
    have h2 := φ.norm_le s p
    have h3 : ((φ.hom ≫ ψ.hom).val.app (op ⊤)).hom s
        = (ψ.hom.val.app (op ⊤)).hom ((φ.hom.val.app (op ⊤)).hom s) := rfl
    rw [h3]
    exact le_trans h1 h2

/-- ★★★★**射は `Γ(L̄)` を `Γ(M̄)` へ写す**——これが原文の射の定義そのものである。 -/
theorem mapsTo_arithSections {mF : GreenMetric X F} {mG : GreenMetric X G}
    (φ : ArithHom mF mG) :
    Set.MapsTo (fun s => (φ.hom.val.app (op ⊤)).hom s) mF.arithSections mG.arithSections := by
  intro s hs p
  exact le_trans (φ.norm_le s p) (hs p)

end ArithHom

/-! ## ★★★★★計量をずらすと `Γ` は増える -/

/-- ★★**計量を `c ≥ 0` だけ大きくすると `Γ` は増える**。

★`Proposition 1.4, (ii)` の『ある正冪が大域切断で生成されるなら `ht ≳ 0`』の
「定数だけずらせば切断が `Γ` に入る」という段の形である。 -/
theorem GreenMetric.arithSections_mono_scale (m : GreenMetric X F) {c : ℝ} (hc : 0 ≤ c) :
    m.arithSections ⊆ (GreenMetric.scale c m).arithSections := by
  intro s hs p
  rw [GreenMetric.norm_scale]
  have h1 : Real.exp (-c) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
  have h2 : 0 ≤ m.norm s p := m.norm_nonneg s p
  calc Real.exp (-c) * m.norm s p ≤ 1 * m.norm s p :=
        mul_le_mul_of_nonneg_right h1 h2
    _ = m.norm s p := one_mul _
    _ ≤ 1 := hs p

/-- ★★★★★**どんな切断も、計量を十分ずらせば `Γ` に入る**（`X^arc` がコンパクトなら）。

原文 (GenEll p.7):
> to a section of L ⊗M over X such that |t|L⊗M ≤1 on Xarc [where we recall that

★★角括弧の中の「`X^arc` はコンパクトだった！」が、この段の**唯一の**根拠である
——`|s|` は連続なのでコンパクト空間上で最大値 `C` を持ち、
`c ≝ log(C+1)` と取れば `exp(-c)·C = C/(C+1) ≤ 1` になる。 -/
theorem GreenMetric.exists_scale_mem_arithSections
    (hcpt : @CompactSpace _ (arcTopology X)) (hne : Nonempty (Spec (CommRingCat.of ℂ) ⟶ X))
    (m : GreenMetric X F) (s : (F.val.obj (op ⊤) : Type)) :
    ∃ c : ℝ, 0 ≤ c ∧ s ∈ (GreenMetric.scale c m).arithSections := by
  letI := arcTopology X
  haveI := hcpt
  haveI := hne
  obtain ⟨p₀, -, hmax⟩ := (isCompact_univ (X := (Spec (CommRingCat.of ℂ) ⟶ X))).exists_isMaxOn
    Set.univ_nonempty (m.norm_continuous s).continuousOn
  set C := m.norm s p₀ with hC
  have hC0 : 0 ≤ C := m.norm_nonneg s p₀
  refine ⟨Real.log (C + 1), Real.log_nonneg (by linarith), fun p => ?_⟩
  rw [GreenMetric.norm_scale, ← Real.log_inv, Real.exp_log (by positivity)]
  have hp : m.norm s p ≤ C := hmax (Set.mem_univ p)
  rw [inv_mul_le_iff₀ (by linarith)]
  linarith

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` はまだ置かない。** `Definition 1.1` を条なしで閉じるには
(i) の算術直線束・射・`Γ` と (ii) の引き戻しを**1 つに組み上げる**段が要る
——素材はすべて `Found/Arakelov` にあるが、組み上げは別の仕事である。 -/

def GreenMetric.arithSections.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(Γ(L̄)——ノルムが 1 以下の大域切断)",
    sectionId := "genell-def-1-1-i" }

def ArithHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(算術直線束の射——ノルムを増やさない層の射)",
    sectionId := "genell-def-1-1-i" }

def GreenMetric.exists_scale_mem_arithSections.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (ii)(証明中の段——計量をずらせば切断は Γ に入る)",
    sectionId := "genell-prop-1-4" }

def ArithHom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "GreenMetric(計量 = 基準計量 × exp(-Green))"
      (.inProject "ABC3" "ABC3.Found.Arakelov.GreenMetric") 3,
    .citation "[ABC3]" "GreenMetric.IsConjCompatible(ι_X と両立する計量)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.GreenMetric.IsConjCompatible") 3,
    .implicitStep
      ("★★逸脱: 原文の射の条件は「**局所的に** |−|_L ≤ 1 の切断が |−|_M ≤ 1 へ移る」" ++
       "であるが、本ファイルは |φ(s)|_M ≤ |s|_L(作用素ノルム ≤ 1)を要求する。" ++
       "★前者を含意する(mapsTo_arithSections)ので下流は弱くならない。" ++
       "★★逆向きは局所自明化とスカラー倍(ArcMetric.smul)から出るはずだが入れていない") 3,
    .implicitStep
      ("★★★Definition 1.1 を条なしで閉じるには、(i) の算術直線束・射・Γ と " ++
       "(ii) の引き戻しを 1 つに組み上げる段が要る。" ++
       "素材(aPicDataWitness・GreenMetric・IsConjCompatible・本ファイル)は" ++
       "すべて Found/Arakelov にある") 3 ]

end ABC3.Found.Arakelov
