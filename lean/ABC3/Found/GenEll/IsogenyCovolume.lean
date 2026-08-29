/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Algebra.Module.ZLattice.Covolume
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★[DelSB616] Théorème 2.4 の**段 3**（`Found`、無条件）

原典: P. Deligne, *Preuve des conjectures de Tate et de Shafarevitch (d'après G. Faltings)*,
Séminaire Bourbaki 616 (1983/84)、物理 p.14。

原文 (DelSB616 p.14):
> F de nombres premiers tels que, pour u : A - B une

## ★★★★★★★★★★★★★★これは何か —— 3 つの葉のうち 1 つが閉じた

`Skeleton/GenEll/IsogenyHeight.lean`（`§9-1016`、第 571）は
`[DelSB616] Théorème 2.4` を 3 段に分けた:

| 段 | 内容 | 原文 | 状態 |
|---|---|---|---|
| 1 | 完全列に付随する `ω` の自然な同型 | `2.1` | ☆未 |
| 2 | `w(H)/𝒪` は `#H` で消える（有限素点側） | `2.2 (b)(c)` | ☆未 |
| **3** | ★**アルキメデスのエルミート構造は `#H` 倍** | `2.2 (a)` | ★★★**本ファイルで閉じた** |

## ★機構 —— mathlib に既にあった

原文 `2.2 (a)` は、楕円曲線の言葉では**周期格子の共体積**の関係である:

    `Λ ⊆ Λ′` が指数 `l` なら  `covol(Λ) = l · covol(Λ′)`

★これは mathlib の **`ZLattice.covolume_div_covolume_eq_relIndex'`** そのものである
（2026-08-29 に `#check` で発見）。★★台帳は「アルキメデスの周期格子の共体積の比較」を
残っている葉として数えていたが、**mathlib に既にあった**。

## ★★高さに効く形

Faltings 高さのアルキメデス局所項は `−(1/2)·log covol` である。したがって

    `−(1/2)·log covol(Λ′) = −(1/2)·log covol(Λ) + (1/2)·log(l)`

——★**同種写像でアルキメデス側はちょうど `(1/2)·log(l)` だけ増える**。
★★`Skeleton` の仮説 `hwarch : warch ≤ log(l)` はこれで満たされる
（`(1/2)log l ≤ log l`、`l ≥ 1`）。

## ☆残る 2 つの葉

☆段 1・段 2 は Néron モデルと `w(H)` の理論を要する（mathlib に無い）。
★段 2 の中身「`E/H` を Weierstrass 曲線として作る」には Vélu の公式が要る。
-/

namespace ABC3.Found.GenEll

open MeasureTheory

/-! ## ★★★★★★★★★★共体積と指数 -/

/-- ★★★★★★★★★★**`Λ ⊆ Λ′` が指数 `l` なら `covol(Λ) = l·covol(Λ′)`**。

原文 (DelSB616 p.14):
> F de nombres premiers tels que, pour u : A - B une

★原文 `2.2 (a)`「`w(H) ⊗ ℂ` のエルミート構造は標準のものの `#H` 倍」の、
楕円曲線の言葉での中身である（`E = ℂ/Λ`、`E′ = ℂ/Λ′`、`H = Λ′/Λ`）。
★★mathlib の `ZLattice.covolume_div_covolume_eq_relIndex'` から出る。 -/
theorem covolume_eq_index_mul (Λ Λ' : Submodule ℤ ℂ) [DiscreteTopology Λ] [IsZLattice ℝ Λ]
    [DiscreteTopology Λ'] [IsZLattice ℝ Λ'] (hle : Λ ≤ Λ') (l : ℕ)
    (hidx : Λ.toAddSubgroup.relIndex Λ'.toAddSubgroup = l) :
    ZLattice.covolume Λ volume = (l : ℝ) * ZLattice.covolume Λ' volume := by
  have h := ZLattice.covolume_div_covolume_eq_relIndex' Λ Λ' hle
  rw [hidx] at h
  have hne : ZLattice.covolume Λ' volume ≠ 0 := ZLattice.covolume_ne_zero Λ' volume
  field_simp at h
  linarith [h]

/-! ## ★★★★★★★★★★★★★★高さに効く形 -/

/-- ★★★★★★★★★★★★★★**同種写像でアルキメデス側はちょうど `(1/2)·log(l)` 増える**。

原文 (DelSB616 p.14):
> F de nombres premiers tels que, pour u : A - B une

★Faltings 高さのアルキメデス局所項は `−(1/2)·log covol` である。
★★これが `Skeleton/GenEll/IsogenyHeight.lean` の仮説 `hwarch` を満たす
（`(1/2)log l ≤ log l`——`covolume_archContrib_le` を見よ）。 -/
theorem archContrib_of_index (Λ Λ' : Submodule ℤ ℂ) [DiscreteTopology Λ] [IsZLattice ℝ Λ]
    [DiscreteTopology Λ'] [IsZLattice ℝ Λ'] (hle : Λ ≤ Λ') (l : ℕ) (hl : 1 ≤ l)
    (hidx : Λ.toAddSubgroup.relIndex Λ'.toAddSubgroup = l) :
    -(1/2) * Real.log (ZLattice.covolume Λ' volume)
      = -(1/2) * Real.log (ZLattice.covolume Λ volume) + (1/2) * Real.log l := by
  have hcov := covolume_eq_index_mul Λ Λ' hle l hidx
  have hl0 : (0:ℝ) < (l:ℝ) := by exact_mod_cast hl
  have hpos' : (0:ℝ) < ZLattice.covolume Λ' volume := ZLattice.covolume_pos Λ' volume
  rw [hcov, Real.log_mul (ne_of_gt hl0) (ne_of_gt hpos')]
  ring

/-- ★★★★★**`Skeleton` の仮説 `hwarch` の形**——`(1/2)·log(l) ≤ log(l)`。 -/
theorem covolume_archContrib_le (l : ℕ) (hl : 1 ≤ l) :
    (1/2) * Real.log l ≤ Real.log l := by
  have hl1 : (1:ℝ) ≤ (l:ℝ) := by exact_mod_cast hl
  have := Real.log_nonneg hl1
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def covolume_eq_index_mul.src : ABC3.Meta.Source :=
  { paper := "DelSB616", pdfPage := 14,
    item := "Théorème 2.4(段 3——共体積は指数倍。原文 2.2 (a))",
    sectionId := "delsb616-thm-2-4" }

def archContrib_of_index.src : ABC3.Meta.Source :=
  { paper := "DelSB616", pdfPage := 14,
    item := "Théorème 2.4(段 3——アルキメデス側はちょうど (1/2)log(l) 増える。★無条件)",
    sectionId := "delsb616-thm-2-4" }

def archContrib_of_index.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]"
      ("ZLattice.covolume_div_covolume_eq_relIndex'" ++
       "(部分格子の共体積の比は相対指数)")
      (.inMathlib "ZLattice.covolume_div_covolume_eq_relIndex'") 2,
    .implicitStep
      ("★★★★★★★★★★到達点(2026-08-29、第 572): " ++
       "[DelSB616] Théorème 2.4 の**3 つの葉のうち段 3 が閉じた**。" ++
       "★台帳は「アルキメデスの周期格子の共体積の比較」を残っている葉として数えていたが、" ++
       "**mathlib に既にあった**(ZLattice.covolume_div_covolume_eq_relIndex')。" ++
       "★★これで Skeleton/GenEll/IsogenyHeight.lean の仮説 hwarch は満たされる" ++
       "——同種写像でアルキメデス側はちょうど (1/2)log(l) だけ増える") 9,
    .folklore
      ("☆残る葉は段 1(原文 2.1、完全列に付随する ω の同型)と" ++
       "段 2(原文 2.2 (b)(c)、w(H)/𝒪 は #H で消える)である。" ++
       "★どちらも Néron モデルと w(H) の理論を要し、mathlib に無い。" ++
       "★★段 2 の中身「E/H を Weierstrass 曲線として作る」には Vélu の公式が要る") 14,
    .implicitStep
      ("☆本ファイルは**格子の水準**で述べている。" ++
       "★楕円曲線 E = ℂ/Λ と部分群 H = Λ′/Λ の対応を付けるには" ++
       "Found/GenEll/Covolume.lean・LatticeFromInvariants.lean 側との接続が要る" ++
       "——それは段 1・段 2 と同じ束に属する") 10 ]

end ABC3.Found.GenEll
