/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Covolume
import ABC3.Found.GenEll.IsogenyCovolume
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★段 3 を `PeriodPair` の言葉で（`Found`、無条件）

原典: P. Deligne, *Preuve des conjectures de Tate et de Shafarevitch (d'après G. Faltings)*,
Séminaire Bourbaki 616 (1983/84)、物理 p.14。

原文 (DelSB616 p.14):
> F de nombres premiers tels que, pour u : A - B une

## ★★★★★★★★★★★★★★★★これは何か

`Found/GenEll/IsogenyCovolume.lean`（`§9-1017`）は
`[DelSB616] Théorème 2.4` の段 3（アルキメデスの共体積は指数倍）を
**mathlib の `ZLattice.covolume`** で取った。

★★本ファイルは**同じことを本プロジェクトの `covol : PeriodPair → ℝ` で取る**
——`archInv L = ‖Δ_L‖·covol(L)⁶`（`§9-353`）が `PeriodPair` の言葉で書かれているので、
そちらに直結する形が要る。

★★★しかも**測度論を経由しない**: `covol` は行列式の絶対値
`|ω₁.re·ω₂.im − ω₁.im·ω₂.re|` として定義されているので、
基底変換行列 `M ∈ M₂(ℤ)` に対する `covol = |det M|·covol′` は**純粋な代数**である。

## ★機構

`Λ ⊆ Λ′` が指数 `l` なら、`Λ` の基底は `Λ′` の基底の**整数係数の一次結合**で書け、
その変換行列の行列式の絶対値が `l` である。したがって

    `covol(Λ) = |det M|·covol(Λ′) = l·covol(Λ′)`

★これを対数に通すと、Faltings 高さのアルキメデス局所項 `−(1/2)·log covol` は
**ちょうど `(1/2)·log(l)` 増える**。

## ★★これがどこに効くか

`§9-1018`（第 573）の測定:

    `12·htFaltOf = −12·log(2π) − (6/d)·Σ_σ log covol_min(E^σ)`

（積公式で `deg∞` が打ち消し合う）。★したがって同種写像評価は
**極小モデルの周期格子の共体積の比較**に完全に還元され、
その解析側が本ファイルである。

☆残るのは**極小モデルへのスケーリング `u_σ` の制御**
——`[DelSB616]` の段 1・段 2（Néron モデルと `w(H)`）である。
-/

namespace ABC3.Found.GenEll

open Complex

/-! ## ★★★★★★★★★★基底変換と共体積 -/

/-- ★★★★★★★★★★**`covol` は基底変換行列の行列式の絶対値だけ動く**。

原文 (DelSB616 p.14):
> F de nombres premiers tels que, pour u : A - B une

★`covol L = |ω₁.re·ω₂.im − ω₁.im·ω₂.re|` の定義に代入して展開するだけである
——★★**測度論は要らない**。 -/
theorem covol_of_int_change (L L' : PeriodPair) (a b c d : ℤ)
    (h₁ : L.ω₁ = (a:ℂ) * L'.ω₁ + (b:ℂ) * L'.ω₂)
    (h₂ : L.ω₂ = (c:ℂ) * L'.ω₁ + (d:ℂ) * L'.ω₂) :
    covol L = |((a * d - b * c : ℤ) : ℝ)| * covol L' := by
  have hre₁ : L.ω₁.re = (a:ℝ) * L'.ω₁.re + (b:ℝ) * L'.ω₂.re := by rw [h₁]; simp
  have him₁ : L.ω₁.im = (a:ℝ) * L'.ω₁.im + (b:ℝ) * L'.ω₂.im := by rw [h₁]; simp
  have hre₂ : L.ω₂.re = (c:ℝ) * L'.ω₁.re + (d:ℝ) * L'.ω₂.re := by rw [h₂]; simp
  have him₂ : L.ω₂.im = (c:ℝ) * L'.ω₁.im + (d:ℝ) * L'.ω₂.im := by rw [h₂]; simp
  rw [covol, covol, hre₁, him₁, hre₂, him₂]
  rw [show ((a:ℝ) * L'.ω₁.re + (b:ℝ) * L'.ω₂.re) * ((c:ℝ) * L'.ω₁.im + (d:ℝ) * L'.ω₂.im)
        - ((a:ℝ) * L'.ω₁.im + (b:ℝ) * L'.ω₂.im) * ((c:ℝ) * L'.ω₁.re + (d:ℝ) * L'.ω₂.re)
      = ((a:ℝ) * d - (b:ℝ) * c) * (L'.ω₁.re * L'.ω₂.im - L'.ω₁.im * L'.ω₂.re) by ring]
  rw [abs_mul]
  congr 1
  push_cast
  ring_nf

/-- ★★★★★★**指数 `l` の部分格子なら `covol` は `l` 倍**。 -/
theorem covol_eq_index_mul_pair (L L' : PeriodPair) (a b c d : ℤ) (l : ℕ)
    (h₁ : L.ω₁ = (a:ℂ) * L'.ω₁ + (b:ℂ) * L'.ω₂)
    (h₂ : L.ω₂ = (c:ℂ) * L'.ω₁ + (d:ℂ) * L'.ω₂)
    (hdet : (a * d - b * c).natAbs = l) :
    covol L = (l : ℝ) * covol L' := by
  rw [covol_of_int_change L L' a b c d h₁ h₂]
  congr 1
  rw [← hdet, ← Int.cast_abs, Int.abs_eq_natAbs]
  simp

/-! ## ★★★★★★★★★★★★★★高さに効く形 -/

/-- ★★★★★★★★★★★★★★**アルキメデス局所項はちょうど `(1/2)·log(l)` 増える**。

原文 (DelSB616 p.14):
> F de nombres premiers tels que, pour u : A - B une

★Faltings 高さのアルキメデス局所項は `−(1/2)·log covol` である
（`§9-1018` の測定: `12·htFaltOf = −12log(2π) − (6/d)Σ_σ log covol_min`）。
★★これが `Skeleton/GenEll/IsogenyHeight.lean` の仮説 `hwarch` を満たす。 -/
theorem archContrib_of_int_change (L L' : PeriodPair) (a b c d : ℤ) (l : ℕ) (hl : 1 ≤ l)
    (h₁ : L.ω₁ = (a:ℂ) * L'.ω₁ + (b:ℂ) * L'.ω₂)
    (h₂ : L.ω₂ = (c:ℂ) * L'.ω₁ + (d:ℂ) * L'.ω₂)
    (hdet : (a * d - b * c).natAbs = l) :
    -(1/2) * Real.log (covol L') = -(1/2) * Real.log (covol L) + (1/2) * Real.log l := by
  have hcov := covol_eq_index_mul_pair L L' a b c d l h₁ h₂ hdet
  have hl0 : (0:ℝ) < (l:ℝ) := by exact_mod_cast hl
  have hpos' : (0:ℝ) < covol L' := covol_pos L'
  rw [hcov, Real.log_mul (ne_of_gt hl0) (ne_of_gt hpos')]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def covol_of_int_change.src : ABC3.Meta.Source :=
  { paper := "DelSB616", pdfPage := 14,
    item := "Théorème 2.4(段 3——covol は基底変換の行列式だけ動く。PeriodPair の言葉で)",
    sectionId := "delsb616-thm-2-4" }

def archContrib_of_int_change.src : ABC3.Meta.Source :=
  { paper := "DelSB616", pdfPage := 14,
    item := "Théorème 2.4(段 3——アルキメデス局所項は (1/2)log(l) 増える。PeriodPair の言葉で。★無条件)",
    sectionId := "delsb616-thm-2-4" }

def archContrib_of_int_change.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "covol(PeriodPair の共体積、§9-353)"
      (.inProject "ABC3" "ABC3.Found.GenEll.covol") 2,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29、第 574): 段 3 を**本プロジェクトの " ++
       "covol : PeriodPair → ℝ の言葉で**取った。" ++
       "★archInv L = ‖Δ_L‖·covol(L)⁶(§9-353)が PeriodPair の言葉なので、" ++
       "mathlib の ZLattice.covolume 版(§9-1017)より直結する。" ++
       "★★しかも**測度論を経由しない**——covol は行列式の絶対値として定義されているので、" ++
       "基底変換に対する挙動は純粋な代数である") 9,
    .implicitStep
      ("☆入力 `h₁`・`h₂`・`hdet` は「Λ ⊆ Λ′ が指数 l」を基底の言葉で書いたものである。" ++
       "★Λ の基底が Λ′ の基底の整数係数一次結合で書けること、" ++
       "およびその行列式の絶対値が指数に等しいことは、" ++
       "有限生成自由アーベル群の Smith 標準形から出る(mathlib に " ++
       "Matrix.SmithNormalForm 系がある)。本ファイルは基底変換を仮定として受けている") 8,
    .folklore
      ("☆残るのは段 1(原文 2.1)と段 2(原文 2.2 (b)(c))——" ++
       "極小モデルへのスケーリング u_σ の制御である。" ++
       "★Néron モデルと w(H) の理論を要し、mathlib に無い") 14 ]

end ABC3.Found.GenEll
