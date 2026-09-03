/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LogDiffTower
import ABC3.Found.GenEll.MinField
import ABC3.Found.GenEll.UPoint
import ABC3.Found.GenEll.CartierPullback

/-!
# LogDiffPoint —— `[GenEll] Definition 1.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-! ## ★1. 点の定義体における `log-diff` -/

/-- ★★**点の定義体における `log-diff`** —— 原文の `deg_F(δ_x)`。

原文 (GenEll p.8):
> determines a well-defined log-different function log-diffX on X(Q).

★★`p.fld` が**最小定義体であるとき**、これが原文の `log-diff_X(x)` である
（`IsMinimalFieldOfDefinition`、`MinField.lean`）。 -/
noncomputable def logDiffAt {X : Scheme.{0}} {D : X.IdealSheafData}
    (p : AlgPointOff X D) : ℝ :=
  letI := p.instField
  letI := p.instNF
  logDiffOfField p.fld

@[simp] theorem logDiffAt_algPointOff (F : Type) [Field F] [NumberField F]
    {X : Scheme.{0}} {D : X.IdealSheafData} (xF : specRingOfIntegers F ⟶ X)
    (h : pullbackIdeal F D xF ≠ 0) :
    logDiffAt (algPointOff F xF h) = logDiffOfField F := rfl

/-! ## ★2. ★★★★★底変換で不変**でない**こと -/

/-- ★★★★★**`log-diff` は底変換で不変ではない** —— 定義体を上げると増えうる。

原文 (GenEll p.8):
> Let x ∈X(F ) ⊆X(Q), where F is a minimal field of definition of x. Then

★★**これが原文の「`F` は最小定義体」という断りの理由そのもの**である。

| 量 | 底変換 | 最小定義体の指定 |
|---|---|---|
| 高さ `ht` | ★**不変**（`htOff_baseChange`） | 要らない |
| `log-diff` | ★★**増えうる**（本定理） | ★**要る** |

★`Definition 1.2`（高さ）が同じ断りを必要としないこととの対比が、この構造を示している。 -/
theorem logDiffAt_le_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    {X : Scheme.{0}} {D : X.IdealSheafData} (xF : specRingOfIntegers F ⟶ X)
    (hJ : pullbackIdeal F D xF ≠ 0)
    (hJ' : pullbackIdeal K D
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) ≠ 0) :
    logDiffAt (algPointOff F xF hJ)
      ≤ logDiffAt (algPointOff K
          (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) hJ') := by
  simpa using logDiffOfField_le F K

/-! ## ★3. 最小定義体で測る -/

/-- ★★★**原文の `log-diff_X(x)`** —— 最小定義体で測った値。

原文 (GenEll p.8):
> Let x ∈X(F ) ⊆X(Q), where F is a minimal field of definition of x. Then

★`F` が最小定義体（`IsMinimalFieldOfDefinition`、`MinField.lean`）であるとき、
`logDiffAt` の値が原文の `deg_F(δ_x)` である、ということを型で明示する。 -/
theorem logDiffAt_of_minimal (F : Type) [Field F] [NumberField F]
    {X : Scheme.{0}} {D : X.IdealSheafData} (xF : specRingOfIntegers F ⟶ X)
    (h : pullbackIdeal F D xF ≠ 0)
    (xgen : Spec (CommRingCat.of F) ⟶ X)
    (_hmin : IsMinimalFieldOfDefinition F xgen) :
    logDiffAt (algPointOff F xF h) = logDiffOfField F := rfl

/-! ## ★4. `Definition 1.5, (iv)` —— `log-cond` も点の関数にする -/

/-- ★★**点の定義体における `log-cond`** —— 原文の `deg_F(f^D_x)`。

原文 (GenEll p.8):
> determines a well-defined log-conductor function log-condD on UX(Q) ⊆X(Q).

★`AlgPointOff X D` は「`x` が `D` を通らない」（`off`）を持つので、
原文の `U_X(ℚ̄)` の点そのものである。 -/
noncomputable def logCondAt {X : Scheme.{0}} {D : X.IdealSheafData}
    (p : AlgPointOff X D) : ℝ :=
  letI := p.instField
  letI := p.instNF
  logCond p.fld D p.map

/-- ★**引き戻しは合成に沿って分解する**（`log-cond` の底変換の第 1 段）。

原文 (GenEll p.8):
> allows one to pull-back the divisor D to Spec(OF ) so as to obtain

★★中身は mathlib の `Scheme.IdealSheafData.comap_comp` 1 本である。

★★★**次の段が mathlib に無い**（2026-08-27 実測）——
「アフィンでイデアル層の `comap` が `Ideal.map`（拡大）に対応する」。
`comap_ideal` / `ideal_comap` / `equivOfIsAffine_comap` / `map_comap` /
`comap_symm_equivOfIsAffine` のいずれの名前でも見つからない。
★これが取れれば `pullbackIdeal K D (bc ≫ xF) = (pullbackIdeal F D xF).map (𝓞_F → 𝓞_K)`
が出て、`log-cond` の底変換が測れる。 -/
theorem comap_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] {X : Scheme.{0}} (D : X.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) :
    D.comap (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF)
      = (D.comap xF).comap (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) :=
  Scheme.IdealSheafData.comap_comp D _ _

/-! ### ★出典の紐付け(`.src`) -/

def comap_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(引き戻しは合成に沿って分解する)",
    sectionId := "genell-def-1-5" }

def comap_baseChange.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.IdealSheafData.comap_comp"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.comap_comp") 8,
    .implicitStep
      ("★★☆次の段が mathlib に無い(2026-08-27 実測): 「アフィンでイデアル層の comap が " ++
       "Ideal.map(拡大)に対応する」。comap_ideal / ideal_comap / equivOfIsAffine_comap / " ++
       "map_comap / comap_symm_equivOfIsAffine のいずれの名前でも見つからない。" ++
       "★これが取れれば pullbackIdeal K D (bc ≫ xF) = (pullbackIdeal F D xF).map で " ++
       "log-cond の底変換が測れる(Σ_{𝔓|𝔭} f_𝔓 ≤ [K:F] により減りうる)") 8 ]

def logCondAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(点の定義体における log-cond)",
    sectionId := "genell-def-1-5" }

def logCondAt_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(負の対照——空因子なら log-cond = 0)",
    sectionId := "genell-def-1-5" }

def logCondAt.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "logCond(引き戻した因子の被約化の正規化次数)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logCond") 8,
    .implicitStep
      ("☆残る段: log-cond_D を U_X(ℚ̄) の関数として完成させるには、" ++
       "log-diff と同じく最小定義体の一意性が要る。★★なお log-cond は log-diff と" ++
       "逆向きに動くはずである —— 分岐があると (D_x)_red の次数は定義体を上げると" ++
       "減りうる(Σ_{𝔓|𝔭} f_𝔓 ≤ [K:F] で、等号は e = 1 のときだけ)。" ++
       "これは Ideal.sum_ramification_inertia で測れる(FrdI の ArithFundamental.lean と同じ機械)") 8 ]

def logDiffAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(点の定義体における log-diff)",
    sectionId := "genell-def-1-5" }

def logDiffAt_le_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(log-diff は底変換で不変でない——最小定義体を要求する理由)",
    sectionId := "genell-def-1-5" }

def logDiffAt_le_baseChange.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "logDiffOfField_le(差積は定義体を上げると増える)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiffOfField_le") 8,
    .implicitStep
      ("★★これが原文「where F is a minimal field of definition of x」の理由である。" ++
       "高さ ht は底変換で不変(htOff_baseChange)なので同じ断りを要しない —— " ++
       "その対比が Definition 1.2 と Definition 1.5 の構造の違いを示す") 8 ]

def logDiffAt_of_minimal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(最小定義体で測った値が原文の log-diff_X(x))",
    sectionId := "genell-def-1-5" }

def logDiffAt_of_minimal.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "IsMinimalFieldOfDefinition(最小定義体の述語)"
      (.inProject "ABC3" "ABC3.Found.GenEll.IsMinimalFieldOfDefinition") 8,
    .implicitStep
      ("☆残る段: log-diff_X を X(ℚ̄) の関数として完成させるには" ++
       "「同じ点の 2 つの最小定義体は同じ log-diff を与える」が要る。" ++
       "中身は「同型な数体は同じ差積次数をもつ」であって、本ファイルは取っていない") 8 ]

end ABC3.Found.GenEll
