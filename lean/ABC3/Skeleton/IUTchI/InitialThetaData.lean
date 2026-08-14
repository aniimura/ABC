import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure

/-!
# [IUTchI] Definition 3.1 — initial Θ-data(**部分的な**転写)

原典: S. Mochizuki, *Inter-universal Teichmüller Theory I: Construction of Hodge
Theaters* [IUTchI]、物理 p.61(全 186 ページ、May 2020。**400 dpi 目視確認 2026-08-15**)。

原文 (IUTchI p.61):
> We shall refer to as initial Θ-data any collection of data

原文 (IUTchI p.61):
> that satisfies the following conditions:

## なぜここを張るか —— 依存グラフの葉だったから

[IUTchIII] Corollary 3.12 の `.needs` は `[IUTchI] Definition 3.1, (b)` を
**辺**として指しているが、その先にスケルトンが無かった。
`tools/check.mjs` の「依存グラフ」印字で **未展開の葉** として現れる状態である。
ここはその葉を1つ展開したもの。

## ★何を転写し、何を転写していないか(目標を下げた旨の明示)

**転写した**のは条件 (a) だけである:

原文 (IUTchI p.61):
> is a number field such that

すなわち「`F` は数体」「`√−1 ∈ F`」「`F̄` は `F` の代数閉包」「`G_F = Gal(F̄/F)`」。
これらは mathlib の語彙(`NumberField`・`IsAlgClosure`)で書ける。

**転写していない**のは条件 (b) 以降である:

原文 (IUTchI p.61):
> is a once-punctured elliptic curve [i.e., a hyperbolic curve of type

`once-punctured elliptic curve`(1点抜き楕円曲線)・`stable reduction`(安定還元)・
`hyperbolic orbicurve`・`field of moduli`・`maximal solvable extension` は、
2026-08-15 の実測で mathlib に無い(探索範囲は `InitialThetaDataFieldPart.needs` に書いた)。
したがって **`InitialThetaDataFieldPart` は Definition 3.1 ではない。**
名前もそう付けてある——(a) の体の部分だけである。

## ★記法の注意(この論文の risk の実例)

条件 (a) の `Gal(F̄/F)` を `pdftotext` は `Gal(F/F)` と出力する
(オーバーバーが落ちる)。**自明群と読める**ので、`.txt` だけからこの条件を
読んではならない。同じ文に装飾ありの `F̄` と装飾なしの `F` が両方あるため、
出力からの復元もできない。詳細は `ResearchPaper/papers.json` の
`IUTchI.notationNotes` と `1_Structured/.../section-3.html` の `#def-3-1-a`。
-/

namespace ABC3.Skeleton.IUTchI

open ABC3.Meta

/-- [IUTchI] Definition 3.1 の **条件 (a) の部分だけ**。

★これは initial Θ-data ではない。(b) 以降(曲線・素点・レベル `l`・`C_K`・`𝕍`・
`𝕍^bad_mod`・`ε`)を**含んでいない**。名前の `FieldPart` はその意味である。 -/
structure InitialThetaDataFieldPart where
  /-- 数体 `F`。 -/
  F : Type
  [instField : Field F]
  [instNumberField : NumberField F]
  /-- 原文 `such that √−1 ∈ F`。 -/
  sqrtNegOne : ∃ i : F, i ^ 2 = -1
  /-- `F̄`(原文のオーバーバー付き `F`)。 -/
  Fbar : Type
  [instFieldBar : Field Fbar]
  [instAlgebra : Algebra F Fbar]
  /-- 原文 `F̄ is an algebraic closure of F`。 -/
  instIsAlgClosure : IsAlgClosure F Fbar

def InitialThetaDataFieldPart.src : Source :=
  { paper := "IUTchI", pdfPage := 61, item := "Definition 3.1, (b)",
    sectionId := "def-3-1-b" }

/-- 原文の Definition 3.1 が要求するもの(G6)。

★これは **下界** である。(a) からは外部依存が出ないが、(b) 以降は
mathlib に無い語彙を要求する。ここではそれを `.citation` として測定つきで記録し、
原文が明示的に別論文へ渡している箇所を `.otherPaper` の**辺**として書く。 -/
def InitialThetaDataFieldPart.needs : List ProofObligation :=
  [ .otherPaper "[AbsTopIII]" "Definition 5.1, (ii)(field of moduli——原文 (b) が引用)" 113,
    .citation "mathlib v4.31.0-rc2" "once-punctured elliptic curve(1点抜き楕円曲線)"
      (.absent <|
        "(S1) 原文側の呼称は “once-punctured elliptic curve” / “hyperbolic curve of type (1,1)”。" ++
        "(S2) 全出現の列挙: `ls Mathlib/AlgebraicGeometry/EllipticCurve/` は 11 項目" ++
        "(Affine/, DivisionPolynomial/, IsomOfJ.lean, Jacobian/, LFunction.lean, " ++
        "ModelsWithJ.lean, NormalForms.lean, Projective/, Reduction.lean, " ++
        "VariableChange.lean, Weierstrass.lean)——いずれも Weierstrass 型の射影平面曲線で、" ++
        "点を抜いた開曲線も『型 (1,1) の双曲的曲線』も無い。" ++
        "`grep -rniE 'once.?punctured|punctured curve|hyperbolic curve' Mathlib/` → 0 件。" ++
        "(S3) Lean 側の語 `Punctured` → 位相の `punctured neighborhood` のみで別物。" ++
        "(S4) 上記のとおり。★隣接物すら無い(下の stable reduction とはここが違う)。")
      61,
    .citation "mathlib v4.31.0-rc2" "stable reduction(安定還元)"
      (.absent <|
        "(S1) 原文側の呼称は “admits stable reduction over all v”。" ++
        "(S2) `grep -rniE 'stable reduction|semistable reduction|IsStableReduction' Mathlib/` " ++
        "→ 0 件。`grep -rniE 'neron|néron' Mathlib/` → 0 件(Néron モデルも無い)。" ++
        "(S3) ★ただし**隣接物はある**: `Mathlib/AlgebraicGeometry/EllipticCurve/Reduction.lean` に " ++
        "`IsIntegral`(整係数 Weierstrass 式)・`IsMinimal`(判別式の付値が最小)・" ++
        "`reduction`(剰余体上の Weierstrass 曲線)・`IsGoodReduction`(判別式の付値が 0)が" ++
        "DVR の分数体上で定義されている。無いのは『良還元でない場合の分類』" ++
        "(乗法的/加法的)と『有限底変換の後に安定になる』という主張の側である。" ++
        "(S4) 上記のとおり。")
      61,
    .implicitStep
      "条件 (b) 以降(hyperbolic orbicurve `C_F`・field of moduli `F_mod`・maximal solvable extension `F_sol`・素点集合 `𝕍`・`𝕍^bad_mod`)を転写していない。ここで作った structure は initial Θ-data ではなく、その (a) 部分である"
      61 ]

end ABC3.Skeleton.IUTchI
