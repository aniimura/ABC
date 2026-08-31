/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.FaltingsWitness
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Found.GaloisRep.Compositum
import ABC3.Meta.Claim

/-!
# `EllModuliData` の対象 —— 半安定な楕円曲線の族（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17–p.23。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★これは何か

`ResearchPaper/ellmoduli-witness-status.json` の `designChoice` で決めた設計を
**Lean の型として固定する**。★これで次のセッションが設計を検討し直さずに済み、
型検査で設計の齟齬も早く見つかる。

## ★★★設計

    Curve    := SSCurve（数体 `L ⊆ ℂ` と、その上の半安定な楕円曲線）
    EllClass := ℂ（`j` 不変量）
    cls E    := E.j

★`EllClass` は**本物の商**でなければならない——`EllClass := Curve` にすると
`northcott`（高さと次数で抑えた類の有限性）が偽になる（同型な Weierstrass モデルが
無限にあるから）。

★★`Curve` を**半安定なものだけ**に制限するのが要点である。そうすると
`EllModuliData` の `torsionExt` 群（`torsionExt` / `cls_torsionExt` /
`degOfDefinition_torsionExt` / `semiStable_torsionExt` / `hasMultRed_torsionExt` /
`primeToLocalHeights_torsionExt`）が `torsionExt := id` で**すべて自明に埋まる**
——これが唯一の mathlib 欠落（Néron–Ogg–Shafarevich、半安定還元の判定）だった。

☆ただしその結果は「半安定な曲線に対する `Lemma 3.5` 等」であって原文の主張そのもの
ではない（原文は `L(E[3],E[5])` への基底変換で一般の場合を半安定へ帰着させる）。
★したがって `.src` は条つきである。

## ★本ファイルで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `SSCurve` | ★数体とその上の半安定な楕円曲線 |
| `SSCurve.j` | ★`j` 不変量を `ℂ` へ |
| `SSCurve.deg` / `SSCurve.htFalt` / `SSCurve.degInf` | ★界面の欄に渡す量 |
| `SSCurve.deg_pos` | ★`degOfDefinition_pos` の中身 |
| `SSCurve.degInf_nonneg` | ★`deg∞ ≥ 0` |
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField

/-! ## ★★★★★★半安定な楕円曲線 -/

/-- ★★★★★★**数体 `L` とその上の半安定な楕円曲線**。

★`EllModuliData` の `Curve` 欄に渡す対象である。 -/
structure SSCurve where
  /-- 定義体（`ℂ` の部分体として取る——`Prop34.lean` の族の形に合わせる）。 -/
  fld : Type
  [isField : Field fld]
  [isNF : NumberField fld]
  /-- 曲線。 -/
  W : WeierstrassCurve fld
  [isEll : W.IsElliptic]
  /-- ★全ての有限素点で半安定。 -/
  ss : ∀ p : HeightOneSpectrum (𝓞 fld), SemistableAt p W
  /-- `j` を `ℂ` へ送る埋め込み。 -/
  emb : fld →+* ℂ

attribute [instance] SSCurve.isField SSCurve.isNF SSCurve.isEll

namespace SSCurve

/-- ★★★★`j` 不変量（`ℂ` の元として）——`cls` 欄。 -/
noncomputable def j (E : SSCurve) : ℂ := E.emb (E.W.j)

/-- ★★★定義体の次数——`degOfDefinition` 欄。 -/
noncomputable def deg (E : SSCurve) : ℕ := Module.finrank ℚ E.fld

/-- ★★★★Faltings 高さ——`faltingsHeight` を作る材料。 -/
noncomputable def htFalt (E : SSCurve) : ℝ := htFaltOf E.fld E.W

/-- ★★★`deg_∞`——`degInf` を作る材料。 -/
noncomputable def degInf (E : SSCurve) : ℝ := degInfOf E.fld E.W

/-! ## ★★★★基本の性質 -/

/-- ★★`degOfDefinition_pos` 欄の中身。 -/
theorem deg_pos (E : SSCurve) : 0 < E.deg := Module.finrank_pos

/-- ★★`degInf_nonneg` 欄の中身。 -/
theorem degInf_nonneg (E : SSCurve) : 0 ≤ E.degInf := degInfOf_nonneg E.W

/-- ★★★★**変数変換で `ht^Falt` は変わらない**（第 722 を `SSCurve` の言葉で）。 -/
theorem htFalt_variableChange (E : SSCurve) (C : VariableChange E.fld)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p (C • E.W))
    (hell : (C • E.W).IsElliptic) :
    htFalt { fld := E.fld, W := C • E.W, isEll := hell, ss := hss, emb := E.emb }
      = E.htFalt := by
  show htFaltOf E.fld (C • E.W) = htFaltOf E.fld E.W
  exact htFaltOf_variableChange E.W C

end SSCurve

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★`ht^Falt` は `j` だけで決まる -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`j` が同じ `SSCurve` は `ht^Falt` も同じ**——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★これが `EllModuliData` の `faltingsHeight : EllClass → ℝ` の well-defined 性である
（`§9-1168`、第 741 の `htFaltOf_congr_j_of_emb` をそのまま `SSCurve` の言葉で）。 -/
theorem htFalt_congr_j (E E' : SSCurve) (hj : E.j = E'.j) : E.htFalt = E'.htFalt :=
  ABC3.Found.GaloisRep.htFaltOf_congr_j_of_emb E.fld E'.fld E.emb E'.emb E.W E'.W E.ss E'.ss hj

open scoped Classical in
/-- ★★★★★★★★★★★★★★**`EllModuliData` の `faltingsHeight` 欄**——`j` の函数として。

★`j` を持つ半安定曲線が無い `j` では `0` と定める（界面はそこでは何も要求しない）。 -/
noncomputable def faltingsHeightJ (j : ℂ) : ℝ :=
  if h : ∃ E : SSCurve, E.j = j then h.choose.htFalt else 0

/-- ★★★★★★★★★★★★★★★★★★**欄の値は曲線の `ht^Falt` に一致する**。 -/
theorem faltingsHeightJ_eq (E : SSCurve) : faltingsHeightJ E.j = E.htFalt := by
  classical
  have h : ∃ E' : SSCurve, E'.j = E.j := ⟨E, rfl⟩
  rw [faltingsHeightJ, dif_pos h]
  exact htFalt_congr_j h.choose E h.choose_spec

/-! ## ★出典の紐付け(`.src`)——★★条つき（半安定に制限した形） -/

def SSCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(EllModuliData の Curve 欄——半安定な楕円曲線の族)",
    sectionId := "genell-prop-3-4" }

def htFalt_congr_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(j が同じ SSCurve は ht^Falt も同じ。★無条件)",
    sectionId := "genell-prop-3-4" }

def faltingsHeightJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(EllModuliData の faltingsHeight 欄——j の函数として)",
    sectionId := "genell-prop-3-4" }

def SSCurve.htFalt_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(SSCurve の上で ht^Falt は変数変換で不変。★無条件)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
