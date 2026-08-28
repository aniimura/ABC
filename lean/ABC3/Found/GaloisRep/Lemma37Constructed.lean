/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.Prop34Chain
import ABC3.Found.GenEll.Sec3Arith
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★`Lemma 3.7, (a)` は**構成した量で**取れる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★★★★★★★★★測定（2026-08-29）—— 第 3 の `≲` は要らなかった

`Skeleton/GenEll/Section3.lean` の `lemma_3_7` は `prop_3_4` を引くが、
★**実際に使っているのは最初の 2 つの `≲` だけ**である:

    `deg∞ ≲ ht∞`  と  `ht∞ ≲ 12(1+ϵ)·ht^Falt`   ⟹  `deg∞ ≤ 14·ht^Falt + A`（`ϵ₀ = 1/6`）
    `ht^Falt ≥ B`（下に有界）

★★★**その 2 つは、構成した `degInfOf`／`htFaltOf` について本ファイルで証明できる**:

| 入力 | 出どころ |
|---|---|
| `deg∞ ≤ 14·ht^Falt + (7/6)·log((2π)¹²M)` | ★`§9-980` の `degInf_le_htFalt` を `ϵ = 1/6` で |
| `ht^Falt ≥ −log((2π)¹²M)/12` | ★**本ファイル**（`deg∞ ≥ 0` と `archSum ≤ d·log((2π)¹²M)`） |
| 数値の核（`C ≔ \|A\| + 100\|B\| + 1`） | ★`§9-991` の `lemma_3_7_numeric` |

★★★★**したがって `Lemma 3.7, (a)` は `Prop 3.4` の第 3 の `≲`（[Silv2] Prop 2.1、
並行セッションの領分）を**待たずに**取れる**——本ファイルの `lemma_3_7_a_constructed`。

## ★★残るのは (b)(c) の例外集合

☆`Lemma 3.7` を**項目として**閉じるには (b)(c) が要る:
`Exc ≔ lcyclicExc ∪ noMultRedExc K_V` を取り、それが `Galois`-finite であること。
★これは `Proposition 1.4, (iv)` の有限性と compactly bounded subset の機械であり、
`Interface/GenEll/EllModuli.lean` に posit されている。

★`.src` は条つき——指標には数えない。
-/

namespace ABC3.Found.GaloisRep

open NumberField WeierstrassCurve ABC3.Found.GenEll

/-! ## ★★★★★★★★`ht^Falt` は下に有界 -/

/-- ★★★★★★★★**`ht^Falt ≥ −log((2π)¹²M)/12`**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★`12·ht^Falt = deg∞ − archSum/d`、`deg∞ ≥ 0`、`archSum/d ≤ log((2π)¹²M)` から出る。
★★`Prop 3.4` の `faltingsHeight_bddBelow`（界面の posit）に対応する内容が、
**構成した `htFaltOf` について証明できた**。 -/
theorem htFalt_bddBelow (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    (M : ℝ) (hM : 1 ≤ M) (hMb : ∀ W : WeierstrassCurve ℂ, curveArchInv W ≤ M) :
    -(Real.log ((2 * Real.pi) ^ 12 * M)) / 12 ≤ htFaltOf L E := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hnn := degInfOf_nonneg (L := L) E
  have hA := archSum_le L E M hM hMb
  have hAd : archSum L E / (Module.finrank ℚ L : ℝ) ≤ Real.log ((2 * Real.pi) ^ 12 * M) := by
    rw [div_le_iff₀ hd]
    linarith [hA]
  have hexp : 12 * htFaltOf L E
      = degInfOf L E - archSum L E / (Module.finrank ℚ L : ℝ) := by
    rw [htFaltOf]
    field_simp
  linarith [hexp, hAd, hnn]

/-! ## ★★★★★★`deg∞ ≤ 14·ht^Falt + A` -/

/-- ★★★★★★**`deg∞ ≤ 14·ht^Falt + (7/6)·log((2π)¹²M)`**
（`Prop 3.4` を `ϵ₀ = 1/6` で使った形）。

★`§9-980` の `degInf_le_htFalt` に `ϵ = 1/6` を入れるだけである
（`12·(1 + 1/6) = 14`——`Skeleton` の `lemma_3_7` が使う係数そのもの）。 -/
theorem degInf_le_fourteen (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    (M : ℝ) (hM : 1 ≤ M) (hMb : ∀ W : WeierstrassCurve ℂ, curveArchInv W ≤ M) :
    degInfOf L E ≤ 14 * htFaltOf L E + (7/6) * Real.log ((2 * Real.pi) ^ 12 * M) := by
  have h := degInf_le_htFalt L E M hM hMb (1/6) (by norm_num)
  norm_num at h ⊢
  linarith [h]

/-! ## ★★★★★★★★★★★★★★★★★★★★`Lemma 3.7, (a)` -/

/-- ★★★★★★★★★★★★★★★★★★★★**`Lemma 3.7, (a)` を構成した量で**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

    `100d·(ht^Falt + C·d^ϵ) ≤ l`  ⟹  `d·deg∞ < l·log(2)`

（`C ≔ |A| + 100|B| + 1`、`A ≔ (7/6)log((2π)¹²M)`、`B ≔ −log((2π)¹²M)/12`）。

★★★**`Prop 3.4` の第 3 の `≲`（[Silv2] Prop 2.1）を待たずに取れる**
——`Skeleton` の `lemma_3_7` が使うのは最初の 2 つの `≲` だけだからである。
★`§9-980`（第 1・第 2 の合成）＋ 本ファイルの下界 ＋ `§9-991`（数値の核）で閉じる。 -/
theorem lemma_3_7_a_constructed (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    (M : ℝ) (hM : 1 ≤ M) (hMb : ∀ W : WeierstrassCurve ℂ, curveArchInv W ≤ M)
    (d : ℝ) (hd : 1 ≤ d) (eps : ℝ) (heps : 0 < eps) (l : ℝ)
    (hlog : (0.69 : ℝ) ≤ Real.log 2)
    (hl : 100 * d * (htFaltOf L E
            + (|(7/6) * Real.log ((2 * Real.pi) ^ 12 * M)|
               + 100 * |(-(Real.log ((2 * Real.pi) ^ 12 * M)) / 12)| + 1) * d ^ eps) ≤ l) :
    d * degInfOf L E < l * Real.log 2 :=
  lemma_3_7_numeric d (degInfOf L E) (htFaltOf L E)
    ((7/6) * Real.log ((2 * Real.pi) ^ 12 * M))
    (-(Real.log ((2 * Real.pi) ^ 12 * M)) / 12) _ l eps hd heps
    (degInf_le_fourteen L E M hM hMb) (htFalt_bddBelow L E M hM hMb) rfl hl hlog

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def htFalt_bddBelow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(ht^Falt は下に有界——構成した量で)",
    sectionId := "genell-lemma-3-7" }

def degInf_le_fourteen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(deg∞ ≤ 14·ht^Falt + A——ϵ₀ = 1/6)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_a_constructed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7, (a)(構成した量で——Prop 3.4 の第 3 の ≲ を待たない)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_a_constructed.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "degInf_le_htFalt(Prop 3.4 の第 1・第 2 の合成、§9-980)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.degInf_le_htFalt") 3,
    .citation "[ABC3]" "lemma_3_7_numeric(数値の核、§9-991)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_7_numeric") 3,
    .citation "[ABC3]" "archSum_le(アルキメデス和の一様上界、第 355)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.archSum_le") 2,
    .implicitStep
      ("★★★★★★★測定(2026-08-29): Skeleton の lemma_3_7 は prop_3_4 を引くが、" ++
       "**実際に使っているのは最初の 2 つの ≲ だけ**である" ++
       "(deg∞ ≲ ht∞ ≲ 12(1+ϵ)ht^Falt から deg∞ ≤ 14·ht^Falt + A、および ht^Falt ≥ B)。" ++
       "★その 2 つは構成した degInfOf／htFaltOf について証明できる" ++
       "——★★**したがって Lemma 3.7, (a) は Prop 3.4 の第 3 の ≲" ++
       "([Silv2] Prop 2.1、並行セッションの領分)を待たずに取れる**") 7,
    .implicitStep
      ("★★Lemma 3.7 を**項目として**閉じるには (b)(c) が要る: " ++
       "Exc ≔ lcyclicExc ∪ noMultRedExc K_V を取り、それが Galois-finite であること。" ++
       "★これは Proposition 1.4, (iv) の有限性と compactly bounded subset の機械であり、" ++
       "Interface/GenEll/EllModuli.lean に posit されている") 6 ]

end ABC3.Found.GaloisRep
