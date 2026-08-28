/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.FaltingsWitness
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★`Proposition 3.4` の鎖のうち**取れる 2 本**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★これは何か —— 原文の 3 本の `≲` を分解する

原文の主張は

    `deg∞ ≲ ht∞ ≲ 12(1 + ϵ)·ht^Falt ≲ (1 + ϵ)·ht∞`

と、それに続く**有限性**（`ht^Falt ≤ C` なる `M_ell(ℚ̄)^{≤d}` の点は有限個）である。

★★`§9-670`（第 357 ブロック）で `htFaltOf` は**構成されている**:

    `12·ht^Falt(E) = deg∞(E) − (1/d)·Σ_σ log( (2π)¹²·‖Δ‖_arch(E^σ) )`

★★★したがって 3 本の `≲` は、**アルキメデス和 `archSum` の評価だけ**の問題になる。

## ★★★★★★測定（2026-08-29）—— どれが取れてどれが残るか

★★**本ファイルが取るのは第 1 と第 2 の `≲` の合成**、すなわち

    **`deg∞ ≲ 12(1 + ϵ)·ht^Falt`**

である（`degInf_le_htFalt`）。★`archSum` の**上界だけ**で出る。

| 原文の段 | 必要なもの | 状態 |
|---|---|---|
| `deg∞ ≲ ht∞ ≲ 12(1+ϵ)·ht^Falt`（合成） | `archSum` の**上界**だけ | ★★**本ファイル** |
| `12(1+ϵ)·ht^Falt ≲ (1+ϵ)·ht∞` | `archSum` の**下界** | ☆**残る**（[Silv2] `Proposition 2.1`） |
| 有限性 | 上の不等式 ＋ `Proposition 1.4, (iv)` | ★半分（`degInf_le_of_htFalt_le`） |

## ★★★★★逸脱（明示）—— `ht∞` を分離していない

★★★**本ファイルは `ht∞` そのものを扱っていない。** 原文の `ht∞` は
`M̄_ell` の無限遠因子に付随する高さ（実質 `h(j)`）であり、`deg∞`（有限素点だけの寄与）
とは別物である。★合成した形（`deg∞ ≲ 12(1+ϵ)ht^Falt`）は原文の鎖の**帰結**であって、
第 2 の `≲` そのものではない。

★★したがって `Proposition 3.4` に残るのは**2 つ**である:

1. **第 3 の `≲`**（`12(1+ϵ)ht^Falt ≲ (1+ϵ)ht∞`）。その中身は
   `−archSum/d ≲ ϵ·ht∞ + C`——すなわち **`log(1/‖Δ‖_arch)` が `ht∞` で抑えられる**こと。
   ★`archNorm` に**一様な正の下界は無い**（`Im τ → ∞` で `‖Δ‖ ~ e^{−2π Im τ} → 0`）ので、
   `archSum` 単独では出ず **`ht∞` との相殺**が要る。
   ★★これが原文の引く [Silv2], `Proposition 2.1`（および [FC] Chapter V, Proposition 4.5 の
   「`the logarithmic singularities at infinity of the metric defined on ωE`」）である。
2. **`ht∞` の同定**——`M̄_ell` の無限遠因子の高さが `h(j)` であること。
   ★`deg∞ ≲ ht∞`（原文の第 1 の `≲`、「定義から直ちに、cf. `Proposition 1.6` の証明」）は
   その同定のもとで「有限素点の寄与 ≤ 全体」になる。

## ★有限性の帰着（本ファイルが取る半分）

`12·ht^Falt = deg∞ − archSum/d` と `archSum/d ≤ log((2π)¹²M)`（`§9-355` の一様上界）から

    `ht^Falt(E) ≤ C  ⟹  deg∞(E) ≤ 12C + log((2π)¹²M)`

★すなわち **`ht^Falt` が有界なら `deg∞` も有界**である（`degInf_le_of_htFalt_le`）。
★★残るのは「`deg∞` が有界な `M_ell(ℚ̄)^{≤d}` の点は有限個」——
これが原文の言う `Proposition 1.4, (iv)` の適用であり、
`Proposition 1.4, (iv)` 自体は `§9-938`／`§9-961` で取ってある。
★★★**残るのは `deg∞` を `M̄_ell` の豊富な直線束の高さと同定する段**である。
-/

namespace ABC3.Found.GaloisRep

open NumberField WeierstrassCurve ABC3.Found.GenEll

/-! ## ★★★★★★★★★★★★第 1・第 2 の `≲` の合成 -/

/-- ★★★★★★★★★★★★**`deg∞ ≲ 12(1+ϵ)·ht^Falt`**（原文の第 1・第 2 の `≲` の合成）。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

    `deg∞(E) ≤ 12(1+ϵ)·ht^Falt(E) + (1+ϵ)·log((2π)¹²M)`

★★機構: `12·ht^Falt = deg∞ − archSum/d` と `archSum ≤ d·log((2π)¹²M)`
（`archSum_le`、第 355 の一様上界）、そして `deg∞ ≥ 0`・`ϵ > 0`。
★★★**`archSum` の上界だけで出る**——下界は要らない。 -/
theorem degInf_le_htFalt (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    (M : ℝ) (hM : 1 ≤ M) (hMb : ∀ W : WeierstrassCurve ℂ, curveArchInv W ≤ M)
    (eps : ℝ) (heps : 0 < eps) :
    degInfOf L E
      ≤ 12 * (1 + eps) * htFaltOf L E + (1 + eps) * Real.log ((2 * Real.pi) ^ 12 * M) := by
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
  nlinarith [hexp, hAd, hnn, heps]

/-! ## ★★★★★★★★★★有限性の帰着 -/

/-- ★★★★★★★★★★**`ht^Falt` が有界なら `deg∞` も有界**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

    `ht^Falt(E) ≤ C  ⟹  deg∞(E) ≤ 12C + log((2π)¹²M)`

★原文の有限性の主張（`ht^Falt([E]) ≤ C` なる `M_ell(ℚ̄)^{≤d}` の点は有限個）は、
「すでに示した不等式 ＋ `Proposition 1.4, (iv)`」から出る、と原文は書く。
★★本補題がその**前半**である——`ht^Falt` の有界性を `deg∞` の有界性に落とす。
★★★残るのは `deg∞` を `M̄_ell` の豊富な直線束の高さと同定して
`Proposition 1.4, (iv)`（`§9-938`／`§9-961`）を当てる段である。 -/
theorem degInf_le_of_htFalt_le (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    (M : ℝ) (hM : 1 ≤ M) (hMb : ∀ W : WeierstrassCurve ℂ, curveArchInv W ≤ M)
    (C : ℝ) (hC : htFaltOf L E ≤ C) :
    degInfOf L E ≤ 12 * C + Real.log ((2 * Real.pi) ^ 12 * M) := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hA := archSum_le L E M hM hMb
  have hAd : archSum L E / (Module.finrank ℚ L : ℝ) ≤ Real.log ((2 * Real.pi) ^ 12 * M) := by
    rw [div_le_iff₀ hd]
    linarith [hA]
  have hexp : 12 * htFaltOf L E
      = degInfOf L E - archSum L E / (Module.finrank ℚ L : ℝ) := by
    rw [htFaltOf]
    field_simp
  linarith [hexp, hAd, hC]

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def degInf_le_htFalt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(第 1・第 2 の ≲ の合成——deg∞ ≲ 12(1+ϵ)·ht^Falt)",
    sectionId := "genell-prop-3-4" }

def degInf_le_of_htFalt_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(有限性の帰着——ht^Falt が有界なら deg∞ も有界)",
    sectionId := "genell-prop-3-4" }

def degInf_le_htFalt.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "archSum_le(アルキメデス和の一様上界、第 355)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.archSum_le") 3,
    .citation "[ABC3]" "htFaltOf(Faltings 高さの式、第 357)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFaltOf") 3,
    .otherPaper "[Silv2]"
      ("Proposition 2.1——★★原文が第 2・第 3 の ≲ の根拠として引く。" ++
       "★測定(2026-08-29): **第 1・第 2 の ≲ の合成(deg∞ ≲ 12(1+ϵ)ht^Falt)は" ++
       "archSum の上界だけで出る**(本ファイル)。" ++
       "★★残るのは**第 3 の ≲**であり、その中身は −archSum/d ≲ ϵ·ht∞ + C" ++
       "——すなわち log(1/‖Δ‖_arch) が ht∞ で抑えられること。" ++
       "★★★archNorm に一様な正の下界は無い(Im τ → ∞ で ‖Δ‖ ~ e^{−2π Im τ} → 0)ので、" ++
       "archSum 単独では出ず ht∞ との相殺が要る" ++
       "——これが『無限遠での対数的特異性』の意味である") 9,
    .otherPaper "[FC]"
      "Chapter V, Proposition 4.5(ωE 上の計量の無限遠での対数的特異性)" 9,
    .implicitStep
      ("★★★★★★★逸脱(2026-08-29): **本ファイルは ht∞ そのものを扱っていない**。" ++
       "原文の ht∞ は M̄_ell の無限遠因子に付随する高さ(実質 h(j))であり、" ++
       "deg∞(有限素点だけの寄与)とは別物である。" ++
       "★取ったのは合成した形(deg∞ ≲ 12(1+ϵ)ht^Falt)であって、" ++
       "第 2 の ≲ そのものではない") 7,
    .implicitStep
      ("★★★★★測定(2026-08-29): Proposition 3.4 に残るのは**2 つ**である。" ++
       "(1) 第 3 の ≲(12(1+ϵ)ht^Falt ≲ (1+ϵ)ht∞)——[Silv2] Prop 2.1。" ++
       "(2) ht∞ の同定(M̄_ell の無限遠因子の高さが h(j) であること)。" ++
       "★有限性の主張は半分(ht^Falt 有界 ⟹ deg∞ 有界)が取れた" ++
       "——残るのは deg∞ を M̄_ell の豊富な直線束の高さと同定して" ++
       "Proposition 1.4, (iv)(§9-938／§9-961)を当てる段である") 7 ]

end ABC3.Found.GaloisRep
