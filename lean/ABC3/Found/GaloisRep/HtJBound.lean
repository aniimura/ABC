/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.FaltingsWitness
import ABC3.Found.GenEll.JArchBound
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★`h(j) ≤ 12(1+ϵ)·ht^Falt + C`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か——`Prop 3.4` の**第 3 の `≲`**

`ResearchPaper/mathlib-gap.json` の `finalTarget20260829` が
「§3・§4 の 7 項目すべてに残る唯一の解析的入力」と測ったもの:

    **`h(j) ≤ 12(1+ϵ)·ht^Falt + C`**

★★★**本ファイルはそれを、有限素点側の 1 本の不等式を仮定として、証明する。**

## ★段取り

| 段 | 中身 | どこ |
|---|---|---|
| 1 | `‖j(τ)‖·(‖Δ‖_Pet(τ))^{1+ϵ}` は有界 | `Found/GenEll/JPeterssonBound.lean`（§9-1000） |
| 2 | 曲線の言葉に翻訳（`τ` を消す） | `Found/GenEll/JArchBound.lean`（§9-1001） |
| 3 | 埋め込みについて足して `[L:ℚ]` で割る | ★本ファイル `htArchJ_add_archSum_le` |
| 4 | `12·ht^Falt = deg∞ − archSum/d` を代入 | ★本ファイル `htArchJ_le` |
| 5 | 有限素点側を足す | ★本ファイル `htJ_le_htFalt` |

★段 4 の結論は

    `h∞(j) ≤ C + 12(1+ϵ)·ht^Falt − (1+ϵ)·deg∞`

であり、★★有限素点側が `h_fin(j) ≤ deg∞` なら足して

    `h(j) ≤ C + 12(1+ϵ)·ht^Falt − ϵ·deg∞ ≤ 12(1+ϵ)·ht^Falt + C`

（`deg∞ ≥ 0` と `ϵ > 0`）。★★★**`ϵ·deg∞` の余裕まで残る。**

## ☆残っているのは有限素点側の 1 本だけ

☆`h_fin(j) ≤ deg∞`——半安定なら `v(j) = v(c₄³) − v(Δ_min)` で
乗法還元では `v(c₄) = 0` なので `v(j) = −v(Δ_min) ≤ 0`、
したがって `log⁺|j|_v = v(Δ_min)·log q_v` で**等号**になる。
★本ファイルではこれを仮定 `hfin : htFinJ ≤ degInfOf L E` として受ける
——★★**受けた仮定は実装ではない**（`check.mjs` B6）。`.src` は条つきにする。
-/

namespace ABC3.Found.GaloisRep

open NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★アルキメデス素点での `h(j)` -/

/-- ★★★★★**`j` の高さのアルキメデス部分** `h∞(j) = (1/d)·Σ_σ log⁺‖σ(j)‖`。

★`d = [L:ℚ]`。★★`Real.log 0 = 0` なので `j = 0` でも `log⁺ = 0` で正しい。 -/
noncomputable def htArchJ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    [E.IsElliptic] : ℝ :=
  (∑ σ : (L →+* ℂ), max (Real.log ‖σ E.j‖) 0) / (Module.finrank ℚ L : ℝ)

/-! ## ★★★★★★★★段 3——埋め込みについて足す -/

/-- ★★★★★★★★**`h∞(j) + (1+ϵ)·archSum/d ≤ C`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`Found/GenEll/JArchBound.lean` の `logPlus_archNorm_le` を
`σ : L →+* ℂ` について足し、埋め込みの個数が `[L:ℚ]` であること
（`NumberField.Embeddings.card`）で割る。 -/
theorem htArchJ_add_archSum_le (eps C : ℝ)
    (hC : ∀ (W : WeierstrassCurve ℂ) [W.IsElliptic],
      max (Real.log ‖W.j‖) 0
        + (1 + eps) * Real.log ((2 * Real.pi) ^ 12 * curveArchInv W) ≤ C)
    (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic] :
    htArchJ L E + (1 + eps) * (archSum L E / (Module.finrank ℚ L : ℝ)) ≤ C := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hsum : ∑ σ : (L →+* ℂ),
      (max (Real.log ‖σ E.j‖) 0
        + (1 + eps) * Real.log ((2 * Real.pi) ^ 12 * archNorm E σ))
      ≤ ∑ _σ : (L →+* ℂ), C :=
    Finset.sum_le_sum (fun σ _ => logPlus_archNorm_le E σ eps C hC)
  rw [Finset.sum_add_distrib, ← Finset.mul_sum] at hsum
  rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, NumberField.Embeddings.card L ℂ] at hsum
  have heq : htArchJ L E + (1 + eps) * (archSum L E / (Module.finrank ℚ L : ℝ))
      = ((∑ σ : (L →+* ℂ), max (Real.log ‖σ E.j‖) 0)
          + (1 + eps) * ∑ σ : (L →+* ℂ),
              Real.log ((2 * Real.pi) ^ 12 * archNorm E σ))
        / (Module.finrank ℚ L : ℝ) := by
    rw [htArchJ, archSum]; ring
  rw [heq, div_le_iff₀ hd]
  linarith

/-! ## ★★★★★★★★★★段 4——Faltings 高さを代入する -/

/-- ★★★★★★★★★★**`h∞(j) ≤ C + 12(1+ϵ)·ht^Falt − (1+ϵ)·deg∞`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`htFaltOf` の定義（`12·ht^Falt = deg∞ − archSum/d`、`§9-670`）を代入するだけである。 -/
theorem htArchJ_le (eps C : ℝ)
    (hC : ∀ (W : WeierstrassCurve ℂ) [W.IsElliptic],
      max (Real.log ‖W.j‖) 0
        + (1 + eps) * Real.log ((2 * Real.pi) ^ 12 * curveArchInv W) ≤ C)
    (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic] :
    htArchJ L E ≤ C + 12 * (1 + eps) * htFaltOf L E - (1 + eps) * degInfOf L E := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have harch : archSum L E / (Module.finrank ℚ L : ℝ) = degInfOf L E - 12 * htFaltOf L E := by
    rw [htFaltOf]; field_simp; ring
  have h := htArchJ_add_archSum_le eps C hC L E
  rw [harch] at h
  linarith

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★段 5——目標 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`h(j) ≤ 12(1+ϵ)·ht^Falt + C`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`h_fin(j) ≤ deg∞`（仮定 `hfin`）と段 4 を足すと `−ϵ·deg∞` が残り、
`deg∞ ≥ 0` なので落とせる。

☆`hfin` は**受けた仮定**である——半安定のときの `v(j) = −v(Δ_min)` が中身であり、
それは代数の段として残る（`.src` は条つき）。 -/
theorem htJ_le_htFalt (eps : ℝ) (heps : 0 < eps) (C : ℝ)
    (hC : ∀ (W : WeierstrassCurve ℂ) [W.IsElliptic],
      max (Real.log ‖W.j‖) 0
        + (1 + eps) * Real.log ((2 * Real.pi) ^ 12 * curveArchInv W) ≤ C)
    (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic]
    (htFinJ : ℝ) (hfin : htFinJ ≤ degInfOf L E) :
    htFinJ + htArchJ L E ≤ 12 * (1 + eps) * htFaltOf L E + C := by
  have h0 : 0 ≤ degInfOf L E := degInfOf_nonneg E
  have h := htArchJ_le eps C hC L E
  nlinarith [h, hfin, h0, heps]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**定数を消した形**——`ϵ > 0` ごとに `C` が取れる。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`C` は `L` にも `E` にも依らない（`Found/GenEll/JArchBound.lean` の
`exists_logPlus_j_bound` が与える普遍定数）——★★これが `Prop 3.4` に要る形である。 -/
theorem exists_htJ_le_htFalt (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic]
      (htFinJ : ℝ), htFinJ ≤ degInfOf L E →
      htFinJ + htArchJ L E ≤ 12 * (1 + eps) * htFaltOf L E + C := by
  obtain ⟨C, hC⟩ := exists_logPlus_j_bound eps heps
  exact ⟨C, fun L _ _ E _ htFinJ hfin => htJ_le_htFalt eps heps C hC L E htFinJ hfin⟩

/-! ## ★出典の紐付け(`.src`)——★★**条つきである。指標には数えない** -/

def htArchJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(j の高さのアルキメデス部分)",
    sectionId := "genell-prop-3-4" }

def htArchJ_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(h∞(j) ≤ C + 12(1+ϵ)ht^Falt − (1+ϵ)deg∞——★無条件)",
    sectionId := "genell-prop-3-4" }

def exists_htJ_le_htFalt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(h(j) ≤ 12(1+ϵ)ht^Falt + C——★有限素点側 h_fin(j) ≤ deg∞ を仮定)",
    sectionId := "genell-prop-3-4" }

def exists_htJ_le_htFalt.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[Silv2]"
      ("Proposition 2.1——★原文が Prop 3.4 の根拠として引く。" ++
       "★★★**アルキメデス側は §9-1000〜1002 で閉じた**") 3,
    .folklore
      ("☆**有限素点側**: `h_fin(j) ≤ deg∞`。半安定なら " ++
       "v(j) = v(c₄³) − v(Δ_min) で、乗法還元では v(c₄) = 0 なので " ++
       "v(j) = −v(Δ_min) ≤ 0、したがって log⁺|j|_v = v(Δ_min)·log q_v で**等号**。" ++
       "★本ファイルは仮定 hfin として受けている——**受けた仮定は実装ではない**") 7,
    .implicitStep
      ("★★★★★★★★★★到達点(2026-08-29、第 553): " ++
       "ResearchPaper/mathlib-gap.json の finalTarget20260829 が" ++
       "「§3・§4 の 7 項目すべてに残る唯一の解析的入力」と測った " ++
       "h(j) ≤ 12(1+ϵ)ht^Falt + C を、**有限素点側 1 本を仮定として証明した**。" ++
       "★段取りは 5 段: (1) ‖j‖·‖Δ‖^{1+ϵ} 有界(§9-1000)、" ++
       "(2) 曲線の言葉へ(§9-1001)、(3) 埋め込みについて和、" ++
       "(4) 12·ht^Falt = deg∞ − archSum/d を代入、(5) 有限素点側を足す。" ++
       "★★段 4 の htArchJ_le は**無条件**である") 9,
    .implicitStep
      ("★次の接続先: Found/GenEll/FiniteFromNorthcott.lean の " ++
       "finite_of_le_of_northcott の仮定 hle は本ファイルが与える。" ++
       "★★仮定 hN(Northcott)は Found/GenEll/NorthcottImage.lean の " ++
       "northcott_of_log_mulHeight_image が与える。" ++
       "☆残るのは ht∞ の同定(M̄_ell の無限遠因子の高さが h(j) であること)") 8 ]

end ABC3.Found.GaloisRep
