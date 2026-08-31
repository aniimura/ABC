/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtJBound
import ABC3.Found.GenEll.JArchLower
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★`Prop 3.4` の 3 本の `≲`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か

原文の鎖

    `deg_∞ ≲ ht_∞ ≲ 12(1+ϵ)·ht^Falt ≲ (1+ϵ)·ht_∞`

を、**`ht_∞ ≔ htJ`（`j` の Weil 高さ）**で、`Found` の具体的な量として証明する。

| `≲` | 主張 | 状態 |
|---|---|---|
| 1 | `deg∞ ≤ htJ` | ★半安定（`deg∞ ≤ h_fin(j)`）のもとで |
| 2 | `htJ ≤ 12(1+ϵ)·ht^Falt + C` | ★★**無条件**（`§9-1003`） |
| 3 | `12(1+ϵ)·ht^Falt ≤ (1+ϵ)·htJ + C` | ★半安定のもとで（`§9-1006` の下からの評価） |

## ★半安定の入り方

`Found/GaloisRep/HtFinJ.lean` の `htFinJ_le_degInfOf` は
**`h_fin(j) ≤ deg∞` を無条件に**与える。★逆向き `deg∞ ≤ h_fin(j)` は
半安定（乗法還元で `v_p(c₄) = 0`）のときに成り立ち、そのとき**等号**になる。

★★本ファイルはそれを仮定 `hss : degInfOf L E ≤ htFinJ L E` として受ける
——★★★**受けた仮定は実装ではない**（`check.mjs` B6）。`.src` は条つきにする。

## ★★逆向きの評価（第 3 の `≲`）

`12·ht^Falt = deg∞ − archSum/d` なので、第 3 の `≲` は
`−archSum/d ≤ h∞(j) + C`、すなわち素点ごとに

    `log⁺|j| + log((2π)¹²‖Δ‖_arch) ≥ −C`

である。★これが `Found/GenEll/JArchLower.lean` の `exists_logPlus_j_lower`
（`max(‖j‖,1)·‖Δ‖_Pet ≥ m > 0`）である。
-/

namespace ABC3.Found.GaloisRep

open NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★埋め込みについて足す（下から） -/

/-- ★★★★★★★★**`h∞(j) + archSum/d ≥ −C`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`Found/GenEll/JArchLower.lean` の `logPlus_archNorm_lower` を
`σ : L →+* ℂ` について足して `[L:ℚ]` で割る。★★**無条件**。 -/
theorem htArchJ_add_archSum_ge (C : ℝ)
    (hC : ∀ (W : WeierstrassCurve ℂ) [W.IsElliptic],
      -C ≤ max (Real.log ‖W.j‖) 0 + Real.log ((2 * Real.pi) ^ 12 * curveArchInv W))
    (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic] :
    -C ≤ htArchJ L E + archSum L E / (Module.finrank ℚ L : ℝ) := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hsum : ∑ _σ : (L →+* ℂ), (-C)
      ≤ ∑ σ : (L →+* ℂ),
        (max (Real.log ‖σ E.j‖) 0
          + Real.log ((2 * Real.pi) ^ 12 * archNorm E σ)) :=
    Finset.sum_le_sum (fun σ _ => logPlus_archNorm_lower E σ C hC)
  rw [Finset.sum_add_distrib] at hsum
  rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, NumberField.Embeddings.card L ℂ] at hsum
  have heq : htArchJ L E + archSum L E / (Module.finrank ℚ L : ℝ)
      = ((∑ σ : (L →+* ℂ), max (Real.log ‖σ E.j‖) 0)
          + ∑ σ : (L →+* ℂ), Real.log ((2 * Real.pi) ^ 12 * archNorm E σ))
        / (Module.finrank ℚ L : ℝ) := by
    rw [htArchJ, archSum]; ring
  rw [heq, le_div_iff₀ hd]
  linarith

/-! ## ★★★★★★★★★★第 3 の `≲` -/

/-- ★★★★★★★★★★**`12·ht^Falt ≤ h(j) + C`**（半安定のもとで）。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`12·ht^Falt = deg∞ − archSum/d` に上の下からの評価を入れる。
☆`hss : deg∞ ≤ h_fin(j)` が半安定の入り口である。 -/
theorem htFalt_le_htJ (C : ℝ)
    (hC : ∀ (W : WeierstrassCurve ℂ) [W.IsElliptic],
      -C ≤ max (Real.log ‖W.j‖) 0 + Real.log ((2 * Real.pi) ^ 12 * curveArchInv W))
    (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : degInfOf L E ≤ htFinJ L E) :
    12 * htFaltOf L E ≤ htJ L E + C := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have harch : archSum L E / (Module.finrank ℚ L : ℝ) = degInfOf L E - 12 * htFaltOf L E := by
    rw [htFaltOf]; field_simp; ring
  have h := htArchJ_add_archSum_ge C hC L E
  rw [harch] at h
  rw [htJ]
  linarith

/-- ★★★★★**定数を消した形**。 -/
theorem exists_htFalt_le_htJ :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic],
      degInfOf L E ≤ htFinJ L E → 12 * htFaltOf L E ≤ htJ L E + C := by
  obtain ⟨C, hC⟩ := exists_logPlus_j_lower
  exact ⟨C, fun L _ _ E _ hss => htFalt_le_htJ C hC L E hss⟩

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★3 本の `≲` -/

theorem htArchJ_nonneg (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    [E.IsElliptic] : 0 ≤ htArchJ L E := by
  rw [htArchJ]
  exact div_nonneg (Finset.sum_nonneg (fun σ _ => le_max_right _ _)) (by positivity)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**原文の鎖**

    `deg_∞ ≲ ht_∞ ≲ 12(1+ϵ)·ht^Falt ≲ (1+ϵ)·ht_∞`

を `ht_∞ ≔ htJ`（`j` の Weil 高さ）で。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★第 2 の `≲` は**無条件**（`§9-1003`）。
☆第 1・第 3 は `hss : deg∞ ≤ h_fin(j)`（半安定）のもとで。
★★逆向き `h_fin(j) ≤ deg∞` は `Found/GaloisRep/HtFinJ.lean` で**無条件**に取れているので、
`hss` があれば `deg∞ = h_fin(j)` である。 -/
theorem prop_3_4_chain (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic],
      degInfOf L E ≤ htFinJ L E →
        degInfOf L E ≤ htJ L E
        ∧ htJ L E ≤ 12 * (1 + eps) * htFaltOf L E + C
        ∧ 12 * (1 + eps) * htFaltOf L E ≤ (1 + eps) * htJ L E + C := by
  obtain ⟨C₁, hC₁⟩ := exists_htJ_le_htFalt' eps heps
  obtain ⟨C₂, hC₂⟩ := exists_htFalt_le_htJ
  refine ⟨max C₁ ((1 + eps) * C₂), fun L _ _ E _ hss => ?_⟩
  have harch := htArchJ_nonneg L E
  have h1 : degInfOf L E ≤ htJ L E := by rw [htJ]; linarith
  have h2 := hC₁ L E
  have h3 := hC₂ L E hss
  refine ⟨h1, le_trans h2 (by simp), ?_⟩
  have h4 : (0:ℝ) < 1 + eps := by linarith
  have h5 : 12 * (1 + eps) * htFaltOf L E = (1 + eps) * (12 * htFaltOf L E) := by ring
  rw [h5]
  have h6 : (1 + eps) * (12 * htFaltOf L E) ≤ (1 + eps) * (htJ L E + C₂) :=
    mul_le_mul_of_nonneg_left h3 h4.le
  have h7 : (1 + eps) * C₂ ≤ max C₁ ((1 + eps) * C₂) := le_max_right _ _
  nlinarith

/-! ## ★出典の紐付け(`.src`)——★★**条つきである。指標には数えない** -/

def htArchJ_add_archSum_ge.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(h∞(j) + archSum/d ≥ −C——★無条件)",
    sectionId := "genell-prop-3-4" }

def exists_htFalt_le_htJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(第 3 の ≲ ——12·ht^Falt ≤ h(j) + C。★半安定 deg∞ ≤ h_fin(j) を仮定)",
    sectionId := "genell-prop-3-4" }

def prop_3_4_chain.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(3 本の ≲ を ht∞ = htJ で。★半安定 deg∞ ≤ h_fin(j) を仮定)",
    sectionId := "genell-prop-3-4" }

def prop_3_4_chain.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore
      ("☆**半安定のとき deg∞ ≤ h_fin(j)**: 乗法還元では v_p(c₄) = 0 なので " ++
       "v_p(j) = −v_p(Δ_min) となり、log⁺|j|_p = v_p(Δ_min)·log q_p で**等号**。" ++
       "★逆向き h_fin(j) ≤ deg∞ は Found/GaloisRep/HtFinJ.lean で**無条件**に取れている。" ++
       "☆本ファイルは hss として受けている——**受けた仮定は実装ではない**") 7,
    .implicitStep
      ("★★★★★★★★★★★★到達点(2026-08-29、第 558): " ++
       "原文の鎖 deg∞ ≲ ht∞ ≲ 12(1+ϵ)ht^Falt ≲ (1+ϵ)ht∞ の**3 本すべて**を、" ++
       "ht∞ = htJ(j の Weil 高さ)として Found の具体的な量で書いた。" ++
       "★第 2 の ≲ は無条件(§9-1003)、第 1・第 3 は半安定のもとで。" ++
       "★★第 3 の ≲ の中身は §9-1006 の max(‖j‖,1)·‖Δ‖_Pet ≥ m > 0 であり、" ++
       "その機構はカスプで E₄ → 1 であること") 9,
    .implicitStep
      ("☆残るのは ht∞ の**同定**である: Interface/GenEll/EllModuli.lean の htInf を " ++
       "htJ と同定し、EllClass(M_ell(ℚ̄) の点)を構成すること。" ++
       "★northcott 欄の中身は Found/GaloisRep/NorthcottHtJ.lean にある") 8 ]

end ABC3.Found.GaloisRep
