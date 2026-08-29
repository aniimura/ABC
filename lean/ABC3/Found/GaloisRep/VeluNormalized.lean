/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtFaltCovolume
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.IsogenyPeriodPair
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★Vélu の正規化はアルキメデス項を消す（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★これは何か

同種写像評価 `ht^Falt(E/H) ≤ ht^Falt(E) + 2·log(l)` に残っていた唯一の入力は
`Found/GaloisRep/IsogenyReduction.lean` の

    `hArch : (l−1)·d·deg∞(E) − (archSum(E′) − archSum(E)) ≤ 24·d·log(l)`

であった。★本ファイルはこれを **Vélu の正規化を使って純粋に有限素点の不等式へ落とす**。

## ★★★★★★機構 —— 共体積表示の差を取る

`§9-1023`（第 579、`twelve_finrank_htFaltOf_eq`）が**無条件で**与える:

    `12·d·ht^Falt(E) = −12·Σᶠ_p neronExp_p(E)·log N(p) − 12·d·log(2π)`
    `                  + 12·Σ_σ log‖u_σ‖ − 6·Σ_σ log(covol P_σ)`

★★`E` と `E′ = E/H` で引き算すると、`log(2π)` の項が消えて

    `12·d·(ht^Falt(E′) − ht^Falt(E))`
    `  = −12·Σᶠ_p [neronExp_p(E′) − neronExp_p(E)]·log N(p)`
    `    + 12·Σ_σ log(‖u′_σ‖/‖u_σ‖) − 6·Σ_σ log(covol P′_σ / covol P_σ)`

## ★★★★★★★★★★★★★★Vélu の正規化がここで効く

`Found/GenEll/Velu.lean` の `velu_omega_gen`（第 591）が示すとおり、
Vélu の写像は **`φ^*(ω_{E′}) = ω_E`** を満たす。すなわち ℂ 上では

| | 帰結 |
|---|---|
| 不変微分が一致する | ★**モデルのスケーリングが同じ**: `‖u′_σ‖ = ‖u_σ‖` |
| `φ` が `z ↦ z` から誘導される | ★**`Λ_σ ⊆ Λ′_σ` で指数 `l`**: `covol P′_σ = covol P_σ / l` |

★★★この 2 つを入れるとアルキメデスの差は

    `12·0 − 6·(−d·log l) = +6·d·log l`

になり、不等式 `12·d·(ht′−ht) ≤ 24·d·log l` は

    **`Σᶠ_p [neronExp_p(E) − neronExp_p(E′)]·log N(p) ≤ (3/2)·d·log(l)`**

★★★★**純粋に有限素点の主張**になる。☆これが [DelSB616] §2 の段 2 である。

## ★★★本ファイルが取るもの / 取らないもの

★**取るもの**: 上の引き算と変形。★★**`hArch` を `hfin` へ置き換えること**。

☆**取らないもの**: (1) `‖u′_σ‖ = ‖u_σ‖` と `covol P′_σ = covol P_σ/l` を
Vélu の写像から**実際に導く**こと（同種写像の構成が要る）、
(2) `hfin` そのもの（[FC] Ch. I, Prop 2.7）。

★★★しかし**アルキメデス側はこれで済んだ**——残るのは `neronExp` の差だけである。
-/

namespace ABC3.Found.GaloisRep

open NumberField IsDedekindDomain WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★★★★★★★★★★★ω-正規化のもとでの高さ評価 -/

/-- ★★★★★★★★★★★★★★★★★★★★**`ω`-正規化された同種写像の高さ評価**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★**アルキメデスの仮定は 2 つだけ**:

* `hu`  : `‖u′_σ‖ = ‖u_σ‖`——★`φ^*(ω′) = ω` なのでモデルのスケーリングが同じ
* `hcov`: `covol P′_σ = covol P_σ / l`——★`Λ_σ ⊆ Λ′_σ` で指数 `l`

★★この 2 つは `Found/GenEll/Velu.lean` の `velu_omega_gen`（第 591）の帰結である。

★★★★**残る `hfin` は純粋に有限素点の不等式**であり、
もはやアルキメデス量（`archSum`・`covol`・`u`）を一切含まない。 -/
theorem htFalt_isogeny_le_of_omega (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ) (hl : 0 < l)
    (P P' : (L →+* ℂ) → PeriodPair) (C C' : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hPC' : ∀ σ, C' σ • (E'.map σ) = latticeCurve (P' σ))
    (hu : ∀ σ, ‖((C' σ).u : ℂ)‖ = ‖((C σ).u : ℂ)‖)
    (hcov : ∀ σ, covol (P' σ) = covol (P σ) / l)
    (hfin : (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
          - (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        ≤ (3 / 2) * (Module.finrank ℚ L : ℝ) * Real.log l) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hl0 : (0:ℝ) < (l:ℝ) := by exact_mod_cast hl
  have h := twelve_finrank_htFaltOf_eq E P C hPC
  have h' := twelve_finrank_htFaltOf_eq E' P' C' hPC'
  -- ★`u` の和は等しい（ω-正規化）
  have hU : (∑ σ : (L →+* ℂ), Real.log ‖((C' σ).u : ℂ)‖)
      = ∑ σ : (L →+* ℂ), Real.log ‖((C σ).u : ℂ)‖ :=
    Finset.sum_congr rfl fun σ _ => by rw [hu σ]
  -- ★`covol` の和は `d·log l` だけ減る（指数 `l` の部分格子）
  have hV : (∑ σ : (L →+* ℂ), Real.log (covol (P' σ)))
      = (∑ σ : (L →+* ℂ), Real.log (covol (P σ)))
        - (Module.finrank ℚ L : ℝ) * Real.log l := by
    have hterm : ∀ σ : (L →+* ℂ), Real.log (covol (P' σ))
        = Real.log (covol (P σ)) - Real.log l := fun σ => by
      rw [hcov σ, Real.log_div (ne_of_gt (covol_pos (P σ))) (ne_of_gt hl0)]
    rw [Finset.sum_congr rfl fun σ _ => hterm σ, Finset.sum_sub_distrib,
      Finset.sum_const, Finset.card_univ, NumberField.Embeddings.card L ℂ, nsmul_eq_mul]
  rw [hU, hV] at h'
  nlinarith [h, h', hfin, hd]

/-! ## ★★★★★★★★★★★★★★★★★★格子の言葉に落とした版 -/

/-- ★★★★★★★★★★★★★★★★★★**`hcov` を基底変換の言葉に落とした版**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`Λ_σ ⊆ Λ′_σ` を「`P σ` の基底が `P′ σ` の基底の整係数一次結合で書け、
行列式の絶対値が `l`」という形で受ける——これは
`Found/GenEll/IsogenyPeriodPair.lean` の `covol_eq_index_mul_pair`（`§9-1019`）
がそのまま使える形である。

★★★すなわち **`hcov` はもはや仮定ではなく、格子の指数から従う**。
☆残るアルキメデスの仮定は `hu`（`φ^*(ω′) = ω` によるスケーリングの一致）だけになった。 -/
theorem htFalt_isogeny_le_of_omega_lattice
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ) (hl : 0 < l)
    (P P' : (L →+* ℂ) → PeriodPair) (C C' : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hPC' : ∀ σ, C' σ • (E'.map σ) = latticeCurve (P' σ))
    (hu : ∀ σ, ‖((C' σ).u : ℂ)‖ = ‖((C σ).u : ℂ)‖)
    (a b c d : (L →+* ℂ) → ℤ)
    (h₁ : ∀ σ, (P σ).ω₁ = ((a σ : ℤ) : ℂ) * (P' σ).ω₁ + ((b σ : ℤ) : ℂ) * (P' σ).ω₂)
    (h₂ : ∀ σ, (P σ).ω₂ = ((c σ : ℤ) : ℂ) * (P' σ).ω₁ + ((d σ : ℤ) : ℂ) * (P' σ).ω₂)
    (hdet : ∀ σ, (a σ * d σ - b σ * c σ).natAbs = l)
    (hfin : (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
          - (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        ≤ (3 / 2) * (Module.finrank ℚ L : ℝ) * Real.log l) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  have hl0 : (0:ℝ) < (l:ℝ) := by exact_mod_cast hl
  refine htFalt_isogeny_le_of_omega E E' l hl P P' C C' hPC hPC' hu ?_ hfin
  intro σ
  have := covol_eq_index_mul_pair (P σ) (P' σ) (a σ) (b σ) (c σ) (d σ) l
    (h₁ σ) (h₂ σ) (hdet σ)
  rw [eq_div_iff (ne_of_gt hl0), this]
  ring

def htFalt_isogeny_le_of_omega_lattice.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(hcov は格子の指数から従う——残るアルキメデスの仮定は hu だけ)",
    sectionId := "genell-prop-3-4" }

/-! ## ★出典の紐付け（`.src`）——★★条つき（`hu`・`hcov`・`hfin` が残っている） -/

def htFalt_isogeny_le_of_omega.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(ω-正規化のもとではアルキメデス項が消える)",
    sectionId := "genell-prop-3-4" }

def htFalt_isogeny_le_of_omega.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆仮定 hu・hcov は Vélu の写像から導くべきものである。" ++
       "★正規化 φ^*(ω′) = ω 自体は Found/GenEll/Velu.lean の velu_omega_gen" ++
       "(第 591)で取っているが、そこから ℂ 上の周期格子への翻訳" ++
       "(同種写像の構成と一様化定理)が残る") 8,
    .implicitStep
      ("☆仮定 hfin: Σᶠ_p [neronExp_p(E) − neronExp_p(E′)]·log N(p) ≤ (3/2)·d·log(l)。" ++
       "★これが [FC] Ch. I, Prop 2.7(極小モデルからの離れ具合)であり、" ++
       "mathlib に Néron モデルが無い(2026-08-29 に #check で確認)") 9,
    .implicitStep
      ("★★★到達点(2026-08-29、第 592): hArch(アルキメデスと有限が混ざった不等式)が" ++
       "hfin(純粋に有限素点の不等式)に置き換わった。" ++
       "★★共体積表示(第 579)の引き算で log(2π) の項も Σ_σ log‖u_σ‖ も消える") 9 ]

end ABC3.Found.GaloisRep
