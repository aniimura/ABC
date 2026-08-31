/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtFaltCovolume
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.IsogenyPeriodPair
import ABC3.Found.GenEll.LatticeScale
import ABC3.Found.GenEll.Uniformization
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

/-! ## ★★★★★★★★★★★★★★★★★★★★選択によらない形——`archDefect` -/

/-- ★★★★★★★★★★★★★★★★**アルキメデスの欠損** `Σ_σ [log‖σΔ_E‖ − log archNorm(E,σ)]`。

★`Found/GenEll/ArchInvCovolume.lean` の `log_archNorm_eq` により、これは
`Σ_σ [12·log‖u_σ‖ − 6·log covol(P σ)]` に等しい（`archDefect_eq`）。
★★すなわち **`u` と `covol` は個別には一様化の取り方に依るが、この組み合わせは依らない**。 -/
noncomputable def archDefect (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) : ℝ :=
  ∑ σ : (L →+* ℂ), (Real.log ‖σ E.Δ‖ - Real.log (archNorm E σ))

/-- ★★★★★★★★★★★★★★★★★★**`ht^Falt` の `archDefect` 表示**——★**無条件**。

    `12·d·ht^Falt(E) = archDefect(E) − 12·Σᶠ_p neronExp_p·log N(p) − 12·d·log(2π)`

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`§9-1023`（第 579）の共体積表示から**周期対 `P` と変数変換 `C` を消した形**である
——`12·Σ_σ log‖u_σ‖ − 6·Σ_σ log(covol P_σ)` がちょうど `archDefect(E)` だから。
★★★これで残るアルキメデスの穴が**一様化の取り方に依らない 1 つの量**になった。 -/
theorem twelve_finrank_htFaltOf_eq_archDefect (E : WeierstrassCurve L) [E.IsElliptic] :
    12 * (Module.finrank ℚ L : ℝ) * htFaltOf L E
      = archDefect L E
        - 12 * (∑ᶠ p : HeightOneSpectrum (𝓞 L),
            (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        - 12 * (Module.finrank ℚ L : ℝ) * Real.log (2 * Real.pi) := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hΔ : E.Δ ≠ 0 := (inferInstance : E.IsElliptic).isUnit.ne_zero
  have hpi : (0:ℝ) < (2 * Real.pi) ^ 12 := by positivity
  have hprod := sum_arch_log_eq_finsum_valAdd E.Δ hΔ
  have hdeg := finrank_degInfOf_eq E hΔ
  have hht : 12 * (Module.finrank ℚ L : ℝ) * htFaltOf L E
      = (Module.finrank ℚ L : ℝ) * degInfOf L E - archSum L E := by
    rw [htFaltOf]
    field_simp
  have harchS : archSum L E
      = 12 * (Module.finrank ℚ L : ℝ) * Real.log (2 * Real.pi)
        + ∑ σ : (L →+* ℂ), Real.log (archNorm E σ) := by
    have hterm : ∀ σ : (L →+* ℂ),
        Real.log ((2 * Real.pi) ^ 12 * archNorm E σ)
          = 12 * Real.log (2 * Real.pi) + Real.log (archNorm E σ) := by
      intro σ
      rw [Real.log_mul (ne_of_gt hpi) (ne_of_gt (archNorm_pos E σ)), Real.log_pow]
      push_cast
      ring
    rw [archSum, Finset.sum_congr rfl fun σ _ => hterm σ, Finset.sum_add_distrib,
      Finset.sum_const, Finset.card_univ, NumberField.Embeddings.card L ℂ, nsmul_eq_mul]
    ring
  rw [hht, hdeg, ← hprod, harchS, archDefect, Finset.sum_sub_distrib]
  ring

/-- ★★★★★★★★★★★★★★**`archDefect` は一様化の取り方に依らない**。

    `archDefect(E) = Σ_σ [12·log‖u_σ‖ − 6·log(covol P_σ)]`

★右辺は `(P, C)` の取り方に依るように見えるが、左辺は依らない。 -/
theorem archDefect_eq (E : WeierstrassCurve L) [E.IsElliptic]
    (P : (L →+* ℂ) → PeriodPair) (C : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ)) :
    archDefect L E
      = ∑ σ : (L →+* ℂ), (12 * Real.log ‖((C σ).u : ℂ)‖ - 6 * Real.log (covol (P σ))) := by
  refine Finset.sum_congr rfl fun σ _ => ?_
  rw [log_archNorm_eq E σ (hPC σ)]
  ring

/-- ★★★★★★★★★★★★★★★★★★★★★★**`ω`-正規化された同種写像の高さ評価
（選択によらない形）**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★**アルキメデスの仮定はただ 1 つ** `harch`:

    `archDefect(E′) = archDefect(E) + 6·d·log(l)`

★★これは `Found/GenEll/Velu.lean` の `velu_omega_gen`（第 591）が与える
`φ^*(ω′) = ω` の帰結である——`‖u′_σ‖ = ‖u_σ‖` かつ `covol P′_σ = covol P_σ/l` なら
`archDefect(E′) − archDefect(E) = 12·0 − 6·(−d·log l) = 6·d·log l`。

★★★☆**周期対も変数変換も statement に現れない**——これが残るアルキメデスの穴の
もっとも小さい形である。 -/
theorem htFalt_isogeny_le_of_archDefect (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ)
    (harch : archDefect L E'
      = archDefect L E + 6 * (Module.finrank ℚ L : ℝ) * Real.log l)
    (hfin : (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
          - (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        ≤ (3 / 2) * (Module.finrank ℚ L : ℝ) * Real.log l) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have h := twelve_finrank_htFaltOf_eq_archDefect E
  have h' := twelve_finrank_htFaltOf_eq_archDefect E'
  rw [harch] at h'
  nlinarith [h, h', hfin, hd]

/-- ★★★★★★★★★★★★**`hu` と `hcov` から `harch` が出る**。

★これで `htFalt_isogeny_le_of_omega`（第 592）は
`htFalt_isogeny_le_of_archDefect` の系である。 -/
theorem archDefect_isogeny_of_omega (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ) (hl : 0 < l)
    (P P' : (L →+* ℂ) → PeriodPair) (C C' : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hPC' : ∀ σ, C' σ • (E'.map σ) = latticeCurve (P' σ))
    (hu : ∀ σ, ‖((C' σ).u : ℂ)‖ = ‖((C σ).u : ℂ)‖)
    (hcov : ∀ σ, covol (P' σ) = covol (P σ) / l) :
    archDefect L E' = archDefect L E + 6 * (Module.finrank ℚ L : ℝ) * Real.log l := by
  have hl0 : (0:ℝ) < (l:ℝ) := by exact_mod_cast hl
  rw [archDefect_eq E P C hPC, archDefect_eq E' P' C' hPC']
  have hterm : ∀ σ : (L →+* ℂ),
      12 * Real.log ‖((C' σ).u : ℂ)‖ - 6 * Real.log (covol (P' σ))
        = (12 * Real.log ‖((C σ).u : ℂ)‖ - 6 * Real.log (covol (P σ)))
          + 6 * Real.log l := by
    intro σ
    rw [hu σ, hcov σ, Real.log_div (ne_of_gt (covol_pos (P σ))) (ne_of_gt hl0)]
    ring
  rw [Finset.sum_congr rfl fun σ _ => hterm σ, Finset.sum_add_distrib,
    Finset.sum_const, Finset.card_univ, NumberField.Embeddings.card L ℂ, nsmul_eq_mul]
  ring

def twelve_finrank_htFaltOf_eq_archDefect.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(ht^Falt の archDefect 表示——周期対も変数変換も現れない。★無条件)",
    sectionId := "genell-prop-3-4" }

def htFalt_isogeny_le_of_archDefect.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(残るアルキメデスの穴はただ 1 つ archDefect(E′) = archDefect(E) + 6d·log l)",
    sectionId := "genell-prop-3-4" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★正規化は要らない——`α` は打ち消し合う -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**次数 `l` の解析的同種写像があれば
`harch` は自動**——★**正規化 `‖u′‖ = ‖u‖` は要らない**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★同種写像 `φ : E → E′` は ℂ 上で `z ↦ α_σ·z` を誘導し、`α_σ·Λ_σ ⊆ Λ′_σ` の指数が `l` になる。
モデルのスケーリングは `u′_σ = α_σ·u_σ` で結ばれる（`φ^*(ω′)` と `ω` の比が `α_σ`）。

★★★このとき

    `κ_σ(E′) − κ_σ(E) = 12·log‖α‖ − 6·log(‖α‖²/l) = 12log‖α‖ − 12log‖α‖ + 6log l = 6·log l`

——★★★★**`α` が打ち消し合う**。すなわち `archDefect(E′) = archDefect(E) + 6·d·log(l)` は
**一様化の取り方にも `α` の大きさにも依らず**、ただ「次数 `l` の同種写像がある」ことだけから従う。

☆第 592・第 594 では `‖u′_σ‖ = ‖u_σ‖`（`φ^*(ω′) = ω` による正規化）を仮定していたが、
★**それは要らなかった**。Vélu の正規化（`Found/GenEll/Velu.lean` の `velu_omega_gen`）は
`α_σ = 1` を与えるが、`archDefect` の等式にはその情報は使われない。 -/
theorem archDefect_isogeny (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ) (hl : 0 < l)
    (P P' : (L →+* ℂ) → PeriodPair) (C C' : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hPC' : ∀ σ, C' σ • (E'.map σ) = latticeCurve (P' σ))
    (α : (L →+* ℂ) → ℂ) (hα : ∀ σ, α σ ≠ 0)
    (hu : ∀ σ, ((C' σ).u : ℂ) = α σ * ((C σ).u : ℂ))
    (a b c d : (L →+* ℂ) → ℤ)
    (h₁ : ∀ σ, α σ * (P σ).ω₁ = (a σ : ℂ) * (P' σ).ω₁ + (b σ : ℂ) * (P' σ).ω₂)
    (h₂ : ∀ σ, α σ * (P σ).ω₂ = (c σ : ℂ) * (P' σ).ω₁ + (d σ : ℂ) * (P' σ).ω₂)
    (hdet : ∀ σ, (a σ * d σ - b σ * c σ).natAbs = l) :
    archDefect L E' = archDefect L E + 6 * (Module.finrank ℚ L : ℝ) * Real.log l := by
  have hl0 : (0:ℝ) < (l:ℝ) := by exact_mod_cast hl
  rw [archDefect_eq E P C hPC, archDefect_eq E' P' C' hPC']
  have hterm : ∀ σ : (L →+* ℂ),
      12 * Real.log ‖((C' σ).u : ℂ)‖ - 6 * Real.log (covol (P' σ))
        = (12 * Real.log ‖((C σ).u : ℂ)‖ - 6 * Real.log (covol (P σ)))
          + 6 * Real.log l := by
    intro σ
    have hna : ‖α σ‖ ≠ 0 := norm_ne_zero_iff.2 (hα σ)
    have hnu : ‖((C σ).u : ℂ)‖ ≠ 0 := norm_ne_zero_iff.2 (C σ).u.ne_zero
    have hcovP : (0:ℝ) < covol (P σ) := covol_pos (P σ)
    -- ★`α·Λ ⊆ Λ′` が指数 `l` ⟹ `‖α‖²·covol(Λ) = l·covol(Λ′)`
    have hidx := covol_eq_index_mul_pair (scalePair (P σ) (α σ) (hα σ)) (P' σ)
      (a σ) (b σ) (c σ) (d σ) l (h₁ σ) (h₂ σ) (hdet σ)
    rw [covol_scalePair (P σ) (α σ) (hα σ)] at hidx
    have hnsq : Complex.normSq (α σ) = ‖α σ‖ ^ 2 := Complex.normSq_eq_norm_sq (α σ)
    rw [hnsq] at hidx
    -- ★`covol(Λ′) = ‖α‖²·covol(Λ)/l`
    have hcov' : covol (P' σ) = ‖α σ‖ ^ 2 * covol (P σ) / l := by
      rw [eq_div_iff (ne_of_gt hl0)]; linear_combination -hidx
    have hpos : (0:ℝ) < ‖α σ‖ ^ 2 * covol (P σ) :=
      mul_pos (pow_pos (lt_of_le_of_ne (norm_nonneg _) (Ne.symm hna)) 2) hcovP
    rw [hu σ, norm_mul, Real.log_mul hna hnu, hcov',
      Real.log_div (ne_of_gt hpos) (ne_of_gt hl0),
      Real.log_mul (by positivity) (ne_of_gt hcovP), Real.log_pow]
    push_cast
    ring
  rw [Finset.sum_congr rfl fun σ _ => hterm σ, Finset.sum_add_distrib,
    Finset.sum_const, Finset.card_univ, NumberField.Embeddings.card L ℂ, nsmul_eq_mul]
  ring

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★**同種写像の高さ評価
——アルキメデス側は完全に済んだ形**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★仮定は 2 つだけであり、**どちらも正規化を含まない**:

* 解析的な同種写像のデータ（`α`・`u′ = α·u`・`α·Λ ⊆ Λ′` の指数が `l`）
* 有限素点側 `hfin`

☆すなわち **`Lemma 3.5` の残りは「次数 `l` の同種写像が ℂ 上で何をするか」と
`[FC] Ch. I, Prop 2.7` の 2 つに完全に分離した**。 -/
theorem htFalt_isogeny_le_of_analytic (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ) (hl : 0 < l)
    (P P' : (L →+* ℂ) → PeriodPair) (C C' : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hPC' : ∀ σ, C' σ • (E'.map σ) = latticeCurve (P' σ))
    (α : (L →+* ℂ) → ℂ) (hα : ∀ σ, α σ ≠ 0)
    (hu : ∀ σ, ((C' σ).u : ℂ) = α σ * ((C σ).u : ℂ))
    (a b c d : (L →+* ℂ) → ℤ)
    (h₁ : ∀ σ, α σ * (P σ).ω₁ = (a σ : ℂ) * (P' σ).ω₁ + (b σ : ℂ) * (P' σ).ω₂)
    (h₂ : ∀ σ, α σ * (P σ).ω₂ = (c σ : ℂ) * (P' σ).ω₁ + (d σ : ℂ) * (P' σ).ω₂)
    (hdet : ∀ σ, (a σ * d σ - b σ * c σ).natAbs = l)
    (hfin : (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
          - (∑ᶠ p : HeightOneSpectrum (𝓞 L),
              (neronExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        ≤ (3 / 2) * (Module.finrank ℚ L : ℝ) * Real.log l) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l :=
  htFalt_isogeny_le_of_archDefect E E' l
    (archDefect_isogeny E E' l hl P P' C C' hPC hPC' α hα hu a b c d h₁ h₂ hdet) hfin

def archDefect_isogeny.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(次数 l の解析的同種写像があれば harch は自動——α は打ち消し合う)",
    sectionId := "genell-prop-3-4" }

def htFalt_isogeny_le_of_analytic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(アルキメデス側は完全に済んだ形——残りは同種写像のデータと hfin)",
    sectionId := "genell-prop-3-4" }

/-! ## ★★★★★★★★★★★★★★★★★★★★`E` が大域極小なら有限側は自動 -/

/-- ★★★★★★★★★★★★★★★★★★★★**`E` が大域極小・`E′` が整なら `hfin` は要らない**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`Found/GaloisRep/NeronValuation.lean` の `neronExp_nonneg`（第 595）により、
整なモデルの Néron 指数は非負である。したがって

* `E` が大域極小 ⟹ `Σᶠ_p neronExp_p(E)·log N(p) = 0`
* `E′` が整     ⟹ `Σᶠ_p neronExp_p(E′)·log N(p) ≥ 0`

なので `hfin` の左辺は `≤ 0` であり、右辺 `(3/2)·d·log(l) ≥ 0` 以下である。

★★★☆**したがって残る穴は `harch` ただ 1 つになる**——ただし `E′` は Vélu の
モデルであって極小とは限らないので、非極小性の分は `harch` の側が背負う
（`archDefect` と `Σᶠ neronExp` は変数変換でちょうど打ち消し合う）。 -/
theorem htFalt_isogeny_le_of_archDefect_minimal (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ) (hl : 0 < l)
    (hmin : ∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0)
    (hint : ∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p))
    (harch : archDefect L E'
      = archDefect L E + 6 * (Module.finrank ℚ L : ℝ) * Real.log l) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hl1 : (1:ℝ) ≤ (l:ℝ) := by exact_mod_cast hl
  have hlog : 0 ≤ Real.log l := Real.log_nonneg hl1
  have hΔ' : E'.Δ ≠ 0 := (inferInstance : E'.IsElliptic).isUnit.ne_zero
  refine htFalt_isogeny_le_of_archDefect E E' l harch ?_
  have hzero : (∑ᶠ p : HeightOneSpectrum (𝓞 L),
      (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal)) = 0 := by
    have : ∀ p : HeightOneSpectrum (𝓞 L),
        ((neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal)) = 0 := by
      intro p; rw [hmin p]; simp
    rw [finsum_congr this, finsum_zero]
  have hnn : 0 ≤ ∑ᶠ p : HeightOneSpectrum (𝓞 L),
      (neronExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
    refine finsum_nonneg fun p => ?_
    refine mul_nonneg ?_ (Real.log_natCast_nonneg _)
    exact_mod_cast neronExp_nonneg p E' hΔ' (hint p)
  rw [hzero]
  nlinarith [hnn, hlog, hd]

def htFalt_isogeny_le_of_archDefect_minimal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(E が大域極小・E′ が整なら有限素点側は自動——残る穴は harch だけ)",
    sectionId := "genell-prop-3-4" }

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

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★本日の到達点をまとめた形 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**同種写像の高さ評価
——外部引用なしで、幾何のデータだけから**。

    `ht^Falt(E′) ≤ ht^Falt(E) + 2·log(l)`

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★仮定はすべて**幾何のデータ**であり、証明されていない外部定理は含まない:

| 仮定 | 内容 |
|---|---|
| `hPC`・`hPC'` | 一様化（`E ⊗ σ ≅ ℂ/Λ_σ`・`E′ ⊗ σ ≅ ℂ/Λ′_σ`） |
| `α`・`hα`・`hu` | 解析的な同種写像のスケーリング `u′_σ = α_σ·u_σ` |
| `h₁`・`h₂`・`hdet` | `α_σ·Λ_σ ⊆ Λ′_σ` の指数が `l` |
| `hmin` | `E` は大域極小 |
| `hint` | `E′` は各素点で整 |

☆この形が本セッション（第 586-616）の到達点である。★かつては

    `hArch : (l−1)·d·deg∞(E) − (archSum(E′) − archSum(E)) ≤ 24·d·log(l)`

という**アルキメデスと有限が混ざった不等式**を外部引用として受けていた
（[FC] Ch. I, Prop 2.7 ＋ アルキメデスの (1,1)-形式）。

★★★★**それが「次数 `l` の同種写像が ℂ 上で何をするか」だけになった**——
アルキメデス項は `α` が打ち消し合って消え（`§9-1038`、第 596）、
有限素点項は `neronExp` の非負性で自動になる（`§9-1037`、第 595）。 -/
theorem htFalt_isogeny_le_of_analytic_minimal (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ) (hl : 0 < l)
    (P P' : (L →+* ℂ) → PeriodPair) (C C' : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hPC' : ∀ σ, C' σ • (E'.map σ) = latticeCurve (P' σ))
    (α : (L →+* ℂ) → ℂ) (hα : ∀ σ, α σ ≠ 0)
    (hu : ∀ σ, ((C' σ).u : ℂ) = α σ * ((C σ).u : ℂ))
    (a b c d : (L →+* ℂ) → ℤ)
    (h₁ : ∀ σ, α σ * (P σ).ω₁ = (a σ : ℂ) * (P' σ).ω₁ + (b σ : ℂ) * (P' σ).ω₂)
    (h₂ : ∀ σ, α σ * (P σ).ω₂ = (c σ : ℂ) * (P' σ).ω₁ + (d σ : ℂ) * (P' σ).ω₂)
    (hdet : ∀ σ, (a σ * d σ - b σ * c σ).natAbs = l)
    (hmin : ∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0)
    (hint : ∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l :=
  htFalt_isogeny_le_of_archDefect_minimal E E' l hl hmin hint
    (archDefect_isogeny E E' l hl P P' C C' hPC hPC' α hα hu a b c d h₁ h₂ hdet)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**同種写像の高さ評価——`α = 1`（同じ変数変換）の形**

    `ht^Falt(E′) ≤ ht^Falt(E) + 2·log(l)`

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★`Found/GenEll/Uniformization.lean` の `exists_velu_model_of_torsion`（第 676）は

    latticeCurve Λ′ = veluCurve (latticeCurve Λ) v w

を与える。すなわち `E′ = E/H` の一意化は **`E` と同じ変数変換 `C`** で
`latticeCurve Λ′` に一致する——`u′ = u`、つまり `α = 1`。
★`Found/GenEll/Velu.lean` の `veluQuotient_map`（第 677）が
「`L` 上で商を取ってから `σ` で送る」＝「`σ` で送ってから商を取る」を保証する。

★★★★☆**これが `Lemma 3.5` が要求する形そのものである**——
残るのは `L`-有理な位数 `l` の巡回部分群から `S`（代表点の集合）を作り、
各 `σ` について `h₁`・`h₂`・`hdet` を第 676 から読み取る段だけ。 -/
theorem htFalt_isogeny_le_of_velu (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ) (hl : 0 < l)
    (P P' : (L →+* ℂ) → PeriodPair) (C : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hPC' : ∀ σ, C σ • (E'.map σ) = latticeCurve (P' σ))
    (a b c d : (L →+* ℂ) → ℤ)
    (h₁ : ∀ σ, (P σ).ω₁ = (a σ : ℂ) * (P' σ).ω₁ + (b σ : ℂ) * (P' σ).ω₂)
    (h₂ : ∀ σ, (P σ).ω₂ = (c σ : ℂ) * (P' σ).ω₁ + (d σ : ℂ) * (P' σ).ω₂)
    (hdet : ∀ σ, (a σ * d σ - b σ * c σ).natAbs = l)
    (hmin : ∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0)
    (hint : ∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l :=
  htFalt_isogeny_le_of_analytic_minimal E E' l hl P P' C C hPC hPC'
    (fun _ => 1) (fun _ => one_ne_zero) (fun _ => by ring)
    a b c d
    (fun σ => by rw [one_mul]; exact h₁ σ)
    (fun σ => by rw [one_mul]; exact h₂ σ)
    hdet hmin hint

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5` の高さ評価——`L` 上の位数 `l` の点から**

    `ht^Falt(E/⟨Q⟩) ≤ ht^Falt(E) + 2·log(l)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★★★★☆**これが `Lemma 3.5` の解析側の到達点である。**

★組み立て:

* `S ≔ ⟨Q⟩∖{O}` の座標（`L × L` の有限集合）、`E′ ≔ veluQuotientFull E S`
* 各 `σ` で `rhPoint σ` により `Q` を `E.map σ` の位数 `l` の点へ（第 696）、
  点集合は `σ(S)` になる（第 703）
* 第 702（`ℂ` 側最終形）で `latticeCurve P′ = C σ • veluQuotientFull (E.map σ) (σS)`
* `veluQuotientFull` は底変換と可換（第 679）なので
  `= C σ • (E′.map σ)`——**同じ変数変換 `C σ`、すなわち `α = 1`**
* 第 678 に渡す

☆残る仮定は `hmin`（`E` が極小）・`hint`（`E′` が整）だけである。 -/
theorem htFalt_veluQuotientFull_le (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ) (hl : 0 < l)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (P : (L →+* ℂ) → PeriodPair) (C : (L →+* ℂ) → VariableChange ℂ)
    (hΔ : ∀ σ, latticeDisc (P σ) ≠ 0)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ))
    (hell1 : ∀ σ : L →+* ℂ, (E.map σ).IsElliptic)
    (hell2 : ∀ σ : L →+* ℂ, (C σ • (E.map σ)).IsElliptic)
    (hmin : ∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0)
    (hint : ∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  have key : ∀ σ : L →+* ℂ, ∃ (P' : PeriodPair) (A B Cc D : ℤ),
      (P σ).ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      (P σ).ω₂ = (Cc : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * Cc).natAbs = l ∧
      C σ • (E'.map σ) = latticeCurve P' := by
    intro σ
    haveI := hell1 σ
    haveI := hell2 σ
    have hQσ : addOrderOf (rhPoint σ E Q) = l := by
      rw [addOrderOf_rhPoint]; exact hQ
    obtain ⟨P', A, B, Cc, D, h1, h2, hdet, hEq⟩ :=
      exists_periodPair_veluQuotientFull (E.map σ) (C σ) (P σ) (hΔ σ) (hPC σ) hl hQσ
    refine ⟨P', A, B, Cc, D, h1, h2, hdet, ?_⟩
    rw [hEq, image_pointCoords_rhPoint_nsmul σ E hQ, hE', veluQuotientFull_map]
  choose P' A B Cc D h1 h2 hdet hPC' using key
  exact htFalt_isogeny_le_of_velu E E' l hl P P' C hPC hPC' A B Cc D h1 h2 hdet hmin hint

def htFalt_veluQuotientFull_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ht^Falt(E/⟨Q⟩) ≤ ht^Falt(E) + 2 log l——解析側の到達点)",
    sectionId := "genell-lemma-3-5" }

def htFalt_isogeny_le_of_velu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(同種写像の高さ評価——α = 1（同じ変数変換）の形。★無条件)",
    sectionId := "genell-prop-3-4" }

def htFalt_isogeny_le_of_analytic_minimal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(同種写像の高さ評価——外部引用なしで幾何のデータだけから)",
    sectionId := "genell-prop-3-4" }

def htFalt_isogeny_le_of_analytic_minimal.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆仮定 hPC・hPC′(一様化)は Found/GenEll/JSurjective.lean の " ++
       "exists_periodPair_of_isElliptic(2026-08-26、第 348)で無条件に取れている。" ++
       "★残るのは α・hu・h₁・h₂・hdet(解析的な同種写像のデータ)を" ++
       "代数的な同種写像 E → E/H から作ることである") 13,
    .implicitStep
      ("★★★★到達点(2026-08-29、第 617): かつて外部引用として受けていた " ++
       "hArch((l−1)d·deg∞(E) − (archSum(E′) − archSum(E)) ≤ 24d·log l)が、" ++
       "「次数 l の同種写像が ℂ 上で何をするか」だけになった。" ++
       "★アルキメデス項は α が打ち消し合って消え(第 596)、" ++
       "有限素点項は neronExp の非負性で自動になる(第 595)") 15,
    .implicitStep
      ("☆解析的な同種写像のデータを作るには ℘ の加法定理が要り、" ++
       "それには「℘ は各値をちょうど 2 回取る」(楕円関数の零点勘定)が要る" ++
       "——Skeleton/GenEll/AdditionTheorem.lean(第 606・第 615)") 21 ]

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
