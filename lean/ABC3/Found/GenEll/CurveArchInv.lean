import ABC3.Found.GenEll.LatticeFromG
import ABC3.Found.GenEll.JSurjective

/-!
# GenEll 第 354 ブロック —— **★★★★★★★曲線のアルキメデス不変量**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★到達点——`archInv` が**曲線の関数**になった

第 349-353 で `archInv = ‖Δ_lat‖·covol⁶` について

* 相似で不変(第 349)
* 束だけで決まる(第 350)
* **`(g₂,g₃)` で決まる**(第 353)

まで来た。★本ブロックはこれを使って、**周期束の選び方に依らない**

> `curveArchInv : WeierstrassCurve ℂ → ℝ`

を建てる。★★一意化(第 348)で `C • W = latticeCurve P` なる `P` は存在し、
2 つの取り方 `(P,C)`・`(P',C')` に対しては `D := C'C⁻¹` の `u` について
`g₂(P') = u⁻⁴g₂(P)`・`g₃(P') = u⁻⁶g₃(P)`——すなわち `P'` と `scalePair P u` が
**同じ `(g₂,g₃)`** を持つので、第 353 で束が一致し、`archInv` が一致する。

## ★★★★★★これは `j` だけの関数である

`curveArchInv` は変数変換で不変(`curveArchInv_variableChange`)、
したがって **`j` だけで決まる**(`curveArchInv_congr_j`)。★実際

    curveArchInv(W) = 4096·π¹²·‖Δ(τ)‖·(Im τ)⁶      (`j(τ) = j(W)` なる τ)

で、右辺は `Δ` の **Petersson ノルム**である(`exists_curveArchInv_eq_petersson`)。
★★これは正しい——Faltings 高さのアルキメデス局所因子は**同型不変量**であり、
モデル依存性は**有限素点側**(`Δ_min`)にしか現れないからである。

## ★★★★界面へ渡す形

数体 `L` と埋め込み `σ : L →+* ℂ` に対し `archNorm E σ := curveArchInv (E.map σ)`。
★これは**正**であり(`archNorm_pos`)、**変数変換で不変**(`archNorm_variableChange`)。
★★★`FaltingsHeightData` に「`ht^Falt` が `archNorm` を通じて曲線に縛られる」条件を
足すための材料がこれでそろった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `latticeCurve_c₆` | ★`c₆ = 216·g₃` |
| `archInv_eq_of_variableChange` | ★★★★★★**周期対の取り方に依らない** |
| `curveArchInv` | ★★★★★★**曲線のアルキメデス不変量** |
| `curveArchInv_pos`・`curveArchInv_variableChange` | ★★★★正・変数変換不変 |
| `curveArchInv_congr_j` | ★★★★★★**`j` だけで決まる** |
| `exists_curveArchInv_eq_petersson` | ★★★★★★★**Petersson ノルムである** |
| `archNorm`・`archNorm_pos`・`archNorm_variableChange` | ★★★★★界面へ渡す形 |
-/

namespace ABC3.Found.GenEll

open Complex PeriodPair WeierstrassCurve UpperHalfPlane ABC3.Found.GaloisRep

/-! ## ★周期束の取り方に依らないこと -/

theorem latticeCurve_c₆ (P : PeriodPair) : (latticeCurve P).c₆ = 216 * P.g₃ := by
  simp only [latticeCurve, WeierstrassCurve.c₆, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆]
  ring

/-- ★★★★★★**`archInv` は周期対の取り方に依らない**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`D := C'C⁻¹` の `u` について `P'` と `scalePair P u` が同じ `(g₂,g₃)` を持つ。 -/
theorem archInv_eq_of_variableChange {W : WeierstrassCurve ℂ} {P P' : PeriodPair}
    {C C' : WeierstrassCurve.VariableChange ℂ}
    (hC : C • W = latticeCurve P) (hC' : C' • W = latticeCurve P') :
    archInv P = archInv P' := by
  set D : WeierstrassCurve.VariableChange ℂ := C' * C⁻¹ with hDdef
  have hDact : latticeCurve P' = D • latticeCurve P := by
    rw [← hC, ← hC', hDdef, mul_smul, inv_smul_smul]
  have hu : ((D.u : ℂ)) ≠ 0 := D.u.ne_zero
  have h4 : P'.g₂ = ((D.u : ℂ) ^ 4)⁻¹ * P.g₂ := by
    have h := congrArg WeierstrassCurve.c₄ hDact
    rw [latticeCurve_c₄, WeierstrassCurve.variableChange_c₄, latticeCurve_c₄] at h
    have hinv : ((D.u⁻¹ : ℂˣ) : ℂ) = ((D.u : ℂ))⁻¹ := by simp
    rw [hinv] at h
    field_simp at h ⊢
    linear_combination h
  have h6 : P'.g₃ = ((D.u : ℂ) ^ 6)⁻¹ * P.g₃ := by
    have h := congrArg WeierstrassCurve.c₆ hDact
    rw [latticeCurve_c₆, WeierstrassCurve.variableChange_c₆, latticeCurve_c₆] at h
    have hinv : ((D.u⁻¹ : ℂˣ) : ℂ) = ((D.u : ℂ))⁻¹ := by simp
    rw [hinv] at h
    field_simp at h ⊢
    linear_combination h
  have hg2 : (scalePair P (D.u : ℂ) hu).g₂ = P'.g₂ := by rw [g₂_scalePair, h4]
  have hg3 : (scalePair P (D.u : ℂ) hu).g₃ = P'.g₃ := by rw [g₃_scalePair, h6]
  have hlat := lattice_eq_of_g₂_g₃_eq hg2 hg3
  rw [← archInv_scalePair P (D.u : ℂ) hu, archInv_congr hlat]

/-! ## ★★★★★★曲線のアルキメデス不変量 -/

open Classical in
/-- ★★★★★★**曲線のアルキメデス不変量**——周期束の取り方に依らない。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★Faltings 高さのアルキメデス局所因子 `‖Δ‖_arch` である。 -/
noncomputable def curveArchInv (W : WeierstrassCurve ℂ) : ℝ :=
  if h : W.IsElliptic then archInv (exists_periodPair_of_isElliptic W h).choose else 0

theorem curveArchInv_eq {W : WeierstrassCurve ℂ} (h : W.IsElliptic) {P : PeriodPair}
    {C : WeierstrassCurve.VariableChange ℂ} (hC : C • W = latticeCurve P) :
    curveArchInv W = archInv P := by
  classical
  rw [curveArchInv, dif_pos h]
  obtain ⟨C₀, hC₀⟩ := (exists_periodPair_of_isElliptic W h).choose_spec
  exact archInv_eq_of_variableChange hC₀ hC

theorem curveArchInv_latticeCurve (P : PeriodPair) : curveArchInv (latticeCurve P) = archInv P :=
  curveArchInv_eq (isElliptic_latticeCurve' P) (C := 1) (by simp)

theorem curveArchInv_pos {W : WeierstrassCurve ℂ} (h : W.IsElliptic) : 0 < curveArchInv W := by
  obtain ⟨P, C, hC⟩ := exists_periodPair_of_isElliptic W h
  rw [curveArchInv_eq h hC]
  exact archInv_pos P

/-- ★★★★★**変数変換で変わらない**。 -/
theorem curveArchInv_variableChange (W : WeierstrassCurve ℂ) (h : W.IsElliptic)
    (C : WeierstrassCurve.VariableChange ℂ) : curveArchInv (C • W) = curveArchInv W := by
  haveI := h
  obtain ⟨P, C₀, hC₀⟩ := exists_periodPair_of_isElliptic W h
  have h2 : (C₀ * C⁻¹) • (C • W) = latticeCurve P := by
    rw [mul_smul, inv_smul_smul, hC₀]
  rw [curveArchInv_eq (inferInstance : (C • W).IsElliptic) h2, curveArchInv_eq h hC₀]

/-! ## ★★★★★★`j` だけで決まること -/

theorem curveArchInv_tauPair (τ : UpperHalfPlane) :
    curveArchInv (latticeCurve (tauPair τ)) = 4096 * Real.pi ^ 12 * peterssonDelta τ := by
  rw [curveArchInv_latticeCurve, archInv_tauPair, peterssonDelta]
  ring

/-- ★★★★★★**`curveArchInv` は `j` だけで決まる**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★アルキメデス局所因子は同型不変量である。 -/
theorem curveArchInv_congr_j {W W' : WeierstrassCurve ℂ} (h : W.IsElliptic) (h' : W'.IsElliptic)
    (hj : W.j = W'.j) : curveArchInv W = curveArchInv W' := by
  haveI := h
  haveI := h'
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_variableChange_of_j_eq W W' hj
  obtain ⟨P, C₀, hC₀⟩ := exists_periodPair_of_isElliptic W' h'
  rw [curveArchInv_eq h' hC₀,
    curveArchInv_eq h (show (C₀ * C) • W = latticeCurve P by rw [mul_smul, hC, hC₀])]

/-- ★★★★★★★**`curveArchInv` は `Δ` の Petersson ノルムである**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`j(τ) = j(W)` なる `τ` を取れば `curveArchInv W = 4096π¹²·‖Δ(τ)‖(Im τ)⁶`。 -/
theorem exists_curveArchInv_eq_petersson (W : WeierstrassCurve ℂ) (h : W.IsElliptic) :
    ∃ τ : UpperHalfPlane, jFun τ = W.j ∧
      curveArchInv W = 4096 * Real.pi ^ 12 * peterssonDelta τ := by
  obtain ⟨τ, hτ⟩ := jFun_surjective W.j
  refine ⟨τ, hτ, ?_⟩
  haveI := h
  haveI : (latticeCurve (tauPair τ)).IsElliptic := isElliptic_latticeCurve' _
  have hj : W.j = (latticeCurve (tauPair τ)).j := by rw [← hτ, jFun_eq_latticeCurve_j]
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_variableChange_of_j_eq W (latticeCurve (tauPair τ)) hj
  rw [curveArchInv_eq h hC, archInv_tauPair, peterssonDelta]
  ring

/-! ## ★★★★★界面へ渡す形 -/

/-- ★★★★★**アルキメデス素点でのノルム**——埋め込み `σ : L →+* ℂ` に沿った `‖Δ‖_arch`。 -/
noncomputable def archNorm {L : Type*} [Field L] (E : WeierstrassCurve L) (σ : L →+* ℂ) : ℝ :=
  curveArchInv (E.map σ)

theorem archNorm_pos {L : Type*} [Field L] (E : WeierstrassCurve L) [E.IsElliptic]
    (σ : L →+* ℂ) : 0 < archNorm E σ :=
  curveArchInv_pos (inferInstance : (E.map σ).IsElliptic)

/-- ★★★★★**変数変換で変わらない**——`ht^Falt` の界面条件に使える形。 -/
theorem archNorm_variableChange {L : Type*} [Field L] (E : WeierstrassCurve L) [E.IsElliptic]
    (C : WeierstrassCurve.VariableChange L) (σ : L →+* ℂ) :
    archNorm (C • E) σ = archNorm E σ := by
  rw [archNorm, archNorm, ← WeierstrassCurve.map_variableChange E C σ]
  exact curveArchInv_variableChange _ (inferInstance : (E.map σ).IsElliptic) _

/-! ## ★出典の紐付け(`.src`) -/

def curveArchInv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def archInv_eq_of_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def archNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
