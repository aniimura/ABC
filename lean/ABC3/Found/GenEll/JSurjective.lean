import ABC3.Found.GenEll.JFunction
import Mathlib.NumberTheory.Modular
import Mathlib.Analysis.Complex.OpenMapping
import Mathlib.AlgebraicGeometry.EllipticCurve.IsomOfJ

/-!
# GenEll 第 348 ブロック —— **★★★★★★★★一意化が閉じた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★到達点

> **任意の楕円曲線 `W/ℂ` は、ある周期束の曲線と変数変換で移り合う**
> (`exists_periodPair_of_isElliptic`)

★これは `Skeleton/GenEll/Uniformization.lean` の `exists_periodPair` そのものであり、
本ブロックでその `sorry` が消える。

## ★★★★★★`j` の全射性——**開かつ閉**

| 段 | 補題 |
|---|---|
| 像は**開** | `isOpen_range_jFun`(`AnalyticOnNhd.is_constant_or_isOpen` + `jC_not_constant`) |
| 像は**閉** | `isClosed_range_jFun`(基本領域に移す + 切り詰めのコンパクト性) |
| `ℂ` は連結 | `isClopen_iff` |

★★★閉であることの証明が本ブロックの核である:
`wₙ → w` で `wₙ = j(τₙ)` とする。★`j` は `SL(2,ℤ)` 不変なので `τₙ` は基本領域 `𝒟` に移せる。
★★`‖wₙ‖` は有界で、`‖j‖ → ∞`(カスプ)だから `im τₙ` も有界。
★★★したがって `τₙ` は**コンパクト**な `𝒟 ∩ {im ≤ M}` に入り、部分列が収束する。
`j` の連続性で `w = j(極限) ∈ 像`。

★★★★★**valence 公式(留数計算)を使わずに済んだ**——
mathlib の `isCompact_truncatedFundamentalDomain` が効いた。

## ★★★★★★最後の 5 行

mathlib の **`WeierstrassCurve.exists_variableChange_of_j_eq`**
(分離閉体上で `j` が等しければ変数変換で移る)に、
第 346 の `latticeCurve_j_tauPair`(`j(ℤ+τℤ) = E₄³/Δ`)と
`isElliptic_latticeCurve'`(格子曲線はつねに楕円曲線)を渡すだけである。

## ★★スケルトンの見積もりとの差

第 331 の `.needs` は本節点を **25-60 ブロック**と見積もっていた
(「上流(mathlib)に入るべき仕事であり、Galois の義務の中では閉じない」)。
★★★実際には **第 332-348 の 17 ブロック**で閉じた。
★★★★差の大半は **(ii) `ℂ/Λ ≅ E(ℂ)` の群同型は要らなかった**ことによる
——`exists_periodPair` は**曲線の同型**であって群同型ではない。
★見積もりが (ii) を要求に数えていたのが過大の原因である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `range_jFun_eq` | ★`j` の像 = `jC` による上半平面の像 |
| `isOpen_range_jFun` | ★★★★★像は開 |
| `isClosed_range_jFun` | ★★★★★★像は閉 |
| `jFun_surjective` | ★★★★★★★**`j : ℍ → ℂ` は全射** |
| `exists_periodPair_of_isElliptic` | ★★★★★★★★**一意化** |
-/

namespace ABC3.Found.GenEll

open Complex Real ModularForm MatrixGroups UpperHalfPlane EisensteinSeries ABC3.Found.GaloisRep
open Filter Topology

/-! ## ★★★★★像は開 -/

theorem range_jFun_eq : Set.range jFun = jC '' UpperHalfPlane.upperHalfPlaneSet := by
  ext w
  constructor
  · rintro ⟨τ, rfl⟩
    exact ⟨(τ : ℂ), τ.2, jC_apply τ.2⟩
  · rintro ⟨z, hz, rfl⟩
    exact ⟨⟨z, hz⟩, (jC_apply hz).symm⟩

/-- ★★★★★**`j` の像は開**——開写像定理。 -/
theorem isOpen_range_jFun : IsOpen (Set.range jFun) := by
  rcases analyticOnNhd_jC.is_constant_or_isOpen ((convex_halfSpace_im_gt 0).isPreconnected)
    with h | h
  · exact absurd h jC_not_constant
  · rw [range_jFun_eq]
    exact h _ subset_rfl UpperHalfPlane.isOpen_upperHalfPlaneSet

/-! ## ★★★★★★像は閉 -/

/-- ★★★★★★**`j` の像は閉**。

★`SL(2,ℤ)` 不変性で基本領域に移し、`‖j‖ → ∞` で `im` を抑え、
`isCompact_truncatedFundamentalDomain` で部分列を取る。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
theorem isClosed_range_jFun : IsClosed (Set.range jFun) := by
  refine IsSeqClosed.isClosed ?_
  intro wn w hmem hlim
  choose τ hτ using hmem
  have hex : ∀ n, ∃ σ : ℍ, σ ∈ ModularGroup.fd ∧ jFun σ = wn n := by
    intro n
    obtain ⟨g, hg⟩ := ModularGroup.exists_smul_mem_fd (τ n)
    exact ⟨g • τ n, hg, by rw [jFun_smul]; exact hτ n⟩
  choose σ hσfd hjσ using hex
  obtain ⟨C, hC⟩ : ∃ C, ∀ n, ‖wn n‖ ≤ C := by
    obtain ⟨C, hC⟩ := (hlim.norm).bddAbove_range
    exact ⟨C, fun n => hC ⟨n, rfl⟩⟩
  obtain ⟨M, hM⟩ := (UpperHalfPlane.atImInfty_mem _).1
    (tendsto_norm_jFun_atImInfty.eventually_gt_atTop C)
  have hin : ∀ n, σ n ∈ ModularGroup.truncatedFundamentalDomain M := by
    intro n
    refine ⟨hσfd n, ?_⟩
    by_contra hlt
    have h2 : C < ‖jFun (σ n)‖ := hM (σ n) (le_of_lt (by linarith [not_le.1 hlt]))
    rw [hjσ n] at h2
    linarith [hC n]
  obtain ⟨s, _, φ, hφ, hconv⟩ :=
    (ModularGroup.isCompact_truncatedFundamentalDomain M).tendsto_subseq hin
  refine ⟨s, ?_⟩
  have h1 : Tendsto (fun n => jFun (σ (φ n))) atTop (𝓝 (jFun s)) :=
    (continuous_jFun.continuousAt).tendsto.comp hconv
  have h2 : Tendsto (fun n => wn (φ n)) atTop (𝓝 w) := hlim.comp hφ.tendsto_atTop
  have h1' : Tendsto (fun n => wn (φ n)) atTop (𝓝 (jFun s)) := by
    simpa [hjσ] using h1
  exact tendsto_nhds_unique h1' h2

/-! ## ★★★★★★★`j` の全射性 -/

/-- ★★★★★★★**`j : ℍ → ℂ` は全射**——像は空でない開かつ閉で、`ℂ` は連結。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
theorem jFun_surjective : Function.Surjective jFun := by
  have hne : (Set.range jFun).Nonempty := ⟨jFun UpperHalfPlane.I, ⟨_, rfl⟩⟩
  rcases isClopen_iff.1 ⟨isClosed_range_jFun, isOpen_range_jFun⟩ with h | h
  · exact absurd h (Set.nonempty_iff_ne_empty.1 hne)
  · exact Set.range_eq_univ.mp h

/-! ## ★★★★★★★★一意化 -/

/-- ★★★★★★★★**複素楕円曲線の一意化**——任意の `W/ℂ` はある周期束の曲線と
変数変換で移り合う。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`j` の全射性(本ブロック)+ `j(ℤ+τℤ) = E₄³/Δ`(第 346)
+ mathlib の `exists_variableChange_of_j_eq`。 -/
theorem exists_periodPair_of_isElliptic (W : WeierstrassCurve ℂ) (hell : W.IsElliptic) :
    ∃ (P : PeriodPair) (C : WeierstrassCurve.VariableChange ℂ), C • W = latticeCurve P := by
  haveI := hell
  obtain ⟨τ, hτ⟩ := jFun_surjective W.j
  haveI : (latticeCurve (tauPair τ)).IsElliptic := isElliptic_latticeCurve' _
  have hj : W.j = (latticeCurve (tauPair τ)).j := by
    rw [← hτ, jFun_eq_latticeCurve_j]
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_variableChange_of_j_eq W (latticeCurve (tauPair τ)) hj
  exact ⟨tauPair τ, C, hC⟩

/-! ## ★出典の紐付け(`.src`) -/

def jFun_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def exists_periodPair_of_isElliptic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def isClosed_range_jFun.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
