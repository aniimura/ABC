/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HyperplaneHeight
import ABC3.Found.GenEll.HyperplanePullback
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★Fubini–Study の Green 関数と無限素点（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★★これは何か —— `§9-867` の仮定 `hgreen` を潰す

`§9-867` は `htArith F D̄ xF = log H(x)/[F:ℚ]` を、仮定

    `hgreen : D.green (archPoint xF v) = log( sup_k v(x_k) / v(x₀) )`

の下で取った。★本ファイルはその仮定を**チャートの形で計算して潰す**。

## ★★★機構 —— 素点と複素点の対応は 2 行

★`archPoint xF v = archSpecMap v ≫ xF` であり、`xF = Spec ψ ≫ chartA i₀` なら
`Spec` の反変性で

    `archPoint xF v = Spec (ψ ≫ archRingHom v) ≫ chartA i₀`

——つまり**同じチャートを通り、環の射が `archRingHom v ∘ ψ` になるだけ**である。

★★あとは `‖archRingHom F v a‖ = v(a)`（`InfinitePlace.norm_embedding_eq`）と
`v(w_k) = v(x_k)/v(x_{i₀})` で、正規化座標の `sup` が

    `sup_k v(x_k)/v(x_{i₀})  /  (v(x₀)/v(x_{i₀}))  =  sup_k v(x_k)/v(x₀)`

に潰れる。★★★`v(x_{i₀})` が**約分される**のが要点で、これが
「Fubini–Study は射影的に不変」ということである。

## ★残っている段（明示）

★★★★本ファイルが取るのは**チャートを固定した形**の Green 関数
（`greenChartOf`）である。`ArithCartier` の欄に入れる**大域的な** `green` を作るには
「別のチャートで測っても同じ値になる」ことが要る——それは上と同じ約分で出るが、
`projChartHom` を 2 つのチャートで比べる段が別に要る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField MvPolynomial HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★無限素点と複素埋め込み -/

/-- ★**`archRingHom` のノルムは素点の値である**。 -/
theorem norm_archRingHom (F : Type) [Field F] [NumberField F] (v : InfinitePlace F) (a : 𝓞 F) :
    ‖archRingHom F v a‖ = v ((a : F)) := by
  show ‖v.embedding ((algebraMap (𝓞 F) F) a)‖ = v ((a : F))
  rw [InfinitePlace.norm_embedding_eq]

/-- ★**正の定数で割る操作は `⨆` と交換する**。 -/
theorem iSup_div_const {ι : Type} (f : ι → ℝ) (c : ℝ) (hc : 0 ≤ c) :
    (⨆ k, f k) / c = ⨆ k, (f k / c) := by
  simp only [div_eq_mul_inv]
  exact Real.iSup_mul_of_nonneg (inv_nonneg.2 hc) f

/-- ★★★★★**正規化座標で測っても Fubini–Study は変わらない**。

    `log( sup_k ‖σ_v(x_k/x_{i₀})‖ / ‖σ_v(x₀/x_{i₀})‖ ) = log( sup_k v(x_k) / v(x₀) )`

★`v(x_{i₀})` が**約分される**のが要点である。 -/
theorem log_iSup_norm_eq (F : Type) [Field F] [NumberField F] {ι : Type}
    (v : InfinitePlace F) (x w : ι → 𝓞 F) (i₀ j : ι)
    (hw : ∀ k, x k = w k * x i₀) (hi0 : x i₀ ≠ 0) :
    Real.log ((⨆ k, ‖archRingHom F v (w k)‖) / ‖archRingHom F v (w j)‖)
      = Real.log ((⨆ k, v ((x k : F))) / v ((x j : F))) := by
  have hi0pos : (0:ℝ) < v ((x i₀ : F)) := (AbsoluteValue.pos_iff _).2 (by simpa using hi0)
  have hwv : ∀ k, v ((w k : F)) = v ((x k : F)) / v ((x i₀ : F)) := by
    intro k
    have hk : v ((x k : F)) = v ((w k : F)) * v ((x i₀ : F)) := by
      rw [hw k]; push_cast; exact map_mul v _ _
    rw [hk, mul_div_assoc, div_self hi0pos.ne', mul_one]
  simp only [norm_archRingHom, hwv]
  rw [← iSup_div_const (fun k => v ((x k : F))) _ hi0pos.le,
    div_div_div_cancel_right₀ hi0pos.ne']

/-! ## ★★複素点はチャートを通ったまま -/

/-- ★★**無限素点が定める複素点は同じチャートを通る**。

★`Spec` の反変性だけである——環の射が `archRingHom v ∘ ψ` になるにすぎない。 -/
theorem archPoint_chart_factor (F : Type) [Field F] [NumberField F] (N : ℕ) (i₀ : Fin (N+1))
    (ψ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i₀))
      ⟶ CommRingCat.of (NumberField.RingOfIntegers F))
    (v : InfinitePlace F) :
    archPoint (Spec.map ψ ≫ chartA N ℤ i₀) v
      = Spec.map (ψ ≫ CommRingCat.ofHom (archRingHom F v)) ≫ chartA N ℤ i₀ := by
  rw [archPoint, archSpecMap, ← Category.assoc, ← Spec.map_comp]

/-! ## ★★★★★チャートを固定した Fubini–Study の Green 関数 -/

/-- ★★★★★**Fubini–Study の Green 関数**（チャートを固定した形）。

    `g(φ) = log( sup_k ‖φ(x_k/x_i)‖ / ‖φ(x₀/x_i)‖ )`

★超平面 `{x₀ = 0}` に対する `−log‖x₀‖_FS` である。 -/
noncomputable def greenChartOf (N : ℕ) (i : Fin (N+1))
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ) : ℝ :=
  Real.log ((⨆ k, ‖φ.hom (projCoord N ℤ i k)‖) / ‖φ.hom (projCoord N ℤ i 0)‖)

/-- ★★★★★★**`§9-867` の仮定 `hgreen` の中身**。

    `greenChartOf N i₀ (ψ ≫ archRingHom v) = log( sup_k v(x_k) / v(x₀) )` -/
theorem greenChartOf_archRingHom (F : Type) [Field F] [NumberField F]
    (N : ℕ) (i₀ : Fin (N+1))
    (ψ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i₀))
      ⟶ CommRingCat.of (NumberField.RingOfIntegers F))
    (v : InfinitePlace F) (x : Fin (N+1) → 𝓞 F)
    (hw : ∀ k, x k = ψ.hom (projCoord N ℤ i₀ k) * x i₀) (hi0 : x i₀ ≠ 0) :
    greenChartOf N i₀ (ψ ≫ CommRingCat.ofHom (archRingHom F v))
      = Real.log ((⨆ k, v ((x k : F))) / v ((x 0 : F))) := by
  have hcomp : ∀ a, (ψ ≫ CommRingCat.ofHom (archRingHom F v)).hom a
      = archRingHom F v (ψ.hom a) := by
    intro a
    rw [CommRingCat.hom_comp, CommRingCat.hom_ofHom]
    rfl
  rw [greenChartOf]
  simp only [hcomp]
  exact log_iSup_norm_eq F v x (fun k => ψ.hom (projCoord N ℤ i₀ k)) i₀ 0 hw hi0

/-! ## ★★★★★★★★★★★★組み立て -/

/-- ★★★★★★★★★★★★**超平面因子の高さは素朴高さである**（Green 関数まで込み）。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `htArith F D̄ (Spec ψ ≫ chartA i₀) = log H(x) / [F : ℚ]`

★`§9-867` の仮定 `hgreen` が `greenChartOf_archRingHom` で潰れた形である。
★★残る仮定は `hdiv`（因子は超平面）と `hgreenChart`
（`D.green` はチャートの上で Fubini–Study である）だけで、
後者は**大域的な `green` を作る段**（チャート独立性）に置き換わる。 -/
theorem htArith_hyperplane_eq_log_mulHeight (F : Type) [Field F] [NumberField F]
    (N : ℕ) (i₀ : Fin (N+1))
    (ψ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i₀))
      ⟶ CommRingCat.of (NumberField.RingOfIntegers F))
    (D : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (hdiv : D.divisor = hyperplaneIdeal N ℤ)
    (x : Fin (N+1) → 𝓞 F) (hx : x ≠ 0)
    (hw : ∀ k, x k = ψ.hom (projCoord N ℤ i₀ k) * x i₀)
    (h0 : x 0 ≠ 0)
    (hgreenChart : ∀ v : InfinitePlace F,
      D.green (archPoint (Spec.map ψ ≫ chartA N ℤ i₀) v)
        = greenChartOf N i₀ (ψ ≫ CommRingCat.ofHom (archRingHom F v))) :
    htArith F D (Spec.map ψ ≫ chartA N ℤ i₀)
      = Real.log (Height.mulHeight (fun k => (x k : F))) / (Module.finrank ℚ F : ℝ) := by
  have hi0 : x i₀ ≠ 0 := by
    intro hc
    exact h0 (by rw [hw 0, hc, mul_zero])
  have hpull : pullbackIdeal F D.divisor (Spec.map ψ ≫ chartA N ℤ i₀)
      = Ideal.span {ψ.hom (projCoord N ℤ i₀ 0)} := by
    rw [hdiv]
    exact pullbackIdeal_hyperplane_point F N i₀ ψ _ rfl
  refine htArith_eq_log_mulHeight F D _ x hx i₀ 0
    (ψ.hom (projCoord N ℤ i₀ 0)) (hw 0) h0 ?_ hpull ?_
  · exact span_range_eq_span_of_dvd F x i₀
      (fun k => ⟨ψ.hom (projCoord N ℤ i₀ k), hw k⟩)
  · intro v
    rw [hgreenChart v, greenChartOf_archRingHom F N i₀ ψ v x hw hi0]

/-! ## ★出典の紐付け(`.src`) -/

def norm_archRingHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(archRingHom のノルムは素点の値である)",
    sectionId := "genell-prop-1-4" }

def log_iSup_norm_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(正規化座標で測っても Fubini–Study は変わらない)",
    sectionId := "genell-prop-1-4" }

def archPoint_chart_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(無限素点が定める複素点は同じチャートを通る)",
    sectionId := "genell-prop-1-4" }

def greenChartOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(Fubini–Study の Green 関数——チャートを固定した形)",
    sectionId := "genell-prop-1-4" }

def htArith_hyperplane_eq_log_mulHeight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(超平面因子の高さは素朴高さである——Green 関数まで込み)",
    sectionId := "genell-prop-1-4" }

def htArith_hyperplane_eq_log_mulHeight.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htArith_eq_log_mulHeight(段 C2c-5、§9-867)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_eq_log_mulHeight") 2,
    .citation "[ABC3]" "pullbackIdeal_hyperplane_point(段 C2c-1、§9-865)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pullbackIdeal_hyperplane_point") 2,
    .citation "[mathlib]" "InfinitePlace.norm_embedding_eq(‖σ_v(a)‖ = v(a))"
      (.inMathlib "NumberField.InfinitePlace.norm_embedding_eq") 2,
    .implicitStep
      ("★素点と複素点の対応は 2 行である: archPoint xF v = archSpecMap v ≫ xF で、" ++
       "xF = Spec ψ ≫ chartA i₀ なら Spec の反変性で " ++
       "archPoint xF v = Spec (ψ ≫ archRingHom v) ≫ chartA i₀" ++
       "——**同じチャートを通り、環の射が archRingHom v ∘ ψ になるだけ**である") 3,
    .implicitStep
      ("★★v(x_{i₀}) が**約分される**のが要点で、これが「Fubini–Study は射影的に不変」" ++
       "ということである(log_iSup_norm_eq)") 3,
    .implicitStep
      ("★★★残るのは**大域的な green を作る段**である: 本ファイルの greenChartOf は" ++
       "チャートを固定した形なので、ArithCartier の欄に入れるには" ++
       "「別のチャートで測っても同じ値」(チャート独立性)が要る。" ++
       "★同じ約分で出るが、projChartHom を 2 つのチャートで比べる段が別に要る") 4 ]

end ABC3.Found.GenEll
