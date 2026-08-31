/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.Def12Complete
import ABC3.Found.Arakelov.UltraCompact
import ABC3.Found.Arakelov.PullbackNorm
import ABC3.Found.Arakelov.SpanPullSec
import ABC3.Found.GenEll.ArchBound
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★[GenEll] Proposition 1.4, (ii) —— **算術直線束の側で**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6–7。

原文 (GenEll p.6):
> by global sections [for instance, if the line bundle LQ is ample], then

## ★★★★★★★★★★これは何か

原文 (ii) の結論は「高さが**下に有界**」である。原文の証明はこう書く:

原文 (GenEll p.7):
> is immediate from the definitions. This completes the proof of assertion (ii).

★その「definitions から直ちに」の中身は、**次数の 2 つの部分の符号**である:

    degFinPre L̄ s = log #(Γ_pre(L)/Γ(X,⊤)·s) ≥ 0     （★つねに非負）
    archDeg  L̄ s = −(Σ_σ log |s|_σ)/[F:ℚ] ≥ −log B   （★|s| ≤ B なら）

★★`B` は **`X^arc` がコンパクトだから**取れる:

原文 (GenEll p.7):
> to a section of L ⊗M over X such that |t|L⊗M ≤1 on Xarc [where we recall that

★★★角括弧の中の「`X^arc` はコンパクトだった！」がこの段の唯一の根拠である。

## ★★★★★★★★★本ファイルが可能になった理由

`§9-802` で **計量の連続性**（`AMetric.IsContinuous`）と
**`|s|_L̄ : X^arc → ℝ` が連続であること**（`AMetric.continuous_norm`）が入った。
★それが無いと「コンパクト上の連続関数は有界」が使えず、`B` が取れなかった。

★★`Found/GenEll/Prop14Proper.lean` は同じ結論を**因子表示**（`ArithCartier`）で
取っていた。本ファイルは**算術直線束**（`AInv` ＝ `Definition 1.1` の設計）の側で取る。

## ★★★★★★★機構（4 段、すべて本セッションの在庫）

| 段 | 内容 |
|---|---|
| `compactSpace_arc` | 固有性 ⟹ `X^arc` コンパクト（`UltraCompact.lean`） |
| `AMetric.continuous_norm` | `|s|` は連続（`§9-802`） |
| `norm_pullSec` | `|f^*s|(q) = |s|(q ≫ f)`（`§9-784`）——★定数が `F` に依らない理由 |
| `degFinPre_nonneg` | 有限部分は非負（`log` の値域） |

★★★定数 `C = log B` は **`F` にも `x_F` にも依らない**——
`norm_pullSec` が引き戻したノルムを元のノルムに戻すからである。

## ★原文との差（記録）

原文の仮定は「`L_ℚ` のある正冪が大域切断で生成される」であり、
そこから**捻り** `M̄`（`Spec ℤ` から来る算術直線束）を掛けて
`|t| ≤ 1` の切断を作る。★本ファイルの仮定は「大域切断 `s` が 1 つある」であって、
捻りの段は**含まない**（`Found/GenEll/VerticalTwist.lean` がその段を因子表示で持つ）。

★★結論の範囲も原文どおり **`s` が消えない点の上**である
——原文は「`F ⊆ X(ℚ̄)` を `s` が非零な点の集合とする」と書く。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace NumberField
open ABC3.Found.GenEll

/-! ## ★★★(1) 2 つの部分の符号 -/

/-- ★★**有限部分はつねに非負である**。

★`degFinPre = log #(…)` であり、`Nat.card` の対数は非負だからである
（無限のときは `Nat.card = 0`、`Real.log 0 = 0` で、やはり非負）。 -/
theorem degFinPre_nonneg {X : Scheme.{0}} (L : AInv X)
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) : 0 ≤ degFinPre L s :=
  Real.log_natCast_nonneg _

theorem log_le_log_of_le_of_one_le {x B : ℝ} (hx : 0 ≤ x) (hB : 1 ≤ B) (h : x ≤ B) :
    Real.log x ≤ Real.log B := by
  rcases eq_or_lt_of_le hx with h0 | h0
  · rw [← h0, Real.log_zero]
    exact Real.log_nonneg hB
  · exact Real.log_le_log h0 h

/-- ★★★**ノルムが `B` 以下ならアルキメデス部分は `−log B` 以上**。 -/
theorem archDeg_ge_of_norm_le (F : Type) [Field F] [NumberField F]
    (L : AMetric (Spec (CommRingCat.of (𝓞 F)))) (s : L.sheaf.obj (op ⊤))
    (B : ℝ) (hB : 1 ≤ B) (h : ∀ σ : F →+* ℂ, L.norm s (embSpecPoint F σ) ≤ B) :
    -(Real.log B) ≤ archDeg F L s := by
  have hn : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  have hsum : (∑ σ : F →+* ℂ, Real.log (L.norm s (embSpecPoint F σ)))
      ≤ (Module.finrank ℚ F : ℝ) * Real.log B := by
    have h1 : (∑ σ : F →+* ℂ, Real.log (L.norm s (embSpecPoint F σ)))
        ≤ (∑ _σ : F →+* ℂ, Real.log B) :=
      Finset.sum_le_sum (fun σ _ =>
        log_le_log_of_le_of_one_le (AMetric.norm_nonneg L s _) hB (h σ))
    have h2 : (∑ _σ : F →+* ℂ, Real.log B)
        = (Fintype.card (F →+* ℂ) : ℝ) * Real.log B := by
      rw [Finset.sum_const, nsmul_eq_mul, Finset.card_univ]
    rw [h2] at h1
    rwa [NumberField.Embeddings.card F ℂ] at h1
  show -(Real.log B) ≤ -(∑ σ : F →+* ℂ,
    Real.log (L.norm s (embSpecPoint F σ))) / (Module.finrank ℚ F : ℝ)
  rw [le_div_iff₀ hn]
  nlinarith [hsum]

/-- ★★アルキメデス部分は、ノルムが `1` 以下なら非負である。 -/
theorem archDeg_nonneg (F : Type) [Field F] [NumberField F]
    (L : AMetric (Spec (CommRingCat.of (𝓞 F)))) (s : L.sheaf.obj (op ⊤))
    (h : ∀ σ : F →+* ℂ, L.norm s (embSpecPoint F σ) ≤ 1) :
    0 ≤ archDeg F L s := by
  have hF : (0 : ℝ) ≤ (Module.finrank ℚ F : ℝ) := by positivity
  show 0 ≤ -(∑ σ : F →+* ℂ, Real.log (L.norm s (embSpecPoint F σ))) / (Module.finrank ℚ F : ℝ)
  refine div_nonneg ?_ hF
  rw [neg_nonneg]
  exact Finset.sum_nonpos (fun σ _ => Real.log_nonpos (AMetric.norm_nonneg L s _) (h σ))

/-! ## ★★★★★(2) `deg_F` の下界 -/

/-- ★★★★★**`deg_F(L̄) ≥ −log B`**（ノルムが `B` 以下の切断があるとき）。 -/
theorem degArithPre_ge_of_norm_le (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) (B : ℝ) (hB : 1 ≤ B)
    (h : ∀ σ : F →+* ℂ, L.carrier.norm s (embSpecPoint F σ) ≤ B) :
    -(Real.log B) ≤ degArithPre F L s := by
  have hF : (0 : ℝ) ≤ (Module.finrank ℚ F : ℝ) := by positivity
  have h1 : 0 ≤ degFinPre L s / (Module.finrank ℚ F : ℝ) :=
    div_nonneg (degFinPre_nonneg L s) hF
  have h2 := archDeg_ge_of_norm_le F L.carrier s B hB h
  show -(Real.log B) ≤ degFinPre L s / (Module.finrank ℚ F : ℝ) + archDeg F L.carrier s
  linarith

/-- ★★★★★**`Γ(L̄)` の切断があれば `deg_F(L̄) ≥ 0`**。 -/
theorem degArithPre_nonneg (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type))
    (h : ∀ σ : F →+* ℂ, L.carrier.norm s (embSpecPoint F σ) ≤ 1) :
    0 ≤ degArithPre F L s := by
  have hF : (0 : ℝ) ≤ (Module.finrank ℚ F : ℝ) := by positivity
  have h1 : 0 ≤ degFinPre L s / (Module.finrank ℚ F : ℝ) :=
    div_nonneg (degFinPre_nonneg L s) hF
  show 0 ≤ degFinPre L s / (Module.finrank ℚ F : ℝ) + archDeg F L.carrier s
  linarith [archDeg_nonneg F L.carrier s h]

theorem degAInv_nonneg_of_mem_gamma (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) (hs : s ≠ 0)
    (hg : s ∈ AMetric.gamma L.carrier) : 0 ≤ degAInv F L := by
  rw [degAInv_eq F L s hs]
  exact degArithPre_nonneg F L s (fun σ => hg (embSpecPoint F σ))

/-! ## ★★★★★★(3) `X^arc` のコンパクト性から一様な `B` -/

/-- ★★★★★★**コンパクトな `X^arc` の上でノルムは一様に上から抑えられる**。

原文 (GenEll p.7):
> to a section of L ⊗M over X such that |t|L⊗M ≤1 on Xarc [where we recall that

★`AMetric.continuous_norm`（`§9-802`）と `compactSpace_arc`（固有性）の合成。 -/
theorem exists_norm_le_of_proper {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (Spec (CommRingCat.of ℂ) ⟶ X)]
    (L : AMetric X) (hL : L.IsContinuous) (s : L.sheaf.obj (op ⊤)) :
    ∃ B : ℝ, 1 ≤ B ∧ ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, L.norm s p ≤ B := by
  letI := arcTopology X
  haveI := compactSpace_arc hval
  obtain ⟨C, hC, hlo⟩ := exists_lower_bound_of_continuous
    (fun p => -(L.norm s p)) (AMetric.continuous_norm hL s).neg
  refine ⟨max 1 C, le_max_left _ _, fun p => ?_⟩
  have h := hlo p
  simp only [neg_le_neg_iff] at h
  exact le_trans h (le_max_right _ _)

/-! ## ★★★★★★★★★★(4) 高さの下界 -/

theorem htMetricU_eq_degAInv {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (M : AInv X) (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetricU F M xF = degAInv F (AInv.pullback xF M) := rfl

theorem norm_pullSecT {X Y : Scheme.{0}} (f : X ⟶ Y) (L : AMetric Y)
    (s : (L.sheaf.obj (op ⊤) : Type)) (q : Spec (CommRingCat.of ℂ) ⟶ X) :
    (AMetricPullback f L).norm (pullSecT f L.sheaf s) q = L.norm s (q ≫ f) :=
  norm_pullSec f L s q

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★**[GenEll] Proposition 1.4, (ii)**（算術直線束の側）——
`X` が固有なら高さは**下に有界**である。

原文 (GenEll p.6):
> by global sections [for instance, if the line bundle LQ is ample], then

    `∃ C ≥ 0, ∀ F, ∀ x_F,  x_F^*s ≠ 0 → −C ≤ ht_M̄(x_F)`

★★定数 `C = log B` は **`F` にも `x_F` にも依らない**
——`norm_pullSec`（`§9-784`）が引き戻したノルムを元のノルムへ戻すからである。
★★★結論の範囲は原文どおり「`s` が消えない点の上」である。 -/
theorem exists_htMetricU_ge {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (Spec (CommRingCat.of ℂ) ⟶ X)]
    (M : AInv X) (hMc : M.carrier.IsContinuous)
    (s : (M.carrier.sheaf.obj (op ⊤) : Type)) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
      (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X),
      pullSecT xF M.carrier.sheaf s ≠ 0 → -C ≤ htMetricU F M xF := by
  obtain ⟨B, hB, hb⟩ := exists_norm_le_of_proper hval M.carrier hMc s
  refine ⟨Real.log B, Real.log_nonneg hB, fun F _ _ xF hne => ?_⟩
  rw [htMetricU_eq_degAInv, degAInv_eq F _ (pullSecT xF M.carrier.sheaf s) hne]
  refine degArithPre_ge_of_norm_le F (AInv.pullback xF M) _ B hB (fun σ => ?_)
  show (AMetricPullback xF M.carrier).norm (pullSecT xF M.carrier.sheaf s)
    (embSpecPoint F σ) ≤ B
  rw [norm_pullSecT]
  exact hb _

/-- ★★★★★★★★**切断が `Γ(M̄)` に入るなら定数は要らない**——`ht_M̄ ≥ 0`。

原文 (GenEll p.7):
> is immediate from the definitions. This completes the proof of assertion (ii).

★これが原文の「definitions から直ちに」そのものである。 -/
theorem htMetricU_nonneg_of_gamma {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (M : AInv X) (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X)
    (s : (M.carrier.sheaf.obj (op ⊤) : Type))
    (hne : pullSecT xF M.carrier.sheaf s ≠ 0)
    (hg : ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, M.carrier.norm s p ≤ 1) :
    0 ≤ htMetricU F M xF := by
  rw [htMetricU_eq_degAInv]
  refine degAInv_nonneg_of_mem_gamma F _ (pullSecT xF M.carrier.sheaf s) hne (fun q => ?_)
  show (AMetricPullback xF M.carrier).norm (pullSecT xF M.carrier.sheaf s) q ≤ 1
  rw [norm_pullSecT]
  exact hg _

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★`X(ℚ̄)` の点の水準での下界。 -/
theorem exists_htMetricAlg_ge {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (Spec (CommRingCat.of ℂ) ⟶ X)]
    (M : AInv X) (hMc : M.carrier.IsContinuous)
    (s : (M.carrier.sheaf.obj (op ⊤) : Type)) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
      (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X),
      pullSecT xF M.carrier.sheaf s ≠ 0 →
      -C ≤ htMetricAlg M (AlgPointAnyClass.mk (algPointAny F xF)) := by
  obtain ⟨C, hC, h⟩ := exists_htMetricU_ge hval M hMc s
  refine ⟨C, hC, fun F _ _ xF hne => ?_⟩
  rw [htMetricAlg_mk]
  exact h F xF hne

/-! ## ★出典の紐付け(`.src`) -/

def degFinPre_nonneg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (ii)(次数の有限部分はつねに非負)",
    sectionId := "genell-prop-1-4" }

def exists_norm_le_of_proper.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (ii)(X^arc のコンパクト性からノルムの一様上界——計量の連続性を使う)",
    sectionId := "genell-prop-1-4" }

def exists_htMetricU_ge.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(算術直線束の側——高さは下に有界。捻りの段は含まない)",
    sectionId := "genell-prop-1-4" }

def htMetricU_nonneg_of_gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (ii)(Γ(M̄) の切断があれば ht ≥ 0——原文の「definitions から直ちに」)",
    sectionId := "genell-prop-1-4" }

def exists_htMetricU_ge.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "compactSpace_arc(固有性 ⟹ X^arc コンパクト)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.compactSpace_arc") 7,
    .citation "[ABC3]" "AMetric.continuous_norm(|s|_L̄ は X^arc 上連続、§9-802)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.continuous_norm") 3,
    .citation "[ABC3]" "norm_pullSec(|f^*s|(q) = |s|(q ≫ f)、§9-784——定数が F に依らない理由)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.norm_pullSec") 3,
    .citation "[ABC3]" "degAInv_eq(deg_F はどの非零切断で測っても同じ、§9-782)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degAInv_eq") 4,
    .implicitStep
      ("★原文の仮定は「L_ℚ のある正冪が大域切断で生成される」であり、そこから " ++
       "Spec ℤ から来る算術直線束 M̄ を掛けて |t| ≤ 1 の切断を作る。" ++
       "★★本ファイルはその**捻りの段を含まない**——因子表示では " ++
       "Found/GenEll/VerticalTwist.lean が同じ段を持つ") 7,
    .implicitStep
      ("★★(iii)(BD-class が L_ℚ の同型類だけに依る)と (iv)(Northcott)は本ファイルに無い。" ++
       "(iv) は ample ⟹ very ample(Serre)を要し、ResearchPaper/mathlib-gap.json の " ++
       "ample-and-projective-embedding の段 E が未着手である") 6 ]

end ABC3.Found.Arakelov
