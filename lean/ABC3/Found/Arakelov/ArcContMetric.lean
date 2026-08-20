import ABC3.Found.Arakelov.ArcSectionCont

/-!
# Arakelov (C3) 第 283 ブロック —— **★★★★★★連続な弧計量が構成できた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★4 法則がすべて揃った

第 246 ブロックの `ArcMetric`(3 法則)に**連続性**を足した `ContArcMetric` を、
**自明化被覆と 1 の分割から実際に構成する**。

| 法則 | 出所 |
|---|---|
| 非負 | 第 272 `gluedNormX_nonneg` |
| スカラー倍 | 第 272 `gluedNormX_smul` |
| `= 0 ↔ v = 0` | 第 272 `gluedNormX_eq_zero_iff` |
| ★★★**連続性** | 第 282 `continuous_gluedNorm_section` |

## ★★1 の分割から出る 2 つの側条件

`gluedNormX_eq_zero_iff` は `i₀`(その点で `ρ > 0`)と**台の有限性**を要求する。
★どちらも mathlib の `PartitionOfUnity` から出る:

| 要求 | mathlib |
|---|---|
| `∃ i₀, 0 < f i₀ p` | ★`PartitionOfUnity.exists_pos` |
| `p ⁻¹ᵁ (U i₀) = ⊤` | ★★`IsSubordinate` + `subset_tsupport` |
| 台の有限性 | ★`LocallyFinite.point_finite` |

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `ContArcMetric` | ★★連続な弧計量 |
| `pou_exists_top` | ★1 の分割から `i₀` と `p ⁻¹ᵁ U i₀ = ⊤` |
| `pou_finite` | ★台の有限性 |
| `gluedArcMetric` | ★★★★★★**構成された連続な弧計量** |

★これで `metric_nonempty` に要るのは**自明化被覆の存在**と**1 の分割の存在**だけになった。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

open scoped Classical

variable {X : Scheme.{0}} {ι : Type}

/-- ★★**連続な弧計量**——第 246 の 3 法則に連続性を足したもの。 -/
structure ContArcMetric (X : Scheme.{0}) (F : X.Modules) extends ArcMetric X F where
  /-- ★★★切断のノルムは `X^arc` 上連続である。 -/
  cont : ∀ (s : (F.val.obj (op ⊤) : Type)),
    @Continuous _ ℝ (arcTopology X) _ (fun p => nrm p (arcEval p F s))

variable (U : ι → X.Opens)
  (f : @PartitionOfUnity ι (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X) Set.univ)
  (hsub : @PartitionOfUnity.IsSubordinate ι (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X)
    Set.univ f (fun i => arcOpenSet (U i)))

/-- ★1 の分割の重みは非負である。 -/
theorem pou_nonneg (i : ι) (p : Spec (CommRingCat.of ℂ) ⟶ X) : 0 ≤ f i p := by
  letI := arcTopology X
  exact f.nonneg i p

include hsub in
/-- ★**1 の分割は各点で「正の重みを持ち、かつその点を含む」添字を与える**。 -/
theorem pou_exists_top (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    ∃ i₀, 0 < f i₀ p ∧ p ⁻¹ᵁ (U i₀) = ⊤ := by
  letI := arcTopology X
  obtain ⟨i₀, hpos⟩ := f.exists_pos (Set.mem_univ p)
  refine ⟨i₀, hpos, ?_⟩
  exact hsub i₀ (subset_tsupport _ (ne_of_gt hpos))

/-- ★**1 の分割の台は各点で有限である**。 -/
theorem pou_finite (F : X.Modules)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) :
    (Function.support (fun i => f i p * extNormX F U e i p w)).Finite := by
  letI := arcTopology X
  refine Set.Finite.subset (f.locallyFinite.point_finite p) ?_
  intro i hi
  exact fun hz => hi (by simp only [hz, zero_mul])

/-- ★★★★★★**構成された連続な弧計量**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★自明化被覆と、それに従属する 1 の分割があれば、
**4 法則をすべて満たす計量が実際に作れる**。 -/
noncomputable def gluedArcMetric (F : X.Modules)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (ho : ∀ i, @IsOpen _ (arcTopology X) (arcOpenSet (U i))) : ContArcMetric X F where
  nrm := fun p w => gluedNormX F U e (fun i p => f i p) p w
  nonneg := fun p w =>
    gluedNormX_nonneg F U e (fun i p => f i p) (pou_nonneg f) p w
  smul := fun p c w => gluedNormX_smul F U e (fun i p => f i p) p c w
  eq_zero_iff := fun p w => by
    obtain ⟨i₀, hpos, htop⟩ := pou_exists_top U f hsub p
    exact gluedNormX_eq_zero_iff F U e (fun i p => f i p) (pou_nonneg f) p
      i₀ hpos htop w (pou_finite U f F e p w)
  cont := fun s => continuous_gluedNorm_section F U e f ho hsub s

/-! ## ★出典の紐付け(`.src`) -/

def gluedArcMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——連続な弧計量の構成)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
