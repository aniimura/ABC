import ABC3.Found.Arakelov.ArcFactorEval
import ABC3.Found.Arakelov.ArcContTransport

/-!
# Arakelov (C3) 第 282 ブロック —— **★★★★★★★切断のノルムは連続である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★仮定 `hcompat` が消えた

第 281 ブロックまでは、貼り合わせたノルムの連続性は

    hcompat : (arcFiberFactor).hom (arcEval (q ≫ ι) F s) = arcEval q (restrict F ι) (…)

を**仮定**したうえでしか出せなかった。★★第 281 でこれが定理になったので、
**仮定なしの連続性**が出る。

## ★★機構

| 段 | 使うもの |
|---|---|
| 1 | 第 277 `xNorm_comp` —— `X` 側のノルムを `V` 側の `localNorm` に落とす |
| 2 | 第 281 `factor_arcEval` —— 評価が制限と両立する |
| 3 | 第 270 `continuous_localNorm` —— `V` 側では連続 |
| 4 | 第 276 `continuousOn_extNormX` —— `X^arc` の開集合上へ移す |
| 5 | 第 273 `continuous_gluedNormX` —— 1 の分割で貼る |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `continuous_xNorm_section` | ★★★★局所での連続性(仮定なし) |
| `continuousOn_extNormX_section` | ★★開集合上での連続性 |
| `continuous_gluedNorm_section` | ★★★★★★★**貼り合わせたノルムの連続性(仮定なし)** |

★これで Interface の `normSection_continuous` が満たせる。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens) (F : X.Modules)
  (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★★★★**局所での連続性**——`hcompat` の仮定なしで出る。 -/
theorem continuous_xNorm_section (s : (F.val.obj (op ⊤) : Type)) :
    @Continuous _ ℝ (arcTopology V.toScheme) _
      (fun q => xNorm V F e (q ≫ V.ι) (comp_preimage_eq_top V q) (arcEval (q ≫ V.ι) F s)) := by
  have h : (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme =>
        xNorm V F e (q ≫ V.ι) (comp_preimage_eq_top V q) (arcEval (q ≫ V.ι) F s))
      = fun q => localNorm V F e q
        (arcEval q (Scheme.Modules.restrict F V.ι) (restrictSection V.ι F s)) := by
    funext q
    rw [xNorm_comp, factor_arcEval]
  rw [h]
  exact continuous_localNorm V F e (restrictSection V.ι F s)

variable {ι : Type}

open scoped Classical

/-- ★★**開集合 `arcOpenSet (U i)` の上での連続性**。 -/
theorem continuousOn_extNormX_section (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i))) (i : ι) (s : (F.val.obj (op ⊤) : Type)) :
    @ContinuousOn _ ℝ (arcTopology X) _
      (fun p => extNormX F U e i p (arcEval p F s)) (arcOpenSet (U i)) :=
  continuousOn_extNormX F U e i (fun p => arcEval p F s)
    (continuous_xNorm_section (U i) F (e i) s)

/-- ★★★★★★★**貼り合わせたノルムの連続性(仮定なし)**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが Interface の `normSection_continuous` そのものである。 -/
theorem continuous_gluedNorm_section (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (f : @PartitionOfUnity ι (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X) Set.univ)
    (ho : ∀ i, @IsOpen _ (arcTopology X) (arcOpenSet (U i)))
    (hsub : @PartitionOfUnity.IsSubordinate ι (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X)
      Set.univ f (fun i => arcOpenSet (U i)))
    (s : (F.val.obj (op ⊤) : Type)) :
    @Continuous _ ℝ (arcTopology X) _
      (fun p => gluedNormX F U e (fun i p => f i p) p (arcEval p F s)) :=
  continuous_gluedNormX F U e f ho hsub (fun p => arcEval p F s)
    (fun i => continuousOn_extNormX_section F U e i s)

/-! ## ★出典の紐付け(`.src`) -/

def continuous_gluedNorm_section.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——切断のノルムが連続であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
