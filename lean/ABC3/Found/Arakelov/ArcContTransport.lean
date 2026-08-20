import ABC3.Found.Arakelov.ArcEmbedding
import ABC3.Found.Arakelov.ArcLiftOpen

/-!
# Arakelov (C3) 第 276 ブロック —— ★★★★**局所連続性を `X^arc` へ移す**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★第 273 の仮定 `hg` が `V^arc` の話に落ちた

    ContinuousOn (fun p => extNormX … i p (φ p)) (arcOpenSet (U i))
      ⟸ Continuous (fun q => xNorm … (q ≫ (U i).ι) … (φ (q ≫ (U i).ι)))

★機構は 2 つ:第 274(`arcOpenSet = range`)と第 275(開埋め込みでの移送)。
★★`extNormX` の `dite` は `comp_preimage_eq_top` で **`dif_pos` に潰れる**。

## ★`liftToOpenOfTop` は `(· ≫ V.ι)` の逆

`liftToOpenOfTop V (q ≫ V.ι) h = q`——`V.ι` が mono なので消去律で出る。
★★これで `xNorm` の中の `lift` が `q` に戻り、第 270 と繋がる。

## ★摩擦 —— `𝟙_` は `open MonoidalCategory` が要る

`𝟙_ (PresheafModulesOn X (U i))` が **`unexpected token '_'`** で落ちた。
★`𝟙` と `_` に分かれて読まれていた。★★**記法が読めないときは `open` を疑う**
——本セッションの「名前が見えない」3 原因目である
(import 非推移 / 名前空間 / ★記法の `open`)。

| 定理 | 内容 |
|---|---|
| `liftToOpenOfTop_comp` | ★`lift` は合成の逆 |
| `continuousOn_extNormX` | ★★★★連続性の移送 |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Topology

variable {X : Scheme.{0}} (V : X.Opens)
/-- ★★`liftToOpenOfTop` は `(· ≫ V.ι)` の逆である。 -/
theorem liftToOpenOfTop_comp (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme)
    (h : (q ≫ V.ι) ⁻¹ᵁ V = ⊤) :
    liftToOpenOfTop V (q ≫ V.ι) h = q :=
  comp_openImmersion_injective V.ι (liftToOpen_fac V (q ≫ V.ι) h)


variable {ι : Type}

open scoped Classical

/-- ★★★★局所連続性を `X^arc` の開部分集合上の連続性へ移す。 -/
theorem continuousOn_extNormX (F : X.Modules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i))) (i : ι)
    (φ : ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, ↥(arcFiber p F))
    (h : @Continuous _ ℝ (arcTopology (U i).toScheme) _
      (fun q => xNorm (U i) F (e i) (q ≫ (U i).ι) (comp_preimage_eq_top (U i) q)
        (φ (q ≫ (U i).ι)))) :
    @ContinuousOn _ ℝ (arcTopology X) _
      (fun p => extNormX F U e i p (φ p)) (arcOpenSet (U i)) := by
  rw [arcOpenSet_eq_range]
  refine continuousOn_range_of_comp (U i) _ ?_
  have he : (fun q : Spec (CommRingCat.of ℂ) ⟶ (U i).toScheme =>
        extNormX F U e i (q ≫ (U i).ι) (φ (q ≫ (U i).ι)))
      = (fun q => xNorm (U i) F (e i) (q ≫ (U i).ι) (comp_preimage_eq_top (U i) q)
        (φ (q ≫ (U i).ι))) := by
    funext q
    simp only [extNormX, dif_pos (comp_preimage_eq_top (U i) q)]
  rw [he]
  exact h


/-! ## ★出典の紐付け(`.src`) -/

def continuousOn_extNormX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——局所連続性を X^arc へ移す)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
