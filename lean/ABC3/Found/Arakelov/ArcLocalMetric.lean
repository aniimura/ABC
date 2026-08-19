import ABC3.Found.Arakelov.ArcEvalCont

/-!
# Arakelov (C3) 第 261 ブロック —— ★★★★★**局所ノルムは連続である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★局所の段が**完成**した

`V` 上の**消えない切断** `g` があれば、正規化ノルム

    genNorm(w) = nrm(w) / nrm(g(p))

は 4 法則(非負・`0 ↔ 0`・`‖c·v‖ = |c|‖v‖`・正規化)を満たし、
★★**切断のノルムが連続**である:

    genNorm (arcEvalOnTop (c·g)) = ‖evalOn c‖        (第 258)
    evalOn は連続                                     (第 260)

★★★これで (C3) の**局所の段はすべて揃った**。

## ★残る 2 つ

| 段 | 内容 |
|---|---|
| 生成切断の非消滅 | ★自明化 `e` から `g` を取り、`arcEvalOnTop g ≠ 0` を示す |
| 貼り合わせ | ★★1 の分割(コンパクト Hausdorff の仮定つき) |

★★非消滅は本ブロックでは**仮定**として受けている——
`e` から `g := e.inv(1)` を取る段で、`p^*` との整合を示す必要がある。

| 定理 | 内容 |
|---|---|
| `evalOn_one` | ★`1` の値は `1` |
| `continuous_genNorm` | ★★★★★**局所ノルムの連続性** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}} (F : X.Modules) (hF : IsLocallyTrivial X F.val) (V : X.Opens)

/-- ★`1` の値は `1`。 -/
theorem evalOn_one (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤) :
    evalOn p V h 1 = 1 := by
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((p.appLE V ⊤ (le_of_eq h.symm)).hom 1) = 1
  rw [map_one, map_one]

/-- ★★★★★生成切断から作った局所ノルムは連続である。 -/
theorem continuous_genNorm (g : (F.val.obj (op V) : Type))
    (hg : ∀ (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme),
      arcEvalOnTop F (q ≫ V.ι) V (comp_preimage_eq_top V q) g ≠ 0)
    (c : ((X.presheaf.obj (op V)) : Type)) :
    @Continuous _ ℝ (arcTopology V.toScheme) _
      (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme =>
        genNorm F hF V g (q ≫ V.ι) (comp_preimage_eq_top V q)
          (arcEvalOnTop F (q ≫ V.ι) V (comp_preimage_eq_top V q) (c • g))) := by
  letI := arcTopology V.toScheme
  have he : (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme =>
        genNorm F hF V g (q ≫ V.ι) (comp_preimage_eq_top V q)
          (arcEvalOnTop F (q ≫ V.ι) V (comp_preimage_eq_top V q) (c • g)))
      = (fun q => ‖evalOn (q ≫ V.ι) V (comp_preimage_eq_top V q) c‖) :=
    funext (fun q => genNorm_arcEvalOnTop F (q ≫ V.ι) hF V g (comp_preimage_eq_top V q) (hg q) c)
  rw [he]
  exact continuous_norm.comp (continuous_evalOn V c)


/-! ## ★出典の紐付け(`.src`) -/

def continuous_genNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——局所ノルムは連続)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
