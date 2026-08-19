import ABC3.Found.Arakelov.ArcTopologyAffine

/-!
# Arakelov (D1) 第 297 ブロック —— **アフィン間の射は弧空間で連続**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★(D1) が要る新しい道具——**弧空間の関手性**

(D1) `APicData` の `pullback` は、可逆層だけでなく**計量も**引き戻さねばならない。
★計量は Green 関数(連続)なので、引き戻しが連続であるためには

    (· ≫ f) : X^arc ⟶ Y^arc  が連続

が要る。★★これは今まで**開埋め込みの場合しか無かった**(第 275)。

## ★★アフィンなら座標計算だけで出る

`Spec` は充満忠実なので `h : Spec A ⟶ Spec B` は環準同型 `φ = Spec.preimage h : B ⟶ A` から来る。
★したがって

    evalAffine B (p ≫ h) b = evalAffine A p (φ b)

——**像の座標 `b` は、もとの座標 `φ b` そのもの**である。
★★積位相への写像は各座標で連続を見ればよいので、これで終わる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `evalAffine_comp` | ★座標の対応 |
| `continuous_comp_affine` | ★★★★**アフィン間の射は弧空間で連続** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

/-- ★**像の座標 `b` は、もとの座標 `φ b`**。 -/
theorem evalAffine_comp {A B : CommRingCat.{0}} (h : Spec A ⟶ Spec B)
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) (b : B) :
    evalAffine B (p ≫ h) b = evalAffine A p ((Spec.preimage h).hom b) := by
  show (Spec.preimage (p ≫ h)).hom b = (Spec.preimage p).hom ((Spec.preimage h).hom b)
  have hc : Spec.preimage (p ≫ h) = Spec.preimage h ≫ Spec.preimage p := by
    apply Spec.map_injective
    rw [Spec.map_comp, Spec.map_preimage, Spec.map_preimage, Spec.map_preimage]
  rw [hc]
  rfl

/-- ★★★★**アフィン間の射は弧空間で連続である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★弧位相は各点収束なので、**座標ごとに見れば射影**である。 -/
theorem continuous_comp_affine {A B : CommRingCat.{0}} (h : Spec A ⟶ Spec B) :
    @Continuous _ _ (arcTopologyAffine A) (arcTopologyAffine B)
      (fun p : Spec (CommRingCat.of ℂ) ⟶ Spec A => p ≫ h) := by
  letI := arcTopologyAffine A
  refine continuous_induced_rng.2 (continuous_pi (fun b => ?_))
  have he : (fun p : Spec (CommRingCat.of ℂ) ⟶ Spec A => evalAffine B (p ≫ h) b)
      = fun p => evalAffine A p ((Spec.preimage h).hom b) :=
    funext (fun p => evalAffine_comp h p b)
  show Continuous (fun p : Spec (CommRingCat.of ℂ) ⟶ Spec A => evalAffine B (p ≫ h) b)
  rw [he]
  exact (continuous_apply _).comp continuous_induced_dom

/-! ## ★出典の紐付け(`.src`) -/

def continuous_comp_affine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィン間の射が弧空間で連続であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
