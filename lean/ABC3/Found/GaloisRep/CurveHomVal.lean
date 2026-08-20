import ABC3.Found.GaloisRep.ClassSum

/-!
# Galois (G5) 第 167 ブロック —— **★★★★★★座標が整なら像はすべて整**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★`hnn` は仮定でなく**定理**だった

第 143 以降、場合 A の仮定として `hnn : ∀ r, w(μ r) ≤ 1` を置いてきた。★これは実は

    w(μ genX) ≤ 1  かつ  w(μ genY) ≤ 1   ⟹   ∀ r, w(μ r) ≤ 1

という**定理**である。座標環は `AdjoinRoot` なので `mk` が全射で、
`AdjoinRoot.lift_mk` により `μ r` は 2 重の `eval₂` になる。★第 144 の
`valuation_eval₂_le_one` を 2 回当てるだけで済む(`valuation_curveHom_le_one`)。

★★これで「点の座標が付値環に入る」ことだけ確かめれば `hnn` が自動で従う。

## ★★★★★引き戻した素イデアルは還元した点の素イデアル(一般形)

第 164 では `μ = [n]^*` の場合に `P' = XYIdeal(n·Q_v)` を示した。★本ブロックで
**任意の `curveHom`** に対して

    pullbackPrime(w_v, curveHom W heq) = XYIdeal( redConst(x), redConst(y) )

が成り立つことを示す。★★これは「点の還元」を素イデアルの言葉に翻訳したもので、
平行移動 `τ_T` による素点の輸送(`e_v` の一定性)でそのまま使える。

## ★★★実測: 環準同型の一意性は既にあった

「座標環からの環準同型は `genX`・`genY`・定数で決まる」は**第 119 ブロックに既にある**
(`Found/GaloisRep/FieldHomExt.lean` の `coordinateRing_hom_ext`)。★これで
「`τ ∘ μ = μ`(`n·T = 0` のとき)」が追加工数なしで言える。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `valuation_curveHom_le_one` | ★★★★★★**座標が整なら像はすべて整**(`hnn` が自動) |
| `pullbackPrime_curveHom` | ★★★★★★★**引き戻した素点 = 還元した点**(一般形) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)

variable [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)

/-! ## ★★★★★★`hnn` は自動である -/

/-- ★★★★★★**座標が付値環に入れば、座標環の像はすべて付値環に入る**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`AdjoinRoot.mk` が全射なので `curveHom W heq r` は 2 重の `eval₂` になる。
★★第 144 の `valuation_eval₂_le_one` を 2 回当てるだけ。 -/
theorem valuation_curveHom_le_one {x y : W.FunctionField}
    (heq : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hx : v.valuation W.FunctionField x ≤ 1) (hy : v.valuation W.FunctionField y ≤ 1)
    (r : W.CoordinateRing) :
    v.valuation W.FunctionField (curveHom W heq r) ≤ 1 := by
  obtain ⟨p, rfl⟩ := AdjoinRoot.mk_surjective r
  rw [curveHom, AdjoinRoot.lift_mk]
  refine valuation_eval₂_le_one (v.valuation W.FunctionField) _ ?_ hy p
  intro cc
  rw [Polynomial.coe_eval₂RingHom]
  exact valuation_eval₂_le_one (v.valuation W.FunctionField) _
    (fun d => valuation_algebraMap_field W v d) hx cc

variable {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

variable [W.IsElliptic]

include hv in
/-- ★★★★★★★**引き戻した素イデアルは還元した点の素イデアルである**(一般形)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 164 の `μ = [n]^*` の場合を任意の `curveHom` に一般化したもの。
★★平行移動 `τ_T` による素点の輸送でそのまま使える。 -/
theorem pullbackPrime_curveHom {x y : W.FunctionField}
    (heq : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hx : v.valuation W.FunctionField x ≤ 1) (hy : v.valuation W.FunctionField y ≤ 1) :
    pullbackPrime (v.valuation W.FunctionField) (curveHom W heq)
        (valuation_curveHom_le_one W v heq hx hy)
      = CoordinateRing.XYIdeal W (redConst W v h hv x)
        (Polynomial.C (redConst W v h hv y)) := by
  have hR : W.Nonsingular (redConst W v h hv x) (redConst W v h hv y) :=
    equation_iff_nonsingular.mp (equation_redConst W v h hv heq hx hy)
  refine ((xyIdeal_isMaximal W hR.1).eq_of_le
    (pullbackPrime_isPrime (v.valuation W.FunctionField) _ _).ne_top ?_).symm
  rw [CoordinateRing.XYIdeal, Ideal.span_le]
  rintro z (rfl | rfl)
  · show v.valuation W.FunctionField
      (curveHom W heq (CoordinateRing.XClass W (redConst W v h hv x))) < 1
    rw [XClass_eq, map_sub, curveHom_genX, curveHom_algebraMap]
    exact redConst_spec W v h hv hx
  · show v.valuation W.FunctionField
      (curveHom W heq (CoordinateRing.YClass W (Polynomial.C (redConst W v h hv y)))) < 1
    rw [YClass_eq, map_sub, curveHom_genY, curveHom_algebraMap]
    exact redConst_spec W v h hv hy

/-! ## ★出典の紐付け(`.src`) -/

def valuation_curveHom_le_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——座標が整なら座標環の像がすべて整であること)",
    sectionId := "genell-thm-3-8" }

def pullbackPrime_curveHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——引き戻した素点が還元した点であること、一般形)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
