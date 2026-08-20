import ABC3.Skeleton.GaloisRep.WeilFunctionField
import ABC3.Found.GaloisRep.TorsionIdealIntegral

/-!
# スケルトン —— **`n` 乗根 `g_P` の取り出し**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★層 3 —— `f_P` と `[n]^*` の合流点

    g_P^n = f_P ∘ [n]                    (`exists_nthRoot_comp_mulByN`)

★両側の材料は **`Found` に揃った**:

| 材料 | 場所 |
|---|---|
| `f_P ∈ F[W]`(`(XYIdeal)^n = (f_P)`) | `Found/GaloisRep/TorsionIdealIntegral.lean`(第 126) |
| `μ : F[W] →+* F(W)`(`[n]` の引き戻し、群法則で固定) | `Found/GaloisRep/GenericNotTorsion.lean`(第 125) |

★★残るのは**因子の計算**である——`div(f_P ∘ [n])` が `n` で割れること、
およびそこから `n` 乗根を取り出すこと。

## ★★★★★ここが Weil 対の**心臓**である

`e_n(P,Q) := g_P(·+Q)/g_P(·)` が定数になるのは、
`g_P(·+Q)^n = f_P([n]·+[n]Q) = f_P([n]·)` (`[n]Q = 0` だから) `= g_P(·)^n`
——すなわち比の `n` 乗が 1 だからである。★**`g_P` の存在なしにはこの議論は立たない。**

## ★★★★★★mathlib の在庫(2026-08-20 再実測)

mathlib は Dedekind 環の**因子機構を完備している**——
`IsDedekindDomain.HeightOneSpectrum`・`count`・`finprod_heightOneSpectrum_factorization`・
adic 付値。★**「因子群が無い」という当初の記録は誤りであった。**

★★欠けているのは **`IsDedekindDomain W.CoordinateRing` の instance だけ**である:

| 条件 | mathlib |
|---|---|
| `IsDomain W.CoordinateRing` | ✅ mathlib |
| `IsNoetherianRing W.CoordinateRing` | ✅ mathlib |
| `Ring.DimensionLEOne W.CoordinateRing` | ✅ **第 127 ブロック** |
| `IsIntegrallyClosed W.CoordinateRing` | ❌ **これだけ** |

★★★これが入れば `div(f)` は**主分数イデアルの素分解**としてそのまま出る。
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F]

/-- ★★★★★**`g_P` の存在**——`f_P ∘ [n]` は関数体の中で `n` 乗根を持つ。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`f_P` は第 126 ブロックで座標環の元として得られている。
★★`μ` は第 118・125 ブロックで**群法則によって固定された形で**得られている。 -/
theorem exists_nthRoot_comp_mulByN (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n) (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    ∃ g : W.FunctionField, g ^ n = μ fP := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_nthRoot_comp_mulByN.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——g_P の n 乗根の取り出し)",
    sectionId := "genell-thm-3-8" }

def exists_nthRoot_comp_mulByN.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1 の証明(g_P の存在)"
      (.absent "mathlib に Weil 対およびその構成要素は 0 件(2026-08-20、`WeilPairing|weil_pairing` で全文検索して 0 件)") 19,
    .implicitStep
      "★★★★★★両側の材料は **Found に揃った**——`f_P ∈ F[W]`(第 126 `xyIdeal_pow_isPrincipal_integral`)と `μ : F[W] →+* F(W)`(第 118・125 `exists_mulByNHom_charZero`、群法則で固定)(0 ブロック)" 19,
    .implicitStep
      "★★★`div(f_P ∘ [n])` が `n` で割れることの因子計算。`[n]` が次数 `n²` であること、および `[n]^{-1}(P)` の各点の重複度を使う(20-50 ブロック)" 19,
    .implicitStep
      "★★★★★★**2026-08-20 の再実測で見積もりが変わった**: mathlib は Dedekind 環の**因子機構を完備している**——`IsDedekindDomain.HeightOneSpectrum`・`count`・`finprod_heightOneSpectrum_factorization`・adic 付値。「因子群が無い」は誤りであり、**欠けているのは `IsDedekindDomain W.CoordinateRing` の instance だけ**である" 19,
    .implicitStep
      "★★★**座標環が Dedekind 環であること**を示す——`IsDomain` ✅・`IsNoetherianRing` ✅(mathlib)、`Ring.DimensionLEOne` ✅(第 127 ブロック、第 116 の整拡大から)。**残るのは `IsIntegrallyClosed` だけ**である(10-25 ブロック)" 19,
    .implicitStep
      "★★`IsIntegrallyClosed` はアフィン曲線の**正則性**である。mathlib の局所判定 `IsIntegrallyClosed.of_localization_maximal` があるので、各極大イデアルでの局所化が DVR であることを示す形になる——非特異点なら `x − x₀` か `y − y₀` が極大イデアルを生成する(10-25 ブロック)" 19,
    .implicitStep
      "★★Dedekind の instance が入れば、`div(f)` は【主分数イデアルの素分解】として mathlib から直接出る。因子群を自分で積む必要は無い(0 ブロック)" 19,
    .implicitStep
      "★因子が `n` で割れる ⟹ `n` 乗根が取れる、の段。Dedekind なら分数イデアル `J` で `J^n = (μ f_P)` なるものが取れ、`J` が単項であることを `toClass` の単射性で示す形になる(10-25 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
