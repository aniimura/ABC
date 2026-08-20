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

★★★★★★**2026-08-20: `IsDedekindDomain W.CoordinateRing` が出た**(第 137 ブロック):

| 条件 | 出所 |
|---|---|
| `IsDomain W.CoordinateRing` | ✅ mathlib |
| `IsNoetherianRing W.CoordinateRing` | ✅ mathlib |
| `Ring.DimensionLEOne W.CoordinateRing` | ✅ 第 127 ブロック |
| `IsIntegrallyClosed W.CoordinateRing` | ✅ **第 137 ブロック(Dedekind の副産物)** |
| **`IsDedekindDomain W.CoordinateRing`** | ✅ **第 137 ブロック** |

★★★したがって `div(f)` は**主分数イデアルの素分解**としてそのまま出る。
★★★★残るのは **`div(f_P ∘ [n])` が `n` で割れることの因子計算**と
**`n` 乗根の取り出し**の 2 つだけである。
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
      "★★★★★★**座標環が Dedekind 環であることは第 137 ブロックで出た**——`IsDomain` ✅・`IsNoetherianRing` ✅(mathlib)、`Ring.DimensionLEOne` ✅(第 127)、そして `IsDedekindDomainDvr`(Noether + 局所化がすべて DVR)を経由して `IsDedekindDomain` ✅。★`IsIntegrallyClosed` は副産物として付いてきた(0 ブロック)" 19,
    .implicitStep
      "★★★★★★**古典的零点定理は要らなかった**(第 134 ブロック)——極大イデアルの `x` 座標は、整拡大(第 116)+ mathlib の `Ideal.IsIntegral.comap_ne_bot` + `F[X]` が PID + `F` が代数閉、の 3 段で取れる。★さらに第 135 で `P = (x − c, y − y₀)` まで出た(`AdjoinRoot.mk` の全射性と `Polynomial.ringHom_ext` だけ。`Q` が極大であることを示す必要すら無い)(0 ブロック)" 19,
    .implicitStep
      "★★★**別解(判別式)**: `CR = F[x][z]`(`z² = Ψ₂Sq(x)`)は **Found に済**(第 129 `genZ_sq`)。`Ψ₂Sq` が squarefree なら整閉になる(0 ブロック)" 19,
    .implicitStep
      "★★★★`Ψ₂Sq` の判別式は **mathlib に既にあった**(第 130 で確認)——`twoTorsionPolynomial_discr : discr = 16Δ`、`twoTorsionPolynomial_discr_ne_zero`。★§9-438 で「無い」と書いたのは `Weierstrass.lean` を見ていなかったためである(0 ブロック)" 19,
    .implicitStep
      "★★`discr ≠ 0` ⟹ `Squarefree Ψ₂Sq` の段。mathlib の在庫(2026-08-20 実測): `Cubic.discr_ne_zero_iff_roots_nodup` ✅、`PerfectField.separable_iff_squarefree` ✅、`Polynomial.Separable.squarefree` ✅、`UniqueFactorizationMonoid.squarefree_iff_nodup_normalizedFactors` ✅。★★★★★★**これは第 131 ブロックで Found に入った**(`squarefree_of_roots_nodup`・`squarefree_Psi2Sq`)。分解 `p = C a·∏(X−r)` と `separable_prod_X_sub_C_iff'` から 1 ブロックで出た(0 ブロック)" 19,
    .implicitStep
      "★★★★★**二次拡大の整閉性を直接示す必要は無かった**——`IsDedekindDomainDvr` 経由で局所化が DVR であることを示す方が短く、整閉性はその副産物として出る(第 136・137)。mathlib に二次拡大の整閉性が無いという実測自体は正しいが、迂回できた(0 ブロック)" 19,
    .implicitStep
      "★★★★★**より短い経路が見つかった**: mathlib の `IsDedekindDomainDvr.isDedekindDomain` は **instance** である(2026-08-20 実測)。`IsDedekindDomainDvr A` = `IsNoetherian A A` + 「零でない素イデアルでの局所化がすべて DVR」なので、**整閉性を経由せず**に Dedekind が出る。★`IsDedekindDomainDvr.isIntegrallyClosed` も instance なので整閉性は副産物として得られる" 19,
    .implicitStep
      "★★★★★★**局所化が DVR であることは第 136 ブロックで出た**——`Ψ₂Sq(c) = 0`(分岐)なら `z`、`≠ 0`(不分岐)なら `x − c` が一意化元。どちらも第 132・133 の分解の片側が局所化で単元になることによる。mathlib の `IsDiscreteValuationRing.TFAE` の項 0 ⟺ 項 4(`DVR ⟺ 極大イデアルが単項`)に流した(0 ブロック)" 19,
    .implicitStep
      "★★★★実測で効いた技法(第 136): `Localization.AtPrime P` の**具体型**で補題を述べると、文の詰め込みだけで 200000 heartbeats を超えた(8.3 秒でタイムアウト)。★`(S : Type) [IsLocalization P.primeCompl S]` と抽象化したら 0.24 秒になった。具体型のインスタンス探索が重い場合、抽象化は単なる一般化以上の効果を持つ(0 ブロック)" 19,
    .implicitStep
      "★★Dedekind の instance が入れば、`div(f)` は【主分数イデアルの素分解】として mathlib から直接出る。因子群を自分で積む必要は無い(0 ブロック)" 19,
    .implicitStep
      "★因子が `n` で割れる ⟹ `n` 乗根が取れる、の段。Dedekind なら分数イデアル `J` で `J^n = (μ f_P)` なるものが取れ、`J` が単項であることを `toClass` の単射性で示す形になる(10-25 ブロック)" 19,
    .implicitStep
      "★★★★**もう一つの道(イデアルの引き戻し)と、その訂正**: 当初「`[n]^* I := μ(I) が生成する CR-部分加群` と定めれば因子群を経由せずに済む」と書いたが、★**これは誤りである**——`FF` の `CR`-加群構造は `algebraMap` によるもので `μ` によるものではないため、`μ('' (f))` の張る部分加群は `(μ f)` にならない。★★正しくは `μ` を体の埋め込みに延ばし、`μ(CR)` の `FF` における整閉包を取る必要があり、それが `CR` と一致するには**結局 `IsIntegrallyClosed CR` が要る**。★★★2 つの道は `IsIntegrallyClosed` で合流する" 19,
    .implicitStep
      "★★★★★**最終段の単元は定数である**——`g^n` と `μ(f_P)` が同じイデアルを生成することから両者は単元倍で一致するが、`isUnit_coordinateRing`(第 128 ブロック、**Found に済**)により単元は定数であり、代数閉体なら `n` 乗根が取れて吸収できる(0 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
