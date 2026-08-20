import ABC3.Skeleton.GaloisRep.WeilDivisor
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
★★`μ` は第 118・125 ブロックで**群法則によって固定された形で**得られている。

## ★★★★★★2026-08-20: この節点自身の `sorry` は消えた

本定理は `WeilDivisor.lean` の 2 節点(D1・D2)と第 139 ブロックから**証明される**。
★`sorry` は D1・D2 に移った。

## ★逸脱の記録——`[IsAlgClosed F]` を足した

第 139 ブロックが「単元(= 定数)の `n` 乗根を取る」ために `F` の代数閉性を使う。
★消費する側((G5) `det_cyclotomic`)は `[IsAlgClosed L]` の下で述べられているので、
後続の証明に影響は出ない。 -/
theorem exists_nthRoot_comp_mulByN [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F)
    [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n) (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    ∃ g : W.FunctionField, g ^ n = μ fP := by
  obtain ⟨J, hJ⟩ := exists_fractionalIdeal_pow h2 W h n hn hP μ hμinj hμF hns hμP hμx hμy fP hfP
  obtain ⟨g, hg⟩ :=
    fractionalIdeal_isPrincipal h2 W h n hn hP μ hμinj hμF hns hμP hμx hμy fP hfP J hJ
  exact exists_nthRoot_of_fractionalIdeal W hn hJ hg

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_nthRoot_comp_mulByN.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——g_P の n 乗根の取り出し)",
    sectionId := "genell-thm-3-8" }

def exists_nthRoot_comp_mulByN.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1 の証明(g_P の存在)"
      (.absent "mathlib に Weil 対およびその構成要素は 0 件(2026-08-20、`WeilPairing|weil_pairing` で全文検索して 0 件)") 19,
    .implicitStep
      "★★★★★★**2026-08-20: 本節点自身の `sorry` は消えた**——`WeilDivisor.lean` の D1(`exists_fractionalIdeal_pow`)・D2(`fractionalIdeal_isPrincipal`)と第 139 ブロック(`exists_nthRoot_of_fractionalIdeal`)から証明される。残る `sorry` は D1・D2 の 2 つだけである(0 ブロック)" 19,
    .implicitStep
      "★★★★★D1 —— `∃ J, J^n = (μ f_P)`。**不分岐性は要らない**——`ord_Q(μ f) = e_Q · ord_{[n]Q}(f_P)` で右の因子が `{n, −n, 0}` だから、`e_Q` が何であっても `n` で割れる(20-50 ブロック、`WeilDivisor.lean` を見よ)" 19,
    .implicitStep
      "★★★★★D2 —— その `J` が単項であること。平行移動不変性から `e_Q` はファイバー上で一定であり、`Σ_{Q ∈ [n]⁻¹(P)} Q = n·P = 0`・`Σ_{T ∈ E[n]} T = 0` により類は 0(10-25 ブロック、`WeilDivisor.lean` を見よ)" 19,
    .implicitStep
      "★★★★**最後の一歩は Found に済**(第 139 `exists_nthRoot_of_fractionalIdeal`)——`(g^n) = (μ f_P)` から単元倍で一致し、第 128 により単元は定数、`F` が代数閉だからその `n` 乗根で吸収できる(0 ブロック)" 19,
    .implicitStep
      "★★★★★★**因子機構の土台は揃った**——`IsDedekindDomain F[W]`(第 137)、`HeightOneSpectrum ↔ 曲線上の点`(第 138)、`(f_P) = XYIdeal_P^n`(第 126)、`μ`(第 118・125)、`τ_T`(第 120・124)、単元 = 定数(第 128)(0 ブロック)" 19,
    .implicitStep
      "★逸脱の記録: 本定理に `[IsAlgClosed F]` を足した。第 139 が単元の `n` 乗根を取るために使う。消費側((G5) `det_cyclotomic`)は `[IsAlgClosed L]` の下で述べられているので後続に影響しない" 19 ]

end ABC3.Skeleton.GaloisRep
