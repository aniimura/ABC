import ABC3.Skeleton.GaloisRep.WeilFunctionField

/-!
# スケルトン —— **`n` 乗根 `g_P` の取り出し**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★層 3 —— `f_P` と `[n]^*` の合流点

第 113 ブロック(`Found`)で `f_P` が、第 118・125 ブロックで `[n]^*` が入った。
★この節点はその**合流点**である:

    g_P^n = f_P ∘ [n]                    (`exists_nthRoot_comp_mulByN`)

★★古典的には「`div(f_P ∘ [n])` が `n` で割れる」という因子の計算であり、
そこから `g_P` の存在が出る。

## ★★★★★2026-08-20 の再定式化

当初は `[n]^*` を **`F(W) →ₐ[F] F(W)`** として受け、`x([n]P) = Φ_n/ΨSq_n` で
固定していた。★しかし**実際に消費されるのは `f_P ∈ F[W]` への作用だけ**である。

★★第 125 ブロックで `μ : F[W] →+* F(W)` が**群法則で固定された形で**取れたので、
そちらを受ける形に改めた:

    n • (生成点) = Point.some (μ (genX W)) (μ (genY W)) h

★★★これは `Φ_n/ΨSq_n` より**強い**固定である(群法則そのもの)。
式の側は臨界路から外れた——`WeilFunctionField.lean` に非消費の節点として残す。

## ★★★★★ここが Weil 対の**心臓**である

`e_n(P,Q) := g_P(·+Q)/g_P(·)` が定数になるのは、
`g_P(·+Q)^n = f_P([n]·+[n]Q) = f_P([n]·)` (`[n]Q = 0` だから) `= g_P(·)^n`
——すなわち比の `n` 乗が 1 だからである。★**`g_P` の存在なしにはこの議論は立たない。**
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F]

/-- ★★★★★**`g_P` の存在**——`f_P ∘ [n]` は関数体の中で `n` 乗根を持つ。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`f_P` は第 113 ブロック(`Found/GaloisRep/TorsionIdeal.lean`)で
`(XYIdeal' P)^n` の生成元として得られている。
★★`μ` は第 118・125 ブロックの `exists_mulByNHom_charZero` が与える
——**群法則で固定されている**(`n • 生成点` の座標)。 -/
theorem exists_nthRoot_comp_mulByN (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n) (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn) :
    ∃ (fP : W.CoordinateRing) (g : W.FunctionField),
      ((((CoordinateRing.XYIdeal' h) ^ n : (FractionalIdeal (nonZeroDivisors W.CoordinateRing)
          W.FunctionField)ˣ) : FractionalIdeal (nonZeroDivisors W.CoordinateRing)
          W.FunctionField) : Submodule W.CoordinateRing W.FunctionField)
        = Submodule.span W.CoordinateRing
            {algebraMap W.CoordinateRing W.FunctionField fP}
      ∧ g ^ n = μ fP := by
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
      "★★★★★★`f_P` の存在は **Found に済**(第 113 `xyIdeal_pow_isPrincipal`)。`[n]` の引き戻し `μ` も **Found に済**(第 118・125 `exists_mulByNHom_charZero`、群法則で固定)。本節点はその合流点である(0 ブロック)" 19,
    .implicitStep
      "★★生成元が `F[W]` の中に取れること——`XYIdeal'` は `XYIdeal`(整イデアル)に等しい(mathlib の `XYIdeal'_eq`)ので、生成元は座標環の元である(3-8 ブロック)" 19,
    .implicitStep
      "★★★`div(f_P ∘ [n])` が `n` で割れることの因子計算。`[n]` が次数 `n²` であること、および `[n]^{-1}(P)` の各点の重複度を使う(20-50 ブロック)" 19,
    .implicitStep
      "★★★★因子が `n` で割れる ⟹ 関数が `n` 乗根を持つ、の段。関数体の乗法群と因子群の完全列(主因子の列)が要る——mathlib には `ClassGroup` はあるが**因子群そのものは無い**(2026-08-20 実測)(15-40 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
