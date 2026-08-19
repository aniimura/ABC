import ABC3.Skeleton.GaloisRep.WeilFunctionField

/-!
# スケルトン —— **`n` 乗根 `g_P` の取り出し**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★層 2 —— `f_P` と `[n]^*` を合わせる

第 113 ブロック(`Found`)で `f_P` が、
`WeilFunctionField.lean` の葉 2 で `[n]^*` が入る。★この節点はその**合流点**である:

    g_P^n = f_P ∘ [n]                    (`exists_nthRoot_comp_mulByN`)

★★古典的には「`div(f_P ∘ [n]) = n²·(∑_{R ∈ E[n]} (P' + R)) − …` が `n` で割れる」
という因子の計算であり、そこから `g_P` の存在が出る。

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
★★`μ` は `WeilFunctionField.lean` の `exists_mulByNPullback` が与える。 -/
theorem exists_nthRoot_comp_mulByN (W : WeierstrassCurve.Affine F) {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n) (hP : n • (Point.some x y h) = 0)
    (μ : W.FunctionField →ₐ[F] W.FunctionField)
    (hμ : μ (coordX W) = polyToFF W (W.Φ (n : ℤ)) / polyToFF W (W.ΨSq (n : ℤ))) :
    ∃ fP g : W.FunctionField,
      ((((CoordinateRing.XYIdeal' h) ^ n : (FractionalIdeal (nonZeroDivisors W.CoordinateRing)
          W.FunctionField)ˣ) : FractionalIdeal (nonZeroDivisors W.CoordinateRing)
          W.FunctionField) : Submodule W.CoordinateRing W.FunctionField)
        = Submodule.span W.CoordinateRing {fP}
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
    .otherPaper "GenEll" "Theorem 3.8(この節点を消費する主張)" 19,
    .implicitStep
      "★`f_P` の存在は **Found に入った**(第 113 ブロック `xyIdeal_pow_isPrincipal`)。当初「因子の層から積む」と見積もっていたが、mathlib が群法則を類群経由で証明しているため `Point.toClass` がそのまま使えた(0 ブロック)" 19,
    .implicitStep
      "★★`div(f_P ∘ [n])` が `n` で割れることの因子計算。`[n]` が次数 `n²` であること、および `[n]^{-1}(P)` の各点が重複度 `n²/#…` で現れることを使う(20-50 ブロック)" 19,
    .implicitStep
      "★★★因子が `n` で割れる ⟹ 関数が `n` 乗根を持つ、の段。関数体の乗法群と因子群の完全列(主因子の列)が要る——mathlib には `ClassGroup` はあるが**因子群そのものは無い**(2026-08-20 実測)(15-40 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
