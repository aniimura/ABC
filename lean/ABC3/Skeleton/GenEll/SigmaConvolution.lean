/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateSeries
import ABC3.Meta.Claim

/-!
# 約数関数の畳み込み恒等式（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★これは何か（2026-08-31、第 822 で発見）

`tateModel_of_quot_mu`（葉 1）の**手順 4（係数の照合）**を進めると、
Vélu の和 `v = ∑_ζ (3X(ζ)² + a₄ − Y(ζ))` の `X(ζ)²` から

    `∑_{m=1}^{n−1} σ₁(m)·σ₁(n−m)`

という **`σ₁` どうしの畳み込み**が出る。一方 `a₄ = −5s₃` なので、
照合の相手は `σ₃` である。★★両者を繋ぐのが**ラマヌジャンの恒等式**

    `12·∑_{m=1}^{n−1} σ₁(m)σ₁(n−m) = 5σ₃(n) − (6n−1)σ₁(n)`

である（`n = 2,…,11` で数値確認済み）。

## ★★★これは既知の古典的事実である

Eisenstein 級数の言葉では `E₂² = E₄ + (E₂ の微分の項)`、あるいは
`E₄ = E₂² − 12·q d/dq E₂` に対応する。★初等的な証明（Liouville の方法、
あるいは `℘` 関数の Laurent 展開）も知られている。

☆mathlib には**無い**（`ArithmeticFunction.sigma` の畳み込みの明示式は
`sigma_one_eq_sigmaOne` 等の基本性質だけ）。

## ★次に何をすれば終わりか

1. 本補題を証明する（Eisenstein 級数か初等的方法）
2. `Found/GaloisRep/VeluMuSum.lean` の `veluVC_zero` に代入する
3. `tateA4(q^l) = tateA4(q) − 5v` を得る
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep

/-- **[GenEll] 葉 1 の手順 4 に要る古典的恒等式**——ラマヌジャンの畳み込み公式。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★`σ₁` どうしの畳み込みが `σ₃` で書ける。これが Tate 曲線の `a₄` の
変換式（`a₄(q^l) = a₄(q) − 5v`）の数論側の核である。

☆数値確認: `n = 2,…,11` で成立（第 822）。 -/
theorem sigma_one_convolution (n : ℕ) (hn : 2 ≤ n) :
    12 * ∑ m ∈ Finset.Ico 1 n, (ArithmeticFunction.sigma 1 m : ℤ)
        * (ArithmeticFunction.sigma 1 (n - m) : ℤ)
      = 5 * (ArithmeticFunction.sigma 3 n : ℤ)
        - (6 * (n : ℤ) - 1) * (ArithmeticFunction.sigma 1 n : ℤ) := by
  sorry

def sigma_one_convolution.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(ラマヌジャンの畳み込み恒等式——σ₁ ⋆ σ₁ と σ₃)",
    sectionId := "genell-lemma-3-2" }

def sigma_one_convolution.needs : List ProofObligation :=
  [ .citation "[mathlib]" "ArithmeticFunction.sigma の加法畳み込みの明示式"
      (.absent "mathlib は Dirichlet 畳み込みしか持たない(2026-08-31)") 8,
    .citation "[mathlib]" "EisensteinSeries.G2_S_transform（E₂ の準モジュラー変換則）"
      (.inMathlib "EisensteinSeries.G2_S_transform") 2,
    .citation "[mathlib]" "EisensteinSeries.E_qExpansion_coeff（E_k の q 展開係数は σ_{k−1}）"
      (.inMathlib "EisensteinSeries.E_qExpansion_coeff") 2,
    .implicitStep
      ("★★★消費側の正確な形（2026-08-31、第 839、l = 5,7,11 ・n < 20 で数値確認）:"
       ++ "★ 3·[∑_{ζ∈μ_l∖{1}} X(ζ,q)²]_n"
       ++ " = l⁴σ₃(n/l)[l∣n] − σ₃(n) + 5(l−1)σ₃(n) + lσ₁(n) − l²σ₁(n/l)[l∣n]。"
       ++ "☆左辺は σ₁ の畳み込みを含み、右辺は σ₃ で書かれている——"
       ++ "この穋ぎにラマヌジャンの恒等式が要る") 2,
    .implicitStep
      ("★★道 B（本プロジェクトの在庫）: Tate 曲線の方程式 tate_equation（証明済み）は"
       ++ "2 変数 (u, q) の恒等式である。u の冪ごとに係数を比べると σ の恒等式が出る") 6,
    .implicitStep
      ("☆再定式化（第 829、n < 30 で数値確認）: 右辺は 1 つの約数和になる:"
       ++ "12·∑_{m<n} σ₁(m)σ₁(n−m) = ∑_{d∣n} (5d³ + d − 6d²·(n/d))。"
       ++ "☆左辺は ab+cd = n なる四つ組にわたる ∑ ac であり、対称性から"
       ++ "∑ac = ∑ad = ∑bc = ∑bd が出る") 2,
    .citation "[ABC3]" "veluV_coeff_of_ne_zero(この結論の消費側、§9-1254)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluV_coeff_of_ne_zero") 2 ]

end ABC3.Skeleton.GenEll
