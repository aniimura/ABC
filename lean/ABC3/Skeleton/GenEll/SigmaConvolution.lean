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
  [ .citation "[mathlib]" "ArithmeticFunction.sigma の畳み込みの明示式"
      (.absent "mathlib を 'sigma.*convolution|sigma_one_mul_sigma_one' で grep して 0 件(2026-08-31)") 8,
    .citation "[mathlib]" "EisensteinSeries.G2_S_transform（E₂ の準モジュラー変換則）"
      (.inMathlib "EisensteinSeries.G2_S_transform") 2,
    .citation "[mathlib]" "EisensteinSeries.E_qExpansion_coeff（E_k の q 展開係数は σ_{k−1}）"
      (.inMathlib "EisensteinSeries.E_qExpansion_coeff") 2,
    .implicitStep
      ("☆道 A（モジュラー形式）: E₄ = E₂² − 12·q dE₂/dq を使う。"
       ++ "mathlib には E₂ の S 変換則（G2_S_transform）と E_k の q 展開があるが、"
       ++ "★E₂ の q 展開（EisensteinSeries/E2/ に QExpansion.lean は無い）と"
       ++ "準モジュラー形式の微分の理論が要る") 8,
    .implicitStep
      ("★★道 B（本プロジェクトの在庫）: Tate 曲線の方程式"
       ++ "Found/GaloisRep/TateEquation.lean の tate_equation（Y² + XY = X³ + a₄X + a₆、証明済み）は"
       ++ "**2 変数 (u, q) の恒等式**である。u の冪ごとに係数を比べると"
       ++ "σ の恒等式が出る——ラマヌジャンの式もその 1 つのはずである。"
       ++ "☆こちらはモジュラー形式の理論を要らない") 6,
    .citation "[ABC3]" "veluVC_zero(この結論の消費側、§9-1241)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluVC_zero") 2 ]

end ABC3.Skeleton.GenEll
