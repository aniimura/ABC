import ABC3.Found.GenEll.LatticeFubini
import ABC3.Found.GaloisRep.TateEquationAnalytic

/-!
# GenEll 第 342 ブロック —— **★★★★★`ζ(6) = π⁶/945` の級数形**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★なぜ `ζ(6)` が要るのか

第 341 で `Σ_{v≠0} (v₁+v₂τ)⁻ᴷ = (Σ_{d>0} d⁻ᴷ)·(Σ_{w 原始} …)` まで来た。
★`Σ_{d>0} d⁻ᴷ = ζ(K)` であり、`g₂³ − 27g₃²` の係数を合わせるのに
**`ζ(4) = π⁴/90` と `ζ(6) = π⁶/945` の両方**が要る:

    (120·ζ(4))³ = (4π⁴/3)³ = 64π¹²/27
    27·(280·ζ(6))² = 27·(8π⁶/27)² = 64π¹²/27

★★この一致がなければ `g₂³ − 27g₃²` は `E₄³ − E₆²` の定数倍にならず、
`Δ ≠ 0` から (i) を出せない。

## ★★★★★在庫の確認(2026-08-26)

★**Bernoulli 数の値と `riemannZeta 6` は既にある**
——`Found/GaloisRep/TateEquationAnalytic.lean` の
`bernoulli'_three` … `bernoulli_six_val`・`riemannZeta_six`(Tate 曲線の解析で積んだもの)。
★★本ブロックはそれを**再利用**し、足りない**級数形**(`HasSum` の実数版)だけを足す
——格子和は `tsum` で書かれているので、`riemannZeta` の形ではなく級数形が要る。

★★★★**別の目的で積んだ在庫がそのまま効いた**——Tate 曲線の `q` 展開のために
Bernoulli 数を計算してあったのが、Eisenstein 級数の正規化で使えた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `hasSum_zeta_six` | ★★★★★**`Σ_{n} 1/n⁶ = π⁶/945`(級数形)** |
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★★`ζ(6)` の級数形 -/

/-- ★★★★★**`Σ_n 1/n⁶ = π⁶/945`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★mathlib は `hasSum_zeta_two`・`hasSum_zeta_four` と一般式 `hasSum_zeta_nat` を持つが
`ζ(6)` の級数形は持たない(2026-08-26 実測)。
★★`bernoulli 6 = 1/42` は `Found/GaloisRep/TateEquationAnalytic.lean` の
`bernoulli_six_val`(Tate 曲線の解析で積んだもの)を使う。 -/
theorem hasSum_zeta_six : HasSum (fun n : ℕ => 1 / (n : ℝ) ^ 6) (Real.pi ^ 6 / 945) := by
  have h := hasSum_zeta_nat (k := 3) (by norm_num)
  norm_num [ABC3.Found.GaloisRep.bernoulli_six_val] at h ⊢
  convert h using 1
  ring

/-! ## ★出典の紐付け(`.src`) -/

def hasSum_zeta_six.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
