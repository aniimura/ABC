/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuHeadDenomFree
import ABC3.Found.GaloisRep.MuDenomFreeSum
import ABC3.Found.GaloisRep.VeluDYDenomFree
import ABC3.Meta.Claim

/-!
# 第 1103 ブロック —— **`Lemma 3.5` を条なしにする残り**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か

`Found/GaloisRep/Lemma35Unconditional.lean` の `Lemma 3.5` は原典に無い仮説を
2 つ置いている（`l ≠ 2`・`d + 1 < l`）。本ファイルは **`d + 1 < l` を外す**残りを
節点に割り、進捗枠を置く。

☆`l ≠ 2` の方は別枠（`l = 2` の Vélu の分岐、20-40 ブロック）。

## ★★★★★★★★積んだ土台（`Found/`、`sorry` 0）

| 第 | 内容 |
|---|---|
| 1092 | `l^k·∑ f_i·inv(1−ζ^i)^k = ∑ f_i·∏_{j≠i}(1−ζ^j)^k` |
| 1097 | 6 種の頭項を分母なしの和に（橋） |
| 1098 | `DX`・`DY`・`D²X` の「対」を分母なしの和に（橋） |
| 1102 | **`E` 版 6 種の定義**（`S(η) = ∑ k·η^k` で無条件）と橋 |
| 1118-1119 | `tateXpairDF`・`tateYpairDF`・`veluV2DF` と橋 |
| 1125 | **DF 形は環準同型と可換**（`map_*DF`。★`a` 側の単元性を要らない） |

★機構は在庫の `one_sub_mul_sum_nsmul`——`(1 − η)·S(η) = −l`、`IsUnit` 不要。

### ★★★★★★★★★★★★第 1125-1128 で節点 3 が閉じた

☆当初は `ℤ[ζ_l] → ℚ(ζ_l)` を万有な環に取る設計だったが、**円分環は要らなかった**。

| 段 | 場所 | 何をするか | 第 |
|---|---|---|---|
| 1 | `A₁ = PowerSeries K`（`K` は体） | `1 − C α` も `(l)` も単元なので在庫の `hu` つきの補題が効く | 1128 |
| 2 | `PowerSeries A → PowerSeries (FractionRing A)` | `PowerSeries.map` の単射性で降ろす | 1128 |
| 3 | `PowerSeries R → R`（`X ↦ q`） | `evalAdicMapHom` で特殊化する | 1126 |

★通り道は第 1125 の `map_*DF`——第 221 の `map_tateXpair` は `IsUnit (1 − a)` を
要求するが、DF 形は要求しない（頭項が `Ring.inverse` を含まない多項式だから）。

★★**`PowerSeries A` は係数環 `A` を取り替えても `(X)`-adic 完備である**——
第 1091 で「商体に移すと `q` の収束が壊れる」と測った行き止まりの抜け道がこれであった。

☆側条件 `hDX ≠ 0` は `A₁` の定数項 `α(1+α)/(1−α)³` を見るだけで出た（第 1128）。

## ★★★★★★★★★★残り 5 節点（進捗枠 **5 / 5**）

| # | 節点 | 内容 | 重み |
|---|---|---|---|
| 1 | `sum_mu_dxpairE_zero` | `∑ DX` の `E` 版（`hu` なし）**★第 1115 で証明済み** | 8 |
| 2 | `sum_mu_d2xpairE` | `∑ D²X` の `E` 版（`hu` なし）**★第 1116 で証明済み** | 10 |
| 3 | `veluV2DF_eq_tateDYpairDF` | Vélu の `v` と `DY` の一致の `E` 版 **★第 1128 で証明済み** | 6 |
| 4 | `c4_velu_tateDF` / `c6_velu_tateDF` | `c₄`・`c₆` の恒等式の `E` 版 **★第 1129-1130 で証明済み** | 10 |
| 5 | `minDeltaExp_eq_mul_at_bad_prime_full_K` | `p ∣ l` の悪い素点で `Δ_min` が `l` 倍 **★第 1131-1138 で証明済み** | 12 |

☆総重み 46。★★★★**5 節点すべてが閉じた（2026-09-01、第 1138）**。

### ★★★★★★★★★★★★節点 5 の設計（第 1131 で確定）

`exists_vw_tate_mu`（`TateIsogeny.lean:1794`、第 1003）は

* `hlu : IsUnit ((l : R))`（＝ `p ∤ l`）
* Vélu の係数 `v`・`w` が **`R` に取れる**こと

を要求する。★`p ∣ l` では**どちらも成り立たない**——`v = ∑ veluV2` の分母に `l` が入る。

☆**しかしそれは障害ではない**。第 1129-1130 の分母なし恒等式を商体 `K` の中で
`l⁶`・`l⁸` で割ると

    `c₄(E_q/H) = l⁴·c₄(E_{q^l})`,   `c₆(E_q/H) = l⁶·c₆(E_{q^l})`,
    `Δ(E_q/H)  = l¹²·Δ(E_{q^l})`

になる（第 1131 の `exists_vw_tate_mu_field`・`exists_vw_tate_mu_field_Delta`）。
★これは「`E_q/H` は `E_{q^l}` を `u = 1/l` で変数変換したもの」という意味であり、
**`Δ_min` は極小モデルの判別式だからそのずれは消える**。

☆`hlu` を消費しているのは次の 3 か所だけである（第 1132 で実測）。

| # | 消費者 | 何のために | 状態 |
|---|---|---|---|
| A | `j_velu_tate_mu_map` | `j(veluCurve) = j(E_{q^l})` | ★第 1132 の `j_eq_of_c4_c6` |
| B | `pointCoords_tatePtPair` | `μ_l` の点の座標を書く | ★第 1133 の `natCast_pow_mul_tateXK` |
| C | `veluQuotientFull_tate_mu` / `isElliptic_veluQuotient_tate_mu` | 商の形と楕円性 | ★第 1134-1135 の `_K` 版 |

☆**3 つとも `hlu` なしの代替が揃い（第 1132-1135）、
それを使って `minDeltaExp` の連鎖も組み直した**（第 1136-1138）:

| 定理 | `hlu` なし版 | 第 |
|---|---|---|
| `j_veluQuot_eq_j_tate_pow` | `j_veluQuot_eq_j_tate_pow_K` | 1136 |
| `minDeltaExp_eq_mul_of_globalVelu'` | `…_K` | 1137 |
| `exists_vw_tate_mu` | `exists_vw_tate_mu_K`（★無条件） | 1138 |
| `minDeltaExp_eq_mul_at_bad_prime` / `_vc` / `_full` | `…_K` | 1138 |

★★**`minDeltaExp_eq_mul_at_bad_prime_full_K` が節点 5 の本体である**。

★**A は閉じた**——`j = Δ⁻¹·c₄³` なので `c₄ ↦ l⁴c₄`・`Δ ↦ l¹²Δ` では `l¹²` が約分される
（`j_of_c4_Delta`・`j_eq_of_c4_c6`、第 1132）。

★**B の核も閉じた**（第 1133）——`K` 水準の座標 `tateXK` は最初から

    `tateXK a w q hq = φ(tateXpairE a w q hq) · φ(1 − a)⁻²`

と**分母を払った形で定義されており、`R` の単元性を要求していない**（要るのは `φ(1−a) ≠ 0`）。
☆`(1 − a)·S(a) = −l` から `(1 − a)²·tateXpairDF = l²·tateXpairE` が `R` の中で出るので、

    `l²·tateXK = φ(tateXpairDF)`,   `l³·tateYK = φ(tateYpairDF)`

が `IsUnit` 抜きで成り立つ。

### ★★★降下の道（第 1128 で確定した最終形）

節点 1-2（和の形）は `ℤ[ζ_l] → ℚ(ζ_l)` の単射性で降ろした（第 1112-1116）。
節点 3（点ごとの恒等式で尾を含むもの）は `PowerSeries` を経由した（第 1125-1128）。
★**尾を含む主張は係数環の取り替えでは降りない**——`FractionRing` は
`I`-adic 完備でないからである。`PowerSeries` はそこを迂回する。
-/

namespace ABC3.Skeleton.GenEll

open Finset ABC3.Meta ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R} [IsAdicComplete I R]

/-- ★★★★★★★★**節点 1**——`∑ DX` の `E` 版（`hu` を受けない）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`sum_mu_dxpair_zero`（`TateODE.lean:144`）の分母なし版である。 -/
theorem sum_mu_dxpairE_zero {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0,
        (tateDXtermE l (ζ ^ i)
          + (l : R) ^ 3 * (tateDXtail (ζ ^ i) q hq - tateDXtail (ζ ^ (l - i)) q hq)) = 0 :=
  ABC3.Found.GaloisRep.sum_mu_dxpairE_zero hl hζ q hq

def sum_mu_dxpairE_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(節点 1——∑ DX の E 版。hu を受けない)",
    sectionId := "genell-lemma-3-5" }

def sum_mu_dxpairE_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "natCast_pow_mul_tateDXterm(橋、第 1102、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.natCast_pow_mul_tateDXterm") 1,
    .implicitStep
      ("★`ℤ[ζ_l]` は ℤ-捻れ自由なので `ℚ(ζ_l)` へ単射。" ++
       "☆`ℚ(ζ_l)` では `1 − ζ^i` が可逆なので既存の `sum_mu_dxpair_zero` が使え、" ++
       "単射性で降ろせば `hu` が消える。") 8 ]

/-- ★★★★★★★★★★**節点 2**——`∑ D²X` の `E` 版（`hu` を受けない）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`sum_mu_d2xpair`（`TateODE.lean:207`）の分母なし版。★`c4_velu_tate` の本体である。 -/
theorem sum_mu_d2xpairE {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    120 * ∑ i ∈ (range l).erase 0,
        (tateD2XtermE l (ζ ^ i)
          + (l : R) ^ 4 * (tateD2Xtail (ζ ^ i) q hq + tateD2Xtail (ζ ^ (l - i)) q hq))
      = (l : R) ^ 4 * (((l : R) ^ 4 - 1)
          + 240 * ((l : R) ^ 4 * evalAdic (sigmaSeries 3) (q ^ l) hql
              - evalAdic (sigmaSeries 3) q hq)) :=
  ABC3.Found.GaloisRep.sum_mu_d2xpairE hl hζ q hq hql

def sum_mu_d2xpairE.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(節点 2——∑ D²X の E 版。c4_velu_tate の本体)",
    sectionId := "genell-lemma-3-5" }

def sum_mu_d2xpairE.needs : List ProofObligation :=
  [ .citation "[ABC3]" "natCast_pow_mul_tateD2Xterm(橋、第 1102、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.natCast_pow_mul_tateD2Xterm") 1,
    .implicitStep "☆節点 1 と同じ降下（`ℤ[ζ_l] → ℚ(ζ_l)` の単射性）。" 10 ]

/-- ★★★★★★★★★★★★**節点 3**——Vélu の `v` と `DY` の一致の `E` 版（第 1128）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`veluV2_eq_tateDYpair`（`TateODE.lean:89`）の分母なし版である。
★★**`IsUnit (1 − α)` も `IsUnit (l)` も取らない**——`p ∣ l` の悪い素点でも成り立つ。 -/
theorem veluV2DF_eq_tateDYpairDF {l : ℕ} (hl0 : (l : R) ≠ 0) {α : R}
    (hα1 : α ≠ 1) (hα0 : α ≠ 0) (hαneg : α + 1 ≠ 0)
    (hpow : α ^ l = 1) (hsum : ∑ k ∈ range l, α ^ k = 0) (q : R) (hq : q ∈ I) :
    veluV2DF l (tateCurveAt q hq) (tateXpairDF l α (α ^ (l - 1) * q) q hq)
        (tateYpairDF l α (α ^ (l - 1) * q) q hq)
      = (l : R) ^ 2 * tateDYpairDF l α (α ^ (l - 1) * q) q hq :=
  ABC3.Found.GaloisRep.veluV2DF_eq_tateDYpairDF hl0 hα1 hα0 hαneg hpow hsum q hq

def veluV2DF_eq_tateDYpairDF.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(節点 3——veluV2 = DY の E 版。hu も IsUnit(l) も取らない)",
    sectionId := "genell-lemma-3-5" }

def veluV2DF_eq_tateDYpairDF.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluV2DF_eq_tateDYpairDF(Found の本体、第 1128、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluV2DF_eq_tateDYpairDF") 1,
    .implicitStep
      ("★★**2026-09-01（第 1128）**——`PowerSeries` を万有な完備環として使った。" ++
       "☆`PowerSeries A` は係数環 `A` を取り替えても `(X)`-adic 完備なので、" ++
       "`PowerSeries (FractionRing A)` で `hu` を得て単射性で降ろせる。" ++
       "★これが第 1091 で測った行き止まり（商体では `q` の収束が壊れる）の抜け道である。") 6 ]

end ABC3.Skeleton.GenEll
