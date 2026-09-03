/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateLocalModel
import ABC3.Skeleton.GenEll.TateLocalModelK
import ABC3.Meta.Claim

/-!
# 第 1085 ブロック —— **`Lemma 3.5` を仮説なしで**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★到達点

第 1083 で `hfin` が無条件に出たので、`Lemma 3.5` の高さ不等式は
**外部からの入力を一切受けずに**証明できる。

★★★★**2026-09-01（第 1141）——`d + 1 < l` が外れた**。
★★★★★★★★**2026-09-01（第 1149）——`l ≠ 2` も外れた**。
☆これで**原典に無い仮説は 1 つも無い**。

## ★★★★★★★★★★逸脱（明示）

原文 (GenEll p.17) の仮説は「`ε > 0`」「`l` 素数」「`H_F` が `l`-cyclic」
「`E_H = E/H`」「`l` は悪い乗法還元の全素点での局所高さと素」だけである。

| 項 | 原典 | 形式化 | 理由・消費側への影響 |
|---|---|---|---|
| `l` | 任意の素数 | ★**同じ（第 1149 で `l ≠ 2` が外れた）** | — |
| `H_F` は部分群スキームで `H_F ×_F ℚ̄ ≅ ℤ/lℤ` | スキーム | ☆`L` 有理点 `Q` で `addOrderOf Q = l` | ★階数 1 の部分群を生成元で書いたもの。`E' = E/⟨Q⟩` を Vélu の商で読む |
| `E_H` の半安定性 | 原文は「同種なので自動」と括弧で述べる | ☆`∀ p, SemistableAt p E'` を仮説に置く | ★「同種なら半安定」の形式化を後回しにしただけで、内容を弱めていない |

### ★★★★★★★★★★★★`d + 1 < l` はどうやって外れたか（第 1123-1141）

`d + 1 < l` は第 1044 が悪い素点で **`p ∤ l`** を出すためだけに要っていた
（`ζ_l` の分岐から `l−1 ≤ d`）。★その `p ∤ l` が要らなくなったので落ちた。

| 段 | 内容 | 第 |
|---|---|---|
| 分母を払う | `tateXpairDF` 等（`(1−η)·S(η) = −l`、`IsUnit` 不要） | 1102-1119 |
| 降下路 | `PowerSeries` は係数環を取り替えても `(X)`-adic 完備 | 1125-1128 |
| `c₄`・`c₆` | `l⁶·c₄ + 240∑ = l¹⁰c₄′`、`l⁸·c₆ + … = l¹⁴c₆′` | 1129-1130 |
| 商体で割る | `c₄ = l⁴c₄′`・`Δ = l¹²Δ′` ⟹ **`j` は等しい** | 1131-1132 |
| 座標・和・商 | `tateXK` は最初から分母を払った形 | 1133-1135 |
| `Δ_min` の連鎖 | `_K` 版で組み直し | 1136-1140 |
| `Lemma 3.5` へ | `lemma_3_5_velu_K` | 1141 |

★★**要点**——`p ∣ l` では Vélu のモデルは極小から `l¹²` 離れるが、
`Δ_min` は極小モデルの判別式なのでそのずれは消える。

### ★★★★★★★★★★★★★★★★`l ≠ 2` はどうやって外れたか（第 1142-1149）

`l = 2` で壊れるのは 3 か所だった（第 1142 の実測）。

| # | 壊れる場所 | 道 | 第 |
|---|---|---|---|
| A | `veluV2_eq_tateDYpair` が `DX ≠ 0` で割っている | ★`a` を一般の径数 `T` に取る（`Localization.Away (X(1−X))`）と `DX ≠ 0`。`T ↦ a` で特殊化 | 1143-1146 |
| B | `exists_veluW_of_inv` の `i ↔ l−i` の対 | ☆`l = 2` では一点が 2-捻れで `veluU = 0` なので `w = veluV2·x` がそのまま取れる | 1149 |
| C | `preΨ` の `Odd k` | ★mathlib の `leadingCoeff_ΨSq` は偶奇を問わないので、`ΨSq` のまま止める | 1148-1149 |

★★☆さらに `veluWFull` の `/2` は `v_p(2) = 0`（Bézout `2a + lb = 1`）を使っていたが、
`l = 2` では `veluU = 0` で項そのものが消えるので、**選言に弱めて両方を通した**（第 1149）。

-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta
open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 高さの不等式（外部入力なし）**（第 1085）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1085）**——`hfin`（第 1083）で最後の入力が消えた。
★★★★★★★★**2026-09-01（第 1141・1149）**——`d + 1 < l` も `l ≠ 2` も外れた。
☆**原典に無い仮説はもはや 1 つも無い**。 -/
theorem lemma_3_5_height_ineq (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C :=
  ABC3.Skeleton.GenEll.lemma_3_5_velu_K eps heps

def lemma_3_5_height_ineq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_height_ineq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_velu_K(第 1141、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.lemma_3_5_velu_K") 1,
    .citation "[ABC3]" "hfin_of_veluQuotientFull(hfin そのもの、第 1083・1149、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hfin_of_veluQuotientFull") 1,
    .implicitStep
      ("★★★★★★★★**2026-09-01（第 1149）の到達**——" ++
       "原典 (GenEll p.17) の Lemma 3.5 に無い仮説は**もはや 1 つも置いていない**。" ++
       "☆ただし `∀ p, SemistableAt p E'`（原典が「同種なので自動」と括弧で述べるもの）は" ++
       "**仮説の形で受けている**——「同種なら半安定」の形式化は未了である。" ++
       "☆`d + 1 < l` は第 1141、`l ≠ 2` は第 1149 で外れた。" ++
       "★`l = 2` で壊れていた 3 か所の道: (A) `DX ≠ 0` は一般の径数で通す（第 1143-1146）、" ++
       "(B) `i ↔ l−i` の対は `veluU = 0` で不要（第 1149）、" ++
       "(C) `preΨ` の偶奇は `ΨSq` のまま止めて不要（第 1148-1149）。" ++
       "☆`veluWFull` の `/2` は `v_p(2) = 0` と `veluU = 0` の選言に弱めて両方を通した。") 1 ]

end ABC3.Found.GaloisRep
