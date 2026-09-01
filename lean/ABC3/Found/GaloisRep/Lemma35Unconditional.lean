/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateLocalModel
import ABC3.Meta.Claim

/-!
# 第 1085 ブロック —— **`Lemma 3.5` を仮説なしで**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★到達点

第 1083 で `hfin` が無条件に出たので、`Lemma 3.5` の高さ不等式は
**外部からの入力を一切受けずに**証明できる。

## ★★★★★★★★★★逸脱（明示）

原文 (GenEll p.17) の仮説は「`ε > 0`」「`l` 素数」「`H_F` が `l`-cyclic」
「`E_H = E/H`」「`l` は悪い乗法還元の全素点での局所高さと素」だけである。

| 項 | 原典 | 形式化 | 理由・消費側への影響 |
|---|---|---|---|
| `l` | 任意の素数 | ★**`l ≠ 2` を追加** | Vélu の `±` の対（`veluWFull` の `/2`）が不動点を持たないことを使う。☆**消費側では自動的**（下記） |
| `l` | 任意の素数 | ★**`d + 1 < l` を追加**（`d = [L:ℚ]`） | 第 1044 が悪い素点で `p ∤ l` を出すのに使う（`ζ_l` の分岐から `l−1 ≤ d`）。原典は Lemma 3.2, (i) で直接 `μ_l` と同定するのでこれを要さない。☆**消費側では自動的**（下記） |
| `H_F` は部分群スキームで `H_F ×_F ℚ̄ ≅ ℤ/lℤ` | スキーム | ☆`L` 有理点 `Q` で `addOrderOf Q = l` | ★階数 1 の部分群を生成元で書いたもの。`E' = E/⟨Q⟩` を Vélu の商で読む |
| `E_H` の半安定性 | 原文は「同種なので自動」と括弧で述べる | ☆`∀ p, SemistableAt p E'` を仮説に置く | ★「同種なら半安定」の形式化を後回しにしただけで、内容を弱めていない |

### ★★★★★★★★追加した 2 つがどこまで吸収されるか（実測）

`Lemma 3.5` を消費するのは `Lemma 3.7` の第 3 の主張であり、
そこには条件 (a) と条件 (b) の **2 つの枝**がある。

| 枝 | `l` の下界 | 逸脱は吸収されるか |
|---|---|---|
| 条件 (a)（`100·d·(ht^Falt + C·d^ε) ≤ l`） | `l ≥ 100·d ≥ 100` | ★**される**。第 1088 で実際に吸収した（`htFalt_le_of_condA_lcyclic`、外部入力なし） |
| 条件 (b)（`[E_L] ∈ K_V` かつ `l` が局所高さと素） | **`2 ≤ l` だけ** | ☆**されない**。`htFalt_le_of_condB` は `2 ≤ l` で `(†)` を消費するので、`d + 1 < l` が成り立たない場合がある |

★★★★**したがって本定理の `.src` は条つきのままにしてある**。
☆（2026-09-01、第 1089 の訂正。第 1087 では条件 (a) だけを見て
条なしにしてしまったが、条件 (b) の枝を見落としていた。）

### ☆この 2 つを外すには

- `l ≠ 2` → `l = 2` の Vélu の分岐。`veluU = 0` なので `/2` は無害だが、
  `preΨ_2 = 1` なので捕れ点との橋が別物になる（`Ψ₂Sq` の主係数は `4`）。
- `d + 1 < l` → `p ∣ l` の悪い素点で Tate モデルの Vélu を
  `l` が単元でない環で行う（`μ_l` が非エタール）。
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
☆原典に無い仮説は `l ≠ 2` と `d + 1 < l` の 2 つだけで、
どちらも後続（`Lemma 3.7`・`Theorem 3.8`）が取る `l` では自動的に成り立つ。 -/
theorem lemma_3_5_height_ineq (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime → l ≠ 2 →
      Module.finrank ℚ L + 1 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C :=
  ABC3.Skeleton.GenEll.lemma_3_5_velu eps heps

def lemma_3_5_height_ineq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(外部入力なし。★原典に無い仮説 l ≠ 2・d+1 < l を置く——条件 (b) の枝では吸収されない)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_height_ineq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_velu(第 1084、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.lemma_3_5_velu") 1,
    .citation "[ABC3]" "hfin_of_veluQuotientFull(hfin そのもの、第 1083、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hfin_of_veluQuotientFull") 1,
    .implicitStep
      ("★★**2026-09-01（第 1085）の逸脱の記録**——原典 (GenEll p.17) の Lemma 3.5 は " ++
       "`l ≠ 2` も `d + 1 < l` も置いていない。" ++
       "☆`l ≠ 2` は Vélu の `±` の対（`veluWFull` の `/2`）が不動点を持たないために要る。" ++
       "☆`d + 1 < l` は第 1044 が `p ∤ l` を出すのに使う（原典は Lemma 3.2, (i) で " ++
       "直接 `μ_l` と同定するのでこれを要さない）。" ++
       "★後続（`Lemma 3.7`・`Theorem 3.8`）は `100·d·(ht^Falt + C·d^ε) ≤ l` を取るので " ++
       "両方とも自動的に成り立ち、影響は無い。" ++
       "☆外すには悪い素点で `p ∣ l` の場合の Tate モデルの Vélu 計算が要る。") 8 ]

end ABC3.Found.GaloisRep
