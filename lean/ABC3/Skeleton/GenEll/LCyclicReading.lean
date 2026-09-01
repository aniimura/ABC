/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Lemma37Full
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Meta.Claim

/-!
# 第 1152 ブロック —— **`l`-巡回部分群の 2 つの読み方が食い違っている**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18–p.20。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★★★★★これは何か——**逸脱が消費側に響いた**

原文の `H_L ⊆ E_L` は **`l`-巡回部分群スキーム**である。本企画はこれを
**2 通りに読んで**しまっており、その 2 つが繋がっていない。

| 読み | 定義 | どこで使うか |
|---|---|---|
| (V) `HasLCyclicVelu` | ☆**`L` 有理点** `Q` で `addOrderOf Q = l`、その Vélu の商が楕円曲線で半安定 | `Lemma 3.5` → `Lemma 3.7`（第 1151） |
| (J) `HasLCyclicJ` | ☆`E[l]` の中の **`Gal`-安定な直線** | `Theorem 3.8`（`imageContainsSL2J_of_alpha'` の `hno`） |

★**含意は片側だけ立つ**: 有理点 `Q` は `⟨Q⟩` を点ごとに固定するので

    (V) ⟹ (J),   したがって   ¬(J) ⟹ ¬(V)

である。☆(J) ⟹ (V) は**偽**——安定な直線の点は `L` 有理とは限らない。

## ★★★★★★★★★★★★食い違いの現れ方

`Lemma 3.7` の第 3 の主張は「`l`-巡回なら `[E] ∈ Exc`」、
`Theorem 3.8` はその対偶「`[E] ∉ Exc` なら `l`-巡回でない」を使う。

| | `Lemma 3.7` が与えるもの | `Theorem 3.8` が要るもの |
|---|---|---|
| 第 1151 の形 | `¬ HasLCyclicVelu` | `¬ HasLCyclicJ` |

★★**向きが合わない**——`¬(J) ⟹ ¬(V)` であって逆ではない。
☆すなわち第 1151 の `lemma_3_7` は `Theorem 3.8` を**まだ養えない**。

★★★これは `CLAUDE.md` の「逸脱」の規約が名指ししている場合である
——『後続の証明に影響が出ないならば読み替えを許容する』。**影響が出た**ので記録する。

## ★★★★★★★★★★★★★★★★節点（進捗枠 **0.7 / 3**）

`Lemma 3.7` の第 3 の主張を **(J) の側で**言い直すのに要るもの:

| # | 節点 | 内容 | 重み |
|---|---|---|---|
| 1 | `veluQuotient_of_stableLine` | ★**核は第 1153 で取れた**（`fixesCoeffs_veluQuotientFull`）。☆残る 1 歩は `L̄^{Gal} = L` で「すべての `σ` で固定 ⟹ `L` の元」と言う段 | 12 → **4** |
| 2 | `lemma_3_5_height_ineq_stableLine` | ☆第 1 で作った `L` 上の商に対して `Lemma 3.5` の高さ不等式を回す（`degInf` の関係は `Lemma 3.2, (ii)` から） | 10 |
| 3 | `lemma_3_7_stableLine` | ☆第 3 の主張を `HasLCyclicJ` の側で述べ直す。★これで `Theorem 3.8` に繋がる | 4 |

☆総重み 26 → **18**（第 1153 で 8 分進んだ）。

### ★★★★★★★★第 1153 で取れたもの

`Found/GenEll/VeluGaloisInvariant.lean`。★いずれも**無条件**である。

| 定理 | 内容 |
|---|---|
| `veluGx_semi` / `veluGy_semi` | 係数を固定する `Φ` に対し `Φ(veluGx(x,y)) = veluGx(Φx,Φy)` |
| `image_semiPair_eq` | `Φ` で保たれる有限集合は `Φ` の像で自分自身（単射だから） |
| `veluVFull_semi` / `veluWFull_semi` | 座標の集合が `Φ` で保たれるなら和は `Φ` で不変 |
| `fixesCoeffs_veluQuotientFull` | ★**商曲線の係数はすべて `Φ` で固定される** |

☆`/2` は `Φ 2 = 2` なので害が無い。★商の `a₁・a₂・a₃` は元のもののまま、
`a₄ = a₄ − 5v`・`a₆ = a₆ − b₂v − 7w` は `v`・`w` の不変性から出る。

### ★★★★節点 1 の道（測定）

`veluQuotientFull` は体 `F` と `S : Finset (F × F)` を受ける。★`F = L̄` で取れば
安定な直線の点の座標は `L̄` にあり、商 `veluCurve E v w` は `L̄` 上で定義される。
☆`H` が `Gal(L̄/L)`-安定なら `S` は `Gal` の作用で保たれるので、

    `v = Σ_{Q ∈ S} veluV2(Q)`,   `2w = Σ_{Q ∈ S} (veluU(Q) + 2 veluV2(Q)·x_Q)`

はいずれも `Gal` 不変、すなわち `L` の元である（`L̄^{Gal} = L`）。
★したがって `veluCurve E v w` は `L` 上の曲線の底変換である。

☆在庫: `Found/GenEll/VeluPointSet.lean` の `veluQuotientFull_image_eq`（座標の像で取った商）と
`Found/GaloisRep/VeluTateMuK.lean` の `veluQuotientFull_points_eq_field` が素材である。

### ☆別の道（採らない理由）

`L(H)`（`H` の点を有理化する体、次数は `l−1` を割る）へ上げる道もあるが、
★`ht^Falt`・`deg∞` が体を変えると動くので、`Lemma 3.7` の数値の帳簿を
書き直すことになる。☆**節点 1 の降下の方が安い**。

## ★★★★★★何が壊れていないか

☆第 1151 の `lemma_3_7`（`.src` は条なし）は**それ自体としては正しい**——
原文の第 1・第 2 の主張は完全に、第 3 の主張は (V) の読みで述べている。
★壊れているのは `Theorem 3.8` へ**繋ぐ**ところだけである。
★★`Lemma 3.5` も (V) の読みなので、そちらも節点 2 で言い直すことになる。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta

/-! ## ★★★★★★★★片側の含意は立つ -/

/-- ★★★★★★★★**`L` 有理な生成元があれば `Gal`-安定な直線がある**（含意の向きの記録）。

☆逆は偽である——安定な直線の点は `L` 有理とは限らない。
★本ファイルはこの向きを**主張として証明しない**（`HasLCyclicJ` は `galRep` の
言葉で書かれており、`HasLCyclicVelu` の `Q` から直線を作るには
`Q` の `Tate` 加群への持ち上げが要る）。☆節点として記録するにとどめる。 -/
def lCyclicReadingFrame.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(l-巡回の 2 つの読みが繋がっていない——Theorem 3.8 へ渡す枠)",
    sectionId := "genell-lemma-3-7" }

def lCyclicReadingFrame.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_7(第 1151、(V) の読みで証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_7") 1,
    .citation "[ABC3]" "imageContainsSL2J_of_alpha'(Theorem 3.8、残る仮説は α だけ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.imageContainsSL2J_of_alpha'") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1152）の測定**——原文の `H_L` を " ++
       "(V)「`L` 有理点 `Q` で位数 `l`」と読んだのは `Lemma 3.5` が Vélu の商を作るためであり、" ++
       "`Theorem 3.8` が要るのは (J)「`Gal`-安定な直線」である。" ++
       "☆(V) ⟹ (J) は立つが逆は偽なので、`Lemma 3.7` の第 3 の主張の対偶が " ++
       "`Theorem 3.8` の `hno` に届かない。" ++
       "★直す道は節点 1（安定な直線による Vélu の商が `Gal` 不変性で `L` に降りること）→ " ++
       "節点 2（`Lemma 3.5` を安定直線の側で述べ直す）→ 節点 3（`Lemma 3.7` の言い直し）。") 26 ]

end ABC3.Skeleton.GenEll
