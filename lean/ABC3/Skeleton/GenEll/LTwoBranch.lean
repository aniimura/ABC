/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateLocalModelK
import ABC3.Meta.Claim

/-!
# 第 1142 ブロック —— **`Lemma 3.5` から `l ≠ 2` を外す枠**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か

第 1141 で `d + 1 < l` が外れ、`Lemma 3.5` の原典に無い仮説は **`l ≠ 2` の 1 つだけ**に
なった。本ファイルはそれを外す残りを節点に割り、進捗枠を置く。

## ★★★★★★★★`l = 2` で何が壊れるか（第 1142 の実測）

`l = 2` では `(range 2).erase 0 = {1}` で `ζ = −1`、`w = −q` である。
★点 `Φ(−1)` は **2-捩れ**なので `P = −P`、すなわち

    `DX = 2Y + X = 0`,   `veluGy = −2Y − X = 0`,   `veluU = veluGy² = 0`

となる。☆壊れるのは次の 3 か所である。

| # | 壊れる場所 | 理由 |
|---|---|---|
| A | `veluV2_eq_tateDYpair`（`TateODE.lean:89`） | `DX ≠ 0` で割っている。`l = 2` では `DX = 0` |
| B | `exists_veluW_of_inv`（`SymmSum.lean:157`） | `i ↔ l−i` の対で偶数にする。`l = 2` では対にならない |
| C | `preΨ_eval_eq_zero_of_addOrderOf_prime` | `Odd k` を使う。`l = 2` では `preΨ₂ = 1` |

☆B は**易しい**——`l = 2` では `veluU = 0` なので `w = v·x` がそのまま `R` に取れる。
★C は `Ψ₂Sq`（主係数 `4`）に取り替える段である。
★★A が本体で、**一般の径数で通す**のが道である（下記）。

## ★★★★★★★★★★★★A の道——一般の径数 `T` で通す

`tate_ode_mul`（`DualTate.lean`）は **`DX·(DY − veluV2) = 0`** を無条件に与える。
`l = 2` では `DX = 0` なので情報が出ない。

☆しかし `DY = veluV2` は **`a` を径数と見れば恒等式**である。★道:

| 段 | 場所 | 何をするか |
|---|---|---|
| 1 | `S ≔ K[T]` を `T` と `1 − T` で局所化 | `1 − T` が単元、`T ≠ 0`、`1 + T ≠ 0` |
| 2 | `PowerSeries S` で `a ≔ C T`・`w ≔ C(T⁻¹)·X` | `DX` の定数項は `T(1+T)/(1−T)³ ≠ 0` なので **`hDX` が出る** |
| 3 | 在庫の `veluV2_eq_tateDYpair` を当てる | `hu`・`hDX` とも成り立つ |
| 4 | `T ↦ −1` で特殊化する | `1 − (−1) = 2 ≠ 0`・`−1 ≠ 0` なので像は定義される |

★これは第 1128 と同じ「万有な環を経由する」型であり、道具（`map_*`・`evalAdicMapHom`）は
すべて第 1125-1127 で建ててある。

## ★★★★★★★★★★残り 4 節点（進捗枠 **2 / 4**）

| # | 節点 | 内容 | 重み |
|---|---|---|---|
| 1 | `veluV2_eq_tateDYpair_any` | `hDX` を取らない `DY = veluV2` **★第 1143-1144 で証明済み** | 10 |
| 2 | `c4_velu_tate_any` / `c6_velu_tate_any` | `hDX` を取らない `c₄`・`c₆` の恒等式 **★第 1145-1146 で証明済み** | 8 |
| 3 | `exists_veluW_two` | `l = 2` では `veluU = 0` なので `w = v·x` | 6 |
| 4 | `psi2Sq_eval_eq_zero` | `l = 2` の捩れ点と多項式の橋（`Ψ₂Sq`） | 8 |

☆総重み 32。★これが閉じれば `Lemma 3.5` は**原典どおりの仮説だけ**になる。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta

/-! ## ★出典の紐付け(`.src`) -/

def lTwoBranchFrame.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l = 2 の枝を外す枠——残る唯一の逸脱)",
    sectionId := "genell-lemma-3-5" }

def lTwoBranchFrame.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_ode_mul(DX·(DY − veluV2) = 0。★無条件、在庫)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tate_ode_mul") 1,
    .citation "[ABC3]" "evalAdicMapHom(S⟦q⟧ →+* R、第 1126、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.evalAdicMapHom") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1142）の測定**——`l = 2` で壊れるのは 3 か所である: " ++
       "(A) `veluV2_eq_tateDYpair` が `DX ≠ 0` で割っている、" ++
       "(B) `exists_veluW_of_inv` の `i ↔ l−i` の対、" ++
       "(C) `preΨ` の `Odd k`。" ++
       "☆B は易しい（`l = 2` では `veluU = 0`）。★A が本体で、" ++
       "`a` を一般の径数 `T` に取れば `DX ≠ 0` になるので在庫がそのまま効き、" ++
       "`T ↦ −1` で特殊化すればよい。") 10 ]

end ABC3.Skeleton.GenEll
