/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateIsogeny.Discriminant

/-!
# Tate 同種（取りまとめ）

★★★★**2026-09-03（第 1457）——147 KB / 81 宣言 / 2,534 行を 6 枚に割った**。

☆割る軸は**ファイル内の見出し**である。論文のセクションでは割れない
——このファイルは [GenEll] §3 の 2 項目しか持たないからである。
★後方参照は 2 件（冒頭の docstring からの言及）だけで、実コードの依存は順方向であった。

☆本ファイルは既存の import を壊さないための取りまとめである。
必要な部分だけを引きたいときは `ABC3.Skeleton.GenEll.TateIsogeny.<名前>` を直接 import する。

| 枚 | 中身 |
|---|---|
| `Target` | 訂正後の目標（`c₄`・`c₆`） |
| `LocalJ` | `ζ` を消す・`j` で受ける捉れ点版・局所の終点 |
| `GlobalVelu` | `j` の段で止める・大域の Vélu の商・悪い素点の局所データ |
| `LocalData` | 極小化した対・`hvw` を作る |
| `Extension` | 一般の局所体に開く・拡大に載せ替える・非分裂の枝 |
| `Discriminant` | 局所の `Δ` の計算・Vélu の商の整性 |
-/
