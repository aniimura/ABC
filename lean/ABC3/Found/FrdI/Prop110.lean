/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop110.Compact

/-!
# Prop110（取りまとめ）

★★★★**2026-09-03（第 1457）——Prop110 を 5 枚に割った**。

☆割る軸は**ファイル内の見出し**である。論文のセクションでは割れない
——大きいファイルは話題が混ざっているのではなく、1 つの項目に Lean が大量に要ったからである。
★入れ子の名前空間・section を跨がない切れ目だけを使い、`.src`／`.needs` は
対応する宣言のある枚へ配ってある（G1 は同じファイルを見るため）。

☆本ファイルは既存の import を壊さないための取りまとめである。
必要な部分だけを引きたいときは `ABC3.Found.FrdI.Prop110.<名前>` を直接 import する。

| 枚 | 中身 |
|---|---|
| `Vocabulary` | §0 の語彙・(i) の一意性と存在 |
| `PreStep` | (ii) の pre-step と Frobenius・`Div` の公式・(iii) |
| `Sharp` | sharp モノイド・(iv) の逆向き・(vi) の前半 |
| `Moreover` | (iii) の moreover・穴を浅くする段 |
| `Compact` | (iv) の In particular・(vi) 後半・Frobenius-compact |
-/
