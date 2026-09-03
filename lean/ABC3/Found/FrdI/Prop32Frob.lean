/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Frob.Assemble

/-!
# Prop32Frob（取りまとめ）

★★★★**2026-09-03（第 1457）——Prop32Frob を 6 枚に割った**。

☆割る軸は**ファイル内の見出し**である。論文のセクションでは割れない
——大きいファイルは話題が混ざっているのではなく、1 つの項目に Lean が大量に要ったからである。
★入れ子の名前空間・section を跨がない切れ目だけを使い、`.src`／`.needs` は
対応する宣言のある枚へ配ってある（G1 は同じファイルを見るため）。

☆本ファイルは既存の import を壊さないための取りまとめである。
必要な部分だけを引きたいときは `ABC3.Found.FrdI.Prop32Frob.<名前>` を直接 import する。

| 枚 | 中身 |
|---|---|
| `Triple` | 3 つ組の添字対象・合成則・恒等射・同型 |
| `Transition` | 根を上げる道具・遷移は同型を保つ・pre-step |
| `Faithful` | (iii)(c) の順逆・3 脚添字 |
| `Functor` | 単系の補題・(vi) 単元を除く忠実性 |
| `Pullback` | `Λ_k` は関手・(iv)(b)・根が違う合成の組み立て |
| `Assemble` | `𝒞` の pull-back・底の比較・出典 |
-/
