---
name: inhabited-two-paths
description: Unique から Inhabited を取る道が 2 本あり、default が一致せず convert/Subsingleton.elim が止まる。補助定義を [Unique ...] で引数化して揃える
metadata:
  type: feedback
---

`ABC3` で `Unique` から `default` を使うとき、**`Unique.instInhabited` と
`Unique.toInhabited` の 2 経路**があり、同じ `Unique` から来ていても
`default` の項が一致しない(2026-08-24、第 107 ブロックで 20 手以上詰まった)。

★症状: `convert` が `HEq` や深い型等式を出し、`Subsingleton.elim` が
`Subsingleton (Type u)` や `Subsingleton (ModuleCat ...)` を要求して止まる。

**Why:** `default` は `Inhabited` インスタンスに依存する項なので、
インスタンスの導出経路が違うと構文的に別物になる。

**How to apply:**
- 補助定義を **`[Unique ι]` で引数化**し、呼ぶ側の goal が使うのと**同じ class** にする。
  `[Inhabited ι]` で引数化すると経路が分かれるので**駄目**。
- 型を固定する補助定義(`restrictSec` 型)を作ると `HEq` が同次の等式になり、
  `congr 1` まで落ちる。そこから先はインスタンスを揃えるだけ。

関連: [[ring-instance-two-paths]](環)、加群(§9-110)、本項(`Inhabited`)で 3 例目。
