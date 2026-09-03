/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop16.Core

/-!
# Prop16（取りまとめ）

★★★★**2026-09-03（第 1457）——134 KB / 99 宣言 / 2,377 行を 3 枚に割った**。

☆`section Dict` が 2,000 行を覆っているので、入れ子の外では 13 KB / 119 KB にしか割れない。
そこで **`section Dict` を各枚で開き直し**、その `variable` を持たせて割ってある。

| 枚 | 中身 |
|---|---|
| `Restrict` | `Φ` の `𝒟'` への制限・(i) の `𝔽_{Φ'}` |
| `Dict` | 辞書——`CfpCat` の Frobenioid 構造の互換性 |
| `Core` | Frobenioid core と (iv) 以降 |
-/
